library(tidyverse)
library(Seurat)
library(qs)
library(argparse)
library(spatstat)

utils <- new.env()
source("scripts/utils/utils.R", local = utils)

plot_utils <- new.env()
source("scripts/utils/plot_utils.R", local = plot_utils)

sf_utils <- new.env()
source("scripts/utils/sf_utils.R", local = sf_utils)


correct_region_names <- function(region_names){
    # Input: Character vector of region names
    #
    # Takes annotated region names and corrects them, by setting all characters to lower case and correcting typos
    #
    # Output: Character vector of corrected region names

    # Set all region names to lower case
    region_names <- region_names %>% tolower()

    # Correct all region names for typos
    region_names <- case_when(
        region_names == 'dorsal_hyothalamus' ~ 'dorsal',
        region_names == 'anterior_hyopthalamus_left' ~ 'anterior_hypothalamus_left',
        region_names == 'interpeduncular' ~ 'interpreduncular',
        region_names == 'lateral_hypothalamus_right-1' ~ 'lateral_hypothalamus_right',
        region_names == 'pontine_grey' ~ 'pontine_gray',
        region_names == 'retrochiasmaric_left' ~ 'retrochiasmatic_left',
        region_names == 'suprachiamsatic' ~ 'suprachiasmatic',
        region_names == 'dorsomedial_hypothalamus' ~ 'dorsomedial',
        region_names == 'dorsal_hypothalamus' ~ 'dorsal',
        TRUE ~ region_names
    )

    names(region_names) <- NULL

    return(region_names)
}

save_density_plot <- function(lambdahat,plot_path){

    # Save the plot as a PNG file
    png(plot_path)

    # Plot the density
    plot(lambdahat)

    # Close the graphics device
    dev.off()
}

get_k_func_cluster_bool <- function(cell_object_ppp, test_distance, n_sim){

    cell_object_K_func <- Kest(cell_object_ppp)
    cell_object_dist <- envelope(cell_object_ppp, Kest, nsim = n_sim, verbose = FALSE)


    r_index <- which(cell_object_dist$r > test_distance)[1]

    cluster_bool <- cell_object_dist$hi[r_index] < cell_object_dist$obs[r_index]

    return(cluster_bool)
}


        

plot_chi2_p_values <- function(chi2_p_values, plot_folder){

    chi2_p_values <- chi2_p_values %>% 
        mutate(inv_p_value = -log(p_value))
    
    chi2_p_values_dot_plot <- chi2_p_values %>%
        ggplot() + geom_point(aes(x = projection, y = section_index, size = inv_p_value, color = inv_p_value > -log(0.05))) + 
        facet_grid(~region) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        scale_color_manual(values = c("TRUE" = "black", "FALSE" = "transparent"))

    ggsave(paste0(plot_folder,"/chi2_p_values_dot_plot.png"),chi2_p_values_dot_plot, width = 12000, height = 2000, units = "px")

}

plot_K_func_list <- function(K_func_list, plot_folder){

    K_func_list <- K_func_list %>% 
        filter(!is.na(clustering))

    K_func_list_dot_plot <- K_func_list %>% 
        ggplot() + geom_point(aes(x = projection, y = section_index, color = clustering)) + 
        facet_grid(~region) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        scale_color_manual(values = c("TRUE" = "black", "FALSE" = "transparent"))

    ggsave(paste0(plot_folder,"/K_func_list_dot_plot.png"),K_func_list_dot_plot, width = 12000, height = 2000, units = "px")

}



calculate_and_plot_intensity <- function(cell_object, region_object, chi2_p_values, K_func_list, plot_folder, scale = 10){
    run_indices <- chi2_p_values$run_index %>% unique()
    section_indices <- chi2_p_values$section_index %>% unique()
    chosen_regions <- cell_object$region_corrected %>% unique()

    for (run_index_ in run_indices){
        for (section_index_ in section_indices){
            projections <- chi2_p_values %>% filter(run_index == run_index_, section_index == section_index_, p_value < 0.05) %>% .[["projection"]] %>% unique()
            
            for (projection_ in projections){
                positive_regions_chi2 <- chi2_p_values %>% filter(run_index == run_index_, section_index == section_index_, projection == projection_, p_value < 0.05) %>%
                    .[["region"]]

                positive_regions_k_func <- K_func_list %>% filter(run_index == run_index_, section_index == section_index_, projection == projection_, clustering == TRUE) %>%
                    .[["region"]]
                
                positive_regions <- intersect(positive_regions_chi2, positive_regions_k_func)

                region_object_sub <- region_object %>% filter(run_index == run_index_, section_index == section_index_, region_corrected %in% chosen_regions)

                cell_object_sub_centroid <- cell_object %>% 
                    filter(run_index == run_index_, section_index == section_index_, single_projection == projection_) %>% 
                    st_centroid()

                plot <- ggplot()
                for (region_ in positive_regions){
                    cell_object_sub_centroid_region <- cell_object_sub_centroid %>% 
                        filter(region_corrected == region_)

                    region_object_sub_region <- region_object_sub %>% 
                        filter(region_corrected == region_)

                    cell_object_ppp <- as.ppp(st_coordinates(cell_object_sub_centroid_region), as.owin(region_object_sub_region))

                    xrange <- cell_object_ppp$window$xrange
                    yrange <- cell_object_ppp$window$yrange

                    scale <- 10

                    density_object <- density(cell_object_ppp, dimyx = c((xrange[2] - xrange[1])/scale, (yrange[2] - yrange[1])/scale)) %>% 
                        as.data.frame()


                    plot <- plot + geom_tile(data = density_object,  mapping = aes(x, y, fill=value))

                }

                custom_colors <- c("blue", "red", "yellow")

                plot <- plot + scale_fill_gradientn(colors = custom_colors) +
                    geom_sf(data = cell_object_sub_centroid, size = 0.2, color = "black") + 
                    geom_sf(data = region_object_sub, fill = NA) + 
                    geom_sf_text(data = region_object_sub, aes(label = region_corrected), color = "black", size = 2, fontface = "bold") +
                    theme_classic() + 
                    theme(
                        axis.line = element_blank(),          # Remove axis lines
                        axis.title = element_blank(),         # Remove axis titles
                        axis.text = element_blank(),          # Remove axis text (numbers)
                        axis.ticks = element_blank()          # Remove axis ticks (optional)
                    )
                
                ggsave(paste0(plot_folder, "/intensity_estimate_",run_index_,"_", section_index_,"_",projection_,".png"), plot)
            }
            

        }

    }
}

create_parser <- function(){
    # Input: NA
    #
    # Creates and outputs parser with arguments that allows interaction via CLI
    #
    # Output: Parser

    # Define parser
    parser <- ArgumentParser(description = "Script to generate spatial object from ROIs and transcripts data file")

    # Add argument to parser
    parser$add_argument("-o","--output-folder-path", type = "character", help = "Relative path to folder to save output")
    parser$add_argument("-p","--plot-folder-path", type = "character", help = "Relative path to folder to save plots")
    parser$add_argument("-sop","--spatial-object-path", type = "character", help="Relative path to spatial object in seurat format")
    parser$add_argument("-rps","--region-paths", nargs="+", type = "character", help="Relative paths to folders with region ROIs in imageJ format")
    parser$add_argument("-rdp","--rotation-df-path", type = "character", help="Relative path to data frame containing rotation angles for spatial sections")

    return(parser)
}



main <- function(){
    # Input: NA
    #
    # Main function that defines parser argument, loads transcripts and ROIs, creates count matrix, creates spatial seurat object, saves the spatial seurat object
    #
    # Output: NA

    # Create parser and parse arguments
    parser <- create_parser()
    args <- parser$parse_args()

    # Retrieve arguments into corresponding variables
    path_to_output_folder <- args$output_folder_path
    path_to_plot_folder <- args$plot_folder_path
    spatial_object_path <- args$spatial_object_path
    region_paths <- args$region_paths
    rotation_df_path <- args$rotation_df_path

    #seurat_to_sf
    # Already projection termination if preprocessed object is used (And Later I will make this contain projection_temrination)
    #Remember to change create_shiny_data accrdongly when projection termination and correct region names have been added in preprocessing

    dir.create(paste0(path_to_output_folder),showWarnings = FALSE)
    dir.create(paste0(path_to_plot_folder),showWarnings = FALSE)

    cat(paste0(utils$current_time()," Saving all output files to ",path_to_output_folder,"\n"))
    cat(paste0(utils$current_time()," Saving all plots to ",path_to_plot_folder,"\n"))

    cat(paste0(utils$current_time()," Loading spatial object from ",spatial_object_path,"\n"))
    spatial_object <- qread(spatial_object_path)

    rotation_df <- read.table(rotation_df_path, sep = "\t", header = TRUE)

    cell_object <- sf_utils$seurat_to_sf(spatial_object)

    cat(paste0(utils$current_time()," Saving cell object to ",path_to_output_folder,"\n"))
    qsave(cell_object, paste0(path_to_output_folder,"/cell_object.qs"))

    #cell_object <- qread("processed_data/cell_object.qs")

    # Remove cells in injection site
    cell_object <- cell_object[!cell_object$in_injection_site,]

    cell_object <- cell_object %>% filter(region_reduced != "Other" & single_projection != "Other")

    cell_object <- sf_utils$rotate_sf_object(cell_object, rotation_df)

    run_indices <- cell_object$run_index %>% unique()
    section_indices <- cell_object$section_index %>% unique()

    run_section_index_pairs <- expand.grid(run_indices, section_indices)

    ####




    region_paths <- c()
    for (run_index in run_indices){
        region_paths <- c(region_paths,paste0(paste0("raw_data/",run_index),"/",list.files(path = paste0("raw_data/",run_index), pattern = "regions")))
    }


    ####

    cat(paste0(utils$current_time()," Loading region object\n"))
    region_object <- sf_utils$load_regions(region_paths, run_section_index_pairs)

    region_object <- sf_utils$rotate_sf_object(region_object, rotation_df)

    region_object$region_corrected <- correct_region_names(region_object$region_name)

    regions <- cell_object$region_corrected %>% na.omit() %>% unique()
    projections <- cell_object %>% .[["single_projection"]] %>% unique()

    region_projection_index_pairs <- expand.grid(regions, projections, run_indices, section_indices)

    chi2_p_values <- data.frame()
    K_func_list <- data.frame()
    
    test_distance <- 100
    n_sim <- 19

    cat(paste0(utils$current_time()," Perform clustering analysis on all region-projection-section pairs with test_distance ", test_distance," and n_sim ", n_sim,":\n"))

    for (row_num in 1:nrow(region_projection_index_pairs)){
        region_ <- region_projection_index_pairs[row_num, "Var1"]
        projection_ <- region_projection_index_pairs[row_num, "Var2"]
        run_index_ <- region_projection_index_pairs[row_num, "Var3"]
        section_index_ <- region_projection_index_pairs[row_num, "Var4"]

        cell_object_sub_centroid <- cell_object %>% 
            filter(run_index == run_index_, section_index == section_index_, single_projection == projection_, region_corrected == region_) %>% 
            st_centroid()

        if (nrow(cell_object_sub_centroid) > 10){
            region_object_sub <- region_object %>%
                filter(run_index == run_index_, section_index == section_index_, region_corrected == region_)

            cell_object_ppp <- as.ppp(st_coordinates(cell_object_sub_centroid), as.owin(region_object_sub))

            # Choose grid size
            cell_object_qq <- quadratcount(cell_object_ppp, nx = 10, ny = 10)

            chi2_p_value <- quadrat.test(cell_object_qq, alternative = "clustered")$p.value %>% 
                data.frame(p_value = ., region = region_, projection = projection_, section_index = section_index_, run_index = run_index_)
            
            chi2_p_values <- rbind(chi2_p_values, chi2_p_value)

            cluster_bool <- get_k_func_cluster_bool(cell_object_ppp, test_distance = test_distance, n_sim = n_sim)

            K_func_list <- rbind(K_func_list, data.frame(clustering = cluster_bool, region = region_, projection = projection_, section_index = section_index_, run_index = run_index_))
        

            #if (chi2_p_value$p_value < 0.05 & cluster_bool){
            #    lambdahat <- density(cell_object_ppp)

            #    save_density_plot(lambdahat, paste0(path_to_plot_folder, "/intensity_estimate_", region_, "_", projection_, "_plot.png"))
                

                #INLA later. With it's own function
            #}
        }

        utils$loop_progress(row_num, nrow(region_projection_index_pairs))

    }

    cat(paste0(utils$current_time()," Calculate and plot density estimate\n"))
    calculate_and_plot_intensity(cell_object, region_object, chi2_p_values, K_func_list, path_to_plot_folder, scale = 10)

    cat(paste0(utils$current_time()," Generate plots for chi2 p-values\n"))
    plot_chi2_p_values(chi2_p_values, path_to_plot_folder)

    cat(paste0(utils$current_time()," Generate plots for K function results\n"))
    plot_K_func_list(K_func_list, path_to_plot_folder)
    
}

main()
