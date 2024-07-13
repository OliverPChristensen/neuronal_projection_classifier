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
    region_names <- region_names %>% tolower()
    region_names <- case_when(
        region_names == 'dorsal_hyothalamus' ~ 'dorsal_hypothalamus',
        region_names == 'anterior_hyopthalamus_left' ~ 'anterior_hypothalamus_left',
        region_names == 'interpeduncular' ~ 'interpreduncular',
        region_names == 'lateral_hypothalamus_right-1' ~ 'lateral_hypothalamus_right',
        region_names == 'pontine_grey' ~ 'pontine_gray',
        region_names == 'retrochiasmaric_left' ~ 'retrochiasmatic_left',
        region_names == 'suprachiamsatic' ~ 'suprachiasmatic',
        TRUE ~ region_names
    )

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

    #seurat_to_sf
    # Already projection termination if preprocessed object is used (And Later I will make this contain projection_temrination)
    #Remember to change create_shiny_data accrdongly when projection termination and correct region names have been added in preprocessing

    dir.create(paste0(path_to_output_folder),showWarnings = FALSE)
    dir.create(paste0(path_to_plot_folder),showWarnings = FALSE)

    cat(paste0(utils$current_time()," Saving all output files to ",path_to_output_folder,"\n"))
    cat(paste0(utils$current_time()," Saving all plots to ",path_to_plot_folder,"\n"))

    cat(paste0(utils$current_time()," Loading spatial object from ",spatial_object_path,"\n"))
    spatial_object <- qread(spatial_object_path)

    cell_object <- sf_utils$seurat_to_sf(spatial_object)

    run_indices <- cell_object$run_index %>% unique()
    section_indices <- cell_object$section_index %>% unique()

    run_section_index_pairs <- expand.grid(run_indices, section_indices)

    cat(paste0(utils$current_time()," Loading region object\n"))
    region_object <- sf_utils$load_regions(region_paths, run_section_index_pairs)

    regions <- cell_object$region_corrected %>% na.omit() %>% unique()
    projections <- cell_object %>% filter(projection != "Other") %>% .[["projection"]] %>% unique()

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
            filter(run_index == run_index_, section_index == section_index_, projection == projection_, region_corrected == region_) %>% 
            st_centroid()

        if (nrow(cell_object_sub_centroid) > 10){
            region_object_sub <- region_object %>%
                filter(run_index == run_index_, section_index == section_index_, region == region_)

            cell_object_ppp <- as.ppp(st_coordinates(cell_object_sub_centroid), as.owin(region_object_sub))

            # Choose grid size
            cell_object_qq <- quadratcount(cell_object_ppp, nx = 10, ny = 10)

            chi2_p_value <- quadrat.test(cell_object_qq, alternative = "clustered")$p.value %>% 
                data.frame(p_value = ., region = region_, projection = projection_, section_index = section_index_)
            
            chi2_p_values <- rbind(chi2_p_values, chi2_p_value)

            cluster_bool <- get_k_func_cluster_bool(cell_object_ppp, test_distance = test_distance, n_sim = n_sim)

            K_func_list <- rbind(K_func_list, data.frame(clustering = cluster_bool, region = region_, projection = projection_, section_index = section_index_))
        

            if (chi2_p_value$p_value < 0.05 & cluster_bool){
                lambdahat <- density(cell_object_ppp)

                save_density_plot(lambdahat, paste0(path_to_plot_folder, "/intensity_estimate_", region_, "_", projection_, "_plot.png"))
                

                #INLA later. With it's own function
            }
        }

        utils$loop_progress(row_num, nrow(region_projection_index_pairs))

    }

    cat(paste0(utils$current_time()," Generate plots for chi2 p-values\n"))
    plot_chi2_p_values(chi2_p_values, path_to_plot_folder)

    cat(paste0(utils$current_time()," Generate plots for K function results\n"))
    plot_K_func_list(K_func_list, path_to_plot_folder)
    
}

main()
