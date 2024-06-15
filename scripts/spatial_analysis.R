library(tidyverse)
library(Seurat)
library(qs)
library(argparse)


current_time <- function(){
    # Input: NA
    #
    # Extract and format current time
    #
    # Output: Current time

    # Extract and format current time
    return(paste0("[",format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),"]"))
}


correct_region_names <- function(spatial_object){
    region_lower <- spatial_object$region %>% tolower()

    region_lower_corrected <- case_when(
        region_lower == 'dorsal_hyothalamus' ~ 'dorsal_hypothalamus',
        region_lower == 'anterior_hyopthalamus_left' ~ 'anterior_hypothalamus_left',
        region_lower == 'interpeduncular' ~ 'interpreduncular',
        region_lower == 'lateral_hypothalamus_right-1' ~ 'lateral_hypothalamus_right',
        region_lower == 'pontine_grey' ~ 'pontine_gray',
        region_lower == 'retrochiasmaric_left' ~ 'retrochiasmatic_left',
        region_lower == 'suprachiamsatic' ~ 'suprachiasmatic',
        TRUE ~ region_lower
    )

    spatial_object$region_corrected <- region_lower_corrected %>% gsub("_right|_left", "", .)

    return(spatial_object)
}

add_projection_termination <- function(spatial_object, threshold){

    spatial_object$projection_termination_bonferroni <- case_when(
        spatial_object$virus1 < threshold & spatial_object$virus2 > threshold  & spatial_object$run_index == "spatial1" ~ "BST",
        spatial_object$virus2 < threshold & spatial_object$virus1 > threshold & spatial_object$run_index == "spatial1" ~ "PAG",
        spatial_object$virus1 < threshold & spatial_object$virus2 < threshold & spatial_object$run_index == "spatial1" ~ "BST+PAG",
        spatial_object$virus1 < threshold & spatial_object$virus2 > threshold  & spatial_object$run_index == "spatial2" ~ "LHA",
        spatial_object$virus2 < threshold & spatial_object$virus1 > threshold & spatial_object$run_index == "spatial2" ~ "PBN",
        spatial_object$virus1 < threshold & spatial_object$virus2 < threshold & spatial_object$run_index == "spatial2" ~ "LHA+PBN",
        spatial_object$virus1 < threshold & spatial_object$virus2 > threshold  & spatial_object$run_index == "spatial3" ~ "PVT",
        spatial_object$virus2 < threshold & spatial_object$virus1 > threshold & spatial_object$run_index == "spatial3" ~ "CEA",
        spatial_object$virus1 < threshold & spatial_object$virus2 < threshold & spatial_object$run_index == "spatial3" ~ "PVT+CEA",
        TRUE ~ "No-virus"
    )

    return(spatial_object)
}

region_umap <- function(spatial_object){
    region_umap <- DimPlot(spatial_object, group.by = 'region_corrected', reduction = 'umap_features', raster = FALSE)

    ggsave("plots/spatial_analysis/region_umap.png", region_umap, width = 5000, height = 3000, units = "px")

    spatial_object$region_corrected_reduced <- case_when(

    )

    region_reduced_umap <- DimPlot(spatial_object, group.by = 'region_corrected_reduced', reduction = 'umap_features', raster = FALSE)

    ggsave("plots/spatial_analysis/region_reduced_umap_umap.png", region_reduced_umap, width = 5000, height = 3000, units = "px")

}

projection_umaps <- function(spatial_object){

    projection_termination_no_double <- case_when(
        grepl("\\+", spatial_object$projection_termination) ~ "double",
        TRUE ~ spatial_object$projection_termination
    )

    names(projection_termination_no_double) <- NULL

    spatial_object$projection_termination_no_double <- projection_termination_no_double

    projection_termination_no_double_umap <- DimPlot(spatial_object, group.by = 'projection_termination_no_double', reduction = 'umap_features', raster = FALSE, order = TRUE)

    ggsave("plots/spatial_analysis/projection_termination_no_double_umap.png", projection_termination_no_double_umap, width = 4000, height = 3000, units = "px")


    spatial_object$projection_all <- case_when(
        spatial_object$sun1 == TRUE ~ TRUE,
        spatial_object$virus1 < 0.05 ~ TRUE,
        spatial_object$virus2 < 0.05 ~ TRUE,
        TRUE ~ FALSE
    )

    projection_all_umap <- DimPlot(spatial_object, group.by = 'projection_all', reduction = 'umap_features', raster = FALSE, order = TRUE)

    ggsave("plots/spatial_analysis/projection_all_umap.png", projection_all_umap, width = 4000, height = 3000, units = "px")

    spatial_object$projection_sun1_vs_virus <- case_when(
        spatial_object$sun1 == TRUE ~ "Sun1",
        spatial_object$virus1 < 0.05 ~ "Virus1",
        spatial_object$virus2 < 0.05 ~ "Virus2",
        spatial_object$sun1 == TRUE & (spatial_object$virus1 < 0.05 | spatial_object$virus2 < 0.05) ~ "Sun1-Virus",
        TRUE ~ "Non-projecting"
    )

    projection_sun1_vs_virus_umap <- DimPlot(spatial_object, group.by = 'projection_sun1_vs_virus', reduction = 'umap_features', raster = FALSE, order = TRUE)

    ggsave("plots/spatial_analysis/projection_sun1_vs_virus_umap.png", projection_sun1_vs_virus_umap, width = 4000, height = 3000, units = "px")

}

area_feature_plot <- function(spatial_object){
    area_umap <- FeaturePlot(spatial_object, feature = "area", reduction = "umap_features", raster = FALSE)

    ggsave("plots/spatial_analysis/area_umap.png", area_umap, width = 4000, height = 3000, units = "px")
}

sun1_virus_across_sections <- function(spatial_object, threshold, plot_path){

    spatial_object$sun1_vs_virus_bonferroni <- case_when(
        spatial_object$sun1 == TRUE & spatial_object$virus1 > threshold & spatial_object$virus2 > threshold ~ "Sun1",
        spatial_object$sun1 == FALSE & (spatial_object$virus1 < threshold | spatial_object$virus2 < threshold) ~ "Virus",
        spatial_object$sun1 == TRUE & (spatial_object$virus1 < threshold | spatial_object$virus2 < threshold) ~ "Sun1+Virus",
        spatial_object$sun1 == FALSE & spatial_object$virus1 > threshold & spatial_object$virus2 > threshold ~ "Non-projecting"
    )

    sun1_vs_virus_counts <- spatial_object@meta.data %>% 
        group_by(run_index, section_index, sun1_vs_virus) %>% 
        summarize(count = n()) %>% mutate(total_count = sum(count),relative_count = count / total_count)
    
    p_values <- sun1_vs_virus_counts %>% 
        select(-c(total_count,relative_count)) %>% 
        pivot_wider(names_from = sun1_vs_virus, values_from = count) %>%
            mutate(p_value = fisher.test(matrix(c(`Non-projecting`,Sun1,`Sun1+Virus`,Virus),nrow = 2))$p.value)


    sun1_vs_virus_bar_plot <- sun1_vs_virus_counts %>%
        ggplot() + geom_col(aes(x = sun1_vs_virus, y = relative_count)) + facet_grid(section_index~run_index) + 
        geom_text(aes(x = sun1_vs_virus, y = relative_count + 0.1, label = count)) + 
        geom_text(data = p_values, aes(x = 2.5, y = Inf, label = format(p_value, scientific = TRUE, digits = 2)), hjust = 0.5, vjust = 1.2, inherit.aes = FALSE, size = 5)


    ggsave(plot_path, sun1_vs_virus_bonferroni_bar_plot)
}

projection_termination_plots <- function(spatial_object, projection_termination_column_name, plot_path){

    projection_terminations <- spatial_object@meta.data %>% 
        group_by(run_index, section_index, !!sym(projection_termination_column_name)) %>%
        summarize(count = n()) %>% 
        filter(!!sym(projection_termination_column_name) != "No-virus")

    projection_terminations <- projection_terminations %>% 
        mutate(!!sym(projection_termination_column_name) := factor(
            !!sym(projection_termination_column_name), 
            levels = c("BST","PAG", "BST+PAG", "LHA", "PBN", "LHA+PBN", "PVT", "CEA", "PVT+CEA")
        ))

    projection_terminations_plot <- projection_terminations %>%
        ggplot() + geom_col(aes(x = !!sym(projection_termination_column_name), y = count)) + facet_grid(section_index~run_index, scales = "free_x") + 
        geom_text(aes(x = !!sym(projection_termination_column_name), y = count + 100, label = count))

    ggsave(plot_path, projection_terminations_plot)
}

spatial_projection_plots <- function(spatial_object,plot_path){

    spatial_object$virus <- case_when(
        spatial_object$virus1 < 0.05 & spatial_object$virus2 > 0.05 ~ "virus1",
        spatial_object$virus1 > 0.05 & spatial_object$virus2 < 0.05 ~ "virus2",
        spatial_object$virus1 > 0.05 & spatial_object$virus2 > 0.05 ~ "No-virus",
        spatial_object$virus1 < 0.05 & spatial_object$virus2 < 0.05 ~ "virus1+virus2",
    )



    levels <- unique(spatial_object$virus)


    num_levels <- spatial_object$projection_termination_no_double %>% unique() %>% length()


    levels <- unique(spatial_object$virus)
    default_colors <- hue_pal()(length(levels))
    names(default_colors) <- levels

    default_colors["No-virus"] <- "#7F7F7F"


    run_indices <- c("spatial1", "spatial2", "spatial3")
    section_indices <- c("A1", "B1", "C1", "D1", "A2", "B2", "C2", "D2")

    for (run_index_number in 1:3){
        for (section_index in section_indices){
            run_index_number_adjusted <- ifelse(run_index_number == 1, "", paste0(".",run_index_number))
            fov <- paste0(section_index,run_index_number_adjusted)

            DefaultBoundary(spatial_object[[fov]]) <- "segmentation"
            spatial_plot <- ImageDimPlot(spatial_object, 
                    fov = fov,
                    axes = F, 
                    dark.background = T, 
                    size = 1.5,
                    border.color = "black", 
                    border.size = 0.1, 
                    cols = default_colors, 
                    group.by = "virus",
                    molecules = c("p147_mCherry_2A_iCre","p36_NLS_Cre"),
                    nmols = 2000000,
                    mols.size = 0.05,
                    mols.alpha = 0.1
            )

            ggsave(paste0("./plots/spatial_projection_plot_",run_indices[i],"_",section_index,".png"),plot = spatial_plot, width = 6000, height = 7500, units = "px")
        }
    }
}

check_background_intensity <- function(spatial_object){

    run_indices <- c("spatial1", "spatial2", "spatial3")
    section_indices <- c("A1", "B1", "C1", "D1", "A2", "B2", "C2", "D2")

    intensities <- spatial_object@meta.data %>% group_by(run_index, section_index, virus1_count, virus2_count, section_area) %>% summarize(total_cell_area = sum(area))
    intensities$total_count_virus1_negative_cells <- 0
    intensities$total_count_virus2_negative_cells <- 0

    for (run_index_number in 1:3){
        for (section_index_number in 1:8){
            negative_cells <- spatial_object@meta.data %>% filter(run_index == run_indices[run_index_number], section_index == section_indices[section_index_number], sun1 == FALSE, virus1 > 0.05, virus2 > 0.05) %>% rownames()
            run_section_pair_boolean <- intensities$run_index == run_indices[run_index_number] & intensities$section_index == section_indices[section_index_number]
            layer <- paste0("counts.",section_index_number,".",run_index_number)

            intensities[run_section_pair_boolean ,"total_count_virus1_negative_cells"] <- spatial_object@assays$raw[layer]["p147-mCherry-2A-iCre",negative_cells] %>% sum()
            intensities[run_section_pair_boolean ,"total_count_virus2_negative_cells"] <- spatial_object@assays$raw[layer]["p36-NLS-Cre",negative_cells] %>% sum()

        }
    }

    intensities <- intensities %>% mutate(
            intensity_background_virus1 = virus1_count/section_area, 
            intensity_background_virus2 = virus2_count/section_area,
            intensity_negative_cells_virus1 = total_count_virus1_negative_cells/total_cell_area,
            intensity_negative_cells_virus2 = total_count_virus2_negative_cells/total_cell_area
        )
    
    intensities$intensity_background_virus1_scaled <- scale_intensity

    intensities_long <- intensities %>% pivot_longer(
            cols = c(intensity_background_virus1,intensity_background_virus2,intensity_negative_cells_virus1,intensity_negative_cells_virus2),
            names_to = "type",
            values_to = "intensity"
        )

    # TODO: Nedskalere background intensity s%>% den er sammenlignelig med cell intensity
    
    check_background_intensity_plot <- intensities_long %>% ggplot() + geom_col(aes(x = type, y = intensity, fill = type)) + facet_grid(section_index~run_index) + coord_flip()

    ggsave("plots/spatial_analysis/check_background_intensity_plot.png", check_background_intensity_plot)
}

fisher_region_cell_property_enrichment <- function(regions,cell_properties){

    region_unique <- unique(na.omit(regions))
    cell_property_unique <- unique(na.omit(cell_properties))

    enrichment_matrix <- matrix(0, nrow = length(region_unique), ncol = length(cell_property_unique))

    rownames(enrichment_matrix) <- sort(region_unique)
    colnames(enrichment_matrix) <- sort(cell_property_unique)

    for (region in region_unique){
        for (cell_property in cell_property_unique){
            regions_test <- case_when(regions == region ~ region, TRUE ~ 'Others')
            cell_properties_test <- case_when(cell_properties == cell_property ~ cell_property, TRUE ~ 'Others')
            
            enrichment_matrix[region,cell_property] <- table(regions_test,cell_properties_test) %>% 
                as.data.frame() %>% 
                arrange(desc(regions_test != "Others"),desc(cell_properties_test != "Others")) %>% 
                pivot_wider(names_from = cell_properties_test, values_from = Freq) %>% 
                select(-regions_test) %>% 
                fisher.test(alternative="greater") %>% 
                .[["p.value"]]
            
        }
    }

    return(enrichment_matrix)
}

plot_region_cell_property_enrichment <- function(enrichment_matrix){
    
    max_value <- -log(min(enrichment_matrix[enrichment_matrix != 0]))*2
    
    #Take inverse and log of enrichment matrix
    enrichment_matrix_inv <- -log(enrichment_matrix) %>% ifelse(is.infinite(.),max_value,.) %>% as.data.frame()

    #Pivot longe r to alow plotting
    enrichment_matrix_inv$regions <- rownames(enrichment_matrix_inv)
    enrichment_matrix_inv_long <- enrichment_matrix_inv %>% pivot_longer(-regions,names_to = 'projection_termination', values_to = 'inv_p_value')

    enrichment_matrix_inv_long$inv_p_value_plot <- ""
    to_add <- enrichment_matrix_inv_long$inv_p_value > -log(0.05/length(enrichment_matrix_inv_long$inv_p_value))
    enrichment_matrix_inv_long$inv_p_value_plot[to_add] <- round(enrichment_matrix_inv_long$inv_p_value[to_add])

    enrichment_matrix_inv_long$projection_termination <- factor(
        enrichment_matrix_inv_long$projection_termination, 
        levels = c("BST","PAG", "BST+PAG", "LHA", "PBN", "LHA+PBN", "PVT", "CEA", "PVT+CEA", "No-virus")
    )

    
    # TODO: Køre med faerre region labels naar jeg faar dem af bernd. Eller lave et plot af hver maaske. 
    #Kan bare køre functioner to gange. en for hver region opløsning. Klart det paeneste og bedste
    region_cell_property_enrichment_plot <- enrichment_matrix_inv_long %>%
        filter(!projection_termination == "No-virus") %>%
        ggplot() + 
        #geom_col(aes(x = projection_termination, y = inv_p_value, fill = regions), position = "dodge") +
        geom_col(aes(x = regions, y = inv_p_value)) + facet_grid(~projection_termination) + coord_flip() + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        geom_text(aes(x = regions, y = inv_p_value + 150, label = inv_p_value_plot))



    ggsave("./plots/enrichment_cell_by_region_corrected.png",region_cell_property_enrichment_plot)

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
    parser$add_argument("-tp","--transcript-path", type = "character", help="Relative path to transcripts file in Resolve format")
    parser$add_argument("-cp","--cell-roi-path", type = "character", help="Relative path to cell ROIs in numpy format")
    parser$add_argument("-ip","--image-roi-path", type = "character", help="Relative path to image ROIs in numpy format. If non given, then area correction is not performed")
    parser$add_argument("-rp","--region-path", type = "character", help="Relative path to folder with region ROIs in imageJ format. If non given, then region annotation is not performed")
    parser$add_argument("-si","--section-index", type="character", help="Section index of section")
    parser$add_argument("-ri","--run-index", type="character", help="Run index of section")

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
    transcript_path <- args$transcript_path
    cell_roi_path <- args$cell_roi_path
    image_roi_path <- args$image_roi_path
    region_path <- args$region_path
    run_index <- args$run_index
    section_index <- args$section_index
    
    # for add_projection_termination
    threshold <- 0.05
    threshold <- 0.05/ncol(spatial_object)

    # For sun1_virus_across_sections
    threshold <- 0.05
    threshold <- 0.05/ncol(spatial_object)

}

main()