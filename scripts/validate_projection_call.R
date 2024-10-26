library(tidyverse)
library(Seurat)
library(qs)
library(argparse)

utils <- new.env()
source("scripts/utils/utils.R", local = utils)

plot_utils <- new.env()
source("scripts/utils/plot_utils.R", local = plot_utils)



plot_sun1_vs_virus <- function(meta_data, plot_folder){
    # Input: data frame with meta data, path to plot folder
    #
    # Summarize meta data into counts of sun1 vs virus for each run and section that can be though of as contingency tables. Then calculates p-values for each contingency table.
    # Then plot the counts and p-values and save the plot in plot folder
    #
    # Output: NA

    meta_data <- meta_data %>% filter(!in_injection_site)

    # Calculate count for sun1 vs virus for each run and section
    sun1_vs_virus_counts <- meta_data %>% 
        group_by(sun1_vs_virus) %>% 
        summarize(count = n()) %>% mutate(total_count = sum(count),relative_count = count / total_count)
    
    # Calculate p-values for the sun1 vs virus counts for each run and section using a fisher test
    p_values <- sun1_vs_virus_counts %>% 
        select(-c(total_count,relative_count)) %>% 
        pivot_wider(names_from = sun1_vs_virus, values_from = count) %>%
            mutate(p_value = fisher.test(matrix(c(Other,Sun1,`Sun1_Virus`,Virus),nrow = 2))$p.value)

    # Plot the counts and p-values for sun1 vs virus across runs and sections
    sun1_vs_virus_bar_plot <- sun1_vs_virus_counts %>%
        ggplot() + geom_col(aes(x = sun1_vs_virus, y = relative_count)) + 
        geom_text(aes(x = sun1_vs_virus, y = relative_count + 0.01, label = count)) + 
        geom_text(data = p_values, aes(x = 2.5, y = Inf, label = format(p_value, scientific = TRUE, digits = 2)), hjust = 0.5, vjust = 1.2, inherit.aes = FALSE, size = 5)

    # Save the plot to plot folder
    ggsave(paste0(plot_folder,"/sun1_vs_virus_bar_plot.png"), sun1_vs_virus_bar_plot)
}

plot_sun1_vs_virus_across_sections <- function(meta_data, plot_folder){
    # Input: data frame with meta data, path to plot folder
    #
    # Summarize meta data into counts of sun1 vs virus for each run and section that can be though of as contingency tables. Then calculates p-values for each contingency table.
    # Then plot the counts and p-values and save the plot in plot folder
    #
    # Output: NA

    meta_data <- meta_data %>% filter(!in_injection_site)

    # Calculate count for sun1 vs virus for each run and section
    sun1_vs_virus_counts <- meta_data %>% 
        group_by(run_index, section_index, sun1_vs_virus) %>% 
        summarize(count = n()) %>% mutate(total_count = sum(count),relative_count = count / total_count)
    
    # Calculate p-values for the sun1 vs virus counts for each run and section using a fisher test
    p_values <- sun1_vs_virus_counts %>% 
        select(-c(total_count,relative_count)) %>% 
        pivot_wider(names_from = sun1_vs_virus, values_from = count) %>%
            mutate(p_value = fisher.test(matrix(c(Other,Sun1,`Sun1_Virus`,Virus),nrow = 2))$p.value)

    # Plot the counts and p-values for sun1 vs virus across runs and sections
    sun1_vs_virus_bar_plot <- sun1_vs_virus_counts %>%
        ggplot() + geom_col(aes(x = sun1_vs_virus, y = relative_count)) + facet_grid(section_index~run_index) + 
        geom_text(aes(x = sun1_vs_virus, y = relative_count + 0.1, label = count)) + 
        geom_text(data = p_values, aes(x = 2.5, y = Inf, label = format(p_value, scientific = TRUE, digits = 2)), hjust = 0.5, vjust = 1.2, inherit.aes = FALSE, size = 5)

    # Save the plot to plot folder
    ggsave(paste0(plot_folder,"/sun1_vs_virus_bar_plot_across_sections.png"), sun1_vs_virus_bar_plot)
}


calculate_sun1_virus_confusion_matrix <- function(sun1_vs_virus){

    sun1_virus_table <- sun1_vs_virus %>% 
        table()
    
    sun1_virus_table <- sun1_virus_table[c("Sun1_Virus", "Sun1", "Virus", "Other")]
    sun1_virus_matrix <- matrix(sun1_virus_table, nrow = 2, ncol = 2)

    fisher_test <- fisher.test(sun1_virus_matrix)
    cat(paste0(utils$current_time()," ",fisher_test,"\n"))

    expected_confusion_matrix <- chisq.test(sun1_virus_matrix)$expected
    cat(paste0(utils$current_time()," ",expected_confusion_matrix,"\n"))

}



check_background_intensity <- function(spatial_object, plot_folder){

    run_indices <- c("spatial1", "spatial2", "spatial3")
    section_indices <- c("A1", "B1", "C1", "D1", "A2", "B2", "C2", "D2")

    intensities <- spatial_object@meta.data %>% 
        group_by(run_index, section_index, virus1_count, virus2_count, section_area) %>% 
        summarize(total_cell_area = sum(area))
    
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
    

    intensities_long <- intensities %>% pivot_longer(
            cols = c(intensity_background_virus1,intensity_background_virus2,intensity_negative_cells_virus1,intensity_negative_cells_virus2),
            names_to = "type",
            values_to = "intensity"
        )

    
    check_background_intensity_plot <- intensities_long %>% 
        ggplot() + geom_col(aes(x = type, y = intensity, fill = type)) + facet_grid(section_index~run_index) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

    ggsave(paste0(plot_folder, "/check_background_intensity_plot.png"), check_background_intensity_plot)
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
    
    dir.create(paste0(path_to_output_folder),showWarnings = FALSE)
    dir.create(paste0(path_to_plot_folder),showWarnings = FALSE)
    
    cat(paste0(utils$current_time()," Saving all output files to ",path_to_output_folder,"\n"))
    cat(paste0(utils$current_time()," Saving all plots to ",path_to_plot_folder,"\n"))

    cat(paste0(utils$current_time()," Loading spatial object from ",spatial_object_path,"\n"))
    spatial_object <- qread(spatial_object_path)

    cat(paste0(utils$current_time()," Generating sun1-vs-virus plot\n"))
    plot_sun1_vs_virus(spatial_object@meta.data, path_to_plot_folder)
    
    cat(paste0(utils$current_time()," Generating sun1-vs-virus plot across sections\n"))
    plot_sun1_vs_virus_across_sections(spatial_object@meta.data, path_to_plot_folder)
    
    cat(paste0(utils$current_time()," Calculating confusion matrix\n"))
    calculate_sun1_virus_confusion_matrix(spatial_object$sun1_vs_virus[!spatial_object$in_injection_site])

    cat(paste0(utils$current_time()," Generating plot to check background intensity\n"))
    check_background_intensity(spatial_object, path_to_plot_folder)

}

main()

