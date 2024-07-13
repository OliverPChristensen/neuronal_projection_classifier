library(argparse)
library(tidyverse)
library(Seurat)
library(Matrix)
library(qs)
library(cowplot)

utils <- new.env()
source("scripts/utils/utils.R", local = utils)

plot_utils <- new.env()
source("scripts/utils/plot_utils.R", local = plot_utils)

sf_utils <- new.env()
source("scripts/utils/sf_utils.R", local = sf_utils)


remove_virus_transcripts <- function(spatial_object){
    # Input: Spatial Seurat object
    #
    # Adds layer in spatial Seurat object without virus transcripts and adds non virus count layer 
    # and non virus total count and feature meta data
    #
    # Output: Spatial Seurat object with non virus count layer and non virus total count and feature meta data

    # Extract count matrix from spatial object
    counts <- spatial_object[['raw']]$counts

    # Remove virus counts from count matrix
    counts_non_virus <- counts[!rownames(counts) %in% c("p147-mCherry-2A-iCre","p36-NLS-Cre"),]

    # Add non virus count matrix to the spatial object
    spatial_object[['raw']]$counts_non_virus <- counts_non_virus

    # Add total counts from the non virus count matrix to the meta data
    spatial_object$nCount_non_virus <- colSums(counts_non_virus)

    # Add total non zero features from the non virus count matrix to the meta data
    spatial_object$nFeature_non_virus <- colSums(counts_non_virus > 0)

    return(spatial_object)
}

quality_control <- function(spatial_object){
    # Input: Spatial Seurat object with nCount_non_virus and area meta data
    #
    # Subset spatial Seurat object based on max density of trancript count and cell area
    #
    # Output: Quality controlled spatial Seurat object

    # Create density plot of counts from spatial meta data
    p <- spatial_object@meta.data %>% ggplot() + geom_density(aes(x = nCount_non_virus))

    # Extract the data from the density plot
    plot_data <- ggplot_build(p)$data

    # Extract y coordinate from density plot corresponding to the density
    density <- plot_data[[1]]$y

    # Extract x coordinate from density plot corresponding to the count
    count <- plot_data[[1]]$x

    # Find the count corresponding to the max density
    max_count <- count[which.max(density)]

    # Find 1% and 99% countiles of area meta data from spatial object
    area_quantiles <- quantile(spatial_object$area, c(0.01,0.99))
    
    # Subset spatial data according to count with max density and 1% and 99% area quantiles
    spatial_object <- subset(spatial_object, subset = nCount_non_virus > max_count & area > area_quantiles[1] & area < area_quantiles[2])

    return(spatial_object)
}

area_normalize_data <- function(spatial_object,scale){
    # Input: Spatial Seurat object with 'raw' Assay containing 'counts_non_virus' layer, and area meta data, scale (int/float)
    #
    # Extract non virus counts matrix from raw Assay from spatial Seurat object, divide by the area of the corresponding cell, multiply by scale argument and log transformize
    #
    # Output: Area normalized spatial Seurat object

    # Area normalize the spatial object by taking the raw counts and dividing by the area of the corresponding cell
    # Then the counts are multiplied by the scale input and log transformed
    spatial_object[['raw']]$norm_data <- log1p(t(t(spatial_object[['raw']]$counts_non_virus) / spatial_object$area)*scale)

    return(spatial_object)
}

load_qc_norm_merge_pr_section <- function(spatial_paths_list){
    # Input: List of paths to spatial objects
    #
    # Loads spatial objects iteratively from list of paths to spatial objects and performs quality control and area normalization on each loaded spatial object separately
    # Then merges all spatial objects within a run, but keeps spatial objects between runs separate
    #
    # Output: List of QC'ed, area normalized spatial objects that are merged across sections within a run

    # Define run indices and section indices
    run_indices <- c("spatial1", "spatial2", "spatial3")
    section_indices <- c("A1", "B1", "C1", "D1", "A2", "B2", "C2", "D2")

    # Define list to put merged spatial objects into
    merged_spatial_objects <- list()

    for (run_index in run_indices){
        
        cat(paste0(utils$current_time()," Loading spatial objects from ",run_index,":\n"))

        # Define list to put spatial objects into
        spatial_objects <- list()

        for (section_index in section_indices){

            # Extract path to spatial object according to current spatial index and run index
            spatial_object_path <- grep(run_index, spatial_paths_list, value = TRUE) %>% grep(section_index, ., value = TRUE)

            # Load spatial object
            cat(paste0(utils$current_time(),"  - Loading spatial object ",section_index," from ",spatial_object_path,"\n"))
            spatial_object <- qread(spatial_object_path)

            # Create Seurat Layer without virus transcripts
            cat(paste0(utils$current_time(),"  - Generating non-virus layer for spatial object ",section_index,"\n"))
            spatial_object <- remove_virus_transcripts(spatial_object)

            # Quality control spatial object
            cat(paste0(utils$current_time(),"  - Quality controlling spatial object ",section_index,"\n"))
            spatial_object <- quality_control(spatial_object)

            # Area normalize spatial object
            cat(paste0(utils$current_time(),"  - Area normalising spatial object ",section_index,"\n"))
            spatial_objects[[section_index]] <- area_normalize_data(spatial_object, scale = 1e4)
        }

        # Merge spatial objects within the current run
        cat(paste0(utils$current_time()," Merging spatial seurat objects from ",run_index,"\n"))
        merged_spatial_objects[[run_index]] <- merge(x = spatial_objects[[1]], y = spatial_objects[-1])
    }

    return(merged_spatial_objects)
}


check_batch_effect_summary <- function(spatial_object_list, plot_folder){
    # Input: List of spatial Seurat objects
    #
    # Extracts meta data from list of spatial Seurat objects. Then generates and saves various summary plots to assess batch effects between sections and runs
    #
    # Output: NA

    # Define emtpy data frame to place spatial object meta data across different spatial objects
    spatial_objects_meta_data <- data.frame()

    # Loop over spatial objects in spatial object list and extract meta data
    for (spatial_object in spatial_object_list){
        spatial_objects_meta_data <- rbind(spatial_objects_meta_data,spatial_object@meta.data)
    }

    # Generate and save density plot of cell areas across sections and runs
    area_x_sections <- spatial_objects_meta_data %>% ggplot() + geom_density(aes(x = area, col = section_index)) + facet_grid(~ run_index)
    ggsave(paste0(plot_folder,"/","spatial_object_check_batch_effect_area_x_sections.png"),area_x_sections)
    
    # Generate and save density plot of cell areas across runs
    area_x_runs <- spatial_objects_meta_data %>% ggplot() + geom_density(aes(x = area, col = run_index))
    ggsave(paste0(plot_folder,"/","spatial_object_check_batch_effect_area_x_runs.png"),area_x_runs)

    # Generate and save density plot of transcript count across sections and runs
    counts_x_sections <- spatial_objects_meta_data %>% ggplot() + geom_density(aes(x = nCount_raw, col = section_index)) + facet_grid(~ run_index)
    ggsave(paste0(plot_folder,"/","spatial_object_check_batch_effect_counts_x_sections.png"),counts_x_sections)

    # Generate and save density plot of transcript count across runs
    counts_x_runs <- spatial_objects_meta_data %>% ggplot() + geom_density(aes(x = nCount_raw, col = run_index))
    ggsave(paste0(plot_folder,"/","spatial_object_check_batch_effect_counts_x_runs.png"),counts_x_runs)

    # Extract and calculate the transcript intensity across sections and runs
    intensity_sections <- spatial_objects_meta_data %>% group_by(section_index, run_index) %>% summarise(section_area = unique(section_area), section_count = unique(section_count)) %>% mutate(section_intensity = (section_count/section_area))
    
    # Generate and save plot of the transcript intensity across sections and runs
    intensity_x_sections <- intensity_sections %>% ggplot() + geom_col(aes(x = section_index,y = section_intensity)) + facet_grid(~ run_index)
    ggsave(paste0(plot_folder,"/","spatial_object_check_batch_effect_intensity_x_sections.png"),intensity_x_sections)

    # Generate and save plot of the transcript intensity across runs
    intensity_x_runs <- intensity_sections %>% group_by(run_index) %>% summarise(run_intensity = mean(section_intensity)) %>% ggplot() + geom_col(aes(x = run_index,y = run_intensity))
    ggsave(paste0(plot_folder,"/","spatial_object_check_batch_effect_intensity_x_runs.png"),intensity_x_runs)
}

pca_plot <- function(spatial_object,index_name,pca,variance_df,first_pc, second_pc, plot_folder){
    # Input: Spatial Seurat object, index name to indicate meta data to use for plotting, PCA object from prcomp,
    # data frame with variance proportions of each principle component, index of first PC to use for plotting, index of second PC to use for plotting
    #
    # Creates data frame with meta data according to index name. Then generates and saves PCA plots with PCs according to index og first PC and second PC
    # with x axis and y axis densities
    #
    # Output: NA

    # Create a data frame with the PCA results and index
    pca_df <- data.frame(
        PC1 = pca$x[, first_pc],
        PC2 = pca$x[, second_pc],
        index = spatial_object@meta.data[,index_name]
    )

    pca_df_shuffled <- pca_df[sample(nrow(pca_df)), ]

    # Plot the specified principle components
    plot <- ggplot(pca_df_shuffled, aes(x = PC1, y = PC2, color = index)) +
        geom_point(size = 1) +
        xlab(paste("PC",first_pc," (", round(variance_df$prop_variance[1] * 100, 2), "% Variance)", sep = "")) +
        ylab(paste("PC",second_pc," (", round(variance_df$prop_variance[2] * 100, 2), "% Variance)", sep = "")) + 
        labs(color = index_name) +
        theme_classic()

    # Create the density plots for x axis of PCA plot
    x_density <- ggplot(pca_df_shuffled, aes(x = PC1, color = index)) +
        geom_density() +
        theme_classic() +
        theme(legend.position = "none") + 
        theme(axis.title.x = element_blank())

    # Create the density plots for y axis of PCA plot
    y_density <- ggplot(pca_df_shuffled, aes(x = PC2, color = index)) +
        geom_density() +
        coord_flip() +
        theme_classic() +
        theme(legend.position = "none") + 
        theme(axis.title.y = element_blank())


    # Create an empty plot with a white background for the top right corner
    empty_plot <- ggplot() + 
        theme_classic()

    # Align the plots
    aligned_plots <- align_plots(x_density, plot, y_density, align = "hv", axis = "tblr")

    # Arrange the plots into a single plot
    final_plot <- plot_grid(
        aligned_plots[[1]], empty_plot,
        aligned_plots[[2]], aligned_plots[[3]],
        ncol = 2, nrow = 2,
        rel_widths = c(3, 1),
        rel_heights = c(1, 3)
        )
    
    # Save the combined PCA plot
    ggsave(paste0(plot_folder,"/","spatial_object_check_batch_effect_pca_plot_",index_name,"_PC",first_pc,"_PC",second_pc,".png"),final_plot)

}

check_batch_effect_pca <- function(spatial_object, plot_folder){
    # Input: Spatial Seurat object
    #
    # Extracts scaled count matrix from spatial Seurat object uses it for PCA. Then calculates forportion of variance for each PC
    # and generates and saves Elbow plot and PCA plots
    #
    # Output: NA

    # Extract the scaled count matrix from the spatial object
    count_matrix <- spatial_object[['raw']]['scale.data']

    # Perform PCA on the scaled count matrix
    pca <- prcomp(t(count_matrix))

    # Calculate proportion of variance in each PC and create a data frame
    variance_df <- data.frame(
        pc = 1:nrow(spatial_object[['raw']]['scale.data']),
        prop_variance = summary(pca)$importance["Proportion of Variance",]
    )

    # Genrate and save the elbow plot
    cat(paste0(utils$current_time(),"  - Generating Elbow plot from PCA\n"))
    elbow_plot <- ggplot(variance_df, aes(x = pc, y = prop_variance)) +
        geom_point() +
        geom_line() +
        xlab("Principal Component") +
        ylab("Variance Explained") +
        ggtitle("Elbow Plot")
    
    ggsave(paste0(plot_folder,"/","spatial_object_check_batch_effect_pca_elbow_plot.png"),elbow_plot)

    # Loop over first 10 PCs
    for (pcs in 1:5){
        
        # Generate PCA plot for the current PCs across sections (Useless?)
        cat(paste0(utils$current_time(),"  - Generating PCA plot of PC ",2*pcs-1," and PC ",2*pcs," across sections\n"))
        pca_plot(spatial_object, 'section_index', pca, variance_df, 2*pcs-1, 2*pcs,plot_folder)
        
        # Generate PCA plot for the current PCs across runs
        cat(paste0(utils$current_time(),"  - Generating PCA plot of PC ",2*pcs-1," and PC ",2*pcs," across runs\n"))
        pca_plot(spatial_object, 'run_index', pca, variance_df, 2*pcs-1, 2*pcs,plot_folder)
    }

}

check_batch_effect_umap <- function(spatial_object,reduction_index, plot_folder){
    # Input: Spatial Seurat object
    #
    # Generates and saves UMAP plots of spatial Seurat object of sections within each run and a UMAP plot across runs
    #
    # Output: NA

    # Generate and save UMAP plot for spatial 1 across sections
    cat(paste0(utils$current_time(),"  - Generating UMAP plot of Spatial 1 data across sections\n"))
    umap_plot_runs <- DimPlot(spatial_object, reduction = reduction_index, group.by = "section_index", raster = FALSE, cells = WhichCells(spatial_object, expression = run_index == "spatial1"), shuffle = TRUE, pt.size = 0.1)
    ggsave(paste0(plot_folder,"/","spatial_object_check_batch_effect_spatial1_",reduction_index,".png"),umap_plot_runs)

    # Generate and save UMAP plot for spatial 2 across sections
    cat(paste0(utils$current_time(),"  - Generating UMAP plot of Spatial 2 data across sections\n"))
    umap_plot_runs <- DimPlot(spatial_object, reduction = reduction_index, group.by = "section_index", raster = FALSE,cells = WhichCells(spatial_object, expression = run_index == "spatial2"), shuffle = TRUE, pt.size = 0.1)
    ggsave(paste0(plot_folder,"/","spatial_object_check_batch_effect_spatial2_",reduction_index,".png"),umap_plot_runs)

    # Generate and save UMAP plot for spatial 3 across sections
    cat(paste0(utils$current_time(),"  - Generating UMAP plot of Spatial 3 data across sections\n"))
    umap_plot_runs <- DimPlot(spatial_object, reduction = reduction_index, group.by = "section_index", raster = FALSE,cells = WhichCells(spatial_object, expression = run_index == "spatial3"), shuffle = TRUE, pt.size = 0.1)
    ggsave(paste0(plot_folder,"/","spatial_object_check_batch_effect_spatial3_",reduction_index,".png"),umap_plot_runs)

    # Generate and save UMAP plot across runs
    cat(paste0(utils$current_time(),"  - Generating UMAP plot of data across runs\n"))
    umap_plot_sections <- DimPlot(spatial_object, reduction = reduction_index, group.by = "run_index", raster = FALSE, shuffle = TRUE, pt.size = 0.1)
    ggsave(paste0(plot_folder,"/","spatial_object_check_batch_effect_runs_",reduction_index,".png"),umap_plot_sections)
}


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
        region_names == 'dorsal_hyothalamus' ~ 'dorsal_hypothalamus',
        region_names == 'anterior_hyopthalamus_left' ~ 'anterior_hypothalamus_left',
        region_names == 'interpeduncular' ~ 'interpreduncular',
        region_names == 'lateral_hypothalamus_right-1' ~ 'lateral_hypothalamus_right',
        region_names == 'pontine_grey' ~ 'pontine_gray',
        region_names == 'retrochiasmaric_left' ~ 'retrochiasmatic_left',
        region_names == 'suprachiamsatic' ~ 'suprachiasmatic',
        region_names == 'dorsomedial_hypothalamus' ~ 'dorsomedial',
        TRUE ~ region_names
    )

    return(region_names)
}

reduce_regions <- function(meta_data){
    # Input: Character vector of region names
    #
    # Reduces region names by extracting arcuate|ventromedial|dorsomedial|lateral_hypothalamus and removing left/right discrepancy
    # All other regions are set to "Other". Also sets lateral_hypothalamus_right from spatial2 to "Other", because this region includes an injection site
    #
    # Output: Character vector of reduced region names

    # Set lateral_hypothalamus_right from spatial2 to "Other" and extracts arcuate|ventromedial|dorsomedial|lateral_hypothalamus with removed left/right suffix
    regions_reduced <- case_when(
        meta_data$region_corrected == "lateral_hypothalamus_right" & meta_data$run_index == "spatial2" ~ "Other",
        str_detect(meta_data$region_corrected, "arcuate|ventromedial|dorsomedial|lateral_hypothalamus") ~ gsub("_right|_left", "", meta_data$region_corrected),
        TRUE ~ "Other"
    )
    
    return(regions_reduced)
}

add_projection <- function(meta_data, virus1, virus2, threshold){
    # Input: meta_data from cell object, numeric vector with corrected p-values for virus1 calls, numeric vector with corrected p-values for virus2 calls
    # Threshold alpha value for statistical test
    #
    # Assigns projection labels to cells in meta_data according to the virus1 and virus2 p-values, run_index (indicates injection site) and threshold.
    # Also excludes lateral_hypothalamus_right from spatial2, because this includes the injection site.
    #
    # Output: Character vector with projection labels

    
    projection <- case_when(
        meta_data$region_corrected == "lateral_hypothalamus_right" & meta_data$run_index == "spatial2" ~ "Other",
        virus1 < threshold & virus2 > threshold & meta_data$run_index == "spatial1" ~ "BST",
        virus2 < threshold & virus1 > threshold & meta_data$run_index == "spatial1" ~ "PAG",
        virus1 < threshold & virus2 < threshold & meta_data$run_index == "spatial1" ~ "BST+PAG",
        virus1 < threshold & virus2 > threshold & meta_data$run_index == "spatial2" ~ "LHA",
        virus2 < threshold & virus1 > threshold & meta_data$run_index == "spatial2" ~ "PBN",
        virus1 < threshold & virus2 < threshold & meta_data$run_index == "spatial2" ~ "LHA+PBN",
        virus1 < threshold & virus2 > threshold & meta_data$run_index == "spatial3" ~ "PVT",
        virus2 < threshold & virus1 > threshold & meta_data$run_index == "spatial3" ~ "CEA",
        virus1 < threshold & virus2 < threshold & meta_data$run_index == "spatial3" ~ "PVT+CEA",
        TRUE ~ "Other"
    )

    projection <- factor(projection, levels = c("BST", "PAG", "BST+PAG","LHA","PBN", "LHA+PBN","PVT","CEA","PVT+CEA", "Other"))

    return(projection)
}

extract_single_projections <- function(meta_data){
    # Input: meta_data from cell object
    #
    # Extracts the single projection labels from the projection meta data and converts double projections to "Other"
    #
    # Output: Character vector with single projections

    single_projection <- case_when(
        meta_data$projection %in% c("BST", "PAG", "LHA", "PBN", "PVT", "CEA") ~ meta_data$projection,
        TRUE ~ "Other"
    )

    return(single_projection)
}

add_sun1_vs_virus <- function(meta_data, virus1, virus2){
    # Input: meta_data from cell object, numeric vector with corrected p-values for virus1 calls, numeric vector with corrected p-values for virus2 calls
    #
    # Combines the sun1 label and virus1 and virus2 p-values to give a label that indicates both sun1 label and virus label
    #
    # Output: Character vector with sun1 and virus labels

    sun1_vs_virus <- case_when(
        meta_data$region_corrected == "lateral_hypothalamus_right" & meta_data$run_index == "spatial2" ~ "Other",
        (meta_data$sun1 == TRUE) & !(virus1 < 0.05 | virus2 < 0.05)~ "Sun1",
        !(meta_data$sun1 == TRUE) & (virus1 < 0.05 | virus2 < 0.05) ~ "Virus",
        (meta_data$sun1 == TRUE) & (virus1 < 0.05 | virus2 < 0.05) ~ "Sun1_Virus",
        TRUE ~ "Other"
    )

    return(sun1_vs_virus)
}



plot_and_save_umap <- function(spatial_object, column, plot_folder, cells = NULL, grey_label = NULL){
    # Input: Spatial Seurat object, character column name to plot in UMAP, plot folder path, list of cell names to include (NULL includes all cells), character name of label from column to plot in grey
    #
    # Generates colors to use for UMAP and then plots cell according to column and saves UMAP
    #
    # Output: NA

    # Generate colors to use for UMAP
    plotting_colors <- plot_utils$generate_label_colors(spatial_object@meta.data[[column]],grey_label)

    # Plot UMAP using Seurat DimPlot function
    umap <- DimPlot(spatial_object, group.by = column, cols = plotting_colors, reduction = 'umap_features', cells = cells, raster = FALSE, shuffle = TRUE)

    # Save UMAP to plot folder
    ggsave(paste0(plot_folder,"/umap_",column,".png"), umap, width = 5000, height = 3000, units = "px")

}


plot_and_save_feature_map <- function(spatial_object, column, plot_folder){
    # Input: Spatial Seurat object, character column name to plot in feature MAP, plot folder path
    #
    # Plots feature MAP of cells with meta data from column. Then saves to feature MAP to plot folder
    #
    # Output: NA

    # Plot feature MAP
    feature_map <- FeaturePlot(spatial_object, feature = column, reduction = "umap_features", raster = FALSE)

    # Save feature MAP in plot folder
    ggsave(paste0(plot_folder,"/umap_",column,".png"), feature_map, width = 4000, height = 3000, units = "px")
}

create_parser <- function() {
    # Input: NA
    #
    # Creates and outputs parser with arguments that allows interaction via CLI
    #
    # Output: Parser

    #Define parser
    parser <- ArgumentParser(description = "Script to generate spatial object from ROIs and transcripts data file")

    #Add argumnent to parser
    parser$add_argument("-o","--output-folder-path", type = "character", default = './processed_data', help = "Relative path to folder to save output")
    parser$add_argument("-p","--plot-folder-path", type = "character", help = "Relative path to folder to save plots")
    parser$add_argument("-sops","--spatial-object-paths", nargs="+", type = "character", help="Characters to be used to specify transcript files in raw data folder")

    return(parser)
}


main <- function(){
    # Input: NA
    #
    # Main function that defines parser argument, loads spatial object, performs QC, area normalization and merges first sections and then runs
    # Also generate plots to check for batch effects between sections and between runs. This is done from both summary information, PCA and UMAP
    # Lastly, the preprocessed spatial object is saved
    #
    # Output: NA

    # Create parser and parse arguments
    parser <- create_parser()
    args <- parser$parse_args()

    # Retrieve arguments into corresponding variables
    path_to_output_folder <- args$output_folder_path
    path_to_plot_folder <- args$plot_folder_path
    spatial_object_paths <- args$spatial_object_paths

    dir.create(paste0(path_to_output_folder),showWarnings = FALSE)
    dir.create(paste0(path_to_plot_folder),showWarnings = FALSE)

    cat(paste0(utils$current_time()," Saving all output files to ",path_to_output_folder,"\n"))
    cat(paste0(utils$current_time()," Saving all plots to ",path_to_plot_folder,"\n"))
    
    # Load, quality control, normalise pr section and merge sections within runs
    spatial_object_list <- load_qc_norm_merge_pr_section(spatial_object_paths)

    # Generate summary plots to check for batch effects between sections and runs
    cat(paste0(utils$current_time()," Generating summary plots to check for batch effects between sections and runs\n"))
    check_batch_effect_summary(spatial_object_list, path_to_plot_folder)

    # Merge spatial objects across runs, because summary plots did not show any batch effect between runs
    cat(paste0(utils$current_time()," Merging spatial seurat objects across runs\n"))
    spatial_object <- merge(x = spatial_object_list[[1]], y = spatial_object_list[-1])

    # Scale data by z-normalising pr gene
    cat(paste0(utils$current_time()," Scaling normalized data from merged spatial seurat object\n"))
    spatial_object <- ScaleData(spatial_object)

    # Run PCA on scaled data
    cat(paste0(utils$current_time()," Running PCA on scaled data from spatial seurat object\n"))
    spatial_object <- RunPCA(spatial_object, features = rownames(spatial_object), approx = FALSE)
    
    # Run UMAP on PCA using first 40 principle components
    cat(paste0(utils$current_time()," Running UMAP on PCA from spatial seurat object\n"))
    spatial_object <- RunUMAP(spatial_object, dims = 1:40, reduction.name = "umap_pca")

    # Run UMAP on features
    cat(paste0(utils$current_time()," Running UMAP on features from spatial seurat object\n"))
    spatial_object <- RunUMAP(spatial_object, features = rownames(spatial_object[['raw']]['scale.data']), reduction.name = "umap_features")

    # Generate PCA plots to check for batch effects between sections and runs
    cat(paste0(utils$current_time()," Generating PCA plots to check for batch effects between sections and runs:\n"))
    check_batch_effect_pca(spatial_object, path_to_plot_folder)
    
    # Generate UMAP plots from PCA to check for batch effects between sections and runs
    cat(paste0(utils$current_time()," Generating UMAP plots from PCA to check for batch effects between sections and runs:\n"))
    check_batch_effect_umap(spatial_object, "umap_pca", path_to_plot_folder)

    # Generate UMAP plots from features to check for batch effects between sections and runs
    cat(paste0(utils$current_time()," Generating UMAP plots from features to check for batch effects between sections and runs:\n"))
    check_batch_effect_umap(spatial_object, "umap_features", path_to_plot_folder)

    # Genrate corrected and reduced region names for spatial object
    cat(paste0(utils$current_time()," Correcting and reducing region names\n"))
    spatial_object$region_corrected <- correct_region_names(spatial_object$region)
    spatial_object$region_reduced <- reduce_regions(spatial_object@meta.data)

    # Generate and save UMAPs with corrected and reduced region names
    cat(paste0(utils$current_time()," Generating UMAP plots for corrected region names and reduced region names\n"))
    plot_and_save_umap(spatial_object, 'region_corrected', path_to_plot_folder)
    plot_and_save_umap(spatial_object, 'region_reduced', path_to_plot_folder)

    #TODO: Remove p-values from LHA spatial2 before running p.adjust? so it does not interfere with fdr correction?

    # Add projection to meta data with no multiple testing correction
    cat(paste0(utils$current_time()," Adding projection to meta data with no multiple testing correction\n"))
    spatial_object$projection <- add_projection(
        spatial_object@meta.data,
        spatial_object$virus1,
        spatial_object$virus2,
        threshold = 0.05
    )

    # Add projection to meta data with FDR correction
    cat(paste0(utils$current_time()," Adding projection to meta data with FDR correction\n"))
    spatial_object$projection_fdr <- add_projection(
        spatial_object@meta.data,
        p.adjust(spatial_object$virus1, method = "fdr"),
        p.adjust(spatial_object$virus2, method = "fdr"),
        threshold = 0.05
    )

    # Add projection to meta data with bonferroni correction
    cat(paste0(utils$current_time()," Adding projection to meta data with Bonferroni correction\n"))
    spatial_object$projection_bonferroni <- add_projection(
        spatial_object@meta.data,
        p.adjust(spatial_object$virus1, method = "bonferroni"),
        p.adjust(spatial_object$virus2, method = "bonferroni"),
        threshold = 0.05
    )

    # Add single projection to meta data where double projections are excluded
    cat(paste0(utils$current_time()," Adding single projections to meta data\n"))
    spatial_object$single_projection <- extract_single_projections(spatial_object@meta.data)

    # Add sun1-vs-virus label to meta data
    cat(paste0(utils$current_time()," Adding sun1-vs-virus to meta data\n"))
    spatial_object$sun1_vs_virus <- add_sun1_vs_virus(spatial_object@meta.data, spatial_object$virus1, spatial_object$virus2)

    # Extract cell names of cells with single projections to be used for UMAP
    single_projection_cells <- spatial_object@meta.data %>% 
        filter(single_projection != "Other") %>% 
        rownames()

    # Extract cell names of cells with sun1 and/or virus label to be used for UMAP
    sun1_vs_virus_cells <- spatial_object@meta.data %>% 
        filter(sun1_vs_virus != "Other") %>% 
        rownames()

    # Generate and UMAP plot with cells that have single projections
    cat(paste0(utils$current_time()," Generating UMAP plots for single projections\n"))
    plot_and_save_umap(spatial_object, 'single_projection', path_to_plot_folder, cells = single_projection_cells)

    # Generate and UMAP plot with cells that have sun1 and/or virus label
    cat(paste0(utils$current_time()," Generating UMAP plots for sun1-vs-virus\n"))
    plot_and_save_umap(spatial_object, 'sun1_vs_virus', path_to_plot_folder, cells = sun1_vs_virus_cells)
    
    # Generate and feature MAP plot of cell areas
    cat(paste0(utils$current_time()," Generating feature map plots for area data\n"))
    plot_and_save_feature_map(spatial_object,'area', path_to_plot_folder)

    # Save preprocessed spatial seurat object
    cat(paste0(utils$current_time()," Saving preprocessed spatial seurat object to ", path_to_output_folder,"/","spatial_object_preprocessed.qs","\n"))
    qsave(spatial_object, paste0(path_to_output_folder,"/","spatial_object_preprocessed.qs"))
}

main()