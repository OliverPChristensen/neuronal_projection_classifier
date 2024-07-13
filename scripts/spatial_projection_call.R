library(reticulate)
library(argparse)
library(sf)
library(tidyverse)
library(Seurat)
library(qs)

utils <- new.env()
source("scripts/utils/utils.R", local = utils)

plot_utils <- new.env()
source("scripts/utils/plot_utils.R", local = plot_utils)

sf_utils <- new.env()
source("scripts/utils/sf_utils.R", local = sf_utils)


sun1_call <- function(sun1_rois, rois_cells){
    # Input: Sun1 ROIs as sfc object with POLYGON geometries, cell ROIs as sfc object with POLYGON geometries
    #
    # Takes sun1 ROIs and removes small ROIs that are likely not to be proper sun1 nuclei. Then Intersects cell ROIs with sun1 ROIs
    # and convert the output list into a boolean vector that indicates sun1 positive cells
    #
    # Output: Boolean vector with cell IDs as names. Indicate sun1 positive cells

    # Extract area from sun1 ROIs
    areas <- sun1_rois %>% st_area()

    # Exclude small sun1 ROIs that are likely not to stem from sun1 positive nuclei
    sun1_rois_size_excluded <- sun1_rois[areas > 500]

    # Extract centroids from size exluded sun1 ROIs
    sun1_centroids_size_excluded <- sun1_rois_size_excluded %>% st_centroid()

    # Generate list of cell ROIs that intersect centroids from sun1 ROIs. These will be called as sun1 positive cells
    sparse_intersects <- st_intersects(rois_cells,sun1_centroids_size_excluded)

    # Add cell IDs to intersects list
    names(sparse_intersects) <- names(rois_cells)
    
    # All index are above -1, so effectively turns the integer vector into a boolean that indicate sun1 positive cells
    sun1_bol <- unlist(sparse_intersects) > -1

    return(sun1_bol)

}

virus_call <- function(spatial_object, transcripts, virus){
    # Input: Seurat object with spatial data, image ROIs as sfc object with POLYGON geometries, transcripts as a sf object with POINT geometries, virus name
    #
    # Calculates total area of image, total virus count in image, combined area of cells and combined count of transcripts in cells.
    # Then calculates the null intensity and then null expectation of virus transcripts outside of cells and uses this to perform poissons test on each cell
    # to test whether a given cell has more transcripts than expected by chance, whcih is then used to make virus call
    #
    # Output: Vector with p-values from poissons test of virus counts and cell IDs as names

    # Combines ROIs of individual tiles in a staining image and extracts the combined area
    area_image <- spatial_object$section_area %>% unique()

    # Calculates the total number of virus transcripts within a section
    virus_image <- sum(transcripts$V4 == virus)

    # Extracts area meta data from spatial object
    area_cells <- spatial_object@meta.data %>% select(area) %>% as.matrix()
    
    # Extracts virus count within each cell in the spatial object. '_' is replaced with '-' in the virus name to comply with seurat format
    virus_cells <- spatial_object@assays$raw$counts[gsub("_", "-", virus),]

    # Calculates the null intensity of the null poissons distribution by calculating the total virus count outside of cells
    # and total area outside of cells and then find the intesity of virus trancripts outside of cells
    null_intensity <- (virus_image - sum(virus_cells))/(area_image - sum(area_cells))

    # Calculates the null expected number af trancripts within each cell in the spatial object using the null intensity of the section and the area of the given cell
    null_expectation <- null_intensity*area_cells
    
    # Makes zero vector to store p-values from poissons test
    virus_test <- rep(0,length(virus_cells))

    # Loops over each cell and calculates the p-value of the number of transcripts in the ith cell compared to the null expectation of the given cell
    for (i in 1:length(virus_cells)){
        virus_test[i] <- poisson.test(virus_cells[i], null_expectation[i], alternative = "greater")$p.value
    }

    # Gives corresponding cell IDs to vector of p-values
    names(virus_test) <- rownames(area_cells)

    return(virus_test)
}

create_parser <- function() {
    # Input: NA
    #
    # Creates and outputs parser with arguments that allows interaction via CLI
    #
    # Output: Parser

    # Define parser
    parser <- ArgumentParser(description = "Script to generate spatial object from ROIs and transcripts data file")

    # Add argument to parser
    parser$add_argument("-o","--output-folder-path", type = "character", help = "Relative path to folder to save output")
    parser$add_argument("-sop","--spatial-object-path", type = "character", help="Relative path to spatial object in seurat format")
    parser$add_argument("-tp","--transcript-path", type = "character", help="Relative path to transcripts file in Resolve format")
    parser$add_argument("-cp","--cell-roi-path", type = "character", help="Relative path to cell ROIs in numpy format")
    parser$add_argument("-sp","--sun1-roi-path", type = "character", help="Relative path to sun1 ROIs in numpy format")
    parser$add_argument("-si","--section-index", type="character", help="Section index of section")
    parser$add_argument("-ri","--run-index", type="character", help="Run index of section")

    return(parser)
}


main <- function(){
    # Input: NA
    #
    # Main function that defines parser argument, loads spatial object, transcripts, cell ROIs, image ROIs and sun1 ROIs for a given section.
    # Then sun1 call and virus call for each of the two viruses is performed and added as meta data in the spatial object. The spatial object is then saved
    #
    # Output: NA

    # Create parser and parse arguments
    parser <- create_parser()
    args <- parser$parse_args()

    # Retrieve arguments into corresponding variables
    path_to_output_folder <- args$output_folder_path
    spatial_object_path <- args$spatial_object_path
    transcript_path <- args$transcript_path
    cell_roi_path <- args$cell_roi_path
    sun1_roi_path <- args$sun1_roi_path
    run_index <- args$run_index
    section_index <- args$section_index

    cat(paste0(utils$current_time()," Loading spatial object from ",spatial_object_path,"\n"))
    # Load spatial seurat object to use for projection calls
    spatial_object <- qread(spatial_object_path)

    # Add empty meta data for projection calls
    spatial_object$sun1 <- FALSE
    spatial_object$virus1 <- NA
    spatial_object$virus2 <- NA

    cat(paste0(utils$current_time()," Generating projection calls for section ",section_index,"\n"))

    # Load transcripts
    cat(paste0(utils$current_time()," Loading transcripts from ",transcript_path,"\n"))
    transcripts <- sf_utils$load_transcripts(transcript_path)

    # Load cell ROIs
    cat(paste0(utils$current_time()," Loading cell ROIs from ",cell_roi_path,"\n"))
    rois_cells <- sf_utils$array_to_sf_rois(cell_roi_path)

    # Load sun1 ROIs
    cat(paste0(utils$current_time()," Loading sun1 ROIs from ", sun1_roi_path,"\n"))
    rois_sun1 <- sf_utils$array_to_sf_rois(sun1_roi_path)

    # Generate sun1 calls
    cat(paste0(utils$current_time()," Generating sun1 call for section ", section_index,"\n"))
    sun1_bol <- sun1_call(rois_sun1, rois_cells)

    # Add sun1 call to meta data
    spatial_object <- AddMetaData(spatial_object,sun1_bol,"sun1")

    # Generate virus calls for virus 1
    cat(paste0(utils$current_time()," Generating virus calls for section ", section_index,"\n"))
    virus1_call <- virus_call(spatial_object, transcripts, "p147_mCherry_2A_iCre")

    # Add virus calls for virus 1 to meta data
    spatial_object <- AddMetaData(spatial_object,virus1_call,"virus1")
    
    # Add total virus count of virus 1
    spatial_object$virus1_count <- sum(transcripts$V4 == "p147_mCherry_2A_iCre")

    # Generate virus calls for virus 2
    virus2_call <- virus_call(spatial_object, transcripts, "p36_NLS_Cre")

    # Add virus calls for virus 2 to meta data
    spatial_object <- AddMetaData(spatial_object,virus2_call,"virus2")

    # Add total virus count of virus 2
    spatial_object$virus2_count <- sum(transcripts$V4 == "p36_NLS_Cre")

    cat(paste0(utils$current_time()," Saving spatial seurat object with projection calls to ", path_to_output_folder,"/",section_index,"_spatial_object_projection_call.qs","\n"))
    dir.create(paste0(path_to_output_folder),showWarnings = FALSE)
    qsave(spatial_object,paste0(path_to_output_folder,"/",section_index,"_spatial_object_projection_call.qs"))
}

main()

