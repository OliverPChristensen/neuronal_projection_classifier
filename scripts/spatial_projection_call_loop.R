library(reticulate)
library(argparse)
library(sf)
library(tidyverse)
library(Seurat)
library(qs)


load_transcripts <- function(transcripts_path){
    # Input: Path to transcript file
    #
    # Loades the transcript file, removes false positive (FP) transripts and creates a sf object from the trasncript dataframe
    #
    # Output: transcripts as a sf object with POINT geometries

    # Read transcripts file and remove false positive transcripts
    transcripts <- read.table(transcripts_path, header = F, sep = "\t") %>% filter(!str_detect(V4, "FP "))

    # Creates sf object from transcript data frame with POINT geometries
    transcript_points <- transcripts %>% st_as_sf(.,coords= c("V1","V2"))
    
    return(transcript_points)
}

array_to_sf_rois <- function(rois_path){
    # Input: Path to rois file saved in numpys npz format
    #
    # Loads rois as numpy arrays and converts them to sfc objects containing rois as POLYGON geometries
    #
    # Output: sfc object containing rois as POLYGON geometries

    # Import numpy module
    np <- import("numpy")

    # Load rois as arrays
    arrays <- np$load(rois_path)

    # Convert rois to a list of POLYGON objects
    rois <- list()
    for (roi in arrays$files){
        roi_coords <- arrays$get(roi)

        # Closing the roi coords
        roi_coords_closed <- rbind(roi_coords,roi_coords[1,])
        rois[[roi]] <- st_polygon(list(roi_coords_closed))
    }

    # Converting list of POLYGONs to sfc object
    rois_sfc <- rois %>% st_sfc()

    # Removing non-valid ROIs
    rois_valid <- rois_sfc[st_is_valid(rois_sfc)]

    return(rois_valid)
}

sun1_call <- function(sun1_rois, rois_cells){
    #
    #
    #
    #
    #


    areas <- sun1_rois %>% st_area()
    sun1_rois_size_excluded <- sun1_rois[areas > 500]
    sun1_centroids_size_excluded <- sun1_rois_size_excluded %>% st_centroid()
    sparse_intersects <- st_intersects(rois_cells,sun1_centroids_size_excluded)

    names(sparse_intersects) <- names(rois_cells) #
    
    #All index are above -1, so effectively turns the integer vector into a boolean
    sun1_bol <- unlist(sparse_intersects) > -1
    return(sun1_bol)

}

section_index_to_seurat_integer <- function(seurat_object, section){
    #
    #
    #
    #
    #

    seurat_object@meta.data %>% filter(section_index == section) %>% rownames() -> seurat_cell_ids
    seurat_section_integer_index <- substr(seurat_cell_ids, nchar(seurat_cell_ids), nchar(seurat_cell_ids)) %>% unique()

    return(seurat_section_integer_index)
}

virus_call <- function(spatial_object, rois_image, transcripts, section, virus){
    #
    #
    #
    #
    #

    area_image <- rois_image %>% st_combine() %>% st_area
    virus_image <- sum(transcripts$V4 == virus)

    area_cells <- spatial_object@meta.data %>% filter(section_index == section) %>% select(area) %>% as.matrix()

    seurat_section_integer_index <- section_index_to_seurat_integer(spatial_object, section)
    
    virus_cells <- spatial_object@assays$raw[paste("counts",seurat_section_integer_index,sep = ".")][gsub("_", "-", virus),]

    null_intensity <- (virus_image - sum(virus_cells))/(area_image - sum(area_cells))

    null_expectation <- null_intensity*area_cells
    
    virus_test <- rep(0,length(virus_cells))

    for (i in 1:length(virus_cells)){
        virus_test[i] <- poisson.test(virus_cells[i], null_expectation[i], alternative = "greater")$p.value
    }

    names(virus_test) <- rownames(area_cells)

    return(virus_test)
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
    parser$add_argument("-fpt","--raw-folder-path", type = "character", default = 'default', help="Relative path to raw data folder")
    parser$add_argument("-fpr","--processed-folder-path", type = "character", default = 'default', help="Relative path to processed data folder")
    parser$add_argument("-o","--output-folder-path", type = "character", default = './processed_data', help = "Relative path to folder to save output")
    parser$add_argument("-fmt","--transcript-file-match", type = "character", default = '', help="Characters to be used to specify transcript files in raw data folder")
    parser$add_argument("-fmc","--cell-roi-file-match", type = "character", default = '', help="Characters to be used to specify cell roi files in processed data folder")
    parser$add_argument("-fmi","--image-roi-file-match", type = "character", help="Characters to be used to specify image roi files in processed data folder")
    parser$add_argument("-fms","--sun1-file-match", type = "character", help="Characters to be used to specify sun1 roi files in processed data folder")
    parser$add_argument("-si","--section-indices", nargs="+", type="character", default=c("A1","B1","C1","D1","A2","B2","C2","D2"), help="List of section indices to process")
    parser$add_argument("-ri","--run-index", type="character", help="List of section indices to process")

    return(parser)
}


main <- function(){
    # Create parser and parse arguments
    parser <- create_parser()
    args <- parser$parse_args()

    # Retrieve arguments into corresponding variables
    path_to_output_folder <- args$output_folder_path
    transcript_file_match <- args$transcript_file_match
    cell_roi_file_match <- args$cell_roi_file_match
    image_roi_file_match <- args$image_roi_file_match
    sun1_file_match <- args$sun1_file_match
    section_indices <- args$section_indices
    run_index <- args$run_index

    if (args$raw_folder_path == "default"){
        raw_folder_path <- paste0("./raw_data/",run_index)
    } else {
        raw_folder_path <- args$raw_folder_path
    }

    if (args$processed_folder_path == "default"){
        processed_folder_path <- paste0("./processed_data/",run_index)
    } else {
        processed_folder_path <- args$processed_folder_path
    }

    # Load list of transcripts files and ROIs files
    transcripts_files <- list.files(path = raw_folder_path) %>% grep(transcript_file_match, ., value = TRUE)
    cell_rois_files <- list.files(path = processed_folder_path) %>% grep(cell_roi_file_match, ., value = TRUE)
    image_rois_files <- list.files(path = processed_folder_path) %>% grep(image_roi_file_match, ., value = TRUE)
    sun1_rois_files <- list.files(path = processed_folder_path) %>% grep(sun1_file_match, ., value = TRUE)

    cat(paste0("Loading spatial object from ",processed_folder_path,"/",run_index,"_seurat_raw.qs","\n"))
    # Load spatial seurat object to use for projection calls
    spatial_object <- qread(paste0(processed_folder_path,"/",run_index,"_seurat_raw.qs"))

    # Add empty meta data for projection calls
    spatial_object$sun1 <- FALSE
    spatial_object$virus1 <- NA
    spatial_object$virus2 <- NA

    t <- 0
    for (section_index in section_indices){
        cat(paste0("Generating projection calls for section ",section_index,"\n"))
        # Load path to transcript file and ROIs files corresponding to current section index
        transcript_path <- paste0(raw_folder_path,"/",grep(section_index, transcripts_files, value = TRUE))
        cell_roi_path <- paste0(processed_folder_path,"/",grep(section_index, cell_rois_files, value = TRUE))
        image_roi_path <- paste0(processed_folder_path,"/",grep(section_index, image_rois_files, value = TRUE))
        sun1_roi_path <- paste0(processed_folder_path,"/",grep(section_index, sun1_rois_files, value = TRUE))

        cat(paste0("Loading transcripts from ",transcript_path,"\n"))
         # Load transcripts
        transcripts <- load_transcripts(transcript_path)

        cat(paste0("Loading cell ROIs from ",cell_roi_path,"\n"))
        # Load cell ROIs
        rois_cells <- array_to_sf_rois(cell_roi_path)

        # Load image ROIs
        cat(paste0("Loading image ROIs from ", image_roi_path,"\n"))
        rois_image <- array_to_sf_rois(image_roi_path)

        # Load sun1 ROIs
        cat(paste0("Loading sun1 ROIs from ", sun1_roi_path,"\n"))
        rois_sun1 <- array_to_sf_rois(sun1_roi_path)

        #sun1_bol <- sun1_call(sun1_rois, spatial_object_sf, section_index)

        # Generate sun1 calls
        cat(paste0("Generating sun1 call for section ", section_index,"\n"))
        sun1_bol <- sun1_call(rois_sun1, rois_cells)
        
        #Convert cell IDs to seurat format
        names(sun1_bol) <- section_index_to_seurat_integer(spatial_object, section_index) %>% paste0(names(sun1_bol),"_",.)

        # Add sun1 call to meta data
        spatial_object <- AddMetaData(spatial_object,sun1_bol,"sun1")

        cat(paste0("Generating virus calls for section ", section_index,"\n"))
        # Generate virus calls for virus 1
        virus1_call <- virus_call(spatial_object, rois_image, transcripts, section_index, "p147_mCherry_2A_iCre")

        # Add virus calls for virus 1 to meta data
        spatial_object <- AddMetaData(spatial_object,virus1_call,"virus1")
        
        # Generate virus calls for virus 2
        virus2_call <- virus_call(spatial_object, rois_image, transcripts, section_index, "p36_NLS_Cre")

        # Add virus calls for virus 2 to meta data
        spatial_object <- AddMetaData(spatial_object,virus2_call,"virus2")

        t <- t + 1
        cat(paste0(t,"/",length(section_indices),": ",section_index," done!","\n"))
    }

    cat(paste0("Saving spatial seurat object with projection calls to ", path_to_output_folder,"/",run_index,"/",run_index,"_seurat_raw_projection_call.qs","\n"))
    dir.create(paste0(path_to_output_folder,"/",run_index),showWarnings = FALSE)
    qsave(spatial_object,paste0(path_to_output_folder,"/",run_index,"/",run_index,"_seurat_raw_projection_call.qs"))
}

main()