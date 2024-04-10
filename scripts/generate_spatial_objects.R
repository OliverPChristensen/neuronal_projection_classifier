library(reticulate)
library(sf)
library(tidyverse)
library(Seurat)
library(qs)
library(argparse)


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

    # Convert the list of polygons to sfc object and return the object
    return(rois %>% st_sfc())
}

create_count_matrix <- function(transcripts, rois){
    # Input: Transcripts as sf object with POINT geometries, ROIs as sfc object with POLYGON geometries
    #
    # Intersects the transcripts and ROIs objects which return a sparse list (sgbp) of individuel transcript-cell pairs.
    # This list of pairs converted to a sparse matrix with counts of each transcript-cell pair
    #
    # Output: Sparse count matrix with cells as columns and genes as rows

    # Intersect transcripts and ROIs objects
    sparse_intersects <- st_intersects(transcripts,rois,sparse = T)

    # Name the transcripts of the intersects list
    names(sparse_intersects) <- transcripts$V4

    # Convert the intersect list to vectors. This process only returns non-empty elements in the intersect list (that is, transcripts that intersects a cell)
    # Then the transcript-cell pairs are counted using table
    sparse_count_matrix <- as(table(names(unlist(sparse_intersects)),unlist(sparse_intersects)),"sparseMatrix")
    
    # Give cells a cell ID
    colnames(sparse_count_matrix) <- paste0("Cell",colnames(sparse_count_matrix))

    return(sparse_count_matrix)
}

create_seurat_segmentations <- function(rois,sparse_count_matrix){
    # Input: ROIs as sfc object with POLYGON geometries, sparse count matrix with transcripts as rows and cells as columns
    #
    # Converts the sfc ROIs object to a data frame with column names that fits expected format by Seurat.
    # Filters cells that do not intersect any transcripts
    #
    # Output: Data frame with rois in Seurat format

    # Retrieve coordinates from sfc object, create cell IDs and filter cells that did not intersect any transcripts
    rois_df <- rois %>% st_coordinates() %>% as.data.frame() %>% mutate(L2 = paste0("Cell",L2)) %>% filter(L2 %in% colnames(sparse_count_matrix))

    # Extract information and change column names according to the Seurat format
    seurat_segmentations <- rois_df[,c("L2","X","Y")]
    colnames(seurat_segmentations) <- c("cell","x","y")

    return(seurat_segmentations)
}

create_seurat_centroids <- function(rois,sparse_count_matrix){
    # Input: ROIs as sfc object with POLYGON geometries, sparse count matrix with transcripts as rows and cells as columns
    #
    # Extracts the centroids from the cell ROIs as data frame with format that fits expected by Seurat.
    #
    # Output: Data frame with centroids of cells in Seurat format

    # Extract centroids from ROIs
    centroids_df <- rois %>% st_centroid() %>% st_coordinates()

    # Add cell IDs and filter cells that did not intersect any transcripts
    seurat_centroids <- cbind(1:nrow(centroids_df),centroids_df) %>% as.data.frame() %>% mutate(V1 = paste0("Cell",V1)) %>% filter(V1 %in% colnames(sparse_count_matrix))

    # Change column names according to the Seurat format
    colnames(seurat_centroids) <- c("cell","x","y")

    return(seurat_centroids)
}

create_seurat_molecules <- function(transcripts){
    # Input: Transcripts as sf object with POINT geometries
    #
    # Converts sf transcript object to a data frame with format that fits expected by Seurat. 
    #
    # Output: Data frame with transcripts in Seurat format

    # Extract coordinates of transcripts from sf object and add trasncript ID
    seurat_molecules <- transcripts %>% st_coordinates() %>% as.data.frame() %>% cbind(.,transcripts$V4)

    # Change column names according to the Seurat format
    colnames(seurat_molecules) <- c("x","y","gene")

    return(seurat_molecules)
}

create_seurat_object <- function(sparse_count_matrix, seurat_segmentations, seurat_centroids, seurat_molecules, section){
    # Input: Sparse count matrix with transcripts as rows and cells as columns, data frame with cell segmetation in Seurat format, 
    #        data frame with cell centroids in Seurat format, data frame with transcripts (molecules) in Seurat format, section index
    #
    # Creates spatial Seurat object using the sparse count matrix and various data frames containing segmentation, centroids and transcripts in Seurat format
    #
    # Output: Spatial Seurat object

    # Create Seurat object from the sparse count matrix
    spatial_obj <- SeuratObject::CreateSeuratObject(counts = sparse_count_matrix, assay = 'raw')

    # Create coords list containing Seurat centroids object and Seurat segmentation object form the Seurat centroids data frame and Seurat segmentation data frame
    coords_list <- list(
       "centroids" = SeuratObject::CreateCentroids(seurat_centroids),
       "segmentation" = SeuratObject::CreateSegmentation(seurat_segmentations)
    )

    # Create Seurat FOV from the coords list and the Seurat transcripts data frame
    coords <- SeuratObject::CreateFOV(
       coords = coords_list,
       type = c("segmentation", "centroids"),
       molecules = seurat_molecules,
       assay = 'raw'
    )

    # Add FOV to the Seurat object effectively making it a spatial Seurat object
    spatial_obj[[section]] <- coords

    return(spatial_obj)
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
    parser$add_argument("-fpt","--transcripts-folder-path", type = "character", default = 'default', help="Relative path to folder with transcripts files")
    parser$add_argument("-fpr","--roi-folder-path", type = "character", default = 'default', help="Relative path to folder with ROIs in numpy format")
    parser$add_argument("-fmt","--transcript-file-match", type = "character", default = '', help="Characters to be used to specify transcript files in transcript folder")
    parser$add_argument("-fmr","--roi-file-match", type = "character", default = '', help="Characters to be used to specify roi files in roi folder")
    parser$add_argument("-o","--output-folder-path", type = "character", default = './processed_data', help = "Relative path to folder to save output")
    parser$add_argument("-si","--section-indices", nargs="+", type="character", default=c("A1","B1","C1","D1","A2","B2","C2","D2"), help="List of section indices to process")
    parser$add_argument("-ri","--run-index", type="character", help="List of section indices to process")

    return(parser)
}

main <- function(){
    # Input: NA
    #
    # Main function that defines parser argument, loads transcripts and ROIs, creates count matrix, creates spatial seurat object, merges seurat objects across sections, saves the merged spatial seurat object
    #
    # Output: NA

    # Create parser and parse arguments
    parser <- create_parser()
    args <- parser$parse_args()

    # Retrieve arguments into corresponding variables
    transcript_file_match <- args$transcript_file_match
    roi_file_match <- args$roi_file_match
    path_to_output_folder <- args$output_folder_path
    section_indices <- args$section_indices
    run_index <- args$run_index

    if (args$transcripts_folder_path == "default"){
        path_to_transcripts_folder <- paste0("./raw_data/",run_index)
    } else {
        path_to_transcripts_folder <- args$transcripts_folder_path
    }

    if (args$roi_folder_path == "default"){
        path_to_roi_folder <- paste0("./processed_data/",run_index)
    } else {
        path_to_roi_folder <- args$roi_folder_path
    }

    # Load list of transcripts files and ROIs files
    transcripts_files <- list.files(path = path_to_transcripts_folder) %>% grep(transcript_file_match, ., value = TRUE)
    rois_files <- list.files(path = path_to_roi_folder) %>% grep(roi_file_match, ., value = TRUE)

    spatial_objects = list()

    t <- 0
    for (section_index in section_indices){

        # Load path to transcript file and ROIs file corresponding to current section index
        transcript_path <- paste0(path_to_transcripts_folder,"/",grep(section_index, transcripts_files, value = TRUE))
        roi_path <- paste0(path_to_roi_folder,"/",grep(section_index, rois_files, value = TRUE))

        # Load transcripts
        print("Loading transcripts")
        transcript_points <- load_transcripts(transcript_path)
        
        # Load ROIs
        print("Loading ROIs")
        rois_polygons <- array_to_sf_rois(roi_path)

        # Create sparse count matrix
        print("Creating count matrix")
        sparse_count_matrix <- create_count_matrix(transcript_points,rois_polygons)

        # Create Seurat segmentation object
        seurat_segmentations <- create_seurat_segmentations(rois_polygons,sparse_count_matrix)

        # Create Seurat centroid object
        seurat_centroids <- create_seurat_centroids(rois_polygons,sparse_count_matrix)

        # Create Seurat molecules object
        seurat_molecules <- create_seurat_molecules(transcript_points)

        # Create spatial Seurat object
        print("Creating seurat object")
        spatial_obj <- create_seurat_object(sparse_count_matrix,seurat_segmentations,seurat_centroids,seurat_molecules,section_index)
        
        # Add section index and run index as meta data
        print("Adding meta data")
        spatial_obj$section_index <- section_index
        spatial_obj$run_index <- run_index

        # Extract area of cells and give cell ID names to the object
        areas <- rois_polygons %>% st_area()
        names(areas) <- names(rois_polygons)
        
        # Add the meta data to the spatial seurat object according to cell IDs
        spatial_obj <- AddMetaData(spatial_obj,areas,"area")

        # Add the spatial object of the current section to the spatial object list
        spatial_objects[[section_index]] <- spatial_obj

        t <- t + 1
        print(paste0(t,"/",length(section_indices),": ",section_index," done!"))
    }

    print("Merging spatial seurat objects")
    merged_spatial_object <- merge(x = spatial_objects[[1]], y = spatial_objects[-1])

    print("Saving merged spatial seurat object")
    dir.create(paste0(path_to_output_folder,"/",run_index),showWarnings = FALSE)
    qsave(merged_spatial_object,paste0(path_to_output_folder,"/",run_index,"/",run_index,"_seurat_raw.qs"))
}

main()