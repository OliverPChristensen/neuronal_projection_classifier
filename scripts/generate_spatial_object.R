library(reticulate)
library(sf)
library(tidyverse)
library(Seurat)
library(qs)
library(argparse)
library(RImageJROI)


current_time <- function(){
    # Input: NA
    #
    # Extract and format current time
    #
    # Output: Current time

    # Extract and format current time
    return(paste0("[",format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),"]"))
}

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
    # Output: sfc object containing ROIs as POLYGON geometries

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

regions_to_sf_rois <- function(regions_folder_path){
    # Input: Path to folder with ImageJ ROI files of regions of a section
    #
    # Loads ImageJ ROI files from folder and converts them to sf POLYGONs in a sf object
    #
    # Output: sf object containing region ROIs as POLYGON geometries
    
    # Get a list of region file names in the regions folder
    region_files <- list.files(path = regions_folder_path)

    # Make empty list and vector to place sf POLYGON and corresponding region name respectively
    region_list <- list()
    region_name <- c()

    # Loop over all files in the regions folder
    for (region_file in region_files){
      
      # Read ImageJ ROI file
      region <- read.ijroi(paste0(regions_folder_path,"/",region_file))

      # Extract coordinates for region ROI and make the coordinates closed
      region_coords <- region$coords
      region_coords_closed <- rbind(region_coords,region_coords[1,])

      # Convert the coordinates into a POLYGON in sf
      region_list[[region_file]] <- st_polygon(list(region_coords_closed))

      # Add the region name for the current region
      region_name <- c(region_name,gsub(".roi","",region_file))
    }

    #Create sf object with region ROIs and their region name
    region_rois <- st_sf(region_name = region_name, geometry = region_list)

    return(region_rois)
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
    transcripts_rois_index <- unlist(sparse_intersects)

    # Convert the cell index to cell IDs retrived from the sfc object
    cell_id <- names(rois)[transcripts_rois_index]

    # Then the transcript-cell pairs are counted using table
    sparse_count_matrix <- as(table(names(transcripts_rois_index),cell_id),"sparseMatrix")

    return(sparse_count_matrix)
}

create_seurat_segmentations <- function(rois,sparse_count_matrix){
    # Input: ROIs as sfc object with POLYGON geometries, sparse count matrix with transcripts as rows and cells as columns
    #
    # Converts the sfc ROIs object to a data frame with column names that fits expected format by Seurat.
    # Filters cells that do not intersect any transcripts
    #
    # Output: Data frame with ROIs in Seurat format

    # Retrieve coordinates from sfc object, transfer cell IDs from sfc object and filter cells that did not intersect any transcripts
    rois_df <- rois %>% st_coordinates() %>% as.data.frame() %>% mutate(L2 = names(rois)[L2]) %>% filter(L2 %in% colnames(sparse_count_matrix))

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
    seurat_centroids <- data.frame(cell = names(rois), centroids_df) %>% filter(cell %in% colnames(sparse_count_matrix))

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

extract_corrected_area <- function(rois_cells, rois_image){
    # Input: Cell ROIs as sfc object with POLYGON geometries, image ROIs as sfc object with POLYGON geometries
    #
    # Finds the part of the cell ROIs that are outside of the image and substract the area from the total area
    # This corrects for the overestimated cell areas caused by grid fill
    #
    # Output: Vector with corrected area and cell IDs as vector names

    # Combine image ROIs into a single sf geometry
    rois_image <- st_combine(rois_image)

    # Make sf object with cell names and geometries from cell ROIs
    rois_cells <- st_sf(name = names(rois_cells), rois_cells)

    # Calculate set difference between cell ROIs and image ROIs
    rois_diff <- st_difference(rois_cells,rois_image)

    # Extract area from the difference between ROIs and name area according to cell ID
    rois_diff %>% st_area -> area_diff
    names(area_diff) <- rois_diff$name

    # Extract area from cell ROIs and name area according to cell ID
    rois_cells %>% st_area -> area_cells
    names(area_cells) <- rois_cells$name

    # Reorder areas of differences using the order in areas of cells
    area_diff <- area_diff[names(area_cells)]

    # Set NAs to zero. This correspands to not substracting any area from cells that was completely contained within the image ROIs
    # and whose area therefore do not need to be corrected
    area_diff[is.na(area_diff)] <- 0

    # Calculate corrected area
    area_corrected <- area_cells - area_diff

    return(area_corrected)
}

region_annotation <- function(rois_cells, rois_regions){
    # Input: Cell ROIs as sfc object with POLYGON geometries, region ROIs as sf object with POLYGON geometries
    #
    # Annotates the anatomical region location of cells by intersecting the region ROIs with the centroids of the cell ROIs
    #
    # Output: Vector with region annotation and cell IDs as vector names

    # Get centroids from cells
    centroid_cells <- rois_cells %>% st_centroid()

    # Get index of regions that intersect a given cell. Cell by region pairs
    sparse_intersects <- st_intersects(centroid_cells,rois_regions)

    # Name the intersects list according to cell IDs from cells and convert the sparse list into a vector
    names(sparse_intersects) <- names(rois_cells)
    cell_by_region_index <- unlist(sparse_intersects)

    # Convert from region index to the corresponding region name
    cell_by_region <- rois_regions$region_name[cell_by_region_index]

    # Give transfer cell IDs to the cell by region vector
    names(cell_by_region) <- names(cell_by_region_index)

    return(cell_by_region)
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
    transcript_path <- args$transcript_path
    cell_roi_path <- args$cell_roi_path
    image_roi_path <- args$image_roi_path
    region_path <- args$region_path
    run_index <- args$run_index
    section_index <- args$section_index
        
    cat(paste0(current_time()," Generating spatial object for section ",section_index,"\n"))

    # Define boolean to indicate if area correction and region annotation should be performed respectively
    is_image = !is.null(image_roi_path)
    cat(paste0(current_time()," Performing area correction: ", is_image,"\n"))
    is_region = !is.null(region_path)
    cat(paste0(current_time()," Performing region annotation: ", is_region,"\n"))

    # Load transcripts
    cat(paste0(current_time()," Loading transcripts from ",transcript_path,"\n"))
    transcript_points <- load_transcripts(transcript_path)
    
    # Load cell ROIs
    cat(paste0(current_time()," Loading cell ROIs from ",cell_roi_path,"\n"))
    rois_cells <- array_to_sf_rois(cell_roi_path)

    if (is_image){
        # Load image ROIs
        cat(paste0(current_time()," Loading image ROIs from ", image_roi_path,"\n"))
        rois_image <- array_to_sf_rois(image_roi_path)
    }
    
    if (is_region){
        # Load image ROIs
        cat(paste0(current_time()," Loading region ROIs from ", region_path,"\n"))
        rois_regions <- regions_to_sf_rois(region_path)
    }
    
    # Create sparse count matrix
    cat(paste0(current_time()," Creating count matrix\n"))
    sparse_count_matrix <- create_count_matrix(transcript_points,rois_cells)

    # Create Seurat segmentation object
    seurat_segmentations <- create_seurat_segmentations(rois_cells,sparse_count_matrix)

    # Create Seurat centroid object
    seurat_centroids <- create_seurat_centroids(rois_cells,sparse_count_matrix)

    # Create Seurat molecules object
    seurat_molecules <- create_seurat_molecules(transcript_points)

    # Create spatial Seurat object
    cat(paste0(current_time()," Creating seurat object\n"))
    spatial_obj <- create_seurat_object(sparse_count_matrix,seurat_segmentations,seurat_centroids,seurat_molecules,section_index)
    
    # Add section index and run index as meta data
    cat(paste0(current_time()," Adding meta data:\n"))
    cat(paste0(current_time(),"  - Section index\n"))
    spatial_obj$section_index <- section_index
    
    cat(paste0(current_time(),"  - Run index\n"))
    spatial_obj$run_index <- run_index

    cat(paste0(current_time(),"  - Section count\n"))
    spatial_obj$section_count <- nrow(transcript_points)

    if (is_image){
        # Extract area of cells and correct for overestimation from grid fill
        cat(paste0(current_time(),"  - Corrected area\n"))
        areas <- extract_corrected_area(rois_cells,rois_image)

        # Calculate total area in the section
        cat(paste0(current_time(),"  - Section area\n"))
        spatial_obj$section_area <- rois_image %>% st_combine() %>% st_area()

    } else {
        # Extract area of cells and give cell ID names to the object
        cat(paste0(current_time(),"  - Area\n"))
        areas <- rois_cells %>% st_area()
        names(areas) <- names(rois_cells)
    }
    
    # Add the area meta data to the spatial seurat object according to cell IDs
    spatial_obj <- AddMetaData(spatial_obj,areas,"area")

    if (is_region){
        # Annotate each cell with a region
        cat(paste0(current_time(),"  - Region annotation\n"))
        regions <- region_annotation(rois_cells, rois_regions)
        # Add the region meta data to the spatial seurat object according to cell IDs
        spatial_obj$region <- NA
        spatial_obj <- AddMetaData(spatial_obj,regions,"region")
    }

    cat(paste0(current_time()," Saving spatial seurat object to ", path_to_output_folder,"/",section_index,"_spatial_object_raw.qs","\n"))
    dir.create(paste0(path_to_output_folder),showWarnings = FALSE)
    qsave(spatial_obj,paste0(path_to_output_folder,"/",section_index,"_spatial_object_raw.qs"))
}

main()