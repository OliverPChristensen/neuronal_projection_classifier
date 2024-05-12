library(reticulate)
library(sf)
library(tidyverse)
library(Seurat)
library(qs)
library(argparse)
library(RImageJROI)


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

    # Checks for valid ROIs. Non valid will not be given an area, but are not valid cells anyways 
    # Maaske skulle jeg aendre saa ikke valid rois fjernes allerede i arrays to sf function
    #rois_cells <- rois_cells[st_is_valid(rois_cells),]
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

    #Get centroids from cells
    centroid_cells <- rois_cells %>% st_centroid()

    #Get index of regions that intersect a given cell. Cell by region pairs
    sparse_intersects <- st_intersects(centroid_cells,rois_regions)

    #Name the intersects list according to cell IDs from cells and convert the sparse list into a vector
    names(sparse_intersects) <- names(rois_cells)
    cell_by_region_index <- unlist(sparse_intersects)

    #Convert from region index to the corresponding region name
    cell_by_region <- rois_regions$region_name[cell_by_region_index]

    #Give transfer cell IDs to the cell by region vector
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

    # Add argumnent to parser
    parser$add_argument("-fpt","--raw-folder-path", type = "character", default = 'default', help="Relative path to raw data folder")
    parser$add_argument("-fpr","--processed-folder-path", type = "character", default = 'default', help="Relative path to processed data folder")
    parser$add_argument("-o","--output-folder-path", type = "character", default = './processed_data', help = "Relative path to folder to save output")
    parser$add_argument("-fmt","--transcript-file-match", type = "character", default = '', help="Characters to be used to specify transcript files in raw data folder")
    parser$add_argument("-fmc","--cell-roi-file-match", type = "character", default = '', help="Characters to be used to specify cell roi files in processed data folder")
    parser$add_argument("-fmi","--image-roi-file-match", type = "character", help="Characters to be used to specify image roi files in processed data folder. If non given, then area correction is not performed. Give empty string to perform area correction with no file match")
    parser$add_argument("-fmr","--region-folder-match", type = "character", help="Characters to be used to specify region folders in raw data folder. If non given, then region annotation is not performed. Give empty string to perform region annotation with no file match")
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
    path_to_output_folder <- args$output_folder_path
    transcript_file_match <- args$transcript_file_match
    cell_roi_file_match <- args$cell_roi_file_match
    image_roi_file_match <- args$image_roi_file_match
    region_folder_match <- args$region_folder_match
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

    #Define boolean to indicate if area correction and region annotation should be performed respectively
    is_image = !is.null(image_roi_file_match)
    cat(paste0("Performing area correction: ", is_image,"\n"))
    is_region = !is.null(region_folder_match)
    cat(paste0("Performing region annotation: ", is_region,"\n"))

    # Load list of transcripts files and ROIs files
    transcripts_files <- list.files(path = raw_folder_path) %>% grep(transcript_file_match, ., value = TRUE)
    cell_rois_files <- list.files(path = processed_folder_path) %>% grep(cell_roi_file_match, ., value = TRUE)

    if (is_image){
        image_rois_files <- list.files(path = processed_folder_path) %>% grep(image_roi_file_match, ., value = TRUE)
    }

    if (is_region){
        region_folders <- list.files(path = raw_folder_path) %>% grep(region_folder_match, ., value = TRUE)
    }

    spatial_objects = list()
    t <- 0

    for (section_index in section_indices){
        
        cat(paste0("Generating spatial object for section ",section_index,"\n"))
        # Load path to transcript file and ROIs file corresponding to current section index
        transcript_path <- paste0(raw_folder_path,"/",grep(section_index, transcripts_files, value = TRUE))
        cell_roi_path <- paste0(processed_folder_path,"/",grep(section_index, cell_rois_files, value = TRUE))

        if (is_image){
            image_roi_path <- paste0(processed_folder_path,"/",grep(section_index, image_rois_files, value = TRUE))
        }

        if (is_region){
            region_path <- paste0(raw_folder_path,"/",grep(section_index, region_folders, value = TRUE))
        }
   
        # Load transcripts
        cat(paste0("Loading transcripts from ",transcript_path,"\n"))
        transcript_points <- load_transcripts(transcript_path)
        
        # Load cell ROIs
        cat(paste0("Loading cell ROIs from ",cell_roi_path,"\n"))
        rois_cells <- array_to_sf_rois(cell_roi_path)

        if (is_image){
            # Load image ROIs
            cat(paste0("Loading image ROIs from ", image_roi_path,"\n"))
            rois_image <- array_to_sf_rois(image_roi_path)
        }
        
        if (is_region){
            # Load image ROIs
            cat(paste0("Loading region ROIs from ", region_path,"\n"))
            rois_regions <- regions_to_sf_rois(region_path)
        }
        
        # Create sparse count matrix
        cat("Creating count matrix\n")
        sparse_count_matrix <- create_count_matrix(transcript_points,rois_cells)

        # Create Seurat segmentation object
        seurat_segmentations <- create_seurat_segmentations(rois_cells,sparse_count_matrix)

        # Create Seurat centroid object
        seurat_centroids <- create_seurat_centroids(rois_cells,sparse_count_matrix)

        # Create Seurat molecules object
        seurat_molecules <- create_seurat_molecules(transcript_points)

        # Create spatial Seurat object
        cat("Creating seurat object\n")
        spatial_obj <- create_seurat_object(sparse_count_matrix,seurat_segmentations,seurat_centroids,seurat_molecules,section_index)
        
        # Add section index and run index as meta data
        cat("Adding meta data:\n")
        cat(" -Section index\n")
        spatial_obj$section_index <- section_index
        cat(" -Run index\n")
        spatial_obj$run_index <- run_index

        if (!is_image){
            # Extract area of cells and give cell ID names to the object
            cat(" -Area\n")
            areas <- rois_cells %>% st_area()
            names(areas) <- names(rois_cells)
        } else {
            # Extract area of cells and correct for overestimation from grid fill
            cat(" -Corrected area\n")
            areas <- extract_corrected_area(rois_cells,rois_image)
        }

        # Add the area meta data to the spatial seurat object according to cell IDs
        spatial_obj <- AddMetaData(spatial_obj,areas,"area")

        if (is_region){
            # Annotate each cell with a region
            cat(" -Region annotation\n")
            regions <- region_annotation(rois_cells, rois_regions)
            # Add the region meta data to the spatial seurat object according to cell IDs
            spatial_obj$region <- NA
            spatial_obj <- AddMetaData(spatial_obj,regions,"region")
        }

        # Add the spatial object of the current section to the spatial object list
        spatial_objects[[section_index]] <- spatial_obj

        t <- t + 1
        cat(paste0(t,"/",length(section_indices),": ",section_index," done!","\n"))
    }

    cat("Merging spatial seurat objects\n")
    merged_spatial_object <- merge(x = spatial_objects[[1]], y = spatial_objects[-1])

    cat(paste0("Saving merged spatial seurat object to ", path_to_output_folder,"/",run_index,"/",run_index,"_seurat_raw.qs","\n"))
    dir.create(paste0(path_to_output_folder,"/",run_index),showWarnings = FALSE)
    qsave(merged_spatial_object,paste0(path_to_output_folder,"/",run_index,"/",run_index,"_seurat_raw.qs"))
}

main()