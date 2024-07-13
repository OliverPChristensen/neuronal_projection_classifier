library(sf)
library(reticulate)
library(tidyverse)
library(RImageJROI)

utils <- new.env()
source("scripts/utils/utils.R", local = utils)

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


seurat_to_sf <- function(spatial_object){
        
    #Loops that fetch the cell segmentations from the Seurat object and put them into a sf geometry
    cell_rois <- list()

    for (fov in names(spatial_object@images)){

        cat(paste0(utils$current_time(), " Converting fov ",fov," to sf\n"))
        
        all_coords <- GetTissueCoordinates(spatial_object@images[[fov]], which = "segmentation")
        
        all_coords_matrix <- as.matrix(all_coords[,c("x","y")])
        
        for (cell_name in unique(all_coords$cell)){
            cell_coords <- all_coords_matrix[all_coords$cell == cell_name,]
            cell_rois[[cell_name]] <- st_polygon(list(cell_coords))
        }

    }

    cell_IDs <- data.frame(cell_IDs = names(cell_rois))
    cell_rois <- st_sf(cell_IDs, geometry = cell_rois)

    meta_data <- spatial_object@meta.data
    meta_data$cell_IDs <- rownames(meta_data)

    spatial_object_sf <- left_join(cell_rois, meta_data, by = "cell_IDs")

    return(spatial_object_sf)
}


rotate_sf_object <- function(sf_object, rotation_df){
    
    run_indices <- sf_object$run_index %>% unique()
    section_indices <- sf_object$section_index %>% unique()

    for (run_index_ in run_indices){
        for (section_index_ in section_indices){

            angle_degrees <- rotation_df %>% 
                filter(run_index == run_index_, section_index == section_index_) %>%
                .[["Rotation"]]

            angle_radians <- angle_degrees * pi / 180

            rotation_matrix <- matrix(
                c(cos(angle_radians), -sin(angle_radians), sin(angle_radians), cos(angle_radians)),
                ncol = 2, 
                byrow = TRUE
            )

            geometry <- sf_object %>% 
                filter(run_index == run_index_, section_index == section_index_) %>%
                st_geometry()

            st_geometry(sf_object[sf_object$run_index == run_index_ & sf_object$section_index == section_index_,]) <- geometry * rotation_matrix
        }
    }

    return(sf_object)
}


load_regions <- function(region_paths, run_section_index_pairs){

    region_sf_object <- data.frame()

    for (row_num in 1:nrow(run_section_index_pairs)){
        run_index <- run_section_index_pairs[row_num, "Var1"]
        section_index <- run_section_index_pairs[row_num, "Var2"]
        cat(paste0(utils$current_time(),"  - Create region object for ",run_index," ",section_index,"\n"))

        region_path <- grep(run_index, region_paths, ignore.case = TRUE, value = TRUE) %>% 
            grep(section_index, ., value = TRUE)

        regions <- sf_utils$regions_to_sf_rois(region_path)

        regions$region <- correct_region_names(regions$region_name)

        regions$run_index <- run_index
        regions$section_index <- section_index

        region_sf_object <- rbind(region_sf_object, regions)
    }

    return(region_sf_object)
}