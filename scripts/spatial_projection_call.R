library(reticulate)
library(sf)
library(tidyverse)
library(Seurat)
library(qs)



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

seurat_to_sf <- function(spatial_object){
    meta_data <- spatial_object@meta.data
        
    #Loops that fetch the cell segmentations from the Seurat object and put them into a sf geometry
    cell_list <- list()
    t <- 1

    for (section_index in names(spatial_object@images)){
        print(paste0("Converting section ",section_index," to sf (",t,"/",length(names(spatial_object@images)),")"))
        
        coords <- GetTissueCoordinates(spatial_object@images[[section_index]], which = "segmentation")
        
        cell_names <- unique(coords$cell)
        
        for (cell_name in cell_names){
            next_coords <- coords[coords$cell == cell_name,]
            cell_list[[cell_name]] <- st_polygon(list(as.matrix(next_coords[,c("x","y")])))
        }

        t <- t + 1
    }

    spatial_object_sf <- st_sf(meta_data, cell_list)

    return(spatial_object_sf)
}

sun1_call <- function(sun1_rois, spatial_object_sf, section_index_sub){
    areas <- sun1_rois %>% st_area()
    sun1_rois_size_excluded <- sun1_rois[areas > 500]

    spatial_object_sf_section <- spatial_object_sf %>% filter(section_index == section_index_sub)
    sparse_intersects <- st_intersects(spatial_object_sf_section,sun1_rois_size_excluded)

    names(sparse_intersects) <- spatial_object_sf_section %>% st_geometry() %>% names() #
    
    #All index are above -1, so effectively turns the integer vector into a boolean
    sun1_bol <- unlist(sparse_intersects) > -1
    return(sun1_bol)

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
    parser$add_argument("-fmr","--roi-file-matches", nargs="+", type = "character", default = c("GFP",".npz"), help="Characters to be used to specify roi files in roi folder")
    parser$add_argument("-o","--output-folder-path", type = "character", default = './processed_data', help = "Relative path to folder to save output")
    parser$add_argument("-si","--section-indices", nargs="+", type="character", default=c("A1","B1","C1","D1","A2","B2","C2","D2"), help="List of section indices to process")
    parser$add_argument("-ri","--run-index", type="character", help="List of section indices to process")

    return(parser)
}


main <- function(){
    # Create parser and parse arguments
    parser <- create_parser()
    args <- parser$parse_args()

    # Retrieve arguments into corresponding variables
    transcript_file_match <- args$transcript_file_match
    roi_file_matches <- args$roi_file_matches #GFP & .npz
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

    #### TEST ####

    transcript_file_match <- "results"
    roi_file_matches <- c("GFP",".npz")
    #section_indices <- c("A1","B1","C1","D1","A2","B2","C2","D2")
    section_index <- c("A1")
    run_index <- "spatial1"

    path_to_transcripts_folder <- paste0("./raw_data/",run_index)
    path_to_roi_folder <- paste0("./processed_data/",run_index)

    qsave(spatial_object_sf, "./processed_data/spatial1/spatial_sf_test.qs")
    spatial_object_sf <- qread("./processed_data/spatial1/spatial_sf_test.qs")
    ##############

    # Load list of transcripts files and ROIs files
    transcripts_files <- list.files(path = path_to_transcripts_folder) %>% grep(transcript_file_match, ., value = TRUE)
    rois_files <- list.files(path = path_to_roi_folder) %>% grep(paste0(roi_file_matches[1],".*",roi_file_matches[2]), ., value = TRUE)
    spatial_object_file <- list.files(path = path_to_roi_folder) %>% grep("seurat", ., value = TRUE)
    spatial_object <- qread(paste0(path_to_roi_folder,"/",spatial_object_file))
    spatial_object_sf <- seurat_to_sf(spatial_object)

    spatial_object$sun1 <- FALSE
    spatial_object$virus <- 0
    for (section_index in section_indices){
        # Load path to transcript file and ROIs file corresponding to current section index
        transcript_path <- paste0(path_to_transcripts_folder,"/",grep(section_index, transcripts_files, value = TRUE))
        roi_path <- paste0(path_to_roi_folder,"/",grep(section_index, rois_files, value = TRUE))

        sun1_rois <- array_to_sf_rois(roi_path)
        sun1_bol <- sun1_call(sun1_rois, spatial_object_sf, section_index)
        spatial_object <- AddMetaData(spatial_object,sun1_bol,"sun1")








    }

}

        area <- sun1_rois %>% st_area() #500
        area %>% as_tibble(x = .) %>% ggplot() + geom_density(aes(x = value))
        bol <- area > 500
        sum(!bol)
        cor_area <- area[bol]
        test <- st_sf(cor_area, geometry = sun1_rois[bol])
        test %>% ggplot() + geom_sf(aes(fill = cor_area))
        ggsave("./plots/spatialcheck.png", width = 6000, height = 7500, units = "px")

test <- "a"
grep(paste0(test,".*b|b.*a"), c("a1b","b1c","c1d"), value = TRUE)
paste0(a,".*",b,"|",b,".*",a)

cell_list <- list()

t <- 1
for (section_index in names(spatial_object@images)){
    coords <- GetTissueCoordinates(spatial_object@images[[section_index]], which = "segmentation")
    
    cell_names <- unique(coords$cell)
    print(t)
    for (cell_name in cell_names){
        next_coords <- coords[coords$cell == cell_name,]
        cell_list[[cell_name]] <- list(as.matrix(next_coords[,c("x","y")]))
    }

    cell_sf <- st_multipolygon(cell_list)
    t <- t + 1
}

#next_coords <- coords[coords$cell == cell_name,]
cell_list[[cell_name]] <- list(as.matrix(coords[,c("x","y")]))
cell_roi <- st_polygon(list(as.matrix(coords[,c("x","y")])))
cell_list[[cell_name]] <- cell_roi