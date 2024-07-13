library(tidyverse)
library(Seurat)
library(sf)
library(qs)
library(argparse)
library(RImageJROI)
library(Matrix)

utils <- new.env()
source("scripts/utils/utils.R", local = utils)

sf_utils <- new.env()
source("scripts/utils/sf_utils.R", local = sf_utils)


create_menu_converter <- function(bregma_tsv_path){

    menu_converter <- read.table(bregma_tsv_path, sep = "\t", header = TRUE)

    menu_converter <- menu_converter %>% 
        mutate(
            run_index = paste0("spatial",Spatial_run),
            section_index = Well, 
            bregma_level = paste0(Bregma," (",Orientation,")")
        )

    menu_converter$projection <- case_when(
        menu_converter$run_index == "spatial1" ~ "BST and PAG",
        menu_converter$run_index == "spatial2" ~ "LHA and PBN",
        menu_converter$run_index == "spatial3" ~ "PVT and CEA"
    )

    return(menu_converter)
}


format_p_value <- function(x) {

  ifelse(abs(x) < 1e-2, format(x, scientific = TRUE, digits = 3), round(x, 2))

}

add_count_matrix <- function(sf_object, spatial_object){

    count_layers <- grep("counts", Layers(spatial_object), value = TRUE)
    count_layers <- count_layers[!grepl("virus",count_layers)]

    total_count_df <- data.frame(transcript_ID = gsub("-","_",rownames(spatial_object)))

    for (count_layer in count_layers){

        count_df <- spatial_object[["raw"]][count_layer] %>% as.data.frame()

        count_df$transcript_ID <- gsub("-","_",rownames(count_df))

        total_count_df <- left_join(total_count_df, count_df, by = "transcript_ID")
    }

    total_count_df <- total_count_df %>% t() %>% as.data.frame()
    colnames(total_count_df) <- total_count_df[rownames(total_count_df) == "transcript_ID",]
    total_count_df$cell_IDs <- rownames(total_count_df)

    sf_object <- left_join(sf_object, total_count_df, by = "cell_IDs")
    
    return(sf_object)
}


prepare_cell_sf_object_centroid <- function(cell_sf_object, spatial_object){

    cell_sf_object_centroid <- cell_sf_object %>% st_centroid()


    cell_sf_object_centroid$virus1 <- format_p_value(cell_sf_object_centroid$virus1)
    cell_sf_object_centroid$virus2 <- format_p_value(cell_sf_object_centroid$virus2)
    cell_sf_object_centroid$area <- round(cell_sf_object_centroid$area)

    cell_sf_object_centroid_selected <- cell_sf_object_centroid %>% select(run_index, section_index, virus1, virus2, area, geometry, cell_IDs) %>% 
        rename(Virus1 = virus1, Virus2 = virus2, "Cell area" = area)

    cell_sf_object_centroid_selected <- add_count_matrix(cell_sf_object_centroid_selected, spatial_object)
    cell_sf_object_centroid_selected <- cell_sf_object_centroid_selected %>% select(-("cell_IDs"))

    return(cell_sf_object_centroid_selected)
}

simplify_cell_object <- function(cell_sf_object){
    
    cell_sf_object <- cell_sf_object %>% 
        st_simplify(dTolerance = 2)

    removed_cell_count <- cell_sf_object %>% filter(st_geometry_type(cell_sf_object) == "MULTIPOLYGON") %>%
        st_drop_geometry() %>%
        group_by(run_index, section_index) %>%
        summarize(counts = n())

    cat(paste0(utils$current_time(), " Number of cells removed after simplify:\n"))
    cat(paste0(utils$current_time(), " ", print(removed_cell_count, n = nrow(removed_cell_count)),"\n"))

    cell_sf_object <- cell_sf_object %>% filter(st_geometry_type(cell_sf_object) != "MULTIPOLYGON")

    return(cell_sf_object)
}



create_transcript_sf_object <- function(transcript_paths, run_index, section_index, transcript_sf_object){

    transcripts_path <- grep(run_index, transcript_paths, ignore.case = TRUE, value = TRUE) %>% 
        grep(section_index, ., value = TRUE)

    transcripts <- sf_utils$load_transcripts(transcripts_path)
    
    transcripts_multipoint <- transcripts %>%
        group_by(V4) %>%
        summarize(geometry = st_combine(geometry))

    transcripts_multipoint$run_index <- run_index
    transcripts_multipoint$section_index <- section_index
    transcript_sf_object <- rbind(transcript_sf_object, transcripts_multipoint)

    return(list(
        transcripts = transcripts,
        transcript_sf_object = transcript_sf_object
    ))
}

create_downsampled_transcript_sf_object <- function(transcripts, run_index, section_index, transcript_sf_object_downsampled){

    transcripts_downsampled <- transcripts %>%
        group_by(V4) %>% 
        slice_sample(n = 10000)

    transcripts_downsampled_multipoint <- transcripts_downsampled %>%
        group_by(V4) %>%
        summarize(geometry = st_combine(geometry))

    transcripts_downsampled_multipoint$run_index <- run_index
    transcripts_downsampled_multipoint$section_index <- section_index

    transcript_sf_object_downsampled <- rbind(transcript_sf_object_downsampled, transcripts_downsampled_multipoint)

    return(transcript_sf_object_downsampled)
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


create_region_sf_object <- function(region_paths, run_index, section_index, region_sf_object){

    region_path <- grep(run_index, region_paths, ignore.case = TRUE, value = TRUE) %>% 
        grep(section_index, ., value = TRUE)

    regions <- sf_utils$regions_to_sf_rois(region_path)

    regions$region_name <- correct_region_names(regions$region_name)

    regions$run_index <- run_index
    regions$section_index <- section_index

    region_sf_object <- rbind(region_sf_object, regions)

    return(list(
        regions = regions,
        region_sf_object = region_sf_object
    ))
}


transcripts_regions_intersects <- function(transcripts, regions){

    sparse_intersects <- st_intersects(transcripts,regions,sparse = TRUE)

    # Name the transcripts of the intersects list
    names(sparse_intersects) <- transcripts$V4
    
    # Convert the intersect list to vectors. This process only returns non-empty elements in the intersect list (that is, transcripts that intersects a cell)
    transcript_region_index <- unlist(sparse_intersects)
    
    # Convert the cell index to cell IDs retrived from the sfc object
    region <- regions$region_name[transcript_region_index]

    transcript_region <- data.frame(transcript = names(transcript_region_index),region)

    transcript_region_counts <- transcript_region %>% 
        group_by(transcript, region) %>% 
        summarize(counts = n())
    
    return(transcript_region_counts)
}


create_transcript_region_count_object <- function(transcripts, regions, run_index, section_index, transcript_region_counts_object){

    transcripts_region_counts <- transcripts_regions_intersects(transcripts, regions)
    
    transcripts_region_counts$run_index <- run_index
    transcripts_region_counts$section_index <- section_index
    
    transcript_region_counts_object <- rbind(transcript_region_counts_object, transcripts_region_counts)

    return(transcript_region_counts_object)
}

create_parser <- function(){
    # Input: NA
    #
    # Creates and outputs parser with arguments that allows interaction via CLI
    #
    # Output: Parser

    # Define parser
    parser <- ArgumentParser(description = "Script to generate R shiny data for Spatial Hypothalamus App")

    # Add argument to parser
    parser$add_argument("-o","--output-folder-path", type = "character", help = "Relative path to folder to save output")
    parser$add_argument("-sop","--spatial-object-path", type = "character", help="Relative path to spatial object")
    parser$add_argument("-tps","--transcript-paths", nargs="+", type = "character", help="Relative paths to transcripts files in Resolve format")
    parser$add_argument("-rps","--region-paths", nargs="+", type = "character", help="Relative paths to folders with region ROIs in imageJ format")
    parser$add_argument("-btp","--bregma-tsv-path", type = "character", help="Relative path to tsv with bregma levels for each section")
    
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
    spatial_object_path <- args$spatial_object_path
    transcript_paths <- args$transcript_paths
    region_paths <- args$region_paths
    bregma_tsv_path <- args$bregma_tsv_path


    cat(paste0(utils$current_time()," Load spatial Seurat object from ",spatial_object_path,"\n"))
    spatial_object <- qread(spatial_object_path)

    cat(paste0(utils$current_time()," Create menu converter\n"))
    menu_converter <- create_menu_converter(bregma_tsv_path)

    cat(paste0(utils$current_time()," Create cell object from spatial seurat object\n"))
    cell_sf_object <- sf_utils$seurat_to_sf(spatial_object)

    cat(paste0(utils$current_time()," Create cell centroid object\n"))
    cell_sf_object_centroid <- prepare_cell_sf_object_centroid(cell_sf_object, spatial_object)

    cat(paste0(utils$current_time()," Simplify cell object\n"))
    cell_sf_object <- simplify_cell_object(cell_sf_object)

    cell_sf_object <- cell_sf_object %>% select(nCount_raw, run_index, section_index, area, region, projection, virus1, virus2, geometry)

    
    run_indices <- c("spatial1", "spatial2", "spatial3")
    section_indices <- c("A1", "B1", "C1", "D1", "A2", "B2", "C2", "D2")

    transcript_object_list <- list(transcript_sf_object = data.frame())
    transcript_sf_object_downsampled <- data.frame()
    region_object_list <- list(region_sf_object = data.frame())
    transcript_region_counts_object <- data.frame()

    for (run_index in run_indices){
        cat(paste0(utils$current_time()," Create shiny objects for ",run_index,":\n"))
        for (section_index in section_indices){
            
            cat(paste0(utils$current_time(),"  - Create transcript object for ",section_index,"\n"))
            transcript_object_list <- create_transcript_sf_object(transcript_paths, run_index, section_index, transcript_object_list$transcript_sf_object)

            cat(paste0(utils$current_time(),"  - Create downsampled transcript object for ",section_index,"\n"))
            transcript_sf_object_downsampled <- create_downsampled_transcript_sf_object(transcript_object_list$transcripts, run_index, section_index, transcript_sf_object_downsampled)

            cat(paste0(utils$current_time(),"  - Create region object for ",section_index,"\n"))
            region_object_list <- create_region_sf_object(region_paths, run_index, section_index, region_object_list$region_sf_object)

            cat(paste0(utils$current_time(),"  - Create transcript-region-count object for ",section_index,"\n"))
            transcript_region_counts_object <- create_transcript_region_count_object(transcript_object_list$transcripts, region_object_list$regions, run_index, section_index, transcript_region_counts_object)
        }
    }


    cat(paste0(utils$current_time()," Rotate sf objects:\n"))

    cat(paste0(utils$current_time(),"  - cell object\n"))
    cell_sf_object <- sf_utils$rotate_sf_object(cell_sf_object, menu_converter)

    cat(paste0(utils$current_time(), "  - cell centroid object\n"))
    cell_sf_object_centroid <- sf_utils$rotate_sf_object(cell_sf_object_centroid, menu_converter)

    cat(paste0(utils$current_time(), "  - transcript object\n"))
    transcript_sf_object <- sf_utils$rotate_sf_object(transcript_object_list$transcript_sf_object, menu_converter)

    cat(paste0(utils$current_time(), "  - downsampled transcript object\n"))
    transcript_sf_object_downsampled <- sf_utils$rotate_sf_object(transcript_sf_object_downsampled, menu_converter)

    cat(paste0(utils$current_time(), "  - region object\n"))
    region_sf_object <- sf_utils$rotate_sf_object(region_object_list$region_sf_object, menu_converter)
    
    
    

    colnames(transcript_sf_object) <- c("transcript","geometry", "run_index", "section_index")
    colnames(transcript_sf_object_downsampled) <- c("transcript","geometry", "run_index", "section_index")
    colnames(region_sf_object) <- c("region", "run_index", "section_index","geometry")


    menu_converter <- menu_converter %>% 
        select(run_index, section_index, projection, bregma_level)


    shiny_data <- list(
        cell_object = cell_sf_object,
        cell_object_centroid = cell_sf_object_centroid,
        transcript_object = transcript_sf_object,
        transcript_object_downsampled = transcript_sf_object_downsampled,
        region_object = region_sf_object,
        transcript_region_table = transcript_region_counts_object,
        menu_converter = menu_converter
    )


    # Save shiny data
    cat(paste0(utils$current_time()," Saving shiny data to ", path_to_output_folder,"/","shiny_data.qs","\n"))
    dir.create(paste0(path_to_output_folder),showWarnings = FALSE)
    qsave(shiny_data, paste0(path_to_output_folder,"/","shiny_data.qs"))
}

main()