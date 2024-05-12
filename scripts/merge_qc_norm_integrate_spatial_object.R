library(argparse)
library(tidyverse)
library(Seurat)
library(qs)

merge_spatial_objects <- function(spatial_paths_list){
    spatial_objects = list()
    t <- 0

    for (i in 1:length(spatial_object_paths)){
        cat(paste0("Loading spatial object ",t,"/",length(spatial_object_paths)," from ",spatial_object_path,"\n"))
        # Load spatial seurat object to use for projection calls
        spatial_object <- qread(spatial_object_path)
    }

    cat("Merging spatial seurat objects\n")
    merged_spatial_object <- merge(x = spatial_objects[[1]], y = spatial_objects[-1])

    return(merged_spatial_object)
}

quality_control <- function(spatial_object){
    area_quantiles <- quantile(spatial_object$area, c(0.01.0,99))
    spatial_object <- subset(spatial_object, subset = nCount_raw > 10 & area > area_quantiles[1] & area < area_quantiles[2])
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
    parser$add_argument("-sop","--spatial-object-paths", nargs="+", type = "character", help="Characters to be used to specify transcript files in raw data folder")

    return(parser)
}


main <- function(){
    # Create parser and parse arguments
    parser <- create_parser()
    args <- parser$parse_args()

    # Retrieve arguments into corresponding variables
    path_to_output_folder <- args$output_folder_path
    spatial_object_paths <- args$spatial_object_paths

    spatial_objects = list()
    t <- 0

    # Load spatial seurat objects from list of paths
    for (spatial_object_path in spatial_object_paths){
        cat(paste0("Loading spatial object ",t,"/",length(spatial_object_paths)," from ",spatial_object_path,"\n"))
        spatial_object <- qread(spatial_object_path)
    }

    cat("Merging spatial seurat objects\n")
    spatial_object <- merge(x = spatial_objects[[1]], y = spatial_objects[-1])
 
    #Quality control
    spatial_object <- quality_control(spatial_object)

    ### Plots of count differences between sections and runs (pga permeabilitetsforskelle, eller forskelle i størrelse af celler?)
    ### Plot area forskelle for at se om det evt forskelle i counts kan skyldes forskelle i area. Evt udregne en mean count pr areal for hver section?
    ### Hvis ikke area og count forskelle saa er area norm klart det bedste. Hvis forskelle i disse saa count normalisering og saa haabe at probe set er represntativt og ikke biased. Vigtigst ift projecting and non projecting neurons
    ### SC tranform vs ikke sc transform? Kommer an paa om den valgte normalisering er effektiv eller om der stadig er effekter der ikke er normallisieret ud

    # Normalisering


    cat(paste0("Saving spatial seurat object with projection calls to ", path_to_output_folder,"/",section_index,"_spatial_object_projection_call.qs","\n"))
    dir.create(paste0(path_to_output_folder),showWarnings = FALSE)
    qsave(spatial_object,paste0(path_to_output_folder,"/",section_index,"_spatial_object_projection_call.qs"))
}

main()