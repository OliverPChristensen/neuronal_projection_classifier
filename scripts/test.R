library(Seurat)
library(qs)
library(tidyverse)
library(sf)
library(RImageJROI)

spatial_object <- qread("./processed_data/spatial1/spatial1_seurat_raw.qs")


DefaultBoundary(spatial_object[["A1"]]) <- "segmentation"

spatial_object@meta.data$area %>% as.numeric %>% quantile(., c(0.1)) -> aq

spatial_object$aq <- spatial_object$area > aq

ImageDimPlot(spatial_object, 
             fov = "A1",
             axes = F, 
             dark.background = T, 
             size = 1.5,
             border.color = "black", 
             border.size = 0.1, 
             cols = "polychrome",
             group.by = "aq",
             #cells = WhichCells(spatial_object, expression = prediction.score.super.major.cell.type > 0.9),
            #molecules = c("Agrp","Glp1r")
)

head(spatial_object)

spatial_object@meta.data %>% ggplot() + geom_density(aes(x = nCount_raw))
spatial_object@meta.data %>% ggplot() + geom_density(aes(x = area))
rownames(spatial_object)
ggsave("./plots/spatialcheck.png", width = 6000, height = 7500, units = "px")


spatial_object@meta.data %>% ggplot() + geom_density(aes(x = area, col = section_index)) 



rownames(spatial_object)

spatial_object@meta.data %>% select(area) %>% sum()

spatial_object@assays


spatial_object@assays$raw$counts.1["p147-mCherry-2A-iCre",] %>% sum()
spatial_object@assays$raw$counts.1["p36-NLS-Cre",] %>% sum()


path_to_regions_rois_folder <- "./raw_data/spatial1"
regions_rois_folder_match <- "regions"
section_indices <- c("A1")
section_index <- "A1"
region_file <- "anterior_hypothalamus_left.roi"

region_rois_to_sf <- function(path_to_regions_rois_folder,regions_rois_folder_match,section_indices){
  sections_folders <- list.files(path = path_to_regions_rois_folder) %>% grep(regions_rois_folder_match, ., value = TRUE)
  
  region_data <- data.frame()
  region_list <- list()
  
  for (section_index in section_indices){

    # Get a list of file names in the folder
    region_files <- list.files(path = paste0(path_to_regions_rois_folder, "/", grep(section_index, sections_folders, value = TRUE)))

    for (region_file in region_files){
        
      roi <- read.ijroi(paste0(path_to_regions_rois_folder,"/",grep(section_index, sections_folders, value = TRUE),"/",region_file))
    
      roi_coords <- roi$coords
      roi_coords_closed <- rbind(roi$coords,roi$coords[1,])
      region_list[[region_file]] <- st_polygon(list(roi_coords_closed))
      region_data <- rbind(region_data,c(gsub(".roi","",region_file),section_index))
    }
    
  }
  colnames(region_data) <- c("region","section")
  
  sf_obj <- st_sf(region_data, geometry = region_list)
  return(sf_obj)
}

section_indices <- c("A1")

regions_sf <- region_rois_to_sf("./raw_data/spatial1","regions", section_indices)


area_region <- regions_sf %>% st_area() %>% sum()

transcripts <- read.table("./raw_data/spatial1/Spatial1_serie2_serie2_A1-1_results.txt", header = F, sep = "\t") %>% filter(!str_detect(V4, "FP "))

virus1_region <- sum(transcripts$V4 == "p147_mCherry_2A_iCre")

virus2_region <- sum(transcripts$V4 == "p36_NLS_Cre")
