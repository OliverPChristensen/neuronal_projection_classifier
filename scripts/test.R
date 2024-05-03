library(Seurat)
library(qs)
library(tidyverse)
library(sf)
library(RImageJROI)
library(progress)
library(terra)
library(reticulate)

spatial_object <- qread("./processed_data/spatial1/spatial1_seurat_raw.qs")


DefaultBoundary(spatial_object[["A1"]]) <- "segmentation"

spatial_object@meta.data$area %>% quantile(., c(0.9)) -> aq

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


area_cells <- spatial_object@meta.data %>% filter(section_index == "A1") %>% select(area) %>% sum()

virus1_cells <- spatial_object@assays$raw$counts.1["p147-mCherry-2A-iCre",] %>% sum()
virus2_cells <- spatial_object@assays$raw$counts.1["p36-NLS-Cre",] %>% sum()


null_intensity1 <- (virus1_region - virus1_cells)/(area_region - area_cells)

cell_areas <- spatial_object@meta.data %>% filter(section_index == "A1") %>% select(area)

null_expectation1 <- null_intensity1*cell_areas$area

virus1_cell <- spatial_object@assays$raw$counts.1["p147-mCherry-2A-iCre",]

poisson.test(virus1_cell[702], null_expectation1[702], alternative = "greater")$p.value

spatial_object$virus_test <- 0


pb <- progress_bar$new(format = "[:bar] :percent %",
                       total = length(virus1_cell),
                       complete = "=",   # Completion bar character
                       incomplete = "-", # Incomplete bar character
                       current = ">",    # Current bar character
                       clear = FALSE,    # If TRUE, clears the bar when finish
                       width = 100)


virus_test <- rep(0,length(virus1_cell))

for (i in 1:length(virus1_cell)){
  virus_test[i] <- poisson.test(virus1_cell[i], null_expectation1[i], alternative = "greater")$p.value
}


sum(virus_test < 0.05/length(spatial_object$virus_test))
length(virus_test)

quantile(spatial_object$virus_test,1:10/10)

sum(spatial_object$virus_test == 0)
tibble(x = virus1_cell) %>% ggplot() + geom_density(aes(x = x)) + xlim(0,10)

sum(virus1_cell == 0)
length(virus1_cell)

which(virus1_cell ==2)
sort(virus1_cell)

1 - ppois(virus1_cell[702]-1,null_expectation1[702])






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



#Sun1

library(raster)


rois <- array_to_sf_rois("./processed_data/spatial1/A1_cell_seg_area_correction_rois.npz")
cells <- array_to_sf_rois("./processed_data/spatial1/A1_cell_seg_cells_rois.npz")

cell_area <- cells %>% st_area()
area <- extract_corrected_area(cells,rois)

cor_area <- cell_area - area
max(cor_area)
min(cor_area)
test <- st_sf(cor_area, geometry = cells)

test %>% ggplot() + geom_sf(aes(fill = cor_area > 0))






rois <- st_combine(rois)
test1 <- st_sf(name = names(cells), cells)
test <- st_difference(test1,rois)
test %>% st_area -> areatest
names(areatest) <- test %>% st_drop_geometry() %>% .[1,1]
cells %>% st_area -> cells_area
names(cells_area) <- names(cells)
areatest <- areatest[names(cells_area)]
areatest[is.na(areatest)] <- 0
corrected_area <- cells_area - areatest

return(corrected_area)


polyT <- raster::raster("./raw_data/spatial1/Panorama_Spatial1_serie2_W0A1_Cy5-class_R11_.tiff")
dapi <- rast("./raw_data/spatial1/Spatial1_serie2_serie2_A1-1_DAPI.tiff")
sum(values(dapi) == 0)
sum(values(polyT) == 0)
length(values(polyT))
length(values(mask))
mask <- polyT != 0

polygontest <- rasterToPolygons(mask, dissolve = TRUE)

plot(rois)
roi <- st_as_sf(polygontest)

test1 <- st_sf(names(cells[1:20]), cells[1:20])
test <- st_difference(test1,rois)
test %>% st_geometry() %>% rownames() 

test2 <- 1:2
names(test2) <- c("a","b")
test2["a"]
test %>% st_area

rownames(test)

plot(test)
plot(cells)
test[1:5]

colnames(test)
test <- terra::extract(img, vect(rois), fun=mean, na.rm=FALSE)

colnames(test) <- c("ID", "meanP")

rois %>% st_area() %>% length()
pix <- st_sf(test, geometry = rois)
pix$bol <- pix$meanP > 200
pix %>% ggplot() + geom_density(aes(x = meanP))
pix %>% ggplot() + geom_sf(aes(fill = bol))
ggsave("./plots/spatialcheck.png", width = 6000, height = 7500, units = "px")

# Example named vectors
vector1 <- c(a = 10, b = 20, c = 30, d = 40)
vector2 <- c(b = 5, c = 15, a = 25)

# Ensure both vectors have the same names
all_names <- union(names(vector1), names(vector2))
vector1 <- vector1[all_names]
vector2 <- vector2[names(vector1)]

# Set missing elements to 0
vector2[is.na(vector2)] <- 0

# Sum the vectors element-wise using matching names
result <- vector1 + vector2
print(result)

test <- c("a", "b", "c")
test[c(1,2,1,3,1,2,3,1,2,3)]





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



regions_to_sf_rois <- function(regions_folder_path){
    #
    #
    #
    #
    #
    
    # Get a list of region file names in the regions folder
    region_files <- list.files(path = regions_folder_path)

    # Make empty list and vector to place sf POLYGON and corresponding region name respectively
    region_list <- list()
    region_name <- c()

    # Loop over all files in the regions folder
    for (region_file in region_files){
      
      #Read ImageJ ROI file
      region <- read.ijroi(paste0(regions_folder_path,"/",region_file))

      #Extract coordinates for region ROI and make the coordinates closed
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

regions_folder_path <- "./raw_data/spatial1/spatial1_A1_regions"

test <- regions_to_sf_rois(regions_folder_path)

test2 <- st_overlaps(test)

cells %>% st_centroid() -> test3
test2 <- st_intersects(test3,test)

names(test2) <- names(cells)

test4 <- unlist(test2)

test %>% st_drop_geometry() %>% .[,1] -> test5

test5 <- test$region_name[test4]
names(test5) <- names(test4)

region_annotation <- function(rois_cells, rois_regions){
    centroid_cells <- rois_cells %>% st_centroid()
    intersects_sparse <- st_intersects(centroid_cells,rois_regions)
    names(intersects_sparse) <- names(rois_cells)
    cell_by_region_index <- unlist(intersects_sparse)
    cell_by_region <- rois_regions$region_name[cell_by_region_index]
    names(cell_by_region) <- names(cell_by_region_index)
    return(cell_by_region)
}

test1 <- region_annotation(cells, test)
names(test1)



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

    # TODO: Change such that the index in sparse_intersect is used to get the names from the rois object
    # names(rois)[unlist(sparse_intersect)]
    # Then table(names(unlist(sparse_intersect)), the above)
    #sparse_count_matrix <- as(table(names(unlist(sparse_intersects)),unlist(sparse_intersects)),"sparseMatrix")
    
    transcripts_rois_index <- unlist(sparse_intersects)
    cell_id <- names(rois)[transcripts_rois_index]

    sparse_count_matrix <- as(table(names(transcripts_rois_index),cell_id),"sparseMatrix")

    # Give cells a cell ID
    #colnames(sparse_count_matrix) <- paste0("Cell",colnames(sparse_count_matrix))

    return(sparse_count_matrix)
}

transcripts <- load_transcripts("./raw_data/spatial1/Spatial1_serie2_serie2_A1-1_results.txt")


test <- create_count_matrix(transcripts, cells)
rois <- cells
colnames(sparse_count_matrix)
rownames(sparse_count_matrix)



create_seurat_segmentations <- function(rois,sparse_count_matrix){
    # Input: ROIs as sfc object with POLYGON geometries, sparse count matrix with transcripts as rows and cells as columns
    #
    # Converts the sfc ROIs object to a data frame with column names that fits expected format by Seurat.
    # Filters cells that do not intersect any transcripts
    #
    # Output: Data frame with rois in Seurat format

    # Retrieve coordinates from sfc object, create cell IDs and filter cells that did not intersect any transcripts
    # TO DO: Set names according to names(rois), so mutate(L2 = names(rois)[L2]), because every coordinate has cell name
    rois_df <- rois %>% st_coordinates() %>% as.data.frame() %>% mutate(L2 = names(rois)[L2]) %>% filter(L2 %in% colnames(sparse_count_matrix))

    # Extract information and change column names according to the Seurat format
    seurat_segmentations <- rois_df[,c("L2","X","Y")]
    colnames(seurat_segmentations) <- c("cell","x","y")

    return(seurat_segmentations)
}


names(rois) %>% as.data.frame() %>% cbind(centroids_df) %>% filter(. %in% colnames(sparse_count_matrix)) -> test

data.frame(cell = names(rois), centroids_df) %>% filter(cell %in% colnames(sparse_count_matrix)) -> test
colnames(test)
class(test$X)
ncol(sparse_count_matrix)



create_seurat_centroids <- function(rois,sparse_count_matrix){
    # Input: ROIs as sfc object with POLYGON geometries, sparse count matrix with transcripts as rows and cells as columns
    #
    # Extracts the centroids from the cell ROIs as data frame with format that fits expected by Seurat.
    #
    # Output: Data frame with centroids of cells in Seurat format

    # Extract centroids from ROIs
    centroids_df <- rois %>% st_centroid() %>% st_coordinates()

    # Add cell IDs and filter cells that did not intersect any transcripts
    # TO DO: Change to that I cbind(names(rois))
    seurat_centroids <- data.frame(cell = names(rois), centroids_df) %>% filter(cell %in% colnames(sparse_count_matrix))

    # Change column names according to the Seurat format
    colnames(seurat_centroids) <- c("cell","x","y")

    return(seurat_centroids)
}


test <- create_seurat_centroids(rois, sparse_count_matrix)


if (logical(0)){
  print(1)
} else {
  print(0)
}
class(logical(0))

length(logical(0))
a <- logical(0)
length(a)
!(is.logical(a) & length(a) == 0)

is.null(NULL)






test <- st_sf(name = names(cells), geometry = cells)
data.frame(name = names(test1), region = test1) -> test2

left_join(test, test2, by = "name") -> test3

test3 %>% ggplot() + geom_sf(aes(fill = region))
ggsave("./plots/spatialcheck.png", width = 6000, height = 7500, units = "px")

plot(rois_cells)
plot(rois_regions)
# join the two using tidyverse funciton mentioned in spatial stat book or others. left join? idk check later