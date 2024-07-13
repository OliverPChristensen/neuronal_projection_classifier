library(Seurat)
library(qs)
library(tidyverse)
library(sf)
library(RImageJROI)
library(progress)
library(terra)
library(reticulate)
library(scales)
library(Matrix)
library(spatstat)


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



cells %>% st_centroid()


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


### TEST ###

transcript_file_match <- "results"
cell_roi_file_match <- "cells"
image_roi_file_match <- "area"
region_folder_match <- "regions"
#section_indices <- c("A1","B1","C1","D1","A2","B2","C2","D2")
section_index <- c("A1")
run_index <- "spatial1"

raw_folder_path <- paste0("./raw_data/",run_index)
processed_folder_path <- paste0("./processed_data/",run_index)
############



# Your vectors of samples and expectations
samples <- c(10, 20, 15, 25)
expectations <- c(12, 18, 17, 23)

# Convert expectations to counts by multiplying by some constant
# For example, if your expectations are rates per some unit of time, you might multiply by the total observation time
# If your expectations are percents, you might multiply by the total sample size
# Adjust this constant according to your specific scenario
constant <- 1000  # Example constant
expected_counts <- expectations * constant

# Perform the Poisson test
test_result <- poisson.test(samples, T = expected_counts)

# Print the result
print(test_result)


sparse_count_matrix <- spatial_object$raw['counts.1']

spatial_obj <- SeuratObject::CreateSeuratObject(counts = sparse_count_matrix, assay = 'raw')




#### TEST ####

transcript_file_match <- "results"
cell_roi_file_match <- "cell_seg_cells"
image_roi_file_match <- "area"
sun1_file_match <- "GFP__cells"
#section_indices <- c("A1","B1","C1","D1","A2","B2","C2","D2")
section_index <- c("A1")
run_index <- "spatial1"

raw_folder_path <- paste0("./raw_data/",run_index)
processed_folder_path <- paste0("./processed_data/",run_index)

spatial_object_file <- "spatial1_seurat_raw.qs"

qsave(spatial_object_sf, "./processed_data/spatial2/spatial_sf_test.qs")
spatial_object_sf <- qread("./processed_data/spatial1/spatial_sf_test.qs")
##############


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

convert_cell_id_seurat <- function(rois, seurat_object, section_index){
        seurat_object@meta.data %>% filter(section_index == section_index) %>% rownames() -> seurat_cell_ids
        
        seurat_section_integer_index <- substr(seurat_cell_ids, nchar(seurat_cell_ids), nchar(seurat_cell_ids)) %>% unique()
        new_cell_ids <- names(rois) %>% paste0(.,"_",seurat_section_integer_index)

        return(new_cell_ids)
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
    rois_sfc <- rois %>% st_sfc()

    rois_valid <- rois_sfc[st_is_valid(rois_sfc)]
    # Convert the list of polygons to sfc object and return the object
    return(rois_valid)
}
rois_path <- "./processed_data/spatial1/A2_cell_seg_cells_rois.npz"

cells <- array_to_sf_rois("./processed_data/spatial1/A1_cell_seg_cells_rois.npz")






spatial_object <- qread("./processed_data/spatial1/spatial1_seurat_raw.qs")


head(spatial_object)


# Visualize QC metrics as a violin plot
VlnPlot(spatial_object, features = c("nCount_raw", "area"), ncol = 2, alpha = 0)

spatial_object@meta.data[spatial_object$area > 15000,]

spatial_object@meta.data$area %>% quantile(., c(0.01,0.99)) -> aq

spatial_object$aq <- spatial_object$area < aq

spatial_object$test <- spatial_object$nCount_raw < 10

DefaultBoundary(spatial_object[["A2"]]) <- "segmentation"
ImageDimPlot(spatial_object, 
             fov = "A2",
             axes = F, 
             dark.background = T, 
             size = 1.5,
             border.color = "black", 
             border.size = 0.1, 
             cols = "polychrome",
             group.by = "test",
             #cells = WhichCells(spatial_object, expression = prediction.score.super.major.cell.type > 0.9),
            #molecules = c("Agrp","Glp1r")
)

ggsave("./plots/spatialcheck.png", width = 6000, height = 7500, units = "px")

table(spatial_object$nCount_raw < 10, spatial_object$area < aq)



spatial_object <- qread("./processed_data/spatial1/C1_spatial_object_projection_call.qs")


sun1 <- spatial_object@meta.data %>% select(sun1) %>% as.matrix()
virus1 <- spatial_object@meta.data %>% select(virus1) %>% as.matrix()
virus1_corrected <- virus1 < 0.05/(length(virus1)*2)
virus2 <- spatial_object@meta.data %>% select(virus2) %>% as.matrix()
virus2_corrected <- virus2 < 0.05/(length(virus2)*2)
virus_corrected <- virus1_corrected | virus2_corrected

table(sun1,virus_corrected)

# Calculate expected counts
chisq.test(table(sun1,virus_corrected))$expected
fisher.test(table(sun1,virus_corrected))



library('tidyverse')
library('Seurat')
library('qs')

neurons <- qread('./raw_data/7plex_neurons_filtered_labelled.qs')

head(neurons)

sum(neurons$infected == FALSE)
sum(neurons$infected == TRUE)

neurons$v147

#CR: cell ranger
#TS: Tap seq
#Lane: pos or neg

unique(neurons$sort)

table(neurons$sort)
neurons@meta.data %>% ggplot() + geom_density(aes(x = as.numeric(sort)))

neurons$lane[10000 + 1:100]
neurons$sort[10000 + 1:100]


ord <- order(unique(neurons@meta.data[,c("area","sort")])[,2])

unique(neurons@meta.data[,c("area","sort")])[ord,]


t <- as_tibble(table(neurons@meta.data[,c("area","sort")]))

t %>% filter(n > 0) %>% arrange(sort,n) -> t

print(t,n = 100)

#CEA High and low prob
#BST High and low prob
#Måske PVH, men ikke så mange gode negativer

neurons@meta.data %>% filter(area == 'CEA', sort == 1) %>% ggplot() + geom_histogram(aes(x = v147))
neurons@meta.data %>% filter(area == 'PVH', sort == 1) %>% ggplot() + geom_histogram(aes(x = v147))






spatial_object1 <- qread("./processed_data/spatial2/C1_spatial_object_raw.qs")
spatial_object2 <- qread("./processed_data/spatial1/A1_spatial_object_raw.qs")
spatial_object3 <- qread("./processed_data/spatial2/A1_spatial_object_raw.qs")


sun1 <- spatial_object@meta.data %>% select(sun1) %>% as.matrix()
virus1 <- spatial_object@meta.data %>% select(virus1) %>% as.matrix()
virus1_corrected <- virus1 < 0.05/(length(virus1)*2)
virus2 <- spatial_object@meta.data %>% select(virus2) %>% as.matrix()
virus2_corrected <- virus2 < 0.05/(length(virus2)*2)
virus_corrected <- virus1_corrected | virus2_corrected

table(sun1,virus_corrected)

# Calculate expected counts
chisq.test(table(sun1,virus_corrected))$expected
fisher.test(table(sun1,virus_corrected))

test <- list("C1" = spatial_object1, "A1" = spatial_object2, "B1" = spatial_object3)
spatial_object <- merge(x = test[[1]], y = test[-1], add.cell.ids = names(test))

colnames(spatial_object)

spatial_object@meta.data %>% ggplot() + geom_density(aes(x = area, col = section_index)) + facet_grid(1 ~ run_index)
spatial_object@meta.data %>% ggplot() + geom_density(aes(x = area, col = run_index))

spatial_object@meta.data %>% ggplot() + geom_density(aes(x = nCount_raw, col = section_index)) + facet_grid(1 ~ run_index)
spatial_object@meta.data %>% ggplot() + geom_density(aes(x = nCount_raw, col = run_index))

ntrans <- spatial_object2$section_count %>% unique()
total_area <- spatial_object2$section_area %>% unique()

ntrans/total_area


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

transcripts <- read.table("raw_data/spatial1/Spatial1_serie2_serie2_A1-1_results.txt", header = F, sep = "\t") %>% filter(!str_detect(V4, "FP "))

nrow(transcripts)
nrow(transcripts_sf)


merge_spatial_objects <- function(spatial_paths_list){
    
    run_indices <- c("spatial1", "spatial2", "spatial3")

    merged_spatial_objects <- list()

    for (run_index in run_indices){
        
        spatial_paths_list_run <- grep(run_index, spatial_paths_list, value = TRUE)
        print(spatial_paths_list_run)
        spatial_objects <- list()

        t <- 0

        for (spatial_object_path in spatial_paths_list_run){
            cat(paste0("Loading spatial object ",t,"/",length(spatial_paths_list_run)," from ",spatial_object_path,"\n"))
            # Load spatial seurat object to use for projection calls
            t <- t + 1
            spatial_objects[[t]] <- qread(spatial_object_path)
        }

        cat("Merging spatial seurat objects\n")
        merged_spatial_objects[[run_index]] <- merge(x = spatial_objects[[1]], y = spatial_objects[-1])
    }

    return(merged_spatial_objects)
}



image_rois_files <- list.files(path = ) %>% grep(image_roi_file_match, ., value = TRUE)



# Define your placeholder
placeholder <- c("A", "B", "C", "D", "E")

run_indices <- c("spatial1", "spatial2", "spatial3")
section_indices <- c("A1", "B1", "C1", "D1", "A2", "B2", "C2", "D2")
combinations <- expand.grid(run_indices, section_indices)
# Create vector of strings
vector_of_strings <- paste0("processed_data/",combinations$Var1,"/",combinations$Var2,"_spatial_object_projection_call.qs")

test <- merge_spatial_objects(vector_of_strings)



test <- spatial_object@meta.data %>% group_by(section_index, run_index) %>% summarise(sum_area = unique(section_area), sum_count = unique(section_count))



test <- list("a" = 2, "b" = 4)

for (t in names(test)){
  print(t)
}



run_indices <- c("spatial1", "spatial2", "spatial3")
section_indices <- c("A1", "B1", "C1", "D1", "A2", "B2", "C2", "D2")
combinations <- expand.grid(run_indices, section_indices)
# Create vector of strings
paths <- paste0("processed_data/",combinations$Var1,"/",combinations$Var2,"_spatial_object_projection_call.qs")


spatial_paths_list_run <- grep("A1", paths, value = TRUE) %>% grep("spatial1", ., value = TRUE)




spatial_object_list <- merge_spatial_objects(paths)

check_section_run_effect(spatial_object_list)


spatial_objects_meta_data <- data.frame()

for (spatial_object in spatial_object_list){
  print(spatial_object)
  spatial_objects_meta_data <- rbind(spatial_objects_meta_data,spatial_object@meta.data)
}

spatial_object@meta.data



p <- spatial_object1@meta.data %>% ggplot() + geom_density(aes(x = nCount_raw)) %>% ggplot_build(.)

density_values <- ggplot_build(spatial_object1@meta.data %>% ggplot() + geom_density(aes(x = nCount_raw)))$data

y <- density_values[[1]]$y
x <- density_values[[1]]$x


x[which.max(y)]


plot(x,y)

test <- density(smoothed_data, bw = 0.5,na.rm = TRUE)

data.frame(x = test$x, y = test$y) %>% ggplot() + geom_line(aes(x = x, y = y)) + xlim(0,100)



RNA.object[["integrated.adt"]] <- CreateAssayObject(data = adt.data )

print(1e6/3)

log1p()


test <- t(t(spatial_object1[['raw']]$counts) / spatial_object1$area)*1e4

test2 <- log1p(t(t(spatial_object1[['raw']]$counts) / spatial_object1$area)*1e4)

log1p(test[1,20])
test2[1,20]

test == test2


nrow(spatial_object1@meta.data)

spatial_object1 <- NormalizeData(object = spatial_object1)

spatial_object1[['raw']]$norm_data2 <- log1p(t(t(spatial_object1[['raw']]$counts) / spatial_object1$area)*1e4)
nrow(spatial_object1@meta.data)

Layers(spatial_object)

test <- ScaleData(spatial_object1, features = rownames(spatial_object1))
ncol(test[['raw']]$counts.2)
test <- RunPCA(test, features = rownames(test))



obj <- spatial_object

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)


ncol(obj[['raw']]$data.1) + ncol(obj[['raw']]$data.2) + ncol(obj[['raw']]$data.3)
ncol(obj[['raw']]$scale.data)

ncol(test[['raw']]$norm_data1)

spatial_object@meta.data %>% filter(section_index == 'C1', run_index == 'spatial2') %>% select(area) %>% as.matrix() %>% as.numeric() -> areas1
spatial_object@meta.data %>% filter(section_index == 'A1', run_index == 'spatial1') %>% select(area) %>% as.matrix() %>% as.numeric() -> areas2
spatial_object@meta.data %>% filter(section_index == 'A1', run_index == 'spatial2') %>% select(area) %>% as.matrix() %>% as.numeric() -> areas3
spatial_object[['raw']]$norm_data1 <- log1p(t(t(spatial_object[['raw']]$counts.1) / areas1)*1e4)
spatial_object[['raw']]$norm_data2 <- log1p(t(t(spatial_object[['raw']]$counts.2) / areas2)*1e4)
spatial_object[['raw']]$norm_data3 <- log1p(t(t(spatial_object[['raw']]$counts.3) / areas3)*1e4)

test <- ScaleData(spatial_object, features = rownames(spatial_object))

class(areas1)
class(spatial_object$area)


test <- Layers(spatial_object)[1]

spatial_object[['raw']][test]
spatial_object[['raw']]$counts.1

test <- list(a = 1, b = 2)




spatial_object1 <- qread("./processed_data/spatial2/C1_spatial_object_raw.qs")
spatial_object2 <- qread("./processed_data/spatial1/A1_spatial_object_raw.qs")
spatial_object3 <- qread("./processed_data/spatial2/A1_spatial_object_raw.qs")

spatial_object1 <- area_normalize_data(spatial_object1,1e4)
spatial_object2 <- area_normalize_data(spatial_object2,1e4)
spatial_object3 <- area_normalize_data(spatial_object3,1e4)

test <- list("C1" = spatial_object1, "A1" = spatial_object2, "B1" = spatial_object3)

spatial_object <- merge(x = test[[1]], y = test[-1])

spatial_object <- ScaleData(spatial_object, features = rownames(spatial_object))

count_matrix <- spatial_object[['raw']]['scale.data']

pca <- prcomp(t(count_matrix))

summary(pca)
pca$x[,1]


spatial_object1 <- qread("processed_data/spatial1/A1_spatial_object_projection_call.qs")



library(Matrix)
# Creating a sparse matrix
sparse_matrix <- Matrix(c(1, 0, 0, 2, 0, 3, 0, 0, 4), nrow = 3, sparse = TRUE)
print(sparse_matrix)
class(sparse_matrix)
t(sparse_matrix)


    ### Plots of count differences between sections and runs (pga permeabilitetsforskelle, eller forskelle i størrelse af celler?)
    ### Plot area forskelle for at se om det evt forskelle i counts kan skyldes forskelle i area. Evt udregne en mean count pr areal for hver section?
    ### Hvis ikke area og count forskelle saa er area norm klart det bedste. Hvis forskelle i disse saa count normalisering og saa haabe at probe set er represntativt og ikke biased. Vigtigst ift projecting and non projecting neurons
    ### SC tranform vs ikke sc transform? Kommer an paa om den valgte normalisering er effektiv eller om der stadig er effekter der ikke er normallisieret ud

    # Normalisering



current_time <- function(){
  return(paste0("[",format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),"]"))
}

current_time()


### TEST ###

    run_indices <- c("spatial1", "spatial2", "spatial3")
    section_indices <- c("A1", "B1", "C1", "D1", "A2", "B2", "C2", "D2")
    combinations <- expand.grid(run_indices, section_indices)
    # Create vector of strings
    spatial_object_paths <- paste0("processed_data/",combinations$Var1,"/",combinations$Var2,"_spatial_object_projection_call.qs")

    path_to_output_folder <- "processed_data"

    ############

quality_control_pr_spatial_object <- function(spatial_object_list){

    for (spatial_object_name in names(spatial_object_list)){
        spatial_object <- spatial_object_list[[spatial_object_name]]
        area_quantiles <- quantile(spatial_object$area, c(0.01,0.99))
        spatial_object <- subset(spatial_object, subset = nCount_raw > 10 & area > area_quantiles[1] & area < area_quantiles[2])

        spatial_object_list[[spatial_object_name]] <- spatial_object
    }

    return(spatial_object_list)
}



spatial_object1 <- qread("processed_data/spatial1/A1_spatial_object_projection_call.qs")

counts <- spatial_object1[['raw']]$counts

rownames(counts)
colnames(counts)

counts_non_virus <- counts[!rownames(counts) %in% c("p147-mCherry-2A-iCre","p36-NLS-Cre"),]


spatial_object1[['raw']]$non_virus <- counts_non_virus

rownames(spatial_object1[['raw']]$non_virus)

rownames(counts_non_virus)
colnames(counts_non_virus)


rowSums(counts_non_virus)
colSums(counts_non_virus)

colSums(counts_non_virus > 0)



t <- tibble(x = 1:100, y = rnorm(100),z = rep(1:10,10))

p <- t %>% ggplot() + geom_point(aes(x = x, y = y, color = as.factor(z)))

p2 <- ggMarginal(p, type="density")

ggsave(paste0("plots","/","test.png"),p2)



# Install and load the necessary packages

install.packages("cowplot")

library(ggplot2)
library(cowplot)

# Create some example data with groups
set.seed(123)
df <- data.frame(
  x = rnorm(300),
  y = rnorm(300),
  group = factor(rep(1:3, each = 100))
)

# Create the base scatter plot with a white background
scatter_plot <- ggplot(df, aes(x = x, y = y, color = group)) +
  geom_point() +
  theme_classic()


# Create the density plots with a white background
x_density <- ggplot(df, aes(x = x, fill = group)) +
  geom_density(alpha = 0.5) +
  theme_classic()


y_density <- ggplot(df, aes(x = y, fill = group)) +
  geom_density(alpha = 0.5) +
  coord_flip() +
  theme_classic()

# Create an empty plot with a white background for the top right corner
empty_plot <- ggplot() + 
  theme_classic()

# Align the plots
aligned_plots <- align_plots(x_density, scatter_plot, y_density, align = "hv", axis = "tblr")

# Arrange the plots into a single plot
final_plot <- plot_grid(
  aligned_plots[[1]], empty_plot,
  aligned_plots[[2]], aligned_plots[[3]],
  ncol = 2, nrow = 2,
  rel_widths = c(3, 1),
  rel_heights = c(1, 3)
)

theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank())



# Plot the specified principle components
plot <- ggplot(df, aes(x = x, y = y, color = group)) +
    geom_point(size = 1) +
    theme_classic()

# Create the density plots for x axis of PCA plot
x_density <- ggplot(df, aes(x = x, color = group)) +
    geom_density() +
    theme_classic() +
    theme(legend.position = "none") +
    theme(axis.title.x = element_blank())

# Create the density plots for y axis of PCA plot
y_density <- ggplot(df, aes(x = y, color = group)) +
    geom_density() +
    coord_flip() +
    theme_classic() +
    theme(legend.position = "none") + 
    theme(axis.title.y = element_blank())



# Create an empty plot with a white background for the top right corner
empty_plot <- ggplot() + 
    theme_classic()

# Align the plots
aligned_plots <- align_plots(x_density, plot, y_density, align = "hv", axis = "tblr")

# Arrange the plots into a single plot
final_plot <- plot_grid(
    aligned_plots[[1]], empty_plot,
    aligned_plots[[2]], aligned_plots[[3]],
    ncol = 2, nrow = 2,
    rel_widths = c(3, 1),
    rel_heights = c(1, 3)
    )


ggsave(paste0("plots","/","test.png"),final_plot)



spatial_object1 <- qread("processed_data/spatial_object_preprocessed.qs")

spatial_object <- qread("processed_data/spatial1/A1_spatial_object_projection_call.qs")

spatial_object@meta.data %>% ggplot() + geom_density(aes(x = nCount_raw),fill = "lightblue") + xlim(0,300)

spatial_object1$projection <- case_when(
  spatial_object1$sun1 == TRUE ~ TRUE,
  spatial_object1$virus1 < 0.05/length(spatial_object1$virus1) ~ TRUE,
  spatial_object1$virus2 < 0.05/length(spatial_object1$virus2) ~ TRUE,
  TRUE ~ FALSE
)

sum(spatial_object1$projection)
length(spatial_object1$projection)

projection_plot <- DimPlot(spatial_object1, group.by = 'projection', reduction = 'umap_features',raster = FALSE, order = TRUE)

ggsave("plots/projection_umap_plot.png", projection_plot, width = 4000, height = 3000, units = "px")

spatial_object1$region_lower <- spatial_object1$region %>% tolower()

spatial_object1$region_lower_corrected <- case_when(
  spatial_object1$region_lower == 'dorsal_hyothalamus' ~ 'dorsal_hypothalamus',
  spatial_object1$region_lower == 'anterior_hyopthalamus_left' ~ 'anterior_hypothalamus_left',
  spatial_object1$region_lower == 'interpeduncular' ~ 'interpreduncular',
  spatial_object1$region_lower == 'lateral_hypothalamus_right-1' ~ 'lateral_hypothalamus_right',
  spatial_object1$region_lower == 'pontine_grey' ~ 'pontine_gray',
  spatial_object1$region_lower == 'retrochiasmaric_left' ~ 'retrochiasmatic_left',
  spatial_object1$region_lower == 'suprachiamsatic' ~ 'suprachiasmatic',
  TRUE ~ spatial_object1$region_lower
)

spatial_object1$region_reduced <- spatial_object1$region_lower_corrected %>% gsub("_right|_left", "", .)

region_plot <- DimPlot(spatial_object1, group.by = 'region_reduced', reduction = 'umap_features', raster = FALSE, shuffle = TRUE) + 
  labs(title="",
       x ="UMAP 1", 
       y = "UMAP 2"
  )


ggsave("plots/region_umap_plot.png", region_plot, width = 5000, height = 3000, units = "px")


DimPlot(spatial_object1, group.by = 'region_reduced', reduction = 'umap_features', raster = FALSE)

library(spdep)

"/processed_data/spatial1/spatial_object_preprocessed.qs"
"./../processed_data/"

"/processed_data/"
"/processed_data/"


spatial_object1$region_lower_corrected %>% unique() %>% sort()

spatial_object1$region_reduced %>% unique() %>% sort()




spatial_object <- qread("processed_data/spatial_object_preprocessed.qs")



correct_region_names <- function(spatial_object){
    region_lower <- spatial_object$region %>% tolower()

    region_lower_corrected <- case_when(
        region_lower == 'dorsal_hyothalamus' ~ 'dorsal_hypothalamus',
        region_lower == 'anterior_hyopthalamus_left' ~ 'anterior_hypothalamus_left',
        region_lower == 'interpeduncular' ~ 'interpreduncular',
        region_lower == 'lateral_hypothalamus_right-1' ~ 'lateral_hypothalamus_right',
        region_lower == 'pontine_grey' ~ 'pontine_gray',
        region_lower == 'retrochiasmaric_left' ~ 'retrochiasmatic_left',
        region_lower == 'suprachiamsatic' ~ 'suprachiasmatic',
        TRUE ~ region_lower
    )

    spatial_object$region_corrected <- region_lower_corrected %>% gsub("_right|_left", "", .)

    return(spatial_object)
}

spatial_object <- correct_region_names(spatial_object)


head(spatial_object)


spatial_object <- add_projection_termination(spatial_object)
head(spatial_object)
spatial_object$projection_termination %>% unique()

spatial_object <- add_projection_termination(spatial_object)
projection_umaps(spatial_object)



levels <- unique(spatial_object$projection_termination)


num_levels <- spatial_object$projection_termination_no_double %>% unique() %>% length()

# Generate colors
default_colors <- hue_pal()(num_levels)
names(default_colors) <- spatial_object$projection_termination_no_double %>% unique()

default_colors["No-virus"] <- "#7F7F7F"

DimPlot(spatial_object, group.by = 'projection_termination_no_double', reduction = 'umap_features', raster = FALSE, order = TRUE) + scale_color_manual(values = default_colors)


length(grep("\\+", spatial_object$projection_termination))
length(spatial_object$projection_termination)
sum(!spatial_object$projection_termination == "No-virus")

spatial_object$projection_termination_no_double %>% unique()

area_umap(spatial_object)

FeaturePlot(spatial_object, feature = "area", reduction = "umap_features", raster = FALSE)

spatial_object$projection_termination

spatial_object$sun1_vs_virus <- case_when(
        spatial_object$sun1 == TRUE & spatial_object$virus1 > 0.05 & spatial_object$virus2 > 0.05 ~ "Sun1",
        spatial_object$sun1 == FALSE & (spatial_object$virus1 < 0.05 | spatial_object$virus2 < 0.05) ~ "Virus",
        spatial_object$sun1 == TRUE & (spatial_object$virus1 < 0.05 | spatial_object$virus2 < 0.05) ~ "Sun1+Virus",
        spatial_object$sun1 == FALSE & spatial_object$virus1 > 0.05 & spatial_object$virus2 > 0.05 ~ "Non-projecting"
    )

spatial_object@meta.data %>% 
  group_by(run_index, section_index, sun1_vs_virus) %>% 
  summarize(count = n()) %>% mutate(total_count = sum(count),relative_count = count / total_count) %>%
  ggplot() + geom_col(aes(x = sun1_vs_virus, y = relative_count)) + facet_grid(section_index~run_index) + 
  geom_text(aes(x = sun1_vs_virus, y = relative_count + 0.1, label = count))


spatial_object@meta.data %>% 
  group_by(run_index, section_index, sun1_vs_virus) %>% 
  summarize(count = n()) %>% mutate(total_count = sum(count),relative_count = count / total_count) -> test


write.table(as.data.frame(test), file = "./processed_data/sun1_vs_virus.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

threshold <- 0.05/ncol(spatial_object)

spatial_object$sun1_vs_virus_bonferroni <- case_when(
    spatial_object$sun1 == TRUE & spatial_object$virus1 > threshold & spatial_object$virus2 > threshold ~ "Sun1",
    spatial_object$sun1 == FALSE & (spatial_object$virus1 < threshold | spatial_object$virus2 < threshold) ~ "Virus",
    spatial_object$sun1 == TRUE & (spatial_object$virus1 < threshold | spatial_object$virus2 < threshold) ~ "Sun1+Virus",
    spatial_object$sun1 == FALSE & spatial_object$virus1 > threshold & spatial_object$virus2 > threshold ~ "Non-projecting"
)

sun1_vs_virus_bonferroni_bar_plot <- spatial_object@meta.data %>% 
        group_by(run_index, section_index, sun1_vs_virus_bonferroni) %>% 
        summarize(count = n()) %>% 
        mutate(total_count = sum(count),relative_count = count / total_count) %>%
        ggplot() + geom_col(aes(x = sun1_vs_virus_bonferroni, y = relative_count)) + facet_grid(section_index~run_index) + 
        geom_text(aes(x = sun1_vs_virus_bonferroni, y = relative_count + 0.1, label = count))
head(spatial_object)


unique(spatial_object$sun1_vs_virus == spatial_object$sun1_vs_virus_bonferroni)


sun1_vs_virus_bonferroni_bar_plot <- spatial_object@meta.data %>% 
        group_by(run_index, section_index, sun1_vs_virus_bonferroni) %>% 
        summarize(count = n()) %>% mutate(p_value = ) 


spatial_object@meta.data %>% 
        group_by(run_index, section_index) %>% select


spatial_object$sun1
spatial_object$virus1 < 0.05 | spatial_object$virus2 < 0.05



spatial_object@meta.data 

fisher.test(table(
  spatial_object$sun1,
  spatial_object$virus1 < 0.05 | spatial_object$virus2 < 0.05
))

sun1_vs_virus_bonferroni_bar_plot <- spatial_object@meta.data %>% 
        group_by(run_index, section_index, sun1_vs_virus_bonferroni) %>% 
        summarize(count = n()) -> test2
test <- rep("a",3*8)

p_values <- data.frame(
  run_index = test2$run_index,
  section_index = test2$section_index,
  label = test
)

sun1_vs_virus_bonferroni_bar_plot <- spatial_object@meta.data %>% 
        group_by(run_index, section_index, sun1_vs_virus_bonferroni) %>% 
        summarize(count = n()) %>% mutate(total_count = sum(count),relative_count = count / total_count) %>%
        ggplot() + geom_col(aes(x = sun1_vs_virus_bonferroni, y = relative_count)) + facet_grid(section_index~run_index) + 
        geom_text(aes(x = sun1_vs_virus_bonferroni, y = relative_count + 0.1, label = count)) + 
        geom_text(data = p_values, aes(x = 2.5, y = Inf, label = format(p_value, scientific = TRUE, digits = 2)), hjust = 0.5, vjust = 1.2, inherit.aes = FALSE, size = 5)


p_values <- sun1_vs_virus_counts %>% select(-c(total_count,relative_count)) %>% pivot_wider(names_from = sun1_vs_virus, values_from = count) %>%
  mutate(p_value = fisher.test(matrix(c(`Non-projecting`,Sun1,`Sun1+Virus`,Virus),nrow = 2))$p.value)


test$p_value

fisher.test(matrix(c(7634,55,194,449),nrow = 2))$p.value


projection_terminations <- spatial_object@meta.data %>% 
  group_by(run_index, section_index, projection_termination) %>%
  summarize(count = n()) %>% 
  filter(!projection_termination == "No-virus")

projection_terminations$projection_termination <- factor(projection_terminations$projection_termination, levels = c("BST","PAG", "BST+PAG", "LHA", "PBN", "LHA+PBN", "PVT", "CEA", "PVT+CEA"))

projection_terminations %>%
  ggplot() + geom_col(aes(x = projection_termination, y = count)) + 
  facet_grid(section_index~run_index, scales = "free_x")

spatial_object@meta.data %>% 
  group_by(run_index, section_index, projection_termination) %>%
  summarize(count = n()) -> test

print(test, n = 100)

head(spatial_object)


spatial_object$virus1_count <- 2386033
spatial_object$virus2_count <- 2386033

rownames(spatial_object@assays$raw[paste0("counts.",section_index_number,".",run_index_number)])
colnames(spatial_object@assays$raw[paste0("counts.",section_index_number,".",run_index_number)])
test <- spatial_object@assays$raw[paste0("counts.",section_index_number,".",run_index_number)]

test[gsub("_", "-", "p147_mCherry_2A_iCre"),cells]

cells


spatial_object <- qread("processed_data/spatial_object_preprocessed.qs")


spatial_object <- correct_region_names(spatial_object)
spatial_object <- add_projection_termination(spatial_object)


table(spatial_object$region_corrected,spatial_object$projection_termination)


regions <- spatial_object$region_corrected
cell_properties <- spatial_object$projection_termination


region_unique <- unique(na.omit(regions))
cell_property_unique <- unique(na.omit(cell_properties))

enrichment_matrix <- matrix(0, nrow = length(region_unique), ncol = length(cell_property_unique))
rownames(enrichment_matrix) <- sort(region_unique)
colnames(enrichment_matrix) <- sort(cell_property_unique)

for (region in region_unique){
  for (cell_property in cell_property_unique){
    regions_test <- case_when(regions == region ~ region, TRUE ~ 'Others')
    cell_properties_test <- case_when(cell_properties == cell_property ~ cell_property, TRUE ~ 'Others')
     
      table(regions_test,cell_properties_test) %>% 
      as.data.frame() %>% 
      arrange(desc(regions_test != "Others"),desc(cell_properties_test != "Others")) %>% 
      pivot_wider(names_from = cell_properties_test, values_from = Freq) %>% 
      select(-regions_test) %>% 
      fisher.test(alternative="greater") %>% 
      .[["p.value"]] -> enrichment_matrix[region,cell_property]
    
  }
}

enrichment_matix <- fisher_region_cell_property_enrichment(spatial_object$region_corrected,spatial_object$projection_termination)


test <- table(regions, cell_properties)

class(test)

enrichment_matrix == 0


    #Dot plot of enrichment fisher test
    region_cell_property_enrichment_plot <- enrichment_matrix_inv_long %>%
        ggplot() + 
        geom_point(aes(x = projection_termination, y = regions, size = inv_p_value, color = factor(inv_p_value  > -log(0.05/nrow(enrichment_matrix_inv_long))))) +
        scale_size(range = c(1,6), name="-log(P)") + 
        scale_color_manual(values = c("TRUE" = "black", "FALSE" = "transparent")) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

sun1_virus_across_sections <- function(spatial_object){
    
    spatial_object$sun1_vs_virus <- case_when(
        spatial_object$sun1 == TRUE & spatial_object$virus1 > 0.05 & spatial_object$virus2 > 0.05 ~ "Sun1",
        spatial_object$sun1 == FALSE & (spatial_object$virus1 < 0.05 | spatial_object$virus2 < 0.05) ~ "Virus",
        spatial_object$sun1 == TRUE & (spatial_object$virus1 < 0.05 | spatial_object$virus2 < 0.05) ~ "Sun1+Virus",
        spatial_object$sun1 == FALSE & spatial_object$virus1 > 0.05 & spatial_object$virus2 > 0.05 ~ "Non-projecting"
    )


    sun1_vs_virus_counts <- spatial_object@meta.data %>% 
        group_by(run_index, section_index, sun1_vs_virus) %>% 
        summarize(count = n()) %>% mutate(total_count = sum(count),relative_count = count / total_count)
    
    p_values <- sun1_vs_virus_counts %>% 
        select(-c(total_count,relative_count)) %>% 
        pivot_wider(names_from = sun1_vs_virus, values_from = count) %>%
            mutate(p_value = fisher.test(matrix(c(`Non-projecting`,Sun1,`Sun1+Virus`,Virus),nrow = 2))$p.value)


    sun1_vs_virus_bar_plot <- sun1_vs_virus_counts %>%
        ggplot() + geom_col(aes(x = sun1_vs_virus, y = relative_count)) + facet_grid(section_index~run_index) + 
        geom_text(aes(x = sun1_vs_virus, y = relative_count + 0.1, label = count)) + 
        geom_text(data = p_values, aes(x = 2.5, y = Inf, label = format(p_value, scientific = TRUE, digits = 2)), hjust = 0.5, vjust = 1.2, inherit.aes = FALSE, size = 5)


    ggsave("plots/spatial_analysis/sun1_vs_virus_bar_plot.png", sun1_vs_virus_bar_plot)

    threshold <- 0.05/ncol(spatial_object)

    spatial_object$sun1_vs_virus_bonferroni <- case_when(
        spatial_object$sun1 == TRUE & spatial_object$virus1 > threshold & spatial_object$virus2 > threshold ~ "Sun1",
        spatial_object$sun1 == FALSE & (spatial_object$virus1 < threshold | spatial_object$virus2 < threshold) ~ "Virus",
        spatial_object$sun1 == TRUE & (spatial_object$virus1 < threshold | spatial_object$virus2 < threshold) ~ "Sun1+Virus",
        spatial_object$sun1 == FALSE & spatial_object$virus1 > threshold & spatial_object$virus2 > threshold ~ "Non-projecting"
    )

    sun1_vs_virus_counts_bonferroni <- spatial_object@meta.data %>% 
        group_by(run_index, section_index, sun1_vs_virus_bonferroni) %>% 
        summarize(count = n()) %>% mutate(total_count = sum(count),relative_count = count / total_count)
    
    p_values_bonferroni <- sun1_vs_virus_counts_bonferroni %>% 
        select(-c(total_count,relative_count)) %>% 
        pivot_wider(names_from = sun1_vs_virus_bonferroni, values_from = count) %>%
            mutate(p_value = fisher.test(matrix(c(`Non-projecting`,Sun1,`Sun1+Virus`,Virus),nrow = 2))$p.value)


    sun1_vs_virus_bonferroni_bar_plot <- sun1_vs_virus_counts_bonferroni %>%
        ggplot() + geom_col(aes(x = sun1_vs_virus_bonferroni, y = relative_count)) + facet_grid(section_index~run_index) + 
        geom_text(aes(x = sun1_vs_virus_bonferroni, y = relative_count + 0.1, label = count)) + 
        geom_text(data = p_values_bonferroni, aes(x = 2.5, y = Inf, label = format(p_value, scientific = TRUE, digits = 2)), hjust = 0.5, vjust = 1.2, inherit.aes = FALSE, size = 5)


    ggsave("plots/spatial_analysis/sun1_vs_virus_bonferroni_bar_plot.png", sun1_vs_virus_bonferroni_bar_plot)

}

projection_termination_plots <- function(spatial_object, projection_termination_column_name){

    projection_terminations <- spatial_object@meta.data %>% 
        group_by(run_index, section_index, !!sym(projection_termination_column_name)) %>%
        summarize(count = n()) %>% 
        filter(!!sym(projection_termination_column_name) != "No-virus")

    projection_terminations <- projection_terminations %>% 
        mutate(!!sym(projection_termination_column_name) := factor(
                !!sym(projection_termination_column_name), 
                levels = c("BST","PAG", "BST+PAG", "LHA", "PBN", "LHA+PBN", "PVT", "CEA", "PVT+CEA")
            ))

    projection_terminations_plot <- projection_terminations %>%
        ggplot() + geom_col(aes(x = !!sym(projection_termination_column_name), y = count)) + facet_grid(section_index~run_index, scales = "free_x") + 
        geom_text(aes(x = !!sym(projection_termination_column_name), y = count + 100, label = count))

    ggsave("plots/spatial_analysis/projection_terminations_plot.png", projection_terminations_plot)

    projection_terminations_bonferroni <- spatial_object@meta.data %>% 
        group_by(run_index, section_index, projection_termination_bonferroni) %>%
        summarize(count = n()) %>% 
        filter(!projection_termination_bonferroni == "No-virus")

    projection_terminations_bonferroni$projection_termination_bonferroni <- factor(
            projection_terminations_bonferroni$projection_termination_bonferroni, 
            levels = c("BST","PAG", "BST+PAG", "LHA", "PBN", "LHA+PBN", "PVT", "CEA", "PVT+CEA")
        )

    projection_terminations_bonferroni_plot <- projection_terminations_bonferroni %>%
        ggplot() + geom_col(aes(x = projection_termination_bonferroni, y = count)) + facet_grid(section_index~run_index, scales = "free_x") + 
        geom_text(aes(x = projection_termination_bonferroni, y = count + 100, label = count))

    ggsave("plots/spatial_analysis/projection_terminations_bonferroni_plot.png", projection_terminations_bonferroni_plot)
}


projection_termination_column_name <- "projection_termination"
projection_terminations$projection_termination


spatial_object$virus1_count <- 2386033
spatial_object$virus2_count <- 2386033





polyT <- raster::raster("./raw_data/spatial1/Panorama_Spatial1_serie2_W0A1_Cy5-class_R11_.tiff")

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

img <- raster::raster("./raw_data/spatial1/Panorama_Spatial1_serie2_W0A1_Cy5-class_R11_.tiff")
test <- terra::extract(img, vect(rois), fun=mean, na.rm=FALSE)


intensities %>% .[["intensity_background_virus1"]]

rate <- intensities %>% filter(run_index == "spatial1", section_index == "A1") %>% .[["intensity_background_virus1"]]*1000

data <- rpois(10000,rate)

data_scaled <- data[data < quantile(data, 0.95)]

mean(data_scaled)

rpois(10,1.1)



rate <- intensities %>% filter(run_index == "spatial1", section_index == "A1") %>% .[["intensity_background_virus1"]]*1000

data <- rpois(1000000,rate)

print(rate/1000)
print(mean(data)/1000)


virus_index <- "virus1"

rates <- intensities %>% .[[paste0("intensity_background_",virus_index)]]*100000000
intensities_scaled <- rep(0,length(rates))


for (i in 1:length(rates)){
  intensities_scaled[i] <- rpois(1000000,rates[i]) %>% .[. < quantile(., 0.95)] %>% mean(.)/1000  
}

intensities[[paste0("intensity_background_",virus_index)]]/intensities_scaled*intensities[[paste0("intensity_negative_cells_",virus_index)]]

print(intensities$intensity_background_virus1/intensities_scaled)

print(intensities_scaled)

intensity_scales <- rpois(10,rate) %>% mean(.)/1000


scale_intensity <- function(intensities, virus_index){
    rates <- intensities %>% .[[paste0("intensity_background_",virus_index)]]*1000
    intensities_scaled <- rep(0,length(rates))


    for (i in 1:length(rates)){
    intensities_scaled[i] <- rpois(1000000,rates[i]) %>% .[. < quantile(., 0.95)] %>% mean(.)/100000
    }

    intensities[[paste0("intensity_background_",virus_index)]]/intensities_scaled*intensities[[paste0("intensity_negative_cells_",virus_index)]]

}

spatial_object$virus <- case_when(
  spatial_object$virus1 < 0.05 & spatial_object$virus2 > 0.05 ~ "virus1",
  spatial_object$virus1 > 0.05 & spatial_object$virus2 < 0.05 ~ "virus2",
  spatial_object$virus1 > 0.05 & spatial_object$virus2 > 0.05 ~ "No-virus",
  spatial_object$virus1 < 0.05 & spatial_object$virus2 < 0.05 ~ "virus1+virus2",
)



levels <- unique(spatial_object$virus)


num_levels <- spatial_object$projection_termination_no_double %>% unique() %>% length()


levels <- unique(spatial_object$virus)
default_colors <- hue_pal()(length(levels))
names(default_colors) <- levels

default_colors["No-virus"] <- "#7F7F7F"





DefaultBoundary(spatial_object[["A1"]]) <- "segmentation"
spatial_plot <- ImageDimPlot(spatial_object, 
             fov = "A1",
             axes = F, 
             dark.background = T, 
             size = 1.5,
             border.color = "black", 
             border.size = 0.1, 
             cols = default_colors, 
             group.by = "virus",
             molecules = c("p147_mCherry_2A_iCre","p36_NLS_Cre"),
             nmols = 2000000,
             mols.size = 0.05,
             mols.alpha = 0.1
)

ggsave("./plots/spatialcheck.png",plot = spatial_plot, width = 6000, height = 7500, units = "px")



section_indices <- c("A1", "B1", "C1", "D1", "A2", "B2", "C2", "D2")

run_indices <- c("spatial1", "spatial2", "spatial3")
for (run_index_number in 1:3){
    for (section_index in section_indices){
      run_index_number_adjusted <- ifelse(run_index_number == 1, "", paste0(".",run_index_number))
      fov <- paste0(section_index,run_index_number_adjusted)

      DefaultBoundary(spatial_object[[fov]]) <- "segmentation"
      spatial_plot <- ImageDimPlot(spatial_object, 
                  fov = fov,
                  axes = F, 
                  dark.background = T, 
                  size = 1.5,
                  border.color = "black", 
                  border.size = 0.1, 
                  cols = default_colors, 
                  group.by = "virus",
                  molecules = c("p147_mCherry_2A_iCre","p36_NLS_Cre"),
                  nmols = 2000000,
                  mols.size = 0.05,
                  mols.alpha = 0.1
      )

      ggsave("./plots/spatialcheck.png",plot = spatial_plot, width = 6000, height = 7500, units = "px")
    }
}





spatial_object <- qread("processed_data/spatial_object_preprocessed.qs")

fov <- "B1"

seurat_to_sf <- function(spatial_object){
        
    #Loops that fetch the cell segmentations from the Seurat object and put them into a sf geometry
    cell_rois <- list()

    for (fov in names(spatial_object@images)){
        print(paste0("Converting fov ",fov," to sf"))
        
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

spatial_object_sf <- seurat_to_sf(spatial_object)


region_paths <- c()
for (run_index in run_indices){
    region_paths <- c(region_paths,paste0(paste0("raw_data/",run_index),"/",list.files(path = paste0("raw_data/",run_index), pattern = "regions")))
}

transcript_paths <- c()
for (run_index in run_indices){
    transcript_paths <- c(transcript_paths,paste0(paste0("raw_data/",run_index),"/",list.files(path = paste0("raw_data/",run_index), pattern = "results")))
}

run_index <- "spatial1"
section_index <- "A1"

regions_list$section_index %>% unique()

transcripts_list

rownames(transcripts)

transcripts[1:1000,] -> test1
test1 %>% select(-c("V3", "V5")) -> test2
test2 %>% st_cast("MULTIPOINT") -> test

test1 %>%
  group_by(V4) %>%
  summarize(geometry = st_combine(geometry)) %>% as.data.frame() -> test2

test2 %>% st_cast("MULTIPOINT") -> test3

test2 %>% as.data.frame() -> test3
test2$V4

test1$V3 %>% unique()

# Example of pivot wider for shiny app

transcript_region_counts_object %>% 
  filter(run_index == "spatial1",section_index == "A1", transcript == "Acvr1c") %>% 
  pivot_wider(names_from = region, values_from = counts) %>%
  select(-c("run_index", "section_index")) -> test

test %>% mutate_all(~replace(., is.na(.), 0)) -> test2


unique(st_geometry_type(cell_sf_object) == "MULTIPOINT")



radioButtons("section", "Section:",
choices = section,
selected = section[1],
inline = TRUE),




create_ggplot <- eventReactive(input$do, {
    
    filtered_data <- filter_data()
    
    projection_termination_unique <- cell_object$projection_termination %>% unique()
    
    # Generate colors
    plotting_colors <- projection_termination_unique %>% length() %>% hue_pal()(.)
    names(plotting_colors) <- projection_termination_unique
    
    plotting_colors["No-virus"] <- "#7F7F7F"
    
    
    plot <- ggplot() + 
        geom_sf(data = filtered_data$cell_object, aes(fill = cell_types, celltype = data.predicted.cell.type, count = data.nCount_Xenium), lwd = 0.1) + 
        scale_fill_manual(values = default_colors)
    
    if (!is.null(filtered_data$transcript_object)){
  
        plot <- plot + 
            geom_sf(data = filtered_data$transcript_object, aes(color = transcript), size = 0.5, shape = 15) + 
            scale_colour_brewer(palette = "Set1") + 
            guides(color = guide_legend(override.aes = list(size = 6)))
      
    }
    
    if (input$regions){
        plot <- plot + 
            geom_sf(data = filtered_data$region_object,fill = NA,color = "white", lwd = 0.1) + 
            geom_sf_text(data = filtered_data$region_object, aes(label = region), color = "white", size = 6, fontface = "bold")
      
    }
    
    plot <- plot + 
        geom_sf(data = filtered_data$cell_object_centroid, aes(text = tooltip), size = 5, alpha = 0)

    plot <- plot + 
        theme_void() + 
        labs(fill = "Selected Transferred Cell Types", color = "Selected Transcripts") + 
        theme(text = element_text(color = "black"),
            legend.margin = margin(10, 10, 10, 10),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 15),
            legend.key.size = unit(2, "lines"),
            panel.background = element_rect(fill = "black", color = NA),
            panel.border = element_rect(color = "#595959", fill = NA, linewidth = 2))

    return(plot)

  })


selectInput("transcript_menu_item", 
label = "Transcripts:",
choices = genes,
#selected = genes[1],
multiple = TRUE),



observe({
    output$plot_ui <- renderUI({
      if (test()) {
        plotlyOutput("spatial_plotly", width = "1200px", height = "800px")
        output$spatial_plot <- NULL 
        output$spatial_plotly <- renderPlotly({
          create_plotly_plot()
        })
      } else {
        plotOutput("spatial_plot", width = "1200px", height = "800px")
        output$spatial_plotly <- NULL  # Clear the interactive plot output
        output$spatial_plot <- renderPlot({
          create_ggplot()
        })
      }
    })
    
    if (test()) {
      output$spatial_plot <- NULL  # Clear the static plot output
      output$spatial_plotly <- renderPlotly({
        #print("test3")
        #print(values2$old_section)
        #print(values2$same_section)
        plot <- plotly_plot()
      })
      
    } else {
      output$spatial_plotly <- NULL  # Clear the interactive plot output
      output$spatial_plot <- renderPlot({
        plot <- plot_func()
        plot
      })
    }
  })




  plot_interactive <- eventReactive(input$render_plot,{
    return(input$interactive)
  })
  
  create_plotly_plot <- eventReactive(input$render_plot, {

    plot <- create_ggplot()
    if (!is.null(zoom_config$previous_zoom_state)){

      plotly_plot <- ggplotly(plot, tooltip = "text") %>%
        layout(
          xaxis = list(range = c(zoom_config$previous_zoom_state$x1, zoom_config$previous_zoom_state$x2)),
          yaxis = list(range = c(zoom_config$previous_zoom_state$y1, zoom_config$previous_zoom_state$y2))
        )

    } else {

      plotly_plot <- ggplotly(plot, tooltip = "text")
    }

    return(plotly_plot)
  })

reset_zoom <- eventReactive(input$render_plot,{

    if (zoom_config$previous_section$run_menu_item != input$run_menu_item & zoom_config$previous_section$section_menu_item != input$section_menu_item){
        zoom_config$previous_zoom_state <- NULL
    }
    
    zoom_config$previous_zoom_state <- list(
      run_menu_item = input$run_menu_item,
      section_menu_item = input$section_menu_item
    )

    return()
})


  if (!is.null(values2$old_section)){
      print("test")
      if (!values2$same_section){
        values3$old_zoom <- NULL
        print("set null")
      }
  }
    
  values2$old_section <- input$section

  test <- list()

test$test2 <- list(
  test3 = "test4"
)

test$test2$test3 == "test4"

TRUE & FALSE




  # Observe the selected category and update the item choices
  observeEvent(input$run_menu_item, {
    # Define the items based on the selected category
    run_menu_items <- switch(input$run_menu_item,
                    "BST and PAG" = c("A1_test1","B1_test1","C1_test1","D1_test1","A2_test1","B2_test1","C2_test1","D2_test1"),
                    "LHA and PBN" = c("A1_test2","B1_test2","C1_test2","D1_test2","A2_test2","B2_test2","C2_test2","D2_test2"),
                    "PVT and CEA" = c("A1_test3","B1_test3","C1_test3","D1_test3","A2_test3","B2_test3","C2_test3","D2_test3")
                    )

    # Update the choices in the item selectizeInput
    updateSelectizeInput(session, "section_menu_item", choices = run_menu_items, server = TRUE)
  })

shiny_data$cell_object %>% colnames()


shiny_data$transcript_object %>% select(transcript, geometry) -> test

test


shiny_data$region_object$test <- "test"

class(shiny_data$region_object)


shiny_data$cell_object_centroid %>% class()


class(transcripts_multipoint)


transcripts_multipoint <- transcripts %>%
    group_by(V4) %>%
    summarize(geometry = st_combine(geometry))



transcripts_multipoint %>% st_geometry() -> geom

transcripts_multipoint %>% st_drop_geometry() %>%
  as.data.frame() -> test

test2 <- st_sf(test, geometry = geom)


class(test2)

transcripts_multipoint %>% st_centroid()

st_coordinates(transcripts)

FALSE & TRUE
bregma_levels[[run_menu_items[1]]][1]

section_conversion[[bregma_levels[[run_menu_items[1]]][1]]]

test <- list("a" = 1, "a" = 2)
test$a

bregma_levels

names(section_conversion) %>% unique() %>% length()



shiny_data$cell_object_centroid
round(2342.4)

test <- read.table("./processed_data/ROI_Bregma_spatial.tsv", sep = "\t", header = TRUE)

#test %>% mutate(bregma_level = paste0(Bregma," (",Orientation,")")) -> test2

test2 %>% mutate(run_index = paste0("spatial",Spatial_run)) -> test2

test2$projection_termination <- case_when(
  test2$run_index == "spatial1" ~ "BST and PAG",
  test2$run_index == "spatial2" ~ "LHA and PBN",
  test2$run_index == "spatial3" ~ "PVT and CEA"
)

test3 <- test2 %>% 
  mutate(section_index = Well, bregma_level = paste0(Bregma," (",Orientation,")")) %>% 
  select(run_index, section_index, projection_termination, bregma_level)




menu_conversion <- read.table("processed_data/ROI_Bregma_spatial.tsv", sep = "\t", header = TRUE)


menu_conversion <- menu_conversion %>% 
    mutate(
        run_index = paste0("spatial",Spatial_run),
        section_index = Well, 
        bregma_level = paste0(Bregma," (",Orientation,")")
    )

menu_conversion$projection_termination <- case_when(
    menu_conversion$run_index == "spatial1" ~ "BST and PAG",
    menu_conversion$run_index == "spatial2" ~ "LHA and PBN",
    menu_conversion$run_index == "spatial3" ~ "PVT and CEA"
)

menu_conversion <- menu_conversion %>% 
    select(run_index, section_index, projection_termination, bregma_level)

menu_conversion %>% filter(Well %in% c("A1", "A2"))

format(1.23423423e-21, scientific = TRUE, digits = 2)
format(1, scientific = TRUE, digits = 2)
format(0, scientific = TRUE, digits = 2)

round(1.23423423e-21)

format_p_value <- function(x) {
  if (abs(x) < 1e-2) {
    return(format(x, scientific = TRUE, digits = 3))
  } else {
    return(round(x, 2))
  }
}

format_p_value(1.2)


FALSE | FALSE | TRUE



shiny_data$cell_object$virus1

x <- shiny_data$cell_object$virus1
ifelse(abs(x) < 1e-2, format(x, scientific = TRUE, digits = 3), round(x, 2))

test <- tibble("xs_s1" = 1:10, "ys_s1" = 1:10, "zs 2" = 1:10)

test %>%
  mutate(total = rowSums(select(.,-c("z")))) -> test

test %>% 
  mutate("A B" = x)
  select(c("A B"), total, everything())

test %>% 
  rename_with(~ str_replace_all(., "_", " "), everything()) %>% 
  rename_all(~ str_to_title(.))

shiny_data


run_conversion <- shiny_data$cell_object$run_index %>% unique() %>% sort() %>% as.list()
section_conversion <- shiny_data$cell_object$section_index %>% unique() %>% sort() %>% rep(.,3) %>% as.list()


transcript_menu_items <- shiny_data$transcript_object$transcript %>% unique() %>% sort()

#I have to make double check for the conversion between index. Becuase doubliqutes in bregma levels
run_menu_items <- c(
    "BST and PAG",
    "LHA and PBN",
    "PVT and CEA"
)

bregma_levels <- list(
    c("-0.47","-0.83","-1.23","-1.43","-1.79","-2.69","0.675","0.2"),
    c("-0.59","-0.95","-1.31","-1.55","-1.91","-2.45","0.3","0.475"),
    c("-0.71","-1.23","-1.43","-1.67","-1.91","-2.45","0.225","0.2")
)

names(bregma_levels) <- run_menu_items
names(run_conversion) <- run_menu_items
names(section_conversion) <- c(bregma_levels[[run_menu_items[1]]],bregma_levels[[run_menu_items[[2]]]],bregma_levels[[run_menu_items[[3]]]])


test <- shiny_data$transcript_region_table %>% 
  filter(run_index == "spatial1", section_index == "A1", transcript == "Agrp")



transcript_region_table <- test %>% 
    ungroup() %>% 
    pivot_wider(names_from = region, values_from = counts) %>%
    #mutate_all(~replace(., is.na(.), 0)) %>%
    select(-c("run_index", "section_index")) %>% 
    mutate(section = rowSums(select(.,-c("transcript")))) %>%
    select(transcript, section, everything()) %>% 
    rename_with(~ str_replace_all(., "_", " "), everything()) %>% 
    rename_all(~ str_to_title(tolower(.)))

is.na(transcript_region_table) %>% sum()



library(ggplot2)

# Example plot
p <- ggplot(mtcars, aes(x = mpg, y = disp, color = factor(cyl))) +
  geom_point() +
  labs(title = "Scatter Plot", x = "Miles per Gallon", y = "Displacement", color = "Cylinders")

# Build the plot object
plot_data <- ggplot_build(p)

# Access legend information
legend_info <- plot_data$layout$guides$colour$title

# Print legend title
print(legend_info)

# Example data
df <- data.frame(
  x = 1:5,
  y = 1:5,
  group = c("Group 1", "Group 2", "Group 1", "Group 2", "Group 1")
)

# Create a dummy plot with only legend
p <- ggplot(df, aes(x = x, y = y, color = group))

# Print the plot
print(p)



library(ggplot2)



# Example data
df <- data.frame(
  x = rep(1,5),
  y = rep(1,5),
  group = c("Group 1", "Group 2", "Group 1", "Group 2", "Group 1")
)

# Create a plot with the required aesthetic to generate the legend
df %>% ggplot() +
  geom_point(aes(x = x, y = y, color = group)) +  # Invisible points to generate the legend
  labs(color = "Group Type") + 
  xlim(0,0) + 
  ylim(0,0) +
  theme_void()




# Extract the legend
extract_legend <- function(plot) {
  ggplotGrob(plot)$grobs[which(sapply(ggplotGrob(plot)$grobs, function(x) x$name) == "guide-box")]
}

# Display only the legend
legend <- extract_legend(p)
grid.newpage()
grid.draw(legend)


library(ggplot2) 
library(grid)
library(gridExtra) 

my_hist <- ggplot(diamonds, aes(clarity, fill = cut)) + 
    geom_bar() 

# Using the cowplot package
legend <- cowplot::get_legend(my_hist)

shiny_data$cell_object$projection_termination


# Create an empty plot and add the legend
legend_plot <- ggplot() +
  theme_void() +
  annotation_custom(legend)

grid.newpage()
grid.draw(legend)



            #uiOutput("plot_ui"),
            #tags$div(style = "margin-top: 200px;"),
            #plotOutput("legend")
            #fluidRow(
                #column(3, uiOutput("plot_ui"),), # Adjust the width as needed
                #column(3, plotOutput("legend")) # Adjust the width as needed
            #)

            div(
                style = "display: flex; align-items: flex-start;",
                    div(
                        style = "flex: 1; margin-right: 0px;",
                        uiOutput("plot_ui")
                    ),
                    div(
                        style = "flex: 1; margin-left: 0px;",
                        plotOutput("legend")
                    )
            )


shiny_data$transcript_region_table %>% group_by(run_index, section_index, transcript) %>%
  summarize(new_counts = sum(counts)) %>% 
  .[["new_counts"]] %>% sort() 

# Downsample til 10K? Og besked om hvilke transcripter der er downsampled og hvad de originalt var
# Hvilket man ogs%>% kan se i summary table, men fint at have det st%>%ende flere steder

shiny_data$transcript_object


test <- tibble(x = 1:1000, y = rep(c("A","B"),500))

test %>%
  group_by(y) %>%
  slice_sample(n = 1000) -> test2

test2$y


# Realtive vs abs count? i tabel
# COunt vs norm count? i hover
# barplot? eller er tabel nok?
# Is lasso interesting if you have the other features?

#(V2) Instead of hover then click and get information on the side
# Dylan wants barplots
#(V2) Relative option.
#(V2) Lasso even with downsamled transcripts because you can then multipy by the down sampling factor that could be stored
# For dvc dylan wants outline of the labelled cells. I think the option to plot individual cell types will help
#Escpeccially because it is only hihly expressed genes that will be downsampled so the multiplication will be a fine approximation



spatial_object <- qread("processed_data/spatial_object_preprocessed.qs")
count_layer <- "counts.1.1"
test <- spatial_object[["raw"]][count_layer]
object.size(test)

dense_matrix <- as.data.frame(t(test))
dense_matrix$cellID <- rownames(dense_matrix)

object.size(dense_matrix)

test2 <- spatial_object@meta.data %>% filter(section_index == "A1", run_index == "spatial1") 
test2$cellID <- rownames(test2)

colnames(test)

left_join(test2, dense_matrix, by = "cellID")

cbind(test2,dense_matrix)

count_matrix <- count_matrices[1]

spatial_object[["raw"]][count_layer]  %>% rownmes()




transcripts_downsampled %>%
  group_by(V4) %>% 
  summarize(counts = n())

shiny_data$menu_converter


count_layer <- count_layers[1]


rotate_geometry <- function(geometry, angle_degrees) {
  # Convert angle to radians
  angle_radians <- angle_degrees * pi / 180
  
  # Create a rotation matrix
  rotation_matrix <- matrix(
    c(cos(angle_radians), -sin(angle_radians),
      sin(angle_radians),  cos(angle_radians)),
    ncol = 2, byrow = TRUE
  )
  
  # Apply the rotation matrix to each geometry
  rotated_geometry <- lapply(geometry, function(geom) {
    coords <- st_coordinates(geom)
    coords_rotated <- coords
    coords_rotated[, 1:2] <- t(rotation_matrix %*% t(coords[, 1:2]))
    
    if (inherits(geom, "POINT")) {
      return(st_point(coords_rotated))
    } else if (inherits(geom, "LINESTRING")) {
      return(st_linestring(coords_rotated))
    } else if (inherits(geom, "POLYGON")) {
      return(st_polygon(list(coords_rotated)))
    } else if (inherits(geom, "MULTIPOINT")) {
      return(st_multipoint(coords_rotated))
    } else if (inherits(geom, "MULTILINESTRING")) {
      return(st_multilinestring(list(coords_rotated)))
    } else if (inherits(geom, "MULTIPOLYGON")) {
      return(st_multipolygon(list(list(coords_rotated))))
    } else {
      stop("Geometry type not supported")
    }
  })
  
  st_sfc(rotated_geometry, crs = st_crs(geometry))
}

# Example sf object with a point, a linestring, and a polygon
df <- data.frame(
  id = 1:3,
  geometry = st_sfc(
    st_point(c(1, 1)), # Point
    st_linestring(matrix(c(3, 3, 4, 4, 5, 3), ncol = 2, byrow = TRUE)), # Line
    st_polygon(list(matrix(c(6, 6, 7, 6, 7, 7, 6, 7, 6, 6), ncol = 2, byrow = TRUE))) # Polygon
  )
)

# Convert the data frame to an sf object
sf_obj <- st_sf(df)



sf_obj_rotated <- rotate_geometry(sf_obj, 100)


class(sf_obj)

m <- matrix(c(1,0.0,0.2,0.8),ncol = 2, byrow = TRUE)
st_geometry(sf_obj)*m -> sf_obj_rot
sf_obj*m
plot(sf_obj)
plot(sf_obj_rot)



angle_degrees <- 20
angle_radians <- angle_degrees * pi / 180
  
# Create a rotation matrix
rotation_matrix <- matrix(
  c(cos(angle_radians), -sin(angle_radians),
    sin(angle_radians),  cos(angle_radians)),
  ncol = 2, byrow = TRUE
)

round(rotation_matrix,5)

shiny_data$cell_object %>% filter(section_index == "A1", run_index == "spatial1") -> test
shiny_data$cell_object %>% filter(section_index == "A1", run_index == "spatial1") -> test2

geom <- test %>% st_geometry()
st_geometry(test2) <- geom*rotation_matrix

p <- ggplot() + geom_sf(test, mapping = aes(fill = "grey50"))
ggsave("data/spatial_check.png", p)

p2 <- ggplot() + geom_sf(test2, mapping = aes(fill = "grey50"))
ggsave("data/spatial_check2.png", p2)

run_index_ <- "spatial1"
section_index_ <- "A1"
spatial_object@meta.data %>% filter(run_index == run_index_, section_index == section_index_) -> test

test$run_index %>% unique()


test %>% 
  group_by(run_index) %>% 
  summarize(counts = n()) -> test2

bregma_tsv_path <- "./processed_data/ROI_Bregma_spatial.tsv"
spatial_object <- qread("processed_data/spatial_object_preprocessed.qs")
shiny_data <- qread("spatial_hypothalamus_app/data/shiny_data.qs")

colnames(shiny_data$cell_object_centroid$Slc24a4.y.y.y.y.y.y.y.y.y.y.y.y)
sf_object <- test
angle_degrees <- 20

rotation_df <- menu_converter

menu_converter <- read.table(bregma_tsv_path, sep = "\t", header = TRUE)


st_coordinates(sf_object) <- round(st_coordinates(test))


shiny_data$transcript_object

shiny_data$cell_object_centroid %>% 
  select(cell_IDs)


cell_sf_object


sf_object <- cell_sf_object

test <- add_count_matrix(cell_sf_object, spatial_object)

object.size(test)
shiny_data$cell_object_centroid
rownames(spatial_object)

colnames(total_count_df)
rownames(total_count_df)

total_count_df$transcript_ID
colnames(total_count_df1)
rownames(total_count_df1)

colnames(total_count_df1) <- total_count_df1[rownames(total_count_df1)== "transcript_ID"]

colnames(total_count_df1)
class(total_count_df1)
total_count_df1[rownames(total_count_df1) == "transcript_ID",]
total_count_df$transcript_ID
nrow(total_count_df1)


unique(test$cell_IDs == "transcript_ID")

shiny_data$cell_object_centroid %>% 
pivot_longer(-c("run_index", "section_index", "tooltip", "geometry"), names_to = "column", values_to = "value") %>%
  mutate(combined = paste(column, value, sep = ": ")) %>%
  select(-c("value")) %>% 
  pivot_wider(names_from = column, values_from = combined) -> test





columns_to_concatenate <- c("tooltip", "Agrp")

test %>% 
rowwise() %>%
  mutate(all_combined = paste(c_across(columns_to_concatenate), collapse = "<br>"))



test$all_combined



# Sample data frame
df <- data.frame(
  A = c(1, 2, 3),
  B = c(4, 5, 6),
  C = c(7, 8, 9)
)

# View the original data frame
print("Original Data Frame:")
print(df)

concant <- c("run_index","section_index")
# Create a new column that concatenates all other columns row-wise
df <- test %>%
  st_drop_geometry() %>% 
  mutate(all_combined = pmap_chr(select(., all_of(concant)), 
                                 ~ paste(concant, c(...), sep = ": ", collapse = "<br>")))

# View the transformed data frame
print("Transformed Data Frame:")
print(df$all_combined)




selected_section_index <- "A1"

selected_run_index <- "spatial1"
concant <- c("tooltip","Agrp")




test$tooltip


shiny_data$cell_object_centroid %>% colnames()

columns_to_concatenate <- c("Asb4", "Agrp")

test %>% 
  st_drop_geometry() %>%
  mutate(all_combined = pmap_chr(select(., all_of(columns_to_concatenate)), 
                                 ~ paste(..., sep = "<br>"))) -> test2

test2$all_combined


# Sample data frame
df <- data.frame(
  A = c(1, 2, 3),
  B = c(4, 5, 6),
  C = c(7, 8, 9)
)

# View the original data frame
print("Original Data Frame:")
print(df)

# Create a new column that concatenates all other columns row-wise
df <- df %>%
  mutate(all_combined = pmap_chr(select(., all_of(c("A","B"))), 
                                 ~ paste(c("A","B"), c(...), sep = ": ", collapse = ", ")))

# View the transformed data frame
print("Transformed Data Frame:")
print(df)




# Sample data frame
df <- data.frame(
  A = c(1, 2, 3),
  B = c(4, 5, 6),
  C = c(7, 8, 9),
  D = c(10, 11, 12)
)

# View the original data frame
print("Original Data Frame:")
print(df)

# List of columns to concatenate
columns_to_concatenate <- c("A", "C")

# Create a new column that concatenates the specified columns row-wise
df <- df %>%
  mutate(all_combined = pmap_chr(select(., all_of(columns_to_concatenate)), 
                                 ~ paste(..., sep = ", ")))

# View the transformed data frame
print("Transformed Data Frame:")
print(df)

# Sample data frame with string data
df <- data.frame(
  A = c("apple<br>sad", "banana: 1", "cherry"),
  B = c("dog", "elephant", "frog"),
  C = c("grape: 1", "honeydew", "iceberg"),
  D = c("jackfruit", "kiwi", "lemon")
)

# View the original data frame
print("Original Data Frame:")
print(df)

# List of columns to concatenate
columns_to_concatenate <- c("A", "C")

# Create a new column that concatenates the specified columns row-wise
df <- df %>%
  mutate(all_combined = pmap_chr(select(., all_of(columns_to_concatenate)), 
                                 ~ paste(..., sep = "<br>")))

# View the transformed data frame
print("Transformed Data Frame:")
print(df)




# Load the necessary libraries
library(dplyr)
library(purrr)
library(sf)  # for handling spatial data frames

# Sample spatial data frame (as an example, replace with your actual data)
# Assuming you have an sf object named `test`
# For demonstration, creating a mock sf data frame similar to your structure
test <- st_as_sf(
  data.frame(
    run_index = rep("spatial1", 10),
    section_index = rep("A1", 10),
    tooltip = rep("Virus1:...", 10),
    geometry = st_sfc(st_point(c(11896.65, -6638.505)), st_point(c(10304.13, -1323.483)),
                      st_point(c(9742.368, 11426.27)), st_point(c(10858.85, 10775.61)),
                      st_point(c(13036.51, 9522.605)), st_point(c(19522.71, 5787.451)),
                      st_point(c(15599.16, 8078.131)), st_point(c(13201.43, 9441.798)),
                      st_point(c(17150.69, 7190.257)), st_point(c(13532.02, 9272.107))),
    Acvr1c = rep("0", 10),
    Agrp = rep("Agrp:...", 10),
    Asb4 = rep("Asb4:...", 10)
  ), crs = NA
)

# View the original data frame
print("Original Data Frame:")
print(test)

# List of columns to concatenate
columns_to_concatenate <- c("Acvr1c", "Agrp", "Asb4")

# Create a new column that concatenates the specified columns row-wise
test <- test %>%
  st_drop_geometry() %>% 
  mutate(all_combined = pmap_chr(select(., all_of(columns_to_concatenate)), 
                                 ~ paste(..., sep = " ")))

# View the transformed data frame
print("Transformed Data Frame:")
print(test)


# Assuming `test` is your spatial data frame (sf object)
# For demonstration, creating a mock sf data frame similar to your structure
test <- st_as_sf(
  data.frame(
    run_index = rep("spatial1", 10),
    section_index = rep("A1", 10),
    tooltip = rep("Virus1:...", 10),
    geometry = st_sfc(st_point(c(11896.65, -6638.505)), st_point(c(10304.13, -1323.483)),
                      st_point(c(9742.368, 11426.27)), st_point(c(10858.85, 10775.61)),
                      st_point(c(13036.51, 9522.605)), st_point(c(19522.71, 5787.451)),
                      st_point(c(15599.16, 8078.131)), st_point(c(13201.43, 9441.798)),
                      st_point(c(17150.69, 7190.257)), st_point(c(13532.02, 9272.107))),
    Acvr1c = rep("0", 10),
    Agrp = rep("Agrp:...", 10),
    Asb4 = rep("Asb4:...", 10)
  ), crs = NA
)



# Sample data frame with string data
test <- data.frame(
  A = c("apple<br>sad", "banana: 1", "cherry"),
  B = c("dog", "elephant", "frog"),
  C = c("grape: 1", "honeydew", "iceberg"),
  D = c("jackfruit", "kiwi", "lemon")
)

# View the original data frame
print("Original Data Frame:")
print(test)

# List of columns to concatenate (excluding `geometry`)
columns_to_concatenate <- c("A", "C")

# Create a new column that concatenates the specified columns row-wise
test <- test %>%
  mutate(all_combined = paste(select(., all_of(columns_to_concatenate)), collapse = "<br>"))

# View the transformed data frame
print("Transformed Data Frame:")
print(test)




# Sample data frame with string data
test <- data.frame(
  A = c("apple<br>sad", "banana: 1", "cherry"),
  B = c("dog", "elephant", "frog"),
  C = c("grape: 1", "honeydew", "iceberg"),
  D = c("jackfruit", "kiwi", "lemon")
)



# View the original data frame
print("Original Data Frame:")
print(test)



columns_to_concatenate <- c("Acvr1c", "Agrp", "Asb4")
# Create a new column that concatenates specified columns row-wise
test <- test %>%
  rowwise() %>%
  mutate(all_combined = paste(c_across(columns_to_concatenate), collapse = "<br>"))

# View the transformed data frame
print("Transformed Data Frame:")
print(test$all_combined)


names(shiny_data)



spatial_object <- qread("processed_data/spatial_object_preprocessed.qs")


spatial_object@meta.data %>%
  rename("Cell area" = area, Virus1 = virus1) %>% colnames()

transcript_region_barplot <- shiny_data$transcript_region_table %>%
                filter(run_index == "spatial1", section_index == "A1", transcript == "Agrp") %>% 
                ggplot() + geom_col(aes(x = transcript, y = counts, fill = region), position = "dodge") + 
                geom_text(aes(x = transcript, y = counts,label = counts), 
  position = position_dodge(width = 0.9), 
  vjust = -0.5, 
  size = 4, 
  color = "black")


ggsave("data/spatial_check2.png", transcript_region_barplot)



transcript_region_barplot <- shiny_data$transcript_region_table %>%
  filter(run_index == "spatial1", section_index == "A1", transcript %in% c("Agrp", "Acvr1c")) %>% 
  ggplot(aes(x = transcript, y = counts, fill = region)) +
  geom_col(position = position_dodge()) +
  geom_text(aes(label = counts), 
  position = position_dodge(width = 0.9), 
  vjust = -0.5, 
  size = 4, 
  color = "black")


uiOutput("dynamic_plot")
    plot_height <- eventReactive(input$render_plot,{
        return(paste0(length(input$transcript_menu_item)*200,"px"))
    })

    output$dynamic_plot <- renderUI({
        plotOutput("transcript_region_barplot", width = "1100px", height = plot_height())
    })



shiny_data$cell_object %>% 
  filter(run_index == "spatial1", section_index == "A1") -> test
cell_region_barplot <- test %>%
            group_by(region,projection_termination) %>% 
            summarize(counts = n()) %>% 
            filter(!is.na(region)) %>% 
                ggplot(aes(x = projection_termination, y = counts, fill = region)) +
                geom_col(position = position_dodge()) +
                geom_text(
                    aes(label = counts), 
                    position = position_dodge(width = 0.9), 
                    vjust = -0.5, 
                    size = 4, 
                    color = "black"
                ) + 
                labs(title = "Projection termination cell count across regions in chosen section") + 
            theme_classic() + 
            theme(
                text = element_text(size = 15, color = "black"),
                plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
                axis.title.x = element_text(size = 17, face = "bold"),
                axis.title.y = element_text(size = 17, face = "bold"),
            ) + coord_flip()

ggsave("data/spatial_check2.png", cell_region_barplot)

names(shiny_data$cell_object_centroid)


test <- c("Agrp","Bace2")


concant <- c("Virus1","Virus2","Cell area",test)

cell_object_centroid <- shiny_data$cell_object_centroid %>% 
    filter(run_index == "spatial1", section_index == "A1")

print(paste(Sys.time(),"tooltip start"))
cell_object_centroid_tooltip <- cell_object_centroid %>% 
    st_drop_geometry() %>% 
    mutate(tooltip = pmap_chr(select(., all_of(concant)), ~ paste(concant, c(...), sep = ": ", collapse = "<br>")))

cell_object_centroid$tooltip <- cell_object_centroid_tooltip$tooltip
print(paste(Sys.time(),"tooltip start"))


shiny_data$cell_object$section_index %>% unique()

shiny_data$cell_object_centroid %>% colnames()


test <- qread("./data/shiny_data.qs")

test$transcript_object$transcript %>% unique() %>% sort()







spatial_object <- qread("processed_data/spatial_object_preprocessed.qs")

spatial_object$region_corrected <- correct_region_names(spatial_object$region)
spatial_object$projection_termination <- add_projection_termination(spatial_object, 0.05)





library(sf)
library(spatstat)
# Prep
shiny_data <- qread("spatial_hypothalamus_app/data/shiny_data_test.qs")

cell_object <- shiny_data$cell_object_centroid
cell_object %>% head()
region_object <- shiny_data$region_object

region_object %>% head()

cell_object_sub <- cell_object %>% 
  filter(run_index == "spatial1", section_index == "A1", Virus1 < 0.05)

region_object_sub <- region_object %>%
  filter(run_index == "spatial1", section_index == "A1") %>% 
  st_combine()

win <- as.owin(region_object_sub)

X <- as.ppp(st_coordinates(cell_object_sub), win)

plot(X)


# CSR 1)
Q <- quadratcount(X, nx = 10, ny = 10)

plot(X)
axis(1)
axis(2)
plot(Q, add = TRUE, cex = 2)


quadrat.test(Q, alternative = "clustered")$p.value

# CSR 2)

#Forskellige metoder + null
K <- Kest(X)
plot(K)


#Null med unvertianty 
E <- envelope(X, Kest, nsim = 19, verbose = FALSE) 
plot(E)

E %>% names()

dist <- 1000

r_index <- which(E$r > dist)[1]
E$hi[r_index] < E$obs[r_index]


# Intensity 1)
lambdahat <- density(X)
attr(lambdahat, "sigma")

# Save the plot as a PNG file
png("plots/density_plot.png")

# Plot the density
plot(lambdahat)

# Close the graphics device
dev.off()

plot(density(X), main = "Default bandwidth")
plot(density(X, sigma = 500))
plot(density(X, sigma = 1000))
plot(density(X, sigma = 2000))


lambdahat$v


# Intensity 2)

book.mesh.dual <- function(mesh) {
  if (mesh$manifold == 'R2') {
    ce <- t(sapply(1:nrow(mesh$graph$tv), function(i) colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
    library(parallel)
    pls <- mclapply(1:mesh$n, function(i) {
      p <- unique(Reduce('rbind', lapply(1:3, function(k) {
        j <- which(mesh$graph$tv[, k] == i)
        if (length(j) > 0)
          return(rbind(
            ce[j, , drop = FALSE],
            cbind(
              mesh$loc[mesh$graph$tv[j, k], 1] + mesh$loc[mesh$graph$tv[j, c(2:3, 1)[k]], 1],
              mesh$loc[mesh$graph$tv[j, k], 2] + mesh$loc[mesh$graph$tv[j, c(2:3, 1)[k]], 2]
            ) / 2
          ))
        else
          return(ce[j, , drop = FALSE])
      })))
      j1 <- which(mesh$segm$bnd$idx[, 1] == i)
      j2 <- which(mesh$segm$bnd$idx[, 2] == i)
      if ((length(j1) > 0) | (length(j2) > 0)) {
        p <- unique(rbind(
          mesh$loc[i, 1:2],
          p,
          mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2] / 2 + mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2] / 2,
          mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2] / 2 + mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2] / 2
        ))
        yy <- p[, 2] - mean(p[, 2]) / 2 - mesh$loc[i, 2] / 2
        xx <- p[, 1] - mean(p[, 1]) / 2 - mesh$loc[i, 1] / 2
      } else {
        yy <- p[, 2] - mesh$loc[i, 2]
        xx <- p[, 1] - mesh$loc[i, 1]
      }
      Polygon(p[order(atan2(yy, xx)), ])
    })
    return(SpatialPolygons(lapply(1:mesh$n, function(i) Polygons(list(pls[[i]]), i))))
  } else {
    stop("It only works for R2!")
  }
}


coo <- cell_object_sub %>%
  st_coordinates()


coo

bbox <- st_bbox(region_object_sub)

# Create a polygon from the bounding box
bbox_polygon <- st_as_sfc(bbox)

# Convert the sfc object to an sf object
bbox_sf <- st_as_sf(bbox_polygon)

# raster grid covering map
grid <- terra::rast(bbox_sf, nrows = 100, ncols = 100)


# coordinates of all cells
xy <- terra::xyFromCell(grid, 1:ncell(grid))

plot(xy)

# transform points to a sf object
dp <- st_as_sf(as.data.frame(xy), coords = c("x", "y"))
# indices points within the map
indicespointswithin <- which(st_intersects(dp, region_object_sub, sparse = FALSE))
# points within the map
dp <- st_filter(dp, region_object_sub)

coop <- st_coordinates(dp)

ggplot() + geom_sf(data = cell_object_sub)
ggplot() + geom_sf(data = region_object_sub)
ggplot() + geom_sf(data = dp)

min(xy[,'x'])

library(INLA)

loc.d <- cbind(st_coordinates(bbox_sf)[, 1], st_coordinates(bbox_sf)[, 2])
mesh <- inla.mesh.2d(loc.domain = loc.d, max.edge = c(500, 1000), offset = c(500, 1000), cutoff = 10)


p <- plot(mesh)
ggsave("plots/plot.png",p)

plot(mesh)
points(coo, col = "red") 
axis(1)
axis(2)

class(mesh)


# Open a PNG graphics device
png("plots/plot.png", width = 8000, height = 6000)

# Create your plot
plot(mesh)

# Close the graphics device
dev.off()

(nv <- mesh$n)
(n <- nrow(coo))

spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)


dmesh <- book.mesh.dual(mesh) 


# Open a PNG graphics device
png("plots/plot.png", width = 8000, height = 6000)

plot(dmesh)
axis(1)
axis(2)

# Close the graphics device
dev.off()

## Lave rplot function som plotter billede til suragate plot og s%>% have vs key board short cut der aabner plot?


# Domain polygon is converted into a SpatialPolygons
domain.polys <- st_polygon(list(loc.d))
domainSP <- domain.polys
# Because the mesh is larger than the study area, we need to
# compute the intersection between each polygon
# in the dual mesh and the study area

dmesh_sf <- st_as_sf(dmesh)

library(rgeos)
w <- sapply(1:length(dmesh), function(i) { 
    if (gIntersects(dmesh[i, ], domainSP))
        return(gArea(gIntersection(dmesh[i, ], domainSP))) 
    else return(0)
})

# Assuming dmesh and domainSP are sf objects
w <- sapply(1:nrow(dmesh_sf), function(i) {
    intersection <- st_intersection(dmesh_sf[i, ], domainSP)
    if (nrow(intersection) > 0) {
        return(st_area(intersection))
    } else {
        return(0)
    }
})

sum(w)
st_area(bbox_polygon)


png("plots/plot.png", width = 8000, height = 6000)

plot(mesh)
points(mesh$loc[which(w > 0), 1:2], col = "black", pch = 20) 
points(mesh$loc[which(w == 0), 1:2], col = "red", pch = 20)

# Close the graphics device
dev.off()



y.pp <- rep(0:1, c(nv, n))
e.pp <- c(w, rep(0, n))


head(cbind(y.pp, e.pp))
tail(cbind(y.pp, e.pp))


# Projection matrix for the integration points (mesh vertices)
A.int <- Diagonal(nv, rep(1, nv))
# Projection matrix for observed points (event locations) 
A.y <- inla.spde.make.A(mesh = mesh, loc = coo)
# Projection matrix for mesh vertices and event locations 
A.pp <- rbind(A.int, A.y)

Ap.pp <- inla.spde.make.A(mesh = mesh, loc = coop)


sum(A.y[1,which(A.y[1,] !=0 )])

# stack for estimation
stk.e.pp <- inla.stack(
  tag = "est.pp", 
  data = list(y = y.pp, e = e.pp), 
  A = list(1, A.pp), 
  effects = list(list(b0 = rep(1, nv + n)), 
  list(s = 1:nv))
)

# stack for prediction stk.p
stk.p.pp <- inla.stack(
  tag = "pred.pp", 
  data = list(y = rep(NA, nrow(coop)), 
  e = rep(0, nrow(coop))), 
  A = list(1, Ap.pp), 
  effects = list(data.frame(b0 = rep(1, nrow(coop))), 
  list(s = 1:nv))
)


# stk.full has stk.e and stk.p
stk.full.pp <- inla.stack(stk.e.pp, stk.p.pp)


formula <- y ~ 0 + b0 + f(s, model = spde)

res <- inla(
  formula, 
  family = 'poisson', 
  data = inla.stack.data(stk.full.pp),
  control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stk.full.pp)), 
  E = inla.stack.data(stk.full.pp)$e,
  verbose=TRUE
)


inla.pardiso.check()



library(raster)

resolution <- 100
r <- raster(bbox_sf, resolution = resolution)
(nrow <- nrow(r))
(ncol <- ncol(r))
nrow*ncol

grid <- rasterToPolygons(r)

r[] <- 0
tab <- table(cellFromXY(r, coo))
r[as.numeric(names(tab))] <- tab



png("plots/plot.png", width = 8000, height = 6000)

plot(r)

# Close the graphics device
dev.off()

grid <- grid[as.vector(t(matrix(1:nrow(grid), nrow = ncol, ncol = nrow))), ]

grid$id <- 1:nrow(grid)
grid$Y <- grid$layer
grid$cellarea <- resolution*resolution


#grid$cov <- extract(rcov, coordinates(grid))

gridmap <- raster::intersect(grid, bbox_sf)
grid <- grid[grid$id %in% gridmap$id, ]


#indNA <- which(is.na(grid$cov))
#indNA

#grid$cov[indNA] <- grid$cov[indNA+1]

grid$id2 <- grid$id

library(INLA)

formula <- Y ~ 1 +
  f(id, model="rw2d", nrow = nrow, ncol = ncol) +
  f(id2, model="iid")



res <- inla(formula, family = "poisson", data = grid@data,
            E = cellarea, control.predictor = list(compute = TRUE), verbose=TRUE)









# R code to prepare data and run the Stan model
library(rstan)

# Example data (replace with your actual data)
set.seed(42)
M <- 3
y <- c(23,12,11)  # Example: observed successes
n <- c(200,123,231)  # Example: random number of trials
theta <- y/n  # Example: samples from a Beta(2, 5) distribution

# Stan model code (as a string)
stan_model_code <- "
data {
  int<lower=0> M;  // number of samples
  vector<lower=0,upper=1>[M] theta;  // parameters from Beta distribution (samples)
  int<lower=0> n[M];  // number of trials for each sample
  int<lower=0> y[M];  // observed successes for each sample
}
parameters {
  real<lower=0> alpha;  // shape parameter of prior for Beta distribution
  real<lower=0> beta;   // rate parameter of prior for Beta distribution
}
model {
  // Priors for alpha and beta
  alpha ~ gamma(0.001, 0.001);
  beta ~ gamma(0.001, 0.001);
  
  // Likelihood
  for (i in 1:M) {
    // Prior for each theta_i from Beta distribution
    theta[i] ~ beta(alpha, beta);
    
    // Likelihood for each sample
    y[i] ~ binomial(n[i], theta[i]);
  }
}
"

# Compile the model
stan_model <- stan_model(model_code = stan_model_code)

# Sample from the posterior
fit <- sampling(stan_model, data = data_to_fit, iter = 10000, warmup = 500, chains = 4)


beta_dist <- marginalise_beta_dist(fit,1000)





# Convert to a data frame for plotting
beta_draws_df <- data.frame(beta_draws = beta_dist)

# Plot the density of the aggregated samples
ggplot(beta_draws_df, aes(x = beta_draws)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  labs(title = "Marginalized Beta Distribution",
       x = "Value",
       y = "Density") + xlim(c(0,0.3))



M <- 2
y <- c(71,70)  # Example: observed successes
n <- c(481,341)  # Example: random number of trials
theta <- y/n

library(tidyverse)


k <- 100

M <- 3
y <- c(23,12,11)*k  # Example: observed successes
n <- c(200,123,231)*k  # Example: random number of trials
theta <- y/n  # Example: samples from a Beta(2, 5) distribution

# Prepare data list for Stan
data_list <- list(M = M, theta = theta, n = n, y = y)


# Sample from the posterior
fit <- sampling(stan_model, data = data_list, iter = 2000, warmup = 500, chains = 4)


beta_dist <- marginalise_beta_dist(fit,100)

# Convert to a data frame for plotting
beta_draws_df2 <- data.frame(beta_draws = beta_dist)
beta_draws_df1$k <- 1
beta_draws_df2$k <- 5
beta_draws_df <- rbind(beta_draws_df1,beta_draws_df2)

# Plot the density of the aggregated samples
beta_draws_df %>% ggplot() +
  geom_density(aes(x = beta_draws, fill = as.factor(k)), alpha = 0.5) +
  labs(title = "Marginalized Beta Distribution",
       x = "Value",
       y = "Density") + xlim(c(0,0.3))


quantile(beta_dist,c(0.025, 0.975))

mean(theta)


# Print summary of the posterior
print(fit)

# Plot posterior distributions if needed
library(bayesplot)
mcmc_hist(fit, pars = c("alpha", "beta"))










# Example: Generate sample data from a beta distribution
set.seed(42)
samples <- rbeta(1000, shape1 = 2.0, shape2 = 5.0)

samples <- theta
# Fit a beta distribution to the sample data using MLE
fit2 <- fitdist(samples, "beta")

# Print the estimated parameters
print(fit2)
library(tidyverse)
fit2$estimate[1]


# Draw samples
set.seed(42)  # For reproducibility

# Convert to a data frame for plotting
beta_draws_df <- data.frame(beta_draws = rbeta(1000, fit2$estimate[1], fit2$estimate[2]))

# Plot the density of the aggregated samples
ggplot(beta_draws_df, aes(x = beta_draws)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  labs(title = "Marginalized Beta Distribution",
       x = "Value",
       y = "Density")




#library(rstan)
library(ggplot2)

# Assuming `fit` is your Stan fit object
posterior_samples <- extract(fit)
alpha_samples <- posterior_samples$alpha
beta_samples <- posterior_samples$beta

# Number of posterior samples
num_posterior_samples <- length(alpha_samples)

# Number of samples to draw from each Beta distribution
num_draws_per_sample <- 1000

# Initialize an empty vector to store the samples
beta_draws <- numeric(num_posterior_samples * num_draws_per_sample)

# Draw samples
set.seed(42)  # For reproducibility
for (i in 1:num_posterior_samples) {
  beta_draws[((i - 1) * num_draws_per_sample + 1):(i * num_draws_per_sample)] <-
    rbeta(num_draws_per_sample, alpha_samples[i], beta_samples[i])
}

# Convert to a data frame for plotting
beta_draws_df <- data.frame(beta_draws = beta_draws)

# Plot the density of the aggregated samples
ggplot(beta_draws_df, aes(x = beta_draws)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  labs(title = "Marginalized Beta Distribution",
       x = "Value",
       y = "Density")



quantile(beta_draws_df$beta_draws,c(0.05, 0.95))


fit$beta





expand.grid(1:10,1:10)

spatial_object <- qread("processed_data/spatial_object_preprocessed.qs")


spatial_object$region_corrected <- correct_region_names(spatial_object$region)
spatial_object$projection_termination <- add_projection_termination(spatial_object,0.05)
meta_data <- filter_data(spatial_object)

#relative_count_df %>% 

relative_count_df$region_corrected %>% unique() -> test1
relative_count_df$projection_termination %>% unique() -> test2

relative_count_df %>% 
  filter(region_corrected == "ventromedial", projection_termination == test2[1]) -> test3


M <- test3 %>% nrow()
theta <- test3$relative_count
n <- test3$total_count
y <- test3$count

row_num <- 3



beta_draws_df = data.frame(beta_draws = beta_dist)

# Plot the density of the aggregated samples
ggplot(beta_draws_df, aes(x = beta_draws)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  labs(title = "Marginalized Beta Distribution",
       x = "Value",
       y = "Density")




meta_data_selected_regions$region_corrected %>% unique()

data_to_fit$theta

fit2 <- fitdist(data_to_fit$theta, "beta")


# Convert to a data frame for plotting
beta_draws_df <- data.frame(beta_draws = rbeta(1000, fit2$estimate[1], fit2$estimate[2]))

# Plot the density of the aggregated samples
ggplot(beta_draws_df, aes(x = beta_draws)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  labs(title = "Marginalized Beta Distribution",
       x = "Value",
       y = "Density")


object.size(sc_object_neuron)

fgf1_spatial_object <- qread()





library(Seurat)
library(qs)
library(tidyverse)


sc_object_neuron <- qread("/projects/cbmr_shared/data/single_cell_fgf1_neuron.qs")
spatial_object <- qread("/projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/raw_data/obj_rctd_merged_01")

probes <- spatial_object %>% rownames()

subset_matrix <- sc_object_neuron[["RNA"]]@counts[probes,] 
sc_object_neuron_probes <- CreateSeuratObject(subset_matrix) 
meta_data <- sc_object_neuron@meta.data
sc_object_neuron_probes <- AddMetaData(object = sc_object_neuron_probes, metadata = meta_data) 


sc_object_neuron_probes

# pre-process dataset (without integration)
sc_object_neuron_probes <- NormalizeData(sc_object_neuron_probes)
sc_object_neuron_probes <- FindVariableFeatures(sc_object_neuron_probes)
sc_object_neuron_probes <- ScaleData(sc_object_neuron_probes)
sc_object_neuron_probes <- RunPCA(sc_object_neuron_probes)


qsave(sc_object_neuron_probes, "/projects/cbmr_shared/data/single_cell_fgf1_neuron_probes.qs")



spatial_object_sub_neuron <- subset(spatial_object, cell_class == "neuron")



pancreas.ref <- IntegrateLayers(object = pancreas.ref, method = CCAIntegration, orig.reduction = "pca",
                                new.reduction = "integrated.cca", verbose = FALSE)

pancreas.anchors <- FindTransferAnchors(reference = pancreas.ref, query = pancreas.query, dims = 1:30,
                                        reference.reduction = "pca")





pancreas.query <- MapQuery(
  anchorset = pancreas.anchors, 
  reference = pancreas.ref, 
  query = pancreas.query, 
  refdata = list(celltype = "celltype"), 
  reference.reduction = "pca", 
  reduction.model = "umap_integrated"
)

p1 <- DimPlot(pancreas.ref, reduction = "umap_integrated", group.by = "celltype", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")





head(beta_dist_quantiles_total)


beta_dist_quantiles <- beta_dist_quantiles_total



beta_dist_quantiles %>% 
filter(projection_termination != "No-virus") %>% 
    pivot_longer(
        cols = c(q0.025, median, q0.975),
        names_to = "quantiles",
        values_to = "value"
    ) %>% 
    filter(region == "lateral_hypothalamus")

for (pair in region_projection_pairs){
  print(pair)
}



quants <- quantile(beta_dist,c(0.025,0.5, 0.975))

beta_dist_quantiles <- data.frame(
    region = region, 
    projection_termination = projection_termination, 
    q0.025 = quants[1], 
    median = quants[2], 
    q0.975 = quants[3]
)

beta_dist_quantiles_total <- rbind(beta_dist_quantiles_total,beta_dist_quantiles)



beta_dist_quantiles %>% 
filter(projection_termination != "No-virus") %>% 
  pivot_longer(
      cols = c(q0.025, median, q0.975),
      names_to = "quantiles",
      values_to = "value"
  ) %>% 
  ggplot() + geom_point(aes(x = region, y = value)) + facet_grid(~projection_termination) + theme(axis.text.x = element_text(angle = 90, hjust = 1))


  plot(1:10,1:10)

  ggsave("plots/plot.png", p, width = 10000, height = 5000, units = "px")



beta_dists <- beta_dists_bayesian
beta_dists %>% colnames()

plot_path <- "./plots/projection_termination_region_plot_bayesian.png"



# Example data
data <- data.frame(
  x = c(1, 2, 3, 4),
  y = c(4, 3, 2, 1)
)

# Plotting the data with a horizontal line at y = 2.5
ggplot(data, aes(x = x, y = y)) +
  geom_point() +  # Plotting data points
  geom_hline(yintercept = 2.5, linetype = "dashed", color = "red", size = 1) +
  labs(title = "Horizontal Line in ggplot", x = "X axis", y = "Y axis") +
  theme_minimal()



projection_terminations <- 1:2




filter_data <- function(spatial_object){
    meta_data <- spatial_object@meta.data %>% 
        filter(!(region_corrected == "lateral_hypothalamus_right" & run_index == "spatial2"))
    
    meta_data$region_corrected <- meta_data$region_corrected %>% 
        gsub("_right|_left", "", .)

    meta_data_selected_regions <- meta_data %>%
        filter(region_corrected %in% c("arcuate", "ventromedial", "dorsomedial", "lateral_hypothalamus"))
    
    return(meta_data_selected_regions)
}




spatial_object <- qread("processed_data/spatial_object_preprocessed.qs")

spatial_object$region_corrected <- correct_region_names(spatial_object$region)


spatial_object$projection <- add_projection(
  spatial_object$virus1, 
  spatial_object$virus2,
  spatial_object$run_index,
  threshold = 0.05
)

spatial_object$region_reduced <- reduce_regions(spatial_object@meta.data)


meta_data$region_corrected[1:100]
meta_data$regions_reduced[1:100]

meta_data %>% filter(run_index == "spatial2") %>% select(region_corrected, regions_reduced) %>% head(1000)

region_umap.png

column <- 'projection'

plotting_colors %>% class()
plotting_colors[3]

spatial_object@meta.data[[column]] %>% class()

spatial_object$projection %>% unique()

cells <- NULL

any(c("A", "B") == "A")

test <- c("A", "B")
test["Others"] <- "C"


spatial_object@meta.data %>% 
    filter(single_projection != "Other") %>% 
    rownames()



plot_projection_umaps <- function(spatial_object, plot_folder){

    projection_no_double <- case_when(
        grepl("\\+", spatial_object$projection) ~ "double",
        TRUE ~ spatial_object$projection
    )

    names(projection_no_double) <- NULL

    spatial_object$projection_no_double <- projection_no_double

    projection_termination_no_double_umap <- DimPlot(spatial_object, group.by = 'projection_no_double', reduction = 'umap_features', raster = FALSE, order = TRUE)

    ggsave(paste0(plot_folder,"/projection_no_double_umap.png"), projection_termination_no_double_umap, width = 4000, height = 3000, units = "px")


    spatial_object$projection_all <- case_when(
        spatial_object$sun1 == TRUE ~ TRUE,
        spatial_object$virus1 < 0.05 ~ TRUE,
        spatial_object$virus2 < 0.05 ~ TRUE,
        TRUE ~ FALSE
    )

    projection_all_umap <- DimPlot(spatial_object, group.by = 'projection_all', reduction = 'umap_features', raster = FALSE, order = TRUE)

    ggsave(paste0(plot_folder,"/projection_all_umap.png"), projection_all_umap, width = 4000, height = 3000, units = "px")

    spatial_object$projection_sun1_vs_virus <- case_when(
        spatial_object$sun1 == TRUE ~ "Sun1",
        spatial_object$virus1 < 0.05 | spatial_object$virus2 < 0.05 ~ "Virus",
        spatial_object$sun1 == TRUE & (spatial_object$virus1 < 0.05 | spatial_object$virus2 < 0.05) ~ "Sun1-Virus",
        TRUE ~ "Non-projecting"
    )

    projection_sun1_vs_virus_umap <- DimPlot(spatial_object, group.by = 'projection_sun1_vs_virus', reduction = 'umap_features', raster = FALSE, order = TRUE)

    ggsave(paste0(plot_folder,"/projection_sun1_vs_virus_umap.png"), projection_sun1_vs_virus_umap, width = 4000, height = 3000, units = "px")

}





selected_regions <- spatial_object$region_corrected %>% 
    unique() %>% 
    grep("arcuate|ventromedial|dorsomedial|lateral_hypo|paraventricular",.,value = TRUE)

#plotting_colors <- generate_label_colors(selected_regions)


# Across section and run porbbaly does not make sense
#spatial_object@meta.data %>% 
    #filter(projection_termination != "No-virus") %>%
    #filter(region_corrected %in% selected_regions) %>% 
    #ggplot() + geom_bar(aes(x = projection_termination, fill = region_corrected), position = 'dodge') + facet_grid(section_index ~ run_index, scales = "free_x") + 
    #scale_fill_manual(values = plotting_colors)
    #theme(axis.text.x = element_text(angle = 90, hjust = 1))



# Does it make sense to look at termination projections across section? Probably should be regions instead. Is done below
projection_termination_plots <- function(spatial_object, projection_termination_column_name, plot_path){

    projection_terminations <- spatial_object@meta.data %>% 
        group_by(run_index, section_index, !!sym(projection_termination_column_name)) %>%
        summarize(count = n()) %>% 
        filter(!!sym(projection_termination_column_name) != "No-virus")

    projection_terminations <- projection_terminations %>% 
        mutate(!!sym(projection_termination_column_name) := factor(
            !!sym(projection_termination_column_name), 
            levels = c("BST","PAG", "BST+PAG", "LHA", "PBN", "LHA+PBN", "PVT", "CEA", "PVT+CEA")
        ))

    projection_terminations_plot <- projection_terminations %>%
        ggplot() + geom_col(aes(x = !!sym(projection_termination_column_name), y = count)) + facet_grid(section_index~run_index, scales = "free_x") + 
        geom_text(aes(x = !!sym(projection_termination_column_name), y = count + 100, label = count))

    ggsave(plot_path, projection_terminations_plot)
}

#Useless now that we have shiny app?
spatial_projection_plots <- function(spatial_object,plot_path){

    spatial_object$virus <- case_when(
        spatial_object$virus1 < 0.05 & spatial_object$virus2 > 0.05 ~ "virus1",
        spatial_object$virus1 > 0.05 & spatial_object$virus2 < 0.05 ~ "virus2",
        spatial_object$virus1 > 0.05 & spatial_object$virus2 > 0.05 ~ "No-virus",
        spatial_object$virus1 < 0.05 & spatial_object$virus2 < 0.05 ~ "virus1+virus2",
    )



    levels <- unique(spatial_object$virus)


    num_levels <- spatial_object$projection_termination_no_double %>% unique() %>% length()


    levels <- unique(spatial_object$virus)
    default_colors <- hue_pal()(length(levels))
    names(default_colors) <- levels

    default_colors["No-virus"] <- "#7F7F7F"


    run_indices <- c("spatial1", "spatial2", "spatial3")
    section_indices <- c("A1", "B1", "C1", "D1", "A2", "B2", "C2", "D2")

    for (run_index_number in 1:3){
        for (section_index in section_indices){
            run_index_number_adjusted <- ifelse(run_index_number == 1, "", paste0(".",run_index_number))
            fov <- paste0(section_index,run_index_number_adjusted)

            DefaultBoundary(spatial_object[[fov]]) <- "segmentation"
            spatial_plot <- ImageDimPlot(spatial_object, 
                    fov = fov,
                    axes = F, 
                    dark.background = T, 
                    size = 1.5,
                    border.color = "black", 
                    border.size = 0.1, 
                    cols = default_colors, 
                    group.by = "virus",
                    molecules = c("p147_mCherry_2A_iCre","p36_NLS_Cre"),
                    nmols = 2000000,
                    mols.size = 0.05,
                    mols.alpha = 0.1
            )

            ggsave(paste0("./plots/spatial_projection_plot_",run_indices[i],"_",section_index,".png"),plot = spatial_plot, width = 6000, height = 7500, units = "px")
        }
    }
}





sun1_vs_virus <- spatial_object$sun1_vs_virus

spatial_object@meta.data %>% 
  filter(sun1 == TRUE, virus1 < 0.05) %>% 
  head(100)


sun1_virus_table["Other"]





spatial_object <- qread("processed_data/spatial_object_preprocessed.qs")

meta_data <- spatial_object@meta.data

sun1_vs_virus <- meta_data$sun1_vs_virus


spatial_object$virus1_count <- 1000000
spatial_object$virus2_count <- 1000000





shiny_data <- qread("spatial_hypothalamus_app/data/shiny_data.qs")

cell_object <- shiny_data$cell_object


cell_object <- cell_object %>% 
  rename(projection = projection_termination)



p <- plot(cell_object_dist)

ggsave("plots/test.png", K_func_list_dot_plot)

png("plots/test.png")

plot(cell_object_dist)
dev.off()


cell_object_dist$r
path_to_plot_folder <- "plots"


cell_object <- qread("processed_data/cell_object_test.qs")



region_paths <- c()
for (run_index in run_indices){
    region_paths <- c(region_paths,paste0(paste0("raw_data/",run_index),"/",list.files(path = paste0("raw_data/",run_index), pattern = "regions")))
}



region_ <- "lateral_amygdaloid_left"
projection_ <- "PVT"
section_index_ <- "A1"

plot_folder <- "plots"


# import functioner (mangler kun create_shiny_data)
# Lave cats, (X)
# Create_shiny_data restructure
# then test snakemake pipeline (cats before because nice with testing)
# Then documentation, then commit



# Total number of iterations
n <- 99

# For loop
for (i in 1:n) {
  # Simulate work (for example, a short delay)
  Sys.sleep(0.1)
  
  # Check progress and print every 10%
  if (i %% (n / 10) == 0) {
    # Calculate the percentage completed
    percent_completed <- (i / n) * 100
    # Create a formatted string with the percentage
    message <- paste0(utils$current_time()," Completed ",percent_completed, "%\n")
    # Print the message to the console
    cat(message)
  }
}

loop_progress(11,100)



spatial_object <- qread("processed_data/spatial_object_preprocessed.qs")
enrichment_matrix <- calculate_fisher_region_cell_property_enrichment(spatial_object$region_corrected, spatial_object$projection)

regions <- spatial_object$region_corrected
cell_properties <- spatial_object$projection

region <- region_unique[1]
cell_property <- cell_property_unique[1]



test <- table(regions_test,cell_properties_test) %>% 
                as.data.frame() %>% 
                arrange(desc(regions_test != "Others"),desc(cell_properties_test != "Others")) %>% 
                pivot_wider(names_from = cell_properties_test, values_from = Freq)


test %>% select(Other)


test2 <- data.frame(x = 1:10, y = 1:10*2)

test2 %>% dplyr::select(-x)

row_num <- 24

data_to_fit$M



3 %% (11 / 10) == 0



# Define the total number of iterations
total_iterations <- 99

# Calculate the step size for 10% of the loop
step_size <- ceiling(total_iterations / 10)

# For loop
for (i in 1:total_iterations) {
  # Perform your operations here
  Sys.sleep(0.001)
  # Print message every 10% of the loop
  if (i %% step_size == 0) {
    cat("Completed", (i / total_iterations) * 100, "%\n")
  }
}


# Define the total number of iterations
total_iterations <- 99

# Initialize progress tracking variables
next_print <- 10
percent_complete <- 0

# For loop
for (i in 1:total_iterations) {
  # Perform your operations here
  
  # Calculate percentage completion
  percent_complete <- (i / total_iterations) * 100
  
  # Print message every 10% of the loop
  if (percent_complete >= next_print) {
    cat("Completed", round(percent_complete), "%\n")
    next_print <- next_print + 10
  }
}


# Define the total number of iterations
total_iterations <- 99

# Calculate the step size for 10% of the loop
step_size <- total_iterations / 10

# For loop
for (i in 1:total_iterations) {
  # Perform your operations here
  
  # Print message every 10% of the loop
  if (i %% step_size == 0 || i == total_iterations) {
    cat("Completed", round((i / total_iterations) * 100), "%\n")
  }
}



# Define the total number of iterations
total_iterations <- 99

# Calculate the iterations that correspond to each 10% mark
percent_marks <- seq(0.1, 1.0, by = 0.1) * total_iterations
percent_marks <- round(percent_marks)

# For loop
for (i in 1:total_iterations) {
  # Perform your operations here
  Sys.sleep(0.001)
  # Print message at the 10% marks
  if (i %in% percent_marks) {
    cat("Completed", round((i / total_iterations) * 100), "%\n")
  }
}

i <- 50



# Define the total number of iterations
total_iterations <- 123123

# Calculate the iterations that correspond to each 10% mark
percent_marks <- ceiling(seq(0.1, 1.0, by = 0.1) * total_iterations)

# For loop
for (i in 1:total_iterations) {
  # Perform your operations here
  
  # Print message at the 10% marks
  if (i %in% percent_marks) {
    cat("Completed", sum(i >= percent_marks)*10, "%\n")
  }
}



loop_progress <- function(current_iteration, total_number_iteration){
    # Calculate the iterations that correspond to each 10% mark
    percent_marks <- ceiling(seq(0.1, 1.0, by = 0.1) * total_number_iteration)

    # Check progress and print every 10%
    if (current_iteration %in% percent_marks) {
        percent_completed <- sum(i >= percent_marks)*10
        return(cat(paste0(utils$current_time(),"  - Completed ",percent_completed, "%\n")))
    }
}



loop_progress <- function(current_iteration, total_number_iteration){
    # Calculate the iterations that correspond to each 10% mark
    percent_marks <- ceiling(seq(0.1, 1.0, by = 0.1) * total_number_iteration)

    # Check progress and print every 10%
    if (current_iteration %in% percent_marks) {
        percent_completed <- sum(current_iteration >= percent_marks)*10
        return(cat(paste0(utils$current_time(),"  - Completed ",percent_completed, "%\n")))
    }
}

loop_progress(3,3)






spatial_object <- qread("processed_data/spatial_object_preprocessed.qs")

spatial_object$region_corrected %>% unique() %>% grep(c("dorsomedial"),., value = TRUE)
