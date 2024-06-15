library(Seurat)
library(qs)
library(tidyverse)
library(sf)
library(RImageJROI)
library(progress)
library(terra)
library(reticulate)
library(scales)


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