library(reticulate)
library(sf)
library(tidyverse)
library(Seurat)
library(qs)


load_transcripts <- function(transcripts_path){
    transcripts <- read.table(transcripts_path, header = F, sep = "\t") %>% filter(!str_detect(V4, "FP "))
    transcript_points <- transcripts %>% st_as_sf(.,coords= c("V1","V2"))
    
    return(transcript_points)
}

array_to_sf_rois <- function(array_path){
    # Load the NumPy module
    np <- import("numpy")

    # Load the .npz file
    arrays <- np$load(array_path)

    rois <- list()
    for (roi in arrays$files){
        roi_coords <- arrays$get(roi)

        #Closing the roi coords
        roi_coords_closed <- rbind(roi_coords,roi_coords[1,])
        rois[[roi]] <- st_polygon(list(roi_coords_closed))
    }

    return(rois %>% st_sfc())
}

create_count_matrix <- function(transcripts, rois){
    sparse_intersects <- st_intersects(transcripts,rois,sparse = T)
    names(sparse_intersects) <- transcripts$V4
    sparse_count_matrix <- as(table(names(unlist(sparse_intersects)),unlist(sparse_intersects)),"sparseMatrix")
    
    colnames(sparse_count_matrix) <- paste0("Cell",colnames(sparse_count_matrix))
    return(sparse_count_matrix)
}

create_seurat_segmentations <- function(rois,sparse_count_matrix){
    rois_df <- rois %>% st_coordinates() %>% as.data.frame() %>% mutate(L2 = paste0("Cell",L2)) %>% filter(L2 %in% colnames(sparse_count_matrix))
    seurat_segmentations <- rois_df[,c("L2","X","Y")]
    colnames(seurat_segmentations) <- c("cell","x","y")

    return(seurat_segmentations)
}

create_seurat_centroids <- function(rois,sparse_count_matrix){
    centroids_df <- rois %>% st_centroid() %>% st_coordinates()
    seurat_centroids <- cbind(1:nrow(centroids_df),centroids_df) %>% as.data.frame() %>% mutate(V1 = paste0("Cell",V1)) %>% filter(V1 %in% colnames(sparse_count_matrix))
    colnames(seurat_centroids) <- c("cell","x","y")

    return(seurat_centroids)
}

create_seurat_molecules <- function(transcripts){
    seurat_molecules <- transcripts %>% st_coordinates() %>% as.data.frame() %>% cbind(.,transcript_points$V4)
    colnames(seurat_molecules) <- c("x","y","gene")

    return(seurat_molecules)
}

create_seurat_object <- function(sparse_count_matrix, seurat_segmentations, seurat_centroids, seurat_molecules, section){
    coords_list <- list(
       "centroids" = SeuratObject::CreateCentroids(seurat_centroids),
       "segmentation" = SeuratObject::CreateSegmentation(seurat_segmentations)
    )

    coords <- SeuratObject::CreateFOV(
       coords = coords_list,
       type = c("segmentation", "centroids"),
       molecules = seurat_molecules,
       assay = 'raw'
    )

    spatial_obj <- SeuratObject::CreateSeuratObject(counts = sparse_count_matrix, assay = 'raw')

    spatial_obj[[section]] <- coords

    return(spatial_obj)
}






sections <- c("A1","B1","C1","D1","A2","B2","C2","D2")
run <- "spatial1"

transcripts_files <- list.files(path = paste0("./raw_data/",run,"/")) %>% grep("results", ., value = TRUE)
rois_files <- list.files(path = paste0("./processed_data/",run,"/")) %>% grep("rois", ., value = TRUE)

spatial_objects = list()

t <- 0
for (section in sections){
    transcript_path <- paste0("./raw_data/",run,"/",grep(section, transcripts_files, value = TRUE))
    roi_path <- paste0("./processed_data/",run,"/",grep(section, rois_files, value = TRUE))

    print("Loading transcripts")
    transcript_points <- load_transcripts(transcript_path)
    
    print("Loading ROIs")
    rois_polygons <- array_to_sf_rois(roi_path)

    print("Generating count matrix")
    sparse_count_matrix <- create_count_matrix(transcript_points,rois_polygons)
    seurat_segmentations <- create_seurat_segmentations(rois_polygons,sparse_count_matrix)
    seurat_centroids <- create_seurat_centroids(rois_polygons,sparse_count_matrix)
    seurat_molecules <- create_seurat_molecules(transcript_points)

    print("Creating seurat object")
    spatial_obj <- create_seurat_object(sparse_count_matrix,seurat_segmentations,seurat_centroids,seurat_molecules,section)
    
    print("Adding meta data")
    spatial_obj$run <- run
    spatial_obj$section <- section

    areas <- rois_polygons %>% st_area()
    names(areas) <- names(rois_polygons)
    
    spatial_obj <- AddMetaData(spatial_obj,areas,"area")

    spatial_objects[[section]] <- spatial_obj

    t <- t + 1
    print(paste0(t,"/",length(sections),": ",section," done!"))
}

print("Merging spatial seurat objects")
merged_spatial_object <- merge(spatial_objects)

print("Saving merged spatial seurat object")
qsave(merged_spatial_object,"/processed_data/",run,"/",run,"_seurat_raw.qs")



