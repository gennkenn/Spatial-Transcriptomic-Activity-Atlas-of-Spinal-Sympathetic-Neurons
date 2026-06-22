## Load Libraries
library(Seurat)
library(future)
plan("multisession", workers = 8)
options(future.globals.maxSize = 256 * 1024^3)

library(ggplot2)
library(dplyr)
library(spacexr)
library(patchwork)
library(Seurat)
library(data.table)
library(sf)
library(geojsonsf)
library(RColorBrewer)

## Load Libraries
library(Seurat)
library(future)
plan("multisession", workers = 8)
library(ggplot2)
library(dplyr)
library(spacexr)
library(patchwork)
library(Seurat)
library(data.table)
library(sf)
library(geojsonsf)
library(RColorBrewer)
library(Seurat)
library(SeuratData)
library(SeuratDisk)

## Redefine ReadXenium()→for 5K panel
ReadXenium <- function (data.dir, outs = c("matrix", "microns"), type = "centroids", 
                        mols.qv.threshold = 20) 
{
  type <- match.arg(arg = type, choices = c("centroids", "segmentations"), 
                    several.ok = TRUE)
  outs <- match.arg(arg = outs, choices = c("matrix", "microns"), 
                    several.ok = TRUE)
  outs <- c(outs, type)
  has_dt <- requireNamespace("data.table", quietly = TRUE) && 
    requireNamespace("R.utils", quietly = TRUE)
  data <- sapply(outs, function(otype) {
    switch(EXPR = otype, matrix = {
      matrix <- suppressWarnings(Read10X(data.dir = file.path(data.dir, 
                                                              "cell_feature_matrix/")))
      matrix
    }, centroids = {
      if (has_dt) {
        cell_info <- as.data.frame(data.table::fread(file.path(data.dir, 
                                                               "cells.csv.gz")))
      } else {
        cell_info <- read.csv(file.path(data.dir, "cells.csv.gz"))
      }
      cell_centroid_df <- data.frame(x = cell_info$x_centroid, 
                                     y = cell_info$y_centroid, cell = cell_info$cell_id, 
                                     stringsAsFactors = FALSE)
      cell_centroid_df
    }, segmentations = {
      if (has_dt) {
        cell_boundaries_df <- as.data.frame(data.table::fread(file.path(data.dir, 
                                                                        "cell_boundaries.csv.gz")))
      } else {
        cell_boundaries_df <- read.csv(file.path(data.dir, 
                                                 "cell_boundaries.csv.gz"), stringsAsFactors = FALSE)
      }
      names(cell_boundaries_df) <- c("cell", "x", "y")
      cell_boundaries_df
    }, microns = {
      
      transcripts <- arrow::read_parquet(file.path(data.dir, "transcripts.parquet"))
      transcripts <- subset(transcripts, qv >= mols.qv.threshold)
      
      df <- data.frame(x = transcripts$x_location, y = transcripts$y_location, 
                       gene = transcripts$feature_name, stringsAsFactors = FALSE)
      df
    }, stop("Unknown Xenium input type: ", otype))
  }, USE.NAMES = TRUE)
  return(data)
}





# assay?の名前が"RNA"になってしまうので、"Xenium"に変更
# num_slicesは入れるROIの数に合わせて適宜変更
prep_FindMarkers <- function(obj, num_slices =10) {
  # Loop through total number slices
  for (i in 1:num_slices) {
    slot(object = obj@assays$SCT@SCTModel.list[[i]], name = "umi.assay") <- "Xenium"
  }
  obj <- obj %>% PrepSCTFindMarkers()
  # Return the modified obj
  return(obj)}

#
PolygonCropSeurat <- function(object, polygons, pyxel_size = 0.2125){
  #' Crop a Seurat Object using a list of Polygons
  #' 
  #' @description This function crops all FOVs included in the list, other FOVs included unchanged
  #' 
  #' The FOVs to be cropped must contain centroids as those are used to determine if the cells are in the polygon or not. 
  #
  #' @param object Seurat object to be cropped
  #' @param polygons a named list of sf POLYGON object, the names must match those of the FOVs to be cropped.
  #' @param pixel_size ratio between the FOV coordinates and the sf Objects coordinates.
  #' @usage PolygonCropSeurat(object polygons, pyxel_size = 0.2125)
  #' @return A new Seurat object, including the retained cells for the cropped FOVs and all cells for the other FOVs.
  #' @details It relies on `PolygonCropFOV` so if a cropped FOV would contain no cells, the uncropped FOV is included instead
  
  # Crop all FOVs for which a polygon has been provided
  message("Start cropping fovs...")
  cropped_fovs <- lapply(names(polygons), function(name){
    fov <- object@images[[name]]
    fov <- PolygonCropFOV(fov, polygons[[name]], pyxel_size)
    message("Done: ", name)
    return(fov)
  })
  names(cropped_fovs) <- names(polygons)
  message("Done cropping fovs...")
  
  # Make a list of all FOVs, cropped and non-cropped
  message("Gathering fovs...")
  untouched_fovs <- object@images[!(names(object@images) %in% names(cropped_fovs))]
  all_fovs <- c(untouched_fovs, cropped_fovs)
  message("Found fovs: ", names(all_fovs))
  
  # Get the cells to keep and make a subset of the object
  message("Making new object fovs...")
  cells_to_keep <- unlist(sapply(all_fovs, Cells))
  new_object <- subset(object, cells = cells_to_keep)
  new_object@images <- all_fovs # Replace the old FOVs
  message("Done cropping Seurat object")
  
  return(new_object)
}


## Load Xenium Data
data.ctrl1 <- ReadXenium(data.dir = "~/output-XETG00312__0043214__XP009_R03_Ctrl_1__20250319__064312", outs = c("matrix", "microns"), type = c("centroids", "segmentations"))
data.ctrl2 <- ReadXenium(data.dir = "~/output-XETG00312__0043214__XP009_R01_Ctrl_2__20250319__064312", outs = c("matrix", "microns"), type = c("centroids", "segmentations"))
data.ctrl3 <- ReadXenium(data.dir = "~/output-XETG00312__0043214__XP009_R02_Ctrl_3__20250319__064312", outs = c("matrix", "microns"), type = c("centroids", "segmentations"))
data.cold1<- ReadXenium(data.dir = "~/output-XETG00312__0043206__XP009_L01_Cold_1__20250319__064312", outs = c("matrix", "microns"), type = c("centroids", "segmentations"))
data.cold2<- ReadXenium(data.dir = "~/output-XETG00312__0043206__XP009_L02_Cold_2__20250319__064312", outs = c("matrix", "microns"), type = c("centroids", "segmentations"))
data.cold3 <- ReadXenium(data.dir = "~/output-XETG00312__0043206__XP009_L03_Cold_3__20250319__064312", outs = c("matrix", "microns"), type = c("centroids", "segmentations"))
data.coldF <-ReadXenium(data.dir = "~/output-XETG00312__0094296__XP015_R01_female_cold__20260220__080935", outs = c("matrix", "microns"), type = c("centroids", "segmentations"))
data.ctrlF <-ReadXenium(data.dir = "~/output-XETG00312__0094293__XP015_L01_female_RT__20260220__080935", outs = c("matrix", "microns"), type = c("centroids", "segmentations"))
data.2DGF <- ReadXenium(data.dir="~/output-XETG00312__0094756__XP016_L01_2DG_female__20260424__073501", outs = c("matrix", "microns"), type = c("centroids", "segmentations"))
data.2DGM <- ReadXenium(data.dir="~/output-XETG00312__0094768__XP016_R01_2DG_male__20260424__073502", outs = c("matrix", "microns"), type = c("centroids", "segmentations"))


# 
segmentations.data.ctrl1 <- list(
  centroids = CreateCentroids(data.ctrl1$centroids),
  segmentation = CreateSegmentation(data.ctrl1$segmentations))
coords.ctrl1 <- CreateFOV(
  coords = segmentations.data.ctrl1, 
  type = c("segmentation", "centroids"), 
  molecules = data.ctrl1$microns, 
  assay = "Xenium")
xenium.obj.ctrl1 <- CreateSeuratObject(
  counts = data.ctrl1$matrix[["Gene Expression"]], 
  assay = "Xenium", project = "Ctrl1")
xenium.obj.ctrl1[["BlankCodeword"]] <- CreateAssayObject(counts = data.ctrl1$matrix[["Unassigned Codeword"]])
xenium.obj.ctrl1[["ControlCodeword"]] <- CreateAssayObject(counts = data.ctrl1$matrix[["Negative Control Codeword"]])
xenium.obj.ctrl1[["ControlProbe"]] <- CreateAssayObject(counts = data.ctrl1$matrix[["Negative Control Probe"]])
xenium.obj.ctrl1[["fov"]] <- coords.ctrl1
xenium.obj.ctrl1$Condition <- 'Ctrl'

xenium.obj.ctrl1 <- subset(xenium.obj.ctrl1, subset = nCount_Xenium > 0)
# xenium.obj.ctrl1 <- subset(xenium.obj.ctrl1,cells = CellID_ctrl1_Slc5a7$Barcode)

print(xenium.obj.ctrl1)

# 
segmentations.data.ctrl2 <- list(
  centroids = CreateCentroids(data.ctrl2$centroids),
  segmentation = CreateSegmentation(data.ctrl2$segmentations))
coords.ctrl2 <- CreateFOV(
  coords = segmentations.data.ctrl2, 
  type = c("segmentation", "centroids"), 
  molecules = data.ctrl2$microns, 
  assay = "Xenium")
xenium.obj.ctrl2 <- CreateSeuratObject(
  counts = data.ctrl2$matrix[["Gene Expression"]], 
  assay = "Xenium", project = "Ctrl2")
xenium.obj.ctrl2[["BlankCodeword"]] <- CreateAssayObject(counts = data.ctrl2$matrix[["Unassigned Codeword"]])
xenium.obj.ctrl2[["ControlCodeword"]] <- CreateAssayObject(counts = data.ctrl2$matrix[["Negative Control Codeword"]])
xenium.obj.ctrl2[["ControlProbe"]] <- CreateAssayObject(counts = data.ctrl2$matrix[["Negative Control Probe"]])
xenium.obj.ctrl2[["fov"]] <- coords.ctrl2
xenium.obj.ctrl2$Condition <- 'Ctrl'
xenium.obj.ctrl2 <- subset(xenium.obj.ctrl2, subset = nCount_Xenium > 0)
# xenium.obj.ctrl2 <- subset(xenium.obj.ctrl2,cells = CellID_ctrl2_Slc5a7$Barcode)

print(xenium.obj.ctrl2)


# 
segmentations.data.ctrl3 <- list(
  centroids = CreateCentroids(data.ctrl3$centroids),
  segmentation = CreateSegmentation(data.ctrl3$segmentations))
coords.ctrl3 <- CreateFOV(
  coords = segmentations.data.ctrl3, 
  type = c("segmentation", "centroids"), 
  molecules = data.ctrl3$microns, 
  assay = "Xenium")
xenium.obj.ctrl3 <- CreateSeuratObject(
  counts = data.ctrl3$matrix[["Gene Expression"]], 
  assay = "Xenium", project = "Ctrl3")
xenium.obj.ctrl3[["BlankCodeword"]] <- CreateAssayObject(counts = data.ctrl3$matrix[["Unassigned Codeword"]])
xenium.obj.ctrl3[["ControlCodeword"]] <- CreateAssayObject(counts = data.ctrl3$matrix[["Negative Control Codeword"]])
xenium.obj.ctrl3[["ControlProbe"]] <- CreateAssayObject(counts = data.ctrl3$matrix[["Negative Control Probe"]])
xenium.obj.ctrl3[["fov"]] <- coords.ctrl3
xenium.obj.ctrl3$Condition <- 'Ctrl'
xenium.obj.ctrl3 <- subset(xenium.obj.ctrl3, subset = nCount_Xenium > 0)
# xenium.obj.ctrl3 <- subset(xenium.obj.ctrl3,cells = CellID_ctrl3_Slc5a7$Barcode)

print(xenium.obj.ctrl3)


# 
segmentations.data.cold1 <- list(
  centroids = CreateCentroids(data.cold1$centroids),
  segmentation = CreateSegmentation(data.cold1$segmentations))
coords.cold1 <- CreateFOV(
  coords = segmentations.data.cold1, 
  type = c("segmentation", "centroids"), 
  molecules = data.cold1$microns, 
  assay = "Xenium")
xenium.obj.cold1 <- CreateSeuratObject(
  counts = data.cold1$matrix[["Gene Expression"]], 
  assay = "Xenium", project = "Cold1")
xenium.obj.cold1[["BlankCodeword"]] <- CreateAssayObject(counts = data.cold1$matrix[["Unassigned Codeword"]])
xenium.obj.cold1[["ControlCodeword"]] <- CreateAssayObject(counts = data.cold1$matrix[["Negative Control Codeword"]])
xenium.obj.cold1[["ControlProbe"]] <- CreateAssayObject(counts = data.cold1$matrix[["Negative Control Probe"]])
xenium.obj.cold1[["fov"]] <- coords.cold1
xenium.obj.cold1$Condition <- 'Cold'
xenium.obj.cold1 <- subset(xenium.obj.cold1, subset = nCount_Xenium > 0)
# xenium.obj.cold1 <- subset(xenium.obj.cold1,cells = CellID_cold1_Slc5a7$Barcode)

print(xenium.obj.cold1)

# 
segmentations.data.cold2 <- list(
  centroids = CreateCentroids(data.cold2$centroids),
  segmentation = CreateSegmentation(data.cold2$segmentations))
coords.cold2 <- CreateFOV(
  coords = segmentations.data.cold2, 
  type = c("segmentation", "centroids"), 
  molecules = data.cold2$microns, 
  assay = "Xenium")
xenium.obj.cold2 <- CreateSeuratObject(
  counts = data.cold2$matrix[["Gene Expression"]], 
  assay = "Xenium", project = "Cold2")
xenium.obj.cold2[["BlankCodeword"]] <- CreateAssayObject(counts = data.cold2$matrix[["Unassigned Codeword"]])
xenium.obj.cold2[["ControlCodeword"]] <- CreateAssayObject(counts = data.cold2$matrix[["Negative Control Codeword"]])
xenium.obj.cold2[["ControlProbe"]] <- CreateAssayObject(counts = data.cold2$matrix[["Negative Control Probe"]])
xenium.obj.cold2[["fov"]] <- coords.cold2
xenium.obj.cold2$Condition <- 'Cold'
xenium.obj.cold2 <- subset(xenium.obj.cold2, subset = nCount_Xenium > 0)
# xenium.obj.cold2 <- subset(xenium.obj.cold2,cells = CellID_cold2_Slc5a7$Barcode)

print(xenium.obj.cold2)


# 
segmentations.data.cold3 <- list(
  centroids = CreateCentroids(data.cold3$centroids),
  segmentation = CreateSegmentation(data.cold3$segmentations))
coords.cold3 <- CreateFOV(
  coords = segmentations.data.cold3, 
  type = c("segmentation", "centroids"), 
  molecules = data.cold3$microns, 
  assay = "Xenium")
xenium.obj.cold3 <- CreateSeuratObject(
  counts = data.cold3$matrix[["Gene Expression"]], 
  assay = "Xenium", project = "Cold3")
xenium.obj.cold3[["BlankCodeword"]] <- CreateAssayObject(counts = data.cold3$matrix[["Unassigned Codeword"]])
xenium.obj.cold3[["ControlCodeword"]] <- CreateAssayObject(counts = data.cold3$matrix[["Negative Control Codeword"]])
xenium.obj.cold3[["ControlProbe"]] <- CreateAssayObject(counts = data.cold3$matrix[["Negative Control Probe"]])
xenium.obj.cold3[["fov"]] <- coords.cold3
xenium.obj.cold3$Condition <- 'Cold'
xenium.obj.cold3 <- subset(xenium.obj.cold3, subset = nCount_Xenium > 0)
# xenium.obj.cold3 <- subset(xenium.obj.cold3, cells = CellID_cold3_Slc5a7$Barcode)
print(xenium.obj.cold3)


#
segmentations.data.ctrlF <- list(
  centroids = CreateCentroids(data.ctrlF$centroids),
  segmentation = CreateSegmentation(data.ctrlF$segmentations))
coords.ctrlF <- CreateFOV(
  coords = segmentations.data.ctrlF, 
  type = c("segmentation", "centroids"), 
  molecules = data.ctrlF$microns, 
  assay = "Xenium")
xenium.obj.ctrlF <- CreateSeuratObject(
  counts = data.ctrlF$matrix[["Gene Expression"]], 
  assay = "Xenium", project = "ctrlF")
xenium.obj.ctrlF[["BlankCodeword"]] <- CreateAssayObject(counts = data.ctrlF$matrix[["Unassigned Codeword"]])
xenium.obj.ctrlF[["ControlCodeword"]] <- CreateAssayObject(counts = data.ctrlF$matrix[["Negative Control Codeword"]])
xenium.obj.ctrlF[["ControlProbe"]] <- CreateAssayObject(counts = data.ctrlF$matrix[["Negative Control Probe"]])
xenium.obj.ctrlF[["fov"]] <- coords.ctrlF
xenium.obj.ctrlF$Condition <- 'Ctrl'
xenium.obj.ctrlF <- subset(xenium.obj.ctrlF, subset = nCount_Xenium > 0)
# xenium.obj.QIH <- subset(xenium.obj.QIH, cells = CellID_QIH_Slc5a7$Barcode)
print(xenium.obj.ctrlF)


segmentations.data.coldF <- list(
  centroids = CreateCentroids(data.coldF$centroids),
  segmentation = CreateSegmentation(data.coldF$segmentations))
coords.coldF <- CreateFOV(
  coords = segmentations.data.coldF, 
  type = c("segmentation", "centroids"), 
  molecules = data.coldF$microns, 
  assay = "Xenium")
xenium.obj.coldF <- CreateSeuratObject(
  counts = data.coldF$matrix[["Gene Expression"]], 
  assay = "Xenium", project = "coldF")
xenium.obj.coldF[["BlankCodeword"]] <- CreateAssayObject(counts = data.coldF$matrix[["Unassigned Codeword"]])
xenium.obj.coldF[["ControlCodeword"]] <- CreateAssayObject(counts = data.coldF$matrix[["Negative Control Codeword"]])
xenium.obj.coldF[["ControlProbe"]] <- CreateAssayObject(counts = data.coldF$matrix[["Negative Control Probe"]])
xenium.obj.coldF[["fov"]] <- coords.coldF
xenium.obj.coldF$Condition <- 'Cold'
xenium.obj.coldF <- subset(xenium.obj.coldF, subset = nCount_Xenium > 0)
# xenium.obj.QIH <- subset(xenium.obj.QIH, cells = CellID_QIH_Slc5a7$Barcode)
print(xenium.obj.coldF)



segmentations.data.2DGF <- list(
  centroids = CreateCentroids(data.2DGF$centroids),
  segmentation = CreateSegmentation(data.2DGF$segmentations))
coords.2DGF <- CreateFOV(
  coords = segmentations.data.2DGF, 
  type = c("segmentation", "centroids"), 
  molecules = data.2DGF$microns, 
  assay = "Xenium")
xenium.obj.2DGF <- CreateSeuratObject(
  counts = data.2DGF$matrix[["Gene Expression"]], 
  assay = "Xenium", project = "2DGF")
xenium.obj.2DGF[["BlankCodeword"]] <- CreateAssayObject(counts = data.2DGF$matrix[["Unassigned Codeword"]])
xenium.obj.2DGF[["ControlCodeword"]] <- CreateAssayObject(counts = data.2DGF$matrix[["Negative Control Codeword"]])
xenium.obj.2DGF[["ControlProbe"]] <- CreateAssayObject(counts = data.2DGF$matrix[["Negative Control Probe"]])
xenium.obj.2DGF[["fov"]] <- coords.2DGF
xenium.obj.2DGF$Condition <- '2DG'
xenium.obj.2DGF <- subset(xenium.obj.2DGF, subset = nCount_Xenium > 0)
# xenium.obj.QIH <- subset(xenium.obj.QIH, cells = CellID_QIH_Slc5a7$Barcode)
print(xenium.obj.2DGF)


segmentations.data.2DGM <- list(
  centroids = CreateCentroids(data.2DGM$centroids),
  segmentation = CreateSegmentation(data.2DGM$segmentations))
coords.2DGM <- CreateFOV(
  coords = segmentations.data.2DGM, 
  type = c("segmentation", "centroids"), 
  molecules = data.2DGM$microns, 
  assay = "Xenium")
xenium.obj.2DGM <- CreateSeuratObject(
  counts = data.2DGM$matrix[["Gene Expression"]], 
  assay = "Xenium", project = "2DGM")
xenium.obj.2DGM[["BlankCodeword"]] <- CreateAssayObject(counts = data.2DGM$matrix[["Unassigned Codeword"]])
xenium.obj.2DGM[["ControlCodeword"]] <- CreateAssayObject(counts = data.2DGM$matrix[["Negative Control Codeword"]])
xenium.obj.2DGM[["ControlProbe"]] <- CreateAssayObject(counts = data.2DGM$matrix[["Negative Control Probe"]])
xenium.obj.2DGM[["fov"]] <- coords.2DGM
xenium.obj.2DGM$Condition <- '2DG'
xenium.obj.2DGM <- subset(xenium.obj.2DGM, subset = nCount_Xenium > 0)
# xenium.obj.QIH <- subset(xenium.obj.QIH, cells = CellID_QIH_Slc5a7$Barcode)
print(xenium.obj.2DGM)


xenium.obj.ctrlM <- merge(x = xenium.obj.ctrl1, y = c(xenium.obj.ctrl2,xenium.obj.ctrl3), add.cell.ids = c("Ctrl1","Ctrl2","Ctrl3"))
xenium.obj.coldM <- merge(x = xenium.obj.cold1, y = c(xenium.obj.cold2,xenium.obj.cold3), add.cell.ids = c("Cold1","Cold2","Cold3"))
saveRDS(xenium.obj.ctrlM,"xenium.obj.ctrlM.rds")
saveRDS(xenium.obj.coldM,"xenium.obj.coldM.rds")

# xenium.obj.ctrlM
xenium.obj.ctrlM <- SCTransform(xenium.obj.ctrlM, assay = "Xenium")
xenium.obj.ctrlM <- RunPCA(xenium.obj.ctrlM, npcs = 50, features = rownames(xenium.obj.ctrlM))
xenium.obj.ctrlM <- RunUMAP(xenium.obj.ctrlM, dims = 1:50)
xenium.obj.ctrlM <- FindNeighbors(xenium.obj.ctrlM, dims = 1:50)
xenium.obj.ctrlM <- FindClusters(xenium.obj.ctrlM, resolution =1.0)
saveRDS(xenium.obj.ctrlM,"xenium.obj.ctrlM.rds")

# xenium.obj.coldM
xenium.obj.coldM <- SCTransform(xenium.obj.coldM, assay = "Xenium")
xenium.obj.coldM <- RunPCA(xenium.obj.coldM, npcs = 50, features = rownames(xenium.obj.coldM))
# merged.obj <- IntegrateLayers(object = merged.obj, method = CCAIntegration, normalization.method = "SCT", orig.reduction = "pca", new.reduction = "integrated.dr",
# verbose = FALSE)
xenium.obj.coldM <- RunUMAP(xenium.obj.coldM, dims = 1:50)
xenium.obj.coldM <- FindNeighbors(xenium.obj.coldM, dims = 1:50)
xenium.obj.coldM <- FindClusters(xenium.obj.coldM, resolution =1.0)
xenium.obj.coldM@active.ident <- xenium.obj.coldM@meta.data$seurat_clusters
saveRDS(xenium.obj.coldM,"xenium.obj.coldM.rds")

## xenium.obj.ctrlF
options(future.globals.maxSize = 256 * 1024 ^ 3) 

xenium.obj.ctrlF <- SCTransform(xenium.obj.ctrlF, assay = "Xenium")
xenium.obj.ctrlF <- RunPCA(xenium.obj.ctrlF, npcs = 50, features = rownames(xenium.obj.ctrlF))

xenium.obj.ctrlF <- RunUMAP(xenium.obj.ctrlF, dims = 1:50)
xenium.obj.ctrlF <- FindNeighbors(xenium.obj.ctrlF, dims = 1:50)
xenium.obj.ctrlF <- FindClusters(xenium.obj.ctrlF, resolution =1.0)
xenium.obj.ctrlF@active.ident <- xenium.obj.ctrlF@meta.data$seurat_clusters
saveRDS(xenium.obj.ctrlF,"xenium.obj.ctrlF.rds")


## xenium.obj.coldF
options(future.globals.maxSize = 256 * 1024 ^ 3) 

xenium.obj.coldF <- SCTransform(xenium.obj.coldF, assay = "Xenium")
xenium.obj.coldF <- RunPCA(xenium.obj.coldF, npcs = 50, features = rownames(xenium.obj.coldF))
xenium.obj.coldF <- RunUMAP(xenium.obj.coldF, dims = 1:50)
xenium.obj.coldF <- FindNeighbors(xenium.obj.coldF, dims = 1:50)
xenium.obj.coldF <- FindClusters(xenium.obj.coldF, resolution =1.0)
xenium.obj.coldF@active.ident <- xenium.obj.coldF@meta.data$seurat_clusters
saveRDS(xenium.obj.coldF,"xenium.obj.coldF.rds")


# xenium.obj.2DGM

xenium.obj.2DGM <- SCTransform(xenium.obj.2DGM, assay = "Xenium")
xenium.obj.2DGM <- RunPCA(xenium.obj.2DGM, npcs = 50, features = rownames(xenium.obj.2DGM))
# merged.obj <- IntegrateLayers(object = merged.obj, method = CCAIntegration, normalization.method = "SCT", orig.reduction = "pca", new.reduction = "integrated.dr",
# verbose = FALSE)
xenium.obj.2DGM <- RunUMAP(xenium.obj.2DGM, dims = 1:50)
xenium.obj.2DGM <- FindNeighbors(xenium.obj.2DGM, dims = 1:50)
xenium.obj.2DGM <- FindClusters(xenium.obj.2DGM, resolution =1.0)
xenium.obj.2DGM@active.ident <- xenium.obj.2DGM@meta.data$seurat_clusters
saveRDS(xenium.obj.2DGM,"xenium.obj.2DGM.rds")


#xenium.obj.2DGF
xenium.obj.2DGF <- SCTransform(xenium.obj.2DGF, assay = "Xenium")
xenium.obj.2DGF <- RunPCA(xenium.obj.2DGF, npcs = 50, features = rownames(xenium.obj.2DGF))
# merged.obj <- IntegrateLayers(object = merged.obj, method = CCAIntegration, normalization.method = "SCT", orig.reduction = "pca", new.reduction = "integrated.dr",
# verbose = FALSE)
xenium.obj.2DGF <- RunUMAP(xenium.obj.2DGF, dims = 1:50)
xenium.obj.2DGF <- FindNeighbors(xenium.obj.2DGF, dims = 1:50)
xenium.obj.2DGF <- FindClusters(xenium.obj.2DGF, resolution =1.0)
xenium.obj.2DGF@active.ident <- xenium.obj.2DGF@meta.data$seurat_clusters
saveRDS(xenium.obj.2DGF,"xenium.obj.2DGF.rds")



# subset cholinergic neurons
xenium.obj.ctrlM_Choline <- subset(xenium.obj.ctrlM, idents = c(17))
xenium.obj.coldM_Choline <- subset(xenium.obj.coldM, idents = c(19))
xenium.obj.ctrlF_Choline <- subset(xenium.obj.ctrlF, idents = c(28,50))
xenium.obj.coldF_Choline <- subset(xenium.obj.coldF, idents = c(22))
xenium.obj.2DGM_Choline <- subset(xenium.obj.2DGM, idents = c(28,36))
xenium.obj.2DGF_Choline <- subset(xenium.obj.2DGF, idents = c(31,50))






# Cholinergic neurons
merged.obj_Choline <- merge(x = xenium.obj.ctrlM_Choline, y = c(xenium.obj.coldM_Choline), add.cell.ids = c("CtrlM","ColdM"))
merged.obj_Choline2 <- merge(x = xenium.obj.ctrlF_Choline, y = c(xenium.obj.coldF_Choline,xenium.obj.2DGM_Choline,xenium.obj.2DGF_Choline), add.cell.ids = c("CtrlF","ColdF","2DGM","2DGF"))
merged.obj_Choline_woQIHLPS <- merge(x = xenium.obj.ctrlM_Choline, y = c(merged.obj_Choline2))
# xenium.SCT.obj <- readRDS("xenium.SCT.woFOV.obj")
merged.obj_Choline_woQIHLPS <- SCTransform(merged.obj_Choline_woQIHLPS, assay = "Xenium")
merged.obj_Choline_woQIHLPS <- RunPCA(merged.obj_Choline_woQIHLPS, npcs = 50, features = rownames(merged.obj_Choline_woQIHLPS))
merged.obj_Choline_woQIHLPS <- IntegrateLayers(object = merged.obj_Choline_woQIHLPS, method = CCAIntegration, normalization.method = "SCT", orig.reduction = "pca", new.reduction = "integrated.dr",verbose = FALSE)
merged.obj_Choline_woQIHLPS <- RunUMAP(merged.obj_Choline_woQIHLPS, reduction = "integrated.dr", dims = 1:50)
merged.obj_Choline_woQIHLPS <- FindNeighbors(merged.obj_Choline_woQIHLPS, reduction = "integrated.dr", dims = 1:50)
merged.obj_Choline_woQIHLPS <- FindClusters(merged.obj_Choline_woQIHLPS, resolution =1.0)
# merged.obj_Choline_all@active.ident <- merged.obj_Choline_all@meta.data$seurat_clusters
saveRDS(merged.obj_Choline_woQIHLPS,"merged.obj_Choline_woQIHLPS.rds")

# metadata: Sex
male_ids <- c("2DGM", "Cold1", "Cold2", "Cold3", "Ctrl1", "Ctrl2", "Ctrl3")
female_ids <- c("2DGF", "ctrlF","coldF")

merged.obj_Choline_woQIHLPS$Sex <- NA

merged.obj_Choline_woQIHLPS$Sex[
  merged.obj_Choline_woQIHLPS$orig.ident %in% male_ids
] <- "Male"

merged.obj_Choline_woQIHLPS$Sex[
  merged.obj_Choline_woQIHLPS$orig.ident %in% female_ids
] <- "Female"

# QC
VlnPlot(merged.obj_Choline_woQIHLPS, features = c("Camk2a", "Rbfox3", "Slc17a6", "Slc17a8", "Gad1", "Gad2","Opalin","Pdgfra","Cd68","Dcn","Aqp4","Pecam1","Sox10","Spag16","Chat","Slc5a7","Slc6a5","Nos1","Pax2","Bcl6"), 
        pt.size = 0, stack = T, flip = T)


# Cholinergic Neurons afterQC
merged.obj_Choline_woQIHLPS_afterQC <- subset(merged.obj_Choline_woQIHLPS,idents = c(0,3,4,5,6,7,8,9,11,16,24,26,28,30),
                                              invert = TRUE)
DimPlot(merged.obj_Choline_woQIHLPS_afterQC)

merged.obj_Choline_woQIHLPS_afterQC <- SCTransform(merged.obj_Choline_woQIHLPS_afterQC, assay = "Xenium")
merged.obj_Choline_woQIHLPS_afterQC <- RunPCA(merged.obj_Choline_woQIHLPS_afterQC, npcs = 50, features = rownames(merged.obj_Choline_woQIHLPS_afterQC))
merged.obj_Choline_woQIHLPS_afterQC <- IntegrateLayers(object = merged.obj_Choline_woQIHLPS_afterQC, method = CCAIntegration, normalization.method = "SCT", orig.reduction = "pca", new.reduction = "integrated.dr", verbose = FALSE)
merged.obj_Choline_woQIHLPS_afterQC <- RunUMAP(merged.obj_Choline_woQIHLPS_afterQC, reduction = "integrated.dr", dims = 1:50)
merged.obj_Choline_woQIHLPS_afterQC <- FindNeighbors(merged.obj_Choline_woQIHLPS_afterQC, reduction = "integrated.dr", dims = 1:50)
merged.obj_Choline_woQIHLPS_afterQC <- FindClusters(merged.obj_Choline_woQIHLPS_afterQC, resolution =1.0)
saveRDS(merged.obj_Choline_woQIHLPS_afterQC,"merged.obj_Choline_woQIHLPS_afterQC.rds")


### labeling
plots <- VlnPlot(merged.obj_Choline_woQIHLPS_afterQC, features = c("Camk2a", "Rbfox3", "Slc17a6", "Slc17a8", "Gad1", "Gad2","Opalin","Pdgfra","Cd68","Dcn","Aqp4","Pecam1","Sox10","Spag16","Chat","Slc5a7","Slc6a5","Nos1","Pax2","Bcl6"), 
                 pt.size = 0, stack = T, flip = T)
wrap_plots(plots = plots, ncol = 1)


SPN_clusters <- c(6,8,11,14,16,17,19,23,26,28,29,31,34,35)
IN_clusters  <- c(24,36,25,18,32,20,5,33)
MN_clusters  <- c(0,1,2,4,7,9,10,12,13,15,21,22,30,38,39,40)
Mixed_clusters <- c(3,27,37)

# active.ident を取得
cluster_ids <- as.numeric(as.character(Idents(merged.obj_Choline_woQIHLPS_afterQC)))

# 初期化
merged.obj_Choline_woQIHLPS_afterQC$Cholinergic_subtype <- NA

# ラベリング
merged.obj_Choline_woQIHLPS_afterQC$Cholinergic_subtype[
  cluster_ids %in% SPN_clusters
] <- "SPN"

merged.obj_Choline_woQIHLPS_afterQC$Cholinergic_subtype[
  cluster_ids %in% IN_clusters
] <- "IN"

merged.obj_Choline_woQIHLPS_afterQC$Cholinergic_subtype[
  cluster_ids %in% MN_clusters
] <- "MN"

merged.obj_Choline_woQIHLPS_afterQC$Cholinergic_subtype[
  cluster_ids %in% Mixed_clusters
] <- "Mixed"

# 確認
table(
  Idents(merged.obj_Choline_woQIHLPS_afterQC),
  merged.obj_Choline_woQIHLPS_afterQC$Cholinergic_subtype
)


# subset for Display (Fig1c)
merged.obj_SPN_woQIHLPS_display <- subset(merged.obj_Choline_woQIHLPS_afterQC,subset = Cholinergic_subtype != "Mixed")
DimPlot(merged.obj_SPN_woQIHLPS_display, label=F, group.by="Cholinergic_subtype")


### subset SPNs

merged.obj_SPN_woQIHLPS <- subset(merged.obj_Choline_woQIHLPS_afterQC,idents = c(6,8,11,14,16,17,19,23,26,28,29,31,34,35))
DimPlot(merged.obj_SPN_woQIHLPS)

merged.obj_SPN_woQIHLPS <- SCTransform(merged.obj_SPN_woQIHLPS, assay = "Xenium")
merged.obj_SPN_woQIHLPS <- RunPCA(merged.obj_SPN_woQIHLPS, npcs = 50, features = rownames(merged.obj_SPN_woQIHLPS))
merged.obj_SPN_woQIHLPS <- IntegrateLayers(object = merged.obj_SPN_woQIHLPS, method = CCAIntegration, normalization.method = "SCT", orig.reduction = "pca", new.reduction = "integrated.dr", verbose = FALSE)
merged.obj_SPN_woQIHLPS <- RunUMAP(merged.obj_SPN_woQIHLPS, reduction = "integrated.dr", dims = 1:50)
merged.obj_SPN_woQIHLPS <- FindNeighbors(merged.obj_SPN_woQIHLPS, reduction = "integrated.dr", dims = 1:50)
merged.obj_SPN_woQIHLPS <- FindClusters(merged.obj_SPN_woQIHLPS, resolution =1.2)
saveRDS(merged.obj_SPN_woQIHLPS,"merged.obj_SPN_woQIHLPS.rds")

# QC
plots <- VlnPlot(merged.obj_SPN_woQIHLPS, features = c("Camk2a", "Rbfox3", "Slc17a6", "Slc17a8", "Gad1", "Gad2","Opalin","Pdgfra","Cd68","Dcn","Aqp4","Pecam1","Sox10","Spag16","Chat","Slc5a7","Slc6a5","Nos1","Pax2","Bcl6"), 
                 pt.size = 0, stack = T, flip = T)
wrap_plots(plots = plots, ncol = 1)

# subset SPNs qualified for 1st QC
merged.obj_SPN_woQIHLPS_afterQC <- subset(merged.obj_SPN_woQIHLPS,idents = c(17,22,26), invert = T)
DimPlot(merged.obj_SPN_woQIHLPS_afterQC)
merged.obj_SPN_woQIHLPS_afterQC[["Xenium"]] <- split(merged.obj_SPN_woQIHLPS_afterQC[["Xenium"]], f = merged.obj_SPN_woQIHLPS_afterQC$orig.ident)

# xenium.SCT.obj <- readRDS("xenium.SCT.woFOV.obj")
merged.obj_SPN_woQIHLPS_afterQC <- SCTransform(merged.obj_SPN_woQIHLPS_afterQC, assay = "Xenium")
merged.obj_SPN_woQIHLPS_afterQC <- RunPCA(merged.obj_SPN_woQIHLPS_afterQC, npcs = 50, features = rownames(merged.obj_SPN_woQIHLPS_afterQC))
merged.obj_SPN_woQIHLPS_afterQC <- IntegrateLayers(object = merged.obj_SPN_woQIHLPS_afterQC, method = CCAIntegration, normalization.method = "SCT", orig.reduction = "pca", new.reduction = "integrated.dr", verbose = FALSE)
merged.obj_SPN_woQIHLPS_afterQC <- RunUMAP(merged.obj_SPN_woQIHLPS_afterQC, reduction = "integrated.dr", dims = 1:50)
merged.obj_SPN_woQIHLPS_afterQC <- FindNeighbors(merged.obj_SPN_woQIHLPS_afterQC, reduction = "integrated.dr", dims = 1:50)
merged.obj_SPN_woQIHLPS_afterQC <- FindClusters(merged.obj_SPN_woQIHLPS_afterQC, resolution =1.2)
# merged.obj_Choline_all_afterQC@active.ident <- merged.obj_Choline_all_afterQC@meta.data$seurat_clusters
saveRDS(merged.obj_SPN_woQIHLPS_afterQC,"merged.obj_SPN_woQIHLPS_afterQC.rds")

# QC
plots <- VlnPlot(merged.obj_SPN_woQIHLPS_afterQC, features = c("Camk2a", "Rbfox3", "Slc17a6", "Slc17a8", "Gad1", "Gad2","Opalin","Pdgfra","Cd68","Dcn","Aqp4","Pecam1","Sox10","Spag16","Chat","Slc5a7","Slc6a5","Nos1","Pax2","Bcl6"), 
                 pt.size = 0, stack = T, flip = T)
wrap_plots(plots = plots, ncol = 1)


### 2nd round of QC of SPN
merged.obj_SPN_woQIHLPS_afterQC2 <- subset(merged.obj_SPN_woQIHLPS_afterQC,idents = c(23), invert = T)# Opalinが存在
DimPlot(merged.obj_SPN_woQIHLPS_afterQC2)
merged.obj_SPN_woQIHLPS_afterQC2[["Xenium"]] <- split(merged.obj_SPN_woQIHLPS_afterQC2[["Xenium"]], f = merged.obj_SPN_woQIHLPS_afterQC2$orig.ident)

merged.obj_SPN_woQIHLPS_afterQC2 <- SCTransform(merged.obj_SPN_woQIHLPS_afterQC2, assay = "Xenium")
merged.obj_SPN_woQIHLPS_afterQC2 <- RunPCA(merged.obj_SPN_woQIHLPS_afterQC2, npcs = 50, features = rownames(merged.obj_SPN_woQIHLPS_afterQC2))
merged.obj_SPN_woQIHLPS_afterQC2 <- IntegrateLayers(object = merged.obj_SPN_woQIHLPS_afterQC2, method = CCAIntegration, normalization.method = "SCT", orig.reduction = "pca", new.reduction = "integrated.dr", verbose = FALSE)
merged.obj_SPN_woQIHLPS_afterQC2 <- RunUMAP(merged.obj_SPN_woQIHLPS_afterQC2, reduction = "integrated.dr", dims = 1:50)
merged.obj_SPN_woQIHLPS_afterQC2 <- FindNeighbors(merged.obj_SPN_woQIHLPS_afterQC2, reduction = "integrated.dr", dims = 1:50)
merged.obj_SPN_woQIHLPS_afterQC2 <- FindClusters(merged.obj_SPN_woQIHLPS_afterQC2, resolution =1.2)
saveRDS(merged.obj_SPN_woQIHLPS_afterQC2,"merged.obj_SPN_woQIHLPS_afterQC2.rds")

# QC
plots <- VlnPlot(merged.obj_SPN_woQIHLPS_afterQC2, features = c("Camk2a", "Rbfox3", "Slc17a6", "Slc17a8", "Gad1", "Gad2","Opalin","Pdgfra","Cd68","Dcn","Aqp4","Pecam1","Sox10","Spag16","Chat","Slc5a7","Slc6a5","Nos1","Pax2","Bcl6"), 
                 pt.size = 0, stack = T, flip = T)
wrap_plots(plots = plots, ncol = 1)


merged.obj_SPN_woQIHLPS_afterQC3 <- subset(merged.obj_SPN_woQIHLPS_afterQC2,idents = c(26), invert = T)
DimPlot(merged.obj_SPN_woQIHLPS_afterQC3)
merged.obj_SPN_woQIHLPS_afterQC3 <- SCTransform(merged.obj_SPN_woQIHLPS_afterQC3, assay = "Xenium")
merged.obj_SPN_woQIHLPS_afterQC3 <- RunPCA(merged.obj_SPN_woQIHLPS_afterQC3, npcs = 50, features = rownames(merged.obj_SPN_woQIHLPS_afterQC3))
merged.obj_SPN_woQIHLPS_afterQC3 <- IntegrateLayers(object = merged.obj_SPN_woQIHLPS_afterQC3, method = CCAIntegration, normalization.method = "SCT", orig.reduction = "pca", new.reduction = "integrated.dr", verbose = FALSE)
merged.obj_SPN_woQIHLPS_afterQC3 <- RunUMAP(merged.obj_SPN_woQIHLPS_afterQC3, reduction = "integrated.dr", dims = 1:50)
merged.obj_SPN_woQIHLPS_afterQC3 <- FindNeighbors(merged.obj_SPN_woQIHLPS_afterQC3, reduction = "integrated.dr", dims = 1:50)
merged.obj_SPN_woQIHLPS_afterQC3 <- FindClusters(merged.obj_SPN_woQIHLPS_afterQC3, resolution =1.2)
# merged.obj_Choline_all_afterQC@active.ident <- merged.obj_Choline_all_afterQC@meta.data$seurat_clusters
saveRDS(merged.obj_SPN_woQIHLPS_afterQC3,"merged.obj_SPN_woQIHLPS_afterQC3.rds")

# Fig1c
DimPlot(merged.obj_SPN_woQIHLPS_afterQC3, label =T)

# Fig1e
markers.to.plot =  c("Chat","Nos1","Ccdc3", "Npffr2", "Trpc3", "Gfra2", "Ccbe1", "Crh", "Nts", "Palld",
                     "Il1r1", "Mamdc2", "Penk", "Phldb2", "Sst", "Calb2", "Reln")
DotPlot(
  merged.obj_SPN_woQIHLPS_afterQC3,
  features = markers.to.plot,
  cols = c("white","tomato2")
) +
  coord_flip()


# Fig1f
prop_cluster <- prop.table(x=table(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$segments,merged.obj_SPN_woQIHLPS_afterQC3@active.ident), margin =2)
# write.csv(prop_cluster,"prop_cluster.csv")
Cluster_name = factor(colnames(prop_cluster), levels=colnames(prop_cluster))    # To order as I want, making a factor object for the cluster
prop_cluster_df <- data.frame(Cluster=Cluster_name,S2_4=prop_cluster[8,],L6S1=prop_cluster[7,],L45=prop_cluster[6,], L1_3=prop_cluster[5,],T11_13=prop_cluster[4,], T8_10=prop_cluster[3,], T5_7=prop_cluster[2,],T1_4=prop_cluster[1,])
prop_cluster_df2 <- tidyr::gather(prop_cluster_df, key=Sample, value=Proportion, -Cluster, factor_key = TRUE)
# ヒートマップ描画
ggplot(prop_cluster_df2, aes(x = Cluster, y = Sample, fill = Proportion)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "tomato2") +
  labs(x = "Clusters", y = "Segments", fill = "Proportion") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

# Fig1g
DimPlot(merged.obj_SPN_woQIHLPS_afterQC3, label =T, group.by="population")


#  Fig1j
DimPlot(
  merged.obj_SPN_woQIHLPS_afterQC3,
  group.by = "Sex",
  label = TRUE,
  shuffle = TRUE,
  pt.size = 0.3,
  cols = c(
    "Female" = "hotpink2",
    "Male"   = "dodgerblue2"
  )
)


# ED Fig1c
DimPlot(merged.obj_SPN_woQIHLPS_afterQC3, label =T, split.by = "ConditionSex")

# ED Fig1d
library(tibble)
#Conditionごとの各細胞腫の割合
df<-table(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$SPN_subtype,merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex)
df <- as.data.frame(df) 
# データを整形（長い形式に変換）
df_long <- df %>%
  rownames_to_column("SPN_subtype") 
# 各Condition内で割合を計算
df_long <- df_long %>%
  group_by(Var2) %>%
  mutate(Fraction = Freq / sum(Freq))
df_long


# 割合の棒グラフを描画
ggplot(df_long, aes(x = Var2, y = Fraction, fill = Var1)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "各Conditionにおける細胞種の割合",
    x = "Condition",
    y = "割合"
  ) +
  theme_minimal()

# ED Fig1e
DimPlot(merged.obj_SPN_woQIHLPS_afterQC3, label =T, split.by = "predicted.celltype")


# ED Fig1f
library(ggplot2)
library(ggalluvial)
meta_df <- merged.obj_SPN_woQIHLPS_afterQC3@meta.data[, c("predicted.celltype", "SPN_subtype")]

# factor → character（安全のため）
meta_df$predicted.celltype <- as.character(meta_df$predicted.celltype)
meta_df$SPN_subtype <- as.character(meta_df$SPN_subtype)

# NA除去
meta_df <- na.omit(meta_df)

# 組み合わせごとに頻度集計
plot_df <- meta_df %>%
  group_by(predicted.celltype, SPN_subtype) %>%
  summarise(value = n(), .groups = "drop")
plot_df_long <- plot_df %>%
  tidyr::pivot_longer(cols = c(predicted.celltype, SPN_subtype),
                      names_to = "axis",
                      values_to = "category")
ggplot(
  links,
  aes(
    axis1 = predicted.celltype,
    axis2 = SPN_subtype,
    y = value
  )
) +
  geom_alluvium(
    aes(fill = predicted.celltype),
    width = 1/12,
    alpha = 1        #  ←これで不透明
  ) +
  geom_stratum(
    width = 1/12,
    fill = "grey80",
    color = "black",
    alpha = 1        #  ← stratums も不透明に
  ) +
  geom_label(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    alpha = 1        #  ← ラベル背景も不透明
  ) +
  scale_x_discrete(limits = c("Predicted", "Subtype"), expand = c(.1, .1)) +
  theme_minimal()

# ED Fig1g


# ED Fig1h
FeaturePlot(merged.obj_SPN_woQIHLPS_afterQC3, features = c("Cartpt"), col=c("aquamarine", "tomato2"),pt.size=0.5)





# Fig2b
FeaturePlot(merged.obj_SPN_woQIHLPS_afterQC3, features = c("Fos"), col=c("aquamarine", "tomato2"),pt.size=0.5, split.by = "Condition")

# Fig2c

FeaturePlot(merged.obj_SPN_woQIHLPS_afterQC3, features = c("Nos1"), col=c("aquamarine", "tomato2"),pt.size=0.5)
FeaturePlot(merged.obj_SPN_woQIHLPS_afterQC3, features = c("Nts"), col=c("aquamarine", "tomato2"),pt.size=0.5)
FeaturePlot(merged.obj_SPN_woQIHLPS_afterQC3, features = c("Crh"), col=c("aquamarine", "tomato2"),pt.size=0.5)
FeaturePlot(merged.obj_SPN_woQIHLPS_afterQC3, features = c("Palld"), col=c("aquamarine", "tomato2"),pt.size=0.5)


# Fig2d_Fos
library(dplyr)
library(ggplot2)

# metadata取得
meta_df <- merged.obj_SPN_woQIHLPS_afterQC3@meta.data

# cluster情報追加
meta_df$Cluster <- Idents(merged.obj_SPN_woQIHLPS_afterQC3)

# Fos発現量取得
DefaultAssay(merged.obj_SPN_woQIHLPS_afterQC3)<-"SCT"
meta_df$Fos <- FetchData(
  merged.obj_SPN_woQIHLPS_afterQC3,
  vars = "Fos"
)[,1]

# 集計
dot_df <- meta_df %>%
  group_by(Cluster, Condition) %>%
  summarise(
    Mean_Fos = mean(Fos, na.rm = TRUE),
    Fos_Positive = mean(Fos > 1.590272, na.rm = TRUE),
    .groups = "drop"
  )

dot_df$Condition <- factor(
  dot_df$Condition,
  levels = c("Ctrl", "Cold","2DG")
)
# DotPlot風描画
ggplot(
  dot_df,
  aes(
    x = Cluster,
    y = Condition
  )
) +
  scale_size(
    range = c(2,12)
  ) +
  geom_point(
    aes(
      size = Fos_Positive,
      color = Mean_Fos
    )
  ) +
  scale_color_gradient(
    low = "aquamarine",
    high = "tomato"
  ) +
  theme_bw() +
  labs(
    size = "Fos+ proportion",
    color = "Mean Fos"
  )


# Fig2d_IEG
library(dplyr)
library(ggplot2)

# metadata取得
meta_df <- merged.obj_SPN_woQIHLPS_afterQC3@meta.data

# Cluster情報追加
meta_df$Cluster <- Idents(merged.obj_SPN_woQIHLPS_afterQC3)

# IEG_score1取得
meta_df$IEG_score1 <- FetchData(
  merged.obj_SPN_woQIHLPS_afterQC3,
  vars = "IEG_score1"
)[,1]

# 集計
dot_df <- meta_df %>%
  group_by(Cluster, Condition) %>%
  summarise(
    Mean_IEG = mean(IEG_score1, na.rm = TRUE),
    IEG_Positive = mean(IEG_score1 > 0.9756509, na.rm = TRUE),
    .groups = "drop"
  )

# DotPlot風描画
ggplot(
  dot_df,
  aes(
    x = Cluster,
    y = Condition
  )
) +
  geom_point(
    aes(
      size = IEG_Positive,
      color = Mean_IEG
    )
  ) +
  scale_color_gradient(
    low = "aquamarine",
    high = "tomato"
  ) +
  scale_size(range = c(2, 12)) +
  theme_bw() +
  labs(
    size = "IEG_score1 positive proportion",
    color = "Mean IEG_score1"
  )



# Fig2e_Fos
dot_df <- meta_df %>%
  group_by(Cluster, ConditionSex) %>%
  summarise(
    Mean_Fos = mean(Fos, na.rm = TRUE),
    Fos_Positive = mean(Fos > 1.590272, na.rm = TRUE),
    .groups = "drop"
  )
dot_df$ConditionSex <- factor(
  dot_df$ConditionSex,
  levels = c("CtrlM", "ColdM","2DGM","CtrlF", "ColdF","2DGF")
)
dot_df_sub <- dot_df %>%
  filter(Cluster %in% c(6, 7, 8))


ggplot(
  dot_df_sub,
  aes(
    x = factor(Cluster),
    y = ConditionSex
  )
) +
  geom_point(
    aes(
      size = Fos_Positive,
      color = Mean_Fos
    )
  ) +
  scale_color_gradient(
    low = "aquamarine",
    high = "tomato"
  ) +
  scale_size(range = c(2, 12)) +
  theme_bw()


# Fig2e_IEG
dot_df <- meta_df %>%
  group_by(Cluster, ConditionSex) %>%
  summarise(
    Mean_IEG = mean(IEG_score1, na.rm = TRUE),
    IEG_Positive = mean(IEG_score1 > 0.9756509, na.rm = TRUE),
    .groups = "drop"
  )
dot_df$ConditionSex <- factor(
  dot_df$ConditionSex,
  levels = c("CtrlM", "ColdM","2DGM","CtrlF", "ColdF","2DGF")
)
dot_df_sub <- dot_df %>%
  filter(Cluster %in% c(6, 7, 8))

ggplot(
  dot_df_sub,
  aes(
    x = factor(Cluster),
    y = ConditionSex
  )
) +
  geom_point(
    aes(
      size = IEG_Positive,
      color = Mean_IEG
    )
  ) +
  scale_color_gradient(
    low = "aquamarine",
    high = "tomato"
  ) +
  scale_size(range = c(2, 12)) +
  theme_bw()


# Fig2f_Fos
library(dplyr)
library(ggplot2)

# metadata取得
meta_df <- merged.obj_SPN_woQIHLPS_afterQC3@meta.data

# cluster情報追加
meta_df$Cluster <- Idents(merged.obj_SPN_woQIHLPS_afterQC3)

# Fos発現量追加
meta_df$Fos <- FetchData(
  merged.obj_SPN_woQIHLPS_afterQC3,
  vars = "Fos"
)[,1]

# 必要なclusterとsegmentだけ抽出
plot_df <- meta_df %>%
  filter(
    Cluster %in% c(6,7,8),
    segments %in% c(
      "T1-4",
      "T5-7",
      "T8-10",
      "T11-13",
      "L1-3"
    )
  )

# ConditionSex → Conditionへ変換
plot_df$Condition <- dplyr::case_when(
  grepl("^Ctrl", plot_df$ConditionSex) ~ "Ctrl",
  grepl("^Cold", plot_df$ConditionSex) ~ "Cold",
  grepl("^2DG",  plot_df$ConditionSex) ~ "2DG"
)

# Fos+割合計算
heat_df <- plot_df %>%
  group_by(Cluster, segments, Condition) %>%
  summarise(
    Fos_Positive = mean(Fos > 1.590272, na.rm = TRUE),
    .groups = "drop"
  )

# 順番指定
heat_df$segments <- factor(
  heat_df$segments,
  levels = rev(c(
    "T1-4",
    "T5-7",
    "T8-10",
    "T11-13",
    "L1-3"
  ))
)
heat_df$Condition <- factor(
  heat_df$Condition,
  levels = c("Ctrl","Cold","2DG")
)

# ヒートマップ
ggplot(
  heat_df,
  aes(
    x = Condition,
    y = segments,
    fill = Fos_Positive
  )
) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "white",
    high = "tomato",
    limits = c(0,1)
  ) +
  facet_wrap(~Cluster) +
  theme_bw() +
  labs(
    fill = "Fos+ ratio",
    x = "Condition",
    y = "Segment"
  )


# Fig2f_IEG
library(dplyr)
library(ggplot2)

IEGs.to.plot =  c("Fos","Arc","Egr1","Npas4","Homer1","Fosb","Nr4a1","Dusp1")
merged.obj_SPN_woQIHLPS_afterQC3 <- AddModuleScore(object = merged.obj_SPN_woQIHLPS_afterQC3,
                                                   features = list(c("Fos","Fosb","Npas4","Nr4a1","Egr1")),
                                                   name = "IEG_score",ctrl = 5)

# metadata取得
meta_df <- merged.obj_SPN_woQIHLPS_afterQC3@meta.data

# cluster情報追加
meta_df$Cluster <- Idents(merged.obj_SPN_woQIHLPS_afterQC3)

# 必要なclusterとsegmentだけ抽出
plot_df <- meta_df %>%
  filter(
    Cluster %in% c(6,7,8),
    segments %in% c(
      "T1-4",
      "T5-7",
      "T8-10",
      "T11-13",
      "L1-3"
    )
  )

# ConditionSex → Conditionへ変換
plot_df$Condition <- dplyr::case_when(
  grepl("^Ctrl", plot_df$ConditionSex) ~ "Ctrl",
  grepl("^Cold", plot_df$ConditionSex) ~ "Cold",
  grepl("^2DG",  plot_df$ConditionSex) ~ "2DG"
)

# IEG+割合計算
heat_df <- plot_df %>%
  group_by(Cluster, segments, Condition) %>%
  summarise(
    IEG_Positive = mean(IEG_score1 > 0.9756509, na.rm = TRUE),
    .groups = "drop"
  )

# 順番指定
heat_df$segments <- factor(
  heat_df$segments,
  levels = rev(c(
    "T1-4",
    "T5-7",
    "T8-10",
    "T11-13",
    "L1-3"
  ))
)

heat_df$Condition <- factor(
  heat_df$Condition,
  levels = c("Ctrl","Cold","2DG")
)

# ヒートマップ
ggplot(
  heat_df,
  aes(
    x = Condition,
    y = segments,
    fill = IEG_Positive
  )
) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "white",
    high = "tomato",
    limits = c(0,1)
  ) +
  facet_wrap(~Cluster) +
  theme_bw() +
  labs(
    fill = "IEG+ ratio",
    x = "Condition",
    y = "Segment"
  )

# Fig2g
# cFos+ SPNの絶対数を棒グラフ化する
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

obj <- merged.obj_SPN_woQIHLPS_afterQC3

DefaultAssay(obj) <- "SCT"

# 表示するsegment
segment_levels <- c(
  "T1-4",
  "T5-7",
  "T8-10",
  "T11-13"
)

# SCT data layerから取得
df <- FetchData(
  object = obj,
  vars = c(
    "Fos",
    "segments",
    "Condition",
    "SPN_subtype"
  ),
  layer = "data"
)

# Fos+の閾値（全細胞で計算）
fos_threshold <- mean(df$Fos, na.rm = TRUE) +
  2 * sd(df$Fos, na.rm = TRUE)

cat("Fos threshold =", fos_threshold, "\n")

# Fos+判定
df <- df %>%
  mutate(
    Fos_positive = Fos > fos_threshold,
    segments = factor(
      segments,
      levels = segment_levels
    )
  )

# Cluster6のみ集計
count_df_cl6 <- df %>%
  filter(
    SPN_subtype == 6,
    Condition %in% c("Cold", "2DG"),
    segments %in% segment_levels,
    !is.na(segments)
  ) %>%
  group_by(Condition, segments) %>%
  summarise(
    Total_cells = n(),
    Fos_positive = sum(Fos_positive),
    .groups = "drop"
  ) %>%
  mutate(
    Fos_negative = Total_cells - Fos_positive
  ) %>%
  select(
    Condition,
    segments,
    Fos_positive,
    Fos_negative
  ) %>%
  pivot_longer(
    cols = c(Fos_positive, Fos_negative),
    names_to = "CellType",
    values_to = "Cell_number"
  )

########################
# Cold
########################
p_cold <- ggplot(
  filter(count_df_cl6, Condition == "Cold"),
  aes(
    x = segments,
    y = Cell_number,
    fill = CellType
  )
) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = c(
      Fos_positive = "tomato2",
      Fos_negative = "grey80"
    ),
    labels = c(
      Fos_positive = "Fos+",
      Fos_negative = "Fos−"
    )
  ) +
  theme_classic() +
  labs(
    title = "Cluster 6 : Cold",
    x = "Segment",
    y = "Cell number",
    fill = NULL
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

p_cold
########################
# 2DG
########################

p_2dg <- ggplot(
  filter(count_df_cl6, Condition == "2DG"),
  aes(
    x = segments,
    y = Cell_number,
    fill = CellType
  )
) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = c(
      Fos_positive = "tomato2",
      Fos_negative = "grey80"
    ),
    labels = c(
      Fos_positive = "Fos+",
      Fos_negative = "Fos−"
    )
  ) +
  theme_classic() +
  labs(
    title = "Cluster 6 : 2DG",
    x = "Segment",
    y = "Cell number",
    fill = NULL
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

p_2dg



# EDEig2a-1
library(EnhancedVolcano)
cold_marker <-FindMarkers(merged.obj_SPN_woQIHLPS_afterQC3, assay = 'SCT', slot = 'data',fc.slot = "data", ident.1 = "Cold" ,ident.2 = "Ctrl", group.by = "Condition", min.pct = 0, logfc.threshold = 0,test.use = "wilcox",pseudocount.use = 1)

EnhancedVolcano(cold_marker,
                lab = rownames(cold_marker),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = bquote(~Log[2]~ 'Fold change'),
                ylab = bquote(~Log[10]~ 'p val'),
                pCutoff = 0.5*10^(-2),
                FCcutoff = 0.5,
                pointSize = 2.0,
                labSize = 4.0,
                colAlpha = 5/5,
                
                legendLabSize = 12,
                legendIconSize = 4.0,
                legendPosition = 'top',
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black')

# EDEig2a-2

twoDG_marker <-FindMarkers(merged.obj_SPN_woQIHLPS_afterQC3, assay = 'SCT', slot = 'data',fc.slot = "data", ident.1 = "2DG" ,ident.2 = "Ctrl", group.by = "Condition", min.pct = 0, logfc.threshold = 0,test.use = "wilcox",pseudocount.use = 1)

EnhancedVolcano(twoDG_marker,
                lab = rownames(twoDG_marker),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = bquote(~Log[2]~ 'Fold change'),
                ylab = bquote(~Log[10]~ 'p val'),
                pCutoff = 0.5*10^(-2),
                FCcutoff = 0.5,
                pointSize = 2.0,
                labSize = 4.0,
                colAlpha = 5/5,
                
                legendLabSize = 12,
                legendIconSize = 4.0,
                legendPosition = 'top',
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black')


# EDFig2b
IEGs.to.plot =  c("Fos","Arc","Egr1","Npas4","Homer1","Fosb","Nr4a1","Dusp1")
merged.obj_SPN_woQIHLPS_afterQC3 <- AddModuleScore(object = merged.obj_SPN_woQIHLPS_afterQC3,
                                                   features = list(c("Fos","Fosb","Npas4","Nr4a1","Egr1")),
                                                   name = "IEG_score",ctrl = 5)
FeaturePlot(merged.obj_SPN_woQIHLPS_afterQC3, features = c("IEG_score1"), col=c("aquamarine", "tomato2"),split.by="Condition",pt.size=0.5)



# EDFig2c_cl6
merged.obj_SPN_woQIHLPS_afterQC3_cl6 <- subset(
  merged.obj_SPN_woQIHLPS_afterQC3,
  idents = 6
)
merged.obj_SPN_woQIHLPS_afterQC3_cl6[["Xenium"]] <- split(merged.obj_SPN_woQIHLPS_afterQC3_cl6[["Xenium"]], f = merged.obj_SPN_woQIHLPS_afterQC3_cl6$orig.ident)
merged.obj_SPN_woQIHLPS_afterQC3_cl6 <- SCTransform(
  merged.obj_SPN_woQIHLPS_afterQC3_cl6,
  assay = "Xenium",
  verbose = FALSE
)
merged.obj_SPN_woQIHLPS_afterQC3_cl6 <-prep_FindMarkers(merged.obj_SPN_woQIHLPS_afterQC3_cl6)
DefaultAssay(merged.obj_SPN_woQIHLPS_afterQC3_cl6) <- "SCT"
cold_marker <-FindMarkers(merged.obj_SPN_woQIHLPS_afterQC3_cl6, assay = 'SCT', slot = 'data',fc.slot = "data", ident.1 = "Cold" ,ident.2 = "Ctrl", group.by = "Condition",subset.ident = 6, min.pct = 0, logfc.threshold = 0,test.use = "wilcox",pseudocount.use = 1)

EnhancedVolcano(cold_marker,
                lab = rownames(cold_marker),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = bquote(~Log[2]~ 'Fold change'),
                ylab = bquote(~Log[10]~ 'p val'),
                pCutoff = 0.5*10^(-2),
                FCcutoff = 0.5,
                pointSize = 2.0,
                labSize = 4.0,
                colAlpha = 5/5,
                
                legendLabSize = 12,
                legendIconSize = 4.0,
                legendPosition = 'top',
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black')

twoDG_marker <-FindMarkers(merged.obj_SPN_woQIHLPS_afterQC3_cl6, assay = 'SCT', slot = 'data',fc.slot = "data", ident.1 = "2DG" ,ident.2 = "Ctrl", group.by = "Condition",subset.ident = 6, min.pct = 0, logfc.threshold = 0,test.use = "wilcox",pseudocount.use = 1)

EnhancedVolcano(twoDG_marker,
                lab = rownames(twoDG_marker),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = bquote(~Log[2]~ 'Fold change'),
                ylab = bquote(~Log[10]~ 'p val'),
                pCutoff = 0.5*10^(-2),
                FCcutoff = 0.5,
                pointSize = 2.0,
                labSize = 4.0,
                colAlpha = 5/5,
                
                legendLabSize = 12,
                legendIconSize = 4.0,
                legendPosition = 'top',
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black')

# EDFig 2c_cl7
merged.obj_SPN_woQIHLPS_afterQC3_cl7 <- subset(
  merged.obj_SPN_woQIHLPS_afterQC3,
  idents = 7
)
merged.obj_SPN_woQIHLPS_afterQC3_cl7[["Xenium"]] <- split(merged.obj_SPN_woQIHLPS_afterQC3_cl7[["Xenium"]], f = merged.obj_SPN_woQIHLPS_afterQC3_cl7$orig.ident)
merged.obj_SPN_woQIHLPS_afterQC3_cl7 <- SCTransform(
  merged.obj_SPN_woQIHLPS_afterQC3_cl7,
  assay = "Xenium",
  verbose = FALSE
)
merged.obj_SPN_woQIHLPS_afterQC3_cl7 <-prep_FindMarkers(merged.obj_SPN_woQIHLPS_afterQC3_cl7)
DefaultAssay(merged.obj_SPN_woQIHLPS_afterQC3_cl7) <- "SCT"
cold_marker <-FindMarkers(merged.obj_SPN_woQIHLPS_afterQC3_cl7, assay = 'SCT', slot = 'data',fc.slot = "data", ident.1 = "Cold" ,ident.2 = "Ctrl", group.by = "Condition",subset.ident =7, min.pct = 0, logfc.threshold = 0,test.use = "wilcox",pseudocount.use = 1)

EnhancedVolcano(cold_marker,
                lab = rownames(cold_marker),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = bquote(~Log[2]~ 'Fold change'),
                ylab = bquote(~Log[10]~ 'p val'),
                pCutoff = 0.5*10^(-2),
                FCcutoff = 0.5,
                pointSize = 2.0,
                labSize = 4.0,
                colAlpha = 5/5,
                
                legendLabSize = 12,
                legendIconSize = 4.0,
                legendPosition = 'top',
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black')

# EDFig 2c_cl8

merged.obj_SPN_woQIHLPS_afterQC3_cl8 <- subset(
  merged.obj_SPN_woQIHLPS_afterQC3,
  idents = 8
)
merged.obj_SPN_woQIHLPS_afterQC3_cl8[["Xenium"]] <- split(merged.obj_SPN_woQIHLPS_afterQC3_cl8[["Xenium"]], f = merged.obj_SPN_woQIHLPS_afterQC3_cl8$orig.ident)
merged.obj_SPN_woQIHLPS_afterQC3_cl8 <- SCTransform(
  merged.obj_SPN_woQIHLPS_afterQC3_cl8,
  assay = "Xenium",
  verbose = FALSE
)
merged.obj_SPN_woQIHLPS_afterQC3_cl8 <-prep_FindMarkers(merged.obj_SPN_woQIHLPS_afterQC3_cl8)
DefaultAssay(merged.obj_SPN_woQIHLPS_afterQC3_cl8) <- "SCT"


twoDG_marker <-FindMarkers(merged.obj_SPN_woQIHLPS_afterQC3_cl8, assay = 'SCT', slot = 'data',fc.slot = "data", ident.1 = "2DG" ,ident.2 = "Ctrl", group.by = "Condition",subset.ident = 8, min.pct = 0, logfc.threshold = 0,test.use = "wilcox",pseudocount.use = 1)

EnhancedVolcano(twoDG_marker,
                lab = rownames(twoDG_marker),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = bquote(~Log[2]~ 'Fold change'),
                ylab = bquote(~Log[10]~ 'p val'),
                pCutoff = 0.5*10^(-2),
                FCcutoff = 0.5,
                pointSize = 2.0,
                labSize = 4.0,
                colAlpha = 5/5,
                
                legendLabSize = 12,
                legendIconSize = 4.0,
                legendPosition = 'top',
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black')



# ED Fig2d
FeaturePlot(merged.obj_SPN_woQIHLPS_afterQC3, features = c("Fos"), col=c("aquamarine", "tomato2"),split.by="ConditionSex",pt.size=0.5)


# EDFig2e_IEGscore
library(dplyr)
library(ggplot2)

# metadata取得
meta_df <- merged.obj_SPN_woQIHLPS_afterQC3@meta.data

# cluster情報追加
meta_df$Cluster <- Idents(merged.obj_SPN_woQIHLPS_afterQC3)

# 必要なclusterとsegmentだけ抽出
plot_df <- meta_df %>%
  filter(
    Cluster %in% c(6,7,8),
    segments %in% c(
      "T1-4",
      "T5-7",
      "T8-10",
      "T11-13",
      "L1-3"
    )
  )

# Condition作成
plot_df$Condition <- dplyr::case_when(
  grepl("^Ctrl", plot_df$ConditionSex) ~ "Ctrl",
  grepl("^Cold", plot_df$ConditionSex) ~ "Cold",
  grepl("^2DG",  plot_df$ConditionSex) ~ "2DG"
)

# sex作成
plot_df$Sex <- dplyr::case_when(
  grepl("F$", plot_df$ConditionSex) ~ "Female",
  grepl("M$", plot_df$ConditionSex) ~ "Male"
)

# IEG+割合計算
heat_df <- plot_df %>%
  group_by(Cluster, Sex, segments, Condition) %>%
  summarise(
    IEG_Positive = mean(
      IEG_score1 > 0.9756509,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# 順番指定
heat_df$segments <- factor(
  heat_df$segments,
  levels = rev(c(
    "T1-4",
    "T5-7",
    "T8-10",
    "T11-13",
    "L1-3"
  ))
)

heat_df$Condition <- factor(
  heat_df$Condition,
  levels = c("Ctrl","Cold","2DG")
)

heat_df$Sex <- factor(
  heat_df$Sex,
  levels = c("Female","Male")
)

# ヒートマップ
ggplot(
  heat_df,
  aes(
    x = Condition,
    y = segments,
    fill = IEG_Positive
  )
) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "white",
    high = "tomato",
    limits = c(0,1)
  ) +
  facet_grid(
    Sex ~ Cluster
  ) +
  theme_bw() +
  labs(
    fill = "IEG+ ratio",
    x = "Condition",
    y = "Segment"
  )


# EDFig2f
library(dplyr)
library(ggplot2)

Idents(merged.obj_SPN_woQIHLPS_afterQC3)<-"SPN_subtype"
meta_df <- merged.obj_SPN_woQIHLPS_afterQC3@meta.data

# Condition作成
meta_df$Condition <- case_when(
  grepl("^Ctrl", meta_df$ConditionSex) ~ "Ctrl",
  grepl("^Cold", meta_df$ConditionSex) ~ "Cold",
  grepl("^2DG",  meta_df$ConditionSex) ~ "2DG"
)

# seurat cluster取得
meta_df$Cluster <- meta_df$SPN_subtype

# 各Conditionにおけるcluster6割合
plot_df <- meta_df %>%
  group_by(Condition) %>%
  summarise(
    Cluster6 = sum(Cluster == 6),
    Total = n(),
    Percent = 100 * Cluster6 / Total,
    .groups = "drop"
  )

plot_df$Condition <- factor(
  plot_df$Condition,
  levels = c("Ctrl", "Cold", "2DG")
)

ggplot(
  plot_df,
  aes(
    x = Condition,
    y = Percent
  )
) +
  geom_col(
    width = 0.7,
    fill = "tomato2"
  ) +
  geom_text(
    aes(
      label = sprintf("%.1f%%", Percent)
    ),
    vjust = -0.4,
    size = 5
  ) +
  theme_bw(base_size = 14) +
  labs(
    x = "",
    y = "Cluster 6 cells (%)"
  ) +
  scale_y_continuous(
    limits = c(0, 20),
    breaks = seq(0, 20, by = 2),
    expand = c(0, 0)
  )
