# Code that assigns a section name and a segment to each cell based on the loaded GeoJSON
PolygonAddSections <- function(object, polygons, pixel_size = 0.2125,sections){
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
  message("Start cropping cells...")
  cropped_cells <- lapply(names(polygons), function(name){
    fov <- object@images[[name]]
    cells_to_keep <- PolygonCropCentroid(fov, polygons[[name]], pixel_size)
    message("Done: ", name)
    return(cells_to_keep)
  })
  combined_cropped_cells <- unlist(cropped_cells)
  object$sections[combined_cropped_cells] <- sections
  message("Done cropping cells...")
  return(object)
}



PolygonCropCentroid <- function(fov, polygon, pixel_size = 0.2125){
  #' Crop an FOV using a polygon
  #' 
  #' @description This function crops an FOV including centroids, sections and molecules.
  #'
  #' The FOV must contain centroids as those are used to determine if the cells are in the polygon or not. 
  #
  #' @param fov fov Object The FOV to be cropped.
  #' @param polygon sf POLYGON object.
  #' @param pixel_size ratio between the FOV coordinates and the sf Objects coordinates.
  #' @usage PolygonCropFOV(polygon, pyxel_size = 0.2125)
  #' @return A new cropped FOV object.
  #' @details If no cells are found the unchanged fov is returned with a Warning.
  #' @note TODO: if no centroids are present try to use segmentation or fail with an Error.
  
  # Use the centroids to identify the cells that should be kept
  message("Start processing centroids...")
  polygon <- st_set_crs(polygon, NA)
  centroids <- st_as_sfc(fov$centroids, forceMulti = FALSE)
  centroids <- st_set_crs(centroids, NA)
  centroids <- st_combine(centroids)
  centroids <- centroids / pixel_size
  centroids <- st_cast(centroids, "POINT")
  contained <- st_contains_properly(polygon, centroids, sparse = FALSE)
  cells_to_keep <- Cells(fov)[contained]
  message("Done processing centroids...")
  
  # If no cells are found within the polygon, return the FOV unchanged with a warning
  if (!any(contained)){
    warning("No cells, found inside the polygon. Returning original object")
  }
  return(cells_to_keep)
}


merged.obj_SPN_woQIHLPS_afterQC2$sections <- NA


# ColdN1
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T1_4_1), pixel_size = 0.2125,sections="ColdM_1_1_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T1_4_2), pixel_size = 0.2125,sections="ColdM_1_1_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T1_4_3), pixel_size = 0.2125,sections="ColdM_1_1_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T1_4_4), pixel_size = 0.2125,sections="ColdM_1_1_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T1_4_5), pixel_size = 0.2125,sections="ColdM_1_1_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T1_4_6), pixel_size = 0.2125,sections="ColdM_1_1_6")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T5_7_1), pixel_size = 0.2125,sections="ColdM_1_2_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T5_7_2), pixel_size = 0.2125,sections="ColdM_1_2_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T5_7_3), pixel_size = 0.2125,sections="ColdM_1_2_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T5_7_4), pixel_size = 0.2125,sections="ColdM_1_2_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T5_7_5), pixel_size = 0.2125,sections="ColdM_1_2_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T5_7_6), pixel_size = 0.2125,sections="ColdM_1_2_6")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T8_10_1), pixel_size = 0.2125,sections="ColdM_1_3_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T8_10_2), pixel_size = 0.2125,sections="ColdM_1_3_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T8_10_3), pixel_size = 0.2125,sections="ColdM_1_3_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T8_10_4), pixel_size = 0.2125,sections="ColdM_1_3_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T8_10_5), pixel_size = 0.2125,sections="ColdM_1_3_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T8_10_6), pixel_size = 0.2125,sections="ColdM_1_3_6")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T11_13_1), pixel_size = 0.2125,sections="ColdM_1_4_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T11_13_2), pixel_size = 0.2125,sections="ColdM_1_4_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T11_13_3), pixel_size = 0.2125,sections="ColdM_1_4_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T11_13_4), pixel_size = 0.2125,sections="ColdM_1_4_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T11_13_5), pixel_size = 0.2125,sections="ColdM_1_4_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T11_13_6), pixel_size = 0.2125,sections="ColdM_1_4_6")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T11_13_7), pixel_size = 0.2125,sections="ColdM_1_4_7")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_T11_13_8), pixel_size = 0.2125,sections="ColdM_1_4_8")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_L1_3_1), pixel_size = 0.2125,sections="ColdM_1_5_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_L1_3_2), pixel_size = 0.2125,sections="ColdM_1_5_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_L1_3_3), pixel_size = 0.2125,sections="ColdM_1_5_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_L1_3_4), pixel_size = 0.2125,sections="ColdM_1_5_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_L1_3_5), pixel_size = 0.2125,sections="ColdM_1_5_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_L1_3_6), pixel_size = 0.2125,sections="ColdM_1_5_6")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_L45_1), pixel_size = 0.2125,sections="ColdM_1_6_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_L45_2), pixel_size = 0.2125,sections="ColdM_1_6_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_L45_3), pixel_size = 0.2125,sections="ColdM_1_6_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_L45_4), pixel_size = 0.2125,sections="ColdM_1_6_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_L45_5), pixel_size = 0.2125,sections="ColdM_1_6_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_L45_6), pixel_size = 0.2125,sections="ColdM_1_6_6")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_L6S1_1), pixel_size = 0.2125,sections="ColdM_1_7_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_L6S1_2), pixel_size = 0.2125,sections="ColdM_1_7_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_L6S1_3), pixel_size = 0.2125,sections="ColdM_1_7_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_L6S1_4), pixel_size = 0.2125,sections="ColdM_1_7_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_L6S1_5), pixel_size = 0.2125,sections="ColdM_1_7_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_L6S1_6), pixel_size = 0.2125,sections="ColdM_1_7_6")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_S2_4_1), pixel_size = 0.2125,sections="ColdM_1_8_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_S2_4_2), pixel_size = 0.2125,sections="ColdM_1_8_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_S2_4_3), pixel_size = 0.2125,sections="ColdM_1_8_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_S2_4_4), pixel_size = 0.2125,sections="ColdM_1_8_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_S2_4_5), pixel_size = 0.2125,sections="ColdM_1_8_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_S2_4_6), pixel_size = 0.2125,sections="ColdM_1_8_6")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2=Polygon_Cold_N1_S2_4_7), pixel_size = 0.2125,sections="ColdM_1_8_7")


# N2
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_T1_4_1), pixel_size = 0.2125,sections="ColdM_2_1_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_T1_4_2), pixel_size = 0.2125,sections="ColdM_2_1_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_T1_4_3), pixel_size = 0.2125,sections="ColdM_2_1_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_T1_4_4), pixel_size = 0.2125,sections="ColdM_2_1_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_T1_4_5), pixel_size = 0.2125,sections="ColdM_2_1_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_T1_4_6), pixel_size = 0.2125,sections="ColdM_2_1_6")
# 
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N3_T5_7_1), pixel_size = 0.2125,sections="ColdM_3_2_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N3_T5_7_2), pixel_size = 0.2125,sections="ColdM_3_2_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N3_T5_7_3), pixel_size = 0.2125,sections="ColdM_3_2_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N3_T5_7_4), pixel_size = 0.2125,sections="ColdM_3_2_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N3_T5_7_5), pixel_size = 0.2125,sections="ColdM_3_2_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N3_T5_7_6), pixel_size = 0.2125,sections="ColdM_3_2_6")


merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_T8_10_1), pixel_size = 0.2125,sections="ColdM_2_3_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_T8_10_2), pixel_size = 0.2125,sections="ColdM_2_3_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_T8_10_3), pixel_size = 0.2125,sections="ColdM_2_3_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_T8_10_4), pixel_size = 0.2125,sections="ColdM_2_3_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_T8_10_5), pixel_size = 0.2125,sections="ColdM_2_3_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_T8_10_6), pixel_size = 0.2125,sections="ColdM_2_3_6")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_T11_13_1), pixel_size = 0.2125,sections="ColdM_2_4_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_T11_13_2), pixel_size = 0.2125,sections="ColdM_2_4_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_T11_13_3), pixel_size = 0.2125,sections="ColdM_2_4_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_T11_13_4), pixel_size = 0.2125,sections="ColdM_2_4_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_T11_13_5), pixel_size = 0.2125,sections="ColdM_2_4_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_T11_13_6), pixel_size = 0.2125,sections="ColdM_2_4_6")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_T11_13_7), pixel_size = 0.2125,sections="ColdM_2_4_7")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_T11_13_8), pixel_size = 0.2125,sections="ColdM_2_4_8")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_L1_3_1), pixel_size = 0.2125,sections="ColdM_2_5_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_L1_3_2), pixel_size = 0.2125,sections="ColdM_2_5_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_L1_3_3), pixel_size = 0.2125,sections="ColdM_2_5_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_L1_3_4), pixel_size = 0.2125,sections="ColdM_2_5_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_L1_3_5), pixel_size = 0.2125,sections="ColdM_2_5_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_L1_3_6), pixel_size = 0.2125,sections="ColdM_2_5_6")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_L45_1), pixel_size = 0.2125,sections="ColdM_2_6_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_L45_2), pixel_size = 0.2125,sections="ColdM_2_6_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_L45_3), pixel_size = 0.2125,sections="ColdM_2_6_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_L45_4), pixel_size = 0.2125,sections="ColdM_2_6_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_L45_5), pixel_size = 0.2125,sections="ColdM_2_6_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_L45_6), pixel_size = 0.2125,sections="ColdM_2_6_6")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_L6S1_1), pixel_size = 0.2125,sections="ColdM_2_7_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_L6S1_2), pixel_size = 0.2125,sections="ColdM_2_7_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_L6S1_3), pixel_size = 0.2125,sections="ColdM_2_7_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_L6S1_4), pixel_size = 0.2125,sections="ColdM_2_7_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_L6S1_5), pixel_size = 0.2125,sections="ColdM_2_7_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_L6S1_6), pixel_size = 0.2125,sections="ColdM_2_7_6")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_S2_4_1), pixel_size = 0.2125,sections="ColdM_2_8_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_S2_4_2), pixel_size = 0.2125,sections="ColdM_2_8_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_S2_4_3), pixel_size = 0.2125,sections="ColdM_2_8_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_S2_4_4), pixel_size = 0.2125,sections="ColdM_2_8_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_S2_4_5), pixel_size = 0.2125,sections="ColdM_2_8_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_S2_4_6), pixel_size = 0.2125,sections="ColdM_2_8_6")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold2=Polygon_Cold_N2_S2_4_7), pixel_size = 0.2125,sections="ColdM_2_8_7")

# N3
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_T1_4_1), pixel_size = 0.2125,sections="ColdM_3_1_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_T1_4_2), pixel_size = 0.2125,sections="ColdM_3_1_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_T1_4_3), pixel_size = 0.2125,sections="ColdM_3_1_3")
# 
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_T5_7_1), pixel_size = 0.2125,sections="3_2_1")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_T5_7_2), pixel_size = 0.2125,sections="3_2_2")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_T5_7_3), pixel_size = 0.2125,sections="3_2_3")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_T5_7_4), pixel_size = 0.2125,sections="3_2_4")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_T5_7_5), pixel_size = 0.2125,sections="3_2_5")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_T5_7_6), pixel_size = 0.2125,sections="3_2_6")


merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_T8_10_1), pixel_size = 0.2125,sections="ColdM_3_3_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_T8_10_2), pixel_size = 0.2125,sections="ColdM_3_3_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_T8_10_3), pixel_size = 0.2125,sections="ColdM_3_3_3")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_T11_13_1), pixel_size = 0.2125,sections="ColdM_3_4_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_T11_13_2), pixel_size = 0.2125,sections="ColdM_3_4_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_T11_13_3), pixel_size = 0.2125,sections="ColdM_3_4_3")


merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_L1_3_1), pixel_size = 0.2125,sections="ColdM_3_5_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_L1_3_2), pixel_size = 0.2125,sections="ColdM_3_5_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_L1_3_3), pixel_size = 0.2125,sections="ColdM_3_5_3")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_L45_1), pixel_size = 0.2125,sections="ColdM_3_6_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_L45_2), pixel_size = 0.2125,sections="ColdM_3_6_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_L45_3), pixel_size = 0.2125,sections="ColdM_3_6_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_L45_4), pixel_size = 0.2125,sections="ColdM_3_6_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_L45_5), pixel_size = 0.2125,sections="ColdM_3_6_5")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_L6S1_1), pixel_size = 0.2125,sections="ColdM_3_7_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_L6S1_2), pixel_size = 0.2125,sections="ColdM_3_7_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_L6S1_3), pixel_size = 0.2125,sections="ColdM_3_7_3")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_S2_4_1), pixel_size = 0.2125,sections="ColdM_3_8_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_S2_4_2), pixel_size = 0.2125,sections="ColdM_3_8_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Cold3=Polygon_Cold_N3_S2_4_3), pixel_size = 0.2125,sections="ColdM_3_8_3")


# CtrlN1
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_T1_4_1), pixel_size = 0.2125,sections="CtrlM_1_1_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_T1_4_2), pixel_size = 0.2125,sections="CtrlM_1_1_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_T1_4_3), pixel_size = 0.2125,sections="CtrlM_1_1_3")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_T5_7_1), pixel_size = 0.2125,sections="CtrlM_1_2_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_T5_7_2), pixel_size = 0.2125,sections="CtrlM_1_2_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_T5_7_3), pixel_size = 0.2125,sections="CtrlM_1_2_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_T5_7_4), pixel_size = 0.2125,sections="CtrlM_1_2_4")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_T8_10_1), pixel_size = 0.2125,sections="CtrlM_1_3_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_T8_10_2), pixel_size = 0.2125,sections="CtrlM_1_3_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_T8_10_3), pixel_size = 0.2125,sections="CtrlM_1_3_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_T8_10_4), pixel_size = 0.2125,sections="CtrlM_1_3_4")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_T11_13_1), pixel_size = 0.2125,sections="CtrlM_1_4_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_T11_13_2), pixel_size = 0.2125,sections="CtrlM_1_4_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_T11_13_3), pixel_size = 0.2125,sections="CtrlM_1_4_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_T11_13_4), pixel_size = 0.2125,sections="CtrlM_1_4_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_T11_13_5), pixel_size = 0.2125,sections="CtrlM_1_4_5")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_L1_3_1), pixel_size = 0.2125,sections="CtrlM_1_5_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_L1_3_2), pixel_size = 0.2125,sections="CtrlM_1_5_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_L1_3_3), pixel_size = 0.2125,sections="CtrlM_1_5_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_L1_3_4), pixel_size = 0.2125,sections="CtrlM_1_5_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_L1_3_5), pixel_size = 0.2125,sections="CtrlM_1_5_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_L1_3_6), pixel_size = 0.2125,sections="CtrlM_1_5_6")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_L45_1), pixel_size = 0.2125,sections="CtrlM_1_6_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_L45_2), pixel_size = 0.2125,sections="CtrlM_1_6_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_L45_3), pixel_size = 0.2125,sections="CtrlM_1_6_3")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_L6S1_1), pixel_size = 0.2125,sections="CtrlM_1_7_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_L6S1_2), pixel_size = 0.2125,sections="CtrlM_1_7_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_L6S1_3), pixel_size = 0.2125,sections="CtrlM_1_7_3")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov=Polygon_Ctrl_N1_S2_4_1), pixel_size = 0.2125,sections="CtrlM_1_8_1")

# N2
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T1_4_1), pixel_size = 0.2125,sections="CtrlM_2_1_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T1_4_2), pixel_size = 0.2125,sections="CtrlM_2_1_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T1_4_3), pixel_size = 0.2125,sections="CtrlM_2_1_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T1_4_4), pixel_size = 0.2125,sections="CtrlM_2_1_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T1_4_5), pixel_size = 0.2125,sections="CtrlM_2_1_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T1_4_6), pixel_size = 0.2125,sections="CtrlM_2_1_6")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T5_7_1), pixel_size = 0.2125,sections="CtrlM_2_2_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T5_7_2), pixel_size = 0.2125,sections="CtrlM_2_2_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T5_7_3), pixel_size = 0.2125,sections="CtrlM_2_2_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T5_7_4), pixel_size = 0.2125,sections="CtrlM_2_2_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T5_7_5), pixel_size = 0.2125,sections="CtrlM_2_2_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T5_7_6), pixel_size = 0.2125,sections="CtrlM_2_2_6")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T8_10_1), pixel_size = 0.2125,sections="CtrlM_2_3_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T8_10_2), pixel_size = 0.2125,sections="CtrlM_2_3_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T8_10_3), pixel_size = 0.2125,sections="CtrlM_2_3_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T8_10_4), pixel_size = 0.2125,sections="CtrlM_2_3_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T8_10_5), pixel_size = 0.2125,sections="CtrlM_2_3_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T8_10_6), pixel_size = 0.2125,sections="CtrlM_2_3_6")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T8_10_7), pixel_size = 0.2125,sections="CtrlM_2_3_7")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T8_10_8), pixel_size = 0.2125,sections="CtrlM_2_3_8")


merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T11_13_1), pixel_size = 0.2125,sections="CtrlM_2_4_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T11_13_2), pixel_size = 0.2125,sections="CtrlM_2_4_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T11_13_3), pixel_size = 0.2125,sections="CtrlM_2_4_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T11_13_4), pixel_size = 0.2125,sections="CtrlM_2_4_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T11_13_5), pixel_size = 0.2125,sections="CtrlM_2_4_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T11_13_6), pixel_size = 0.2125,sections="CtrlM_2_4_6")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T11_13_7), pixel_size = 0.2125,sections="CtrlM_2_4_7")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_T11_13_8), pixel_size = 0.2125,sections="CtrlM_2_4_8")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_L1_3_1), pixel_size = 0.2125,sections="CtrlM_2_5_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_L1_3_2), pixel_size = 0.2125,sections="CtrlM_2_5_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_L1_3_3), pixel_size = 0.2125,sections="CtrlM_2_5_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_L1_3_4), pixel_size = 0.2125,sections="CtrlM_2_5_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_L1_3_5), pixel_size = 0.2125,sections="CtrlM_2_5_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_L1_3_6), pixel_size = 0.2125,sections="CtrlM_2_5_6")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_L1_3_7), pixel_size = 0.2125,sections="CtrlM_2_5_7")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_L1_3_8), pixel_size = 0.2125,sections="CtrlM_2_5_8")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_L45_1), pixel_size = 0.2125,sections="CtrlM_2_6_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_L45_2), pixel_size = 0.2125,sections="CtrlM_2_6_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_L45_3), pixel_size = 0.2125,sections="CtrlM_2_6_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_L45_4), pixel_size = 0.2125,sections="CtrlM_2_6_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_L45_5), pixel_size = 0.2125,sections="CtrlM_2_6_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_L45_6), pixel_size = 0.2125,sections="CtrlM_2_6_6")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_L6S1_1), pixel_size = 0.2125,sections="CtrlM_2_7_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_L6S1_2), pixel_size = 0.2125,sections="CtrlM_2_7_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_L6S1_3), pixel_size = 0.2125,sections="CtrlM_2_7_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_L6S1_4), pixel_size = 0.2125,sections="CtrlM_2_7_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_L6S1_5), pixel_size = 0.2125,sections="CtrlM_2_7_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_L6S1_6), pixel_size = 0.2125,sections="CtrlM_2_7_6")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_S2_4_1), pixel_size = 0.2125,sections="CtrlM_2_8_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_S2_4_2), pixel_size = 0.2125,sections="CtrlM_2_8_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_S2_4_3), pixel_size = 0.2125,sections="CtrlM_2_8_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_S2_4_4), pixel_size = 0.2125,sections="CtrlM_2_8_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl2=Polygon_Ctrl_N2_S2_4_5), pixel_size = 0.2125,sections="CtrlM_2_8_5")


# N3
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T1_4_1), pixel_size = 0.2125,sections="CtrlM_3_1_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T1_4_2), pixel_size = 0.2125,sections="CtrlM_3_1_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T1_4_3), pixel_size = 0.2125,sections="CtrlM_3_1_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T1_4_4), pixel_size = 0.2125,sections="CtrlM_3_1_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T1_4_5), pixel_size = 0.2125,sections="CtrlM_3_1_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T1_4_6), pixel_size = 0.2125,sections="CtrlM_3_1_6")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T5_7_1), pixel_size = 0.2125,sections="CtrlM_3_2_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T5_7_2), pixel_size = 0.2125,sections="CtrlM_3_2_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T5_7_3), pixel_size = 0.2125,sections="CtrlM_3_2_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T5_7_4), pixel_size = 0.2125,sections="CtrlM_3_2_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T5_7_5), pixel_size = 0.2125,sections="CtrlM_3_2_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T5_7_6), pixel_size = 0.2125,sections="CtrlM_3_2_6")


merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T8_10_1), pixel_size = 0.2125,sections="CtrlM_3_3_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T8_10_2), pixel_size = 0.2125,sections="CtrlM_3_3_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T8_10_3), pixel_size = 0.2125,sections="CtrlM_3_3_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T8_10_4), pixel_size = 0.2125,sections="CtrlM_3_3_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T8_10_5), pixel_size = 0.2125,sections="CtrlM_3_3_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T8_10_6), pixel_size = 0.2125,sections="CtrlM_3_3_6")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T11_13_1), pixel_size = 0.2125,sections="CtrlM_3_4_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T11_13_2), pixel_size = 0.2125,sections="CtrlM_3_4_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T11_13_3), pixel_size = 0.2125,sections="CtrlM_3_4_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T11_13_4), pixel_size = 0.2125,sections="CtrlM_3_4_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T11_13_5), pixel_size = 0.2125,sections="CtrlM_3_4_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T11_13_6), pixel_size = 0.2125,sections="CtrlM_3_4_6")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_T11_13_7), pixel_size = 0.2125,sections="CtrlM_3_4_7")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L1_3_1), pixel_size = 0.2125,sections="CtrlM_3_5_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L1_3_2), pixel_size = 0.2125,sections="CtrlM_3_5_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L1_3_3), pixel_size = 0.2125,sections="CtrlM_3_5_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L1_3_4), pixel_size = 0.2125,sections="CtrlM_3_5_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L1_3_5), pixel_size = 0.2125,sections="CtrlM_3_5_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L1_3_6), pixel_size = 0.2125,sections="CtrlM_3_5_6")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L1_3_7), pixel_size = 0.2125,sections="CtrlM_3_5_7")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L1_3_8), pixel_size = 0.2125,sections="CtrlM_3_5_8")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L1_3_9), pixel_size = 0.2125,sections="CtrlM_3_5_9")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L1_3_10), pixel_size = 0.2125,sections="CtrlM_3_5_10")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L1_3_11), pixel_size = 0.2125,sections="CtrlM_3_5_11")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L1_3_12), pixel_size = 0.2125,sections="CtrlM_3_5_12")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L1_3_13), pixel_size = 0.2125,sections="CtrlM_3_5_13")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L1_3_14), pixel_size = 0.2125,sections="CtrlM_3_5_14")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L1_3_15), pixel_size = 0.2125,sections="CtrlM_3_5_15")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L1_3_16), pixel_size = 0.2125,sections="CtrlM_3_5_16")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L45_1), pixel_size = 0.2125,sections="CtrlM_3_6_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L45_2), pixel_size = 0.2125,sections="CtrlM_3_6_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L45_3), pixel_size = 0.2125,sections="CtrlM_3_6_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L45_4), pixel_size = 0.2125,sections="CtrlM_3_6_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L45_5), pixel_size = 0.2125,sections="CtrlM_3_6_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L45_6), pixel_size = 0.2125,sections="CtrlM_3_6_6")

merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L6S1_1), pixel_size = 0.2125,sections="CtrlM_3_7_1")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L6S1_2), pixel_size = 0.2125,sections="CtrlM_3_7_2")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L6S1_3), pixel_size = 0.2125,sections="CtrlM_3_7_3")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L6S1_4), pixel_size = 0.2125,sections="CtrlM_3_7_4")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L6S1_5), pixel_size = 0.2125,sections="CtrlM_3_7_5")
merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.Ctrl3=Polygon_Ctrl_N3_L6S1_6), pixel_size = 0.2125,sections="CtrlM_3_7_6")



# # CtrlF
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_T1_4_1), pixel_size = 0.2125,sections="1_1_1")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_T1_4_2), pixel_size = 0.2125,sections="1_1_2")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_T1_4_3), pixel_size = 0.2125,sections="1_1_3")
# 
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_T5_7_1), pixel_size = 0.2125,sections="1_2_1")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_T5_7_2), pixel_size = 0.2125,sections="1_2_2")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_T5_7_3), pixel_size = 0.2125,sections="1_2_3")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_T5_7_4), pixel_size = 0.2125,sections="1_2_4")
# 
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_T8_10_1), pixel_size = 0.2125,sections="1_3_1")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_T8_10_2), pixel_size = 0.2125,sections="1_3_2")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_T8_10_3), pixel_size = 0.2125,sections="1_3_3")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_T8_10_4), pixel_size = 0.2125,sections="1_3_4")
# 
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_T11_13_1), pixel_size = 0.2125,sections="1_4_1")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_T11_13_2), pixel_size = 0.2125,sections="1_4_2")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_T11_13_3), pixel_size = 0.2125,sections="1_4_3")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_T11_13_4), pixel_size = 0.2125,sections="1_4_4")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_T11_13_5), pixel_size = 0.2125,sections="1_4_5")
# 
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_L1_3_1), pixel_size = 0.2125,sections="1_5_1")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_L1_3_2), pixel_size = 0.2125,sections="1_5_2")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_L1_3_3), pixel_size = 0.2125,sections="1_5_3")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_L1_3_4), pixel_size = 0.2125,sections="1_5_4")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_L1_3_5), pixel_size = 0.2125,sections="1_5_5")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_L1_3_6), pixel_size = 0.2125,sections="1_5_6")
# 
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_L45_1), pixel_size = 0.2125,sections="1_6_1")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_L45_2), pixel_size = 0.2125,sections="1_6_2")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_L45_3), pixel_size = 0.2125,sections="1_6_3")
# 
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_L6S1_1), pixel_size = 0.2125,sections="1_7_1")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_L6S1_2), pixel_size = 0.2125,sections="1_7_2")
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_L6S1_3), pixel_size = 0.2125,sections="1_7_3")
# 
# merged.obj_SPN_woQIHLPS_afterQC2<-PolygonAddSections(merged.obj_SPN_woQIHLPS_afterQC2, polygons = list(fov.2.3=Polygon_CtrlF_N1_S2_4_1), pixel_size = 0.2125,sections="1_8_1")
# 



# CtrlF
# ============================================
# PolygonAddSections を自動化
# ============================================

region_map <- c(
  "1" = "T1_4",
  "2" = "T5_7",
  "3" = "T8_10",
  "4" = "T11_13",
  "5" = "L1_3",
  "6" = "L45",
  "7" = "L6S1",
  "8" = "S2_4"
)

# geojson ファイル一覧
files <- list.files(
  "./Xenium_SC/ROI_geojson/ctrlF",
  pattern = "\\.geojson$",
  full.names = FALSE
)

# ファイル名順に並べる
files <- sort(files)

for (f in files) {
  
  # 例: 3-4-7.geojson
  fname <- tools::file_path_sans_ext(f)
  
  # c("3","4","7")
  parts <- strsplit(fname, "-")[[1]]
  
  animal_num <- parts[1]
  region_num <- parts[2]
  roi_num <- parts[3]
  
  region_name <- region_map[region_num]
  
  # Polygon object 名
  polygon_name <- paste0(
    "Polygon_CtrlF_N",
    animal_num,
    "_",
    region_name,
    "_",
    roi_num
  )
  
  # sections 名
  section_name <- paste0("CtrlF_",
    animal_num, "_",
    region_num, "_",
    roi_num
  )
  
  # polygon object を取得
  polygon_obj <- get(polygon_name)
  
  # N1 / N2 / N3 に応じて fov 名変更
  fov_name <- switch(
    animal_num,
    "1" = "fov.2.3",
    "2" = "fov.2.3",
    "3" = "fov.2.3"
  )
  
  # list(fov=...) を動的生成
  polygon_list <- setNames(
    list(polygon_obj),
    fov_name
  )
  
  # PolygonAddSections
  merged.obj_SPN_woQIHLPS_afterQC2 <- PolygonAddSections(
    merged.obj_SPN_woQIHLPS_afterQC2,
    polygons = polygon_list,
    pixel_size = 0.2125,
    sections = section_name
  )
  
  cat("Added:", section_name, "\n")
}



# ColdF
# ============================================
# PolygonAddSections を自動化
# ============================================

region_map <- c(
  "1" = "T1_4",
  "2" = "T5_7",
  "3" = "T8_10",
  "4" = "T11_13",
  "5" = "L1_3",
  "6" = "L45",
  "7" = "L6S1",
  "8" = "S2_4"
)

# geojson ファイル一覧
files <- list.files(
  "./Xenium_SC/ROI_geojson/ColdF",
  pattern = "\\.geojson$",
  full.names = FALSE
)

# ファイル名順に並べる
files <- sort(files)

for (f in files) {
  
  # 例: 3-4-7.geojson
  fname <- tools::file_path_sans_ext(f)
  
  # c("3","4","7")
  parts <- strsplit(fname, "-")[[1]]
  
  animal_num <- parts[1]
  region_num <- parts[2]
  roi_num <- parts[3]
  
  region_name <- region_map[region_num]
  
  # Polygon object 名
  polygon_name <- paste0(
    "Polygon_ColdF_N",
    animal_num,
    "_",
    region_name,
    "_",
    roi_num
  )
  
  # sections 名
  section_name <- paste0("ColdF_",
    animal_num, "_",
    region_num, "_",
    roi_num
  )
  
  # polygon object を取得
  polygon_obj <- get(polygon_name)
  
  # N1 / N2 / N3 に応じて fov 名変更
  fov_name <- switch(
    animal_num,
    "1" = "fov.coldF",
    "2" = "fov.coldF",
    "3" = "fov.coldF"
  )
  
  # list(fov=...) を動的生成
  polygon_list <- setNames(
    list(polygon_obj),
    fov_name
  )
  
  # PolygonAddSections
  merged.obj_SPN_woQIHLPS_afterQC2 <- PolygonAddSections(
    merged.obj_SPN_woQIHLPS_afterQC2,
    polygons = polygon_list,
    pixel_size = 0.2125,
    sections = section_name
  )
  
  cat("Added:", section_name, "\n")
}



# 2DGF
# ============================================
# PolygonAddSections を自動化
# ============================================

region_map <- c(
  "1" = "T1_4",
  "2" = "T5_7",
  "3" = "T8_10",
  "4" = "T11_13",
  "5" = "L1_3",
  "6" = "L45",
  "7" = "L6S1",
  "8" = "S2_4"
)

# geojson ファイル一覧
files <- list.files(
  "./Xenium_SC/ROI_geojson/2DGF",
  pattern = "\\.geojson$",
  full.names = FALSE
)

# ファイル名順に並べる
files <- sort(files)

for (f in files) {
  
  # 例: 3-4-7.geojson
  fname <- tools::file_path_sans_ext(f)
  
  # c("3","4","7")
  parts <- strsplit(fname, "-")[[1]]
  
  animal_num <- parts[1]
  region_num <- parts[2]
  roi_num <- parts[3]
  
  region_name <- region_map[region_num]
  
  # Polygon object 名
  polygon_name <- paste0(
    "Polygon_2DGF_N",
    animal_num,
    "_",
    region_name,
    "_",
    roi_num
  )
  
  # sections 名
  section_name <- paste0("2DGF_",
    animal_num, "_",
    region_num, "_",
    roi_num
  )
  
  # polygon object を取得
  polygon_obj <- get(polygon_name)
  
  # N1 / N2 / N3 に応じて fov 名変更
  fov_name <- switch(
    animal_num,
    "1" = "fov.2DGF",
    "2" = "fov.2DGF",
    "3" = "fov.2DGF"
  )
  
  # list(fov=...) を動的生成
  polygon_list <- setNames(
    list(polygon_obj),
    fov_name
  )
  
  # PolygonAddSections
  merged.obj_SPN_woQIHLPS_afterQC2 <- PolygonAddSections(
    merged.obj_SPN_woQIHLPS_afterQC2,
    polygons = polygon_list,
    pixel_size = 0.2125,
    sections = section_name
  )
  
  cat("Added:", section_name, "\n")
}

# 2DGM
# ============================================
# PolygonAddSections を自動化
# ============================================

region_map <- c(
  "1" = "T1_4",
  "2" = "T5_7",
  "3" = "T8_10",
  "4" = "T11_13",
  "5" = "L1_3",
  "6" = "L45",
  "7" = "L6S1",
  "8" = "S2_4"
)

# geojson ファイル一覧
files <- list.files(
  "./Xenium_SC/ROI_geojson/2DGM",
  pattern = "\\.geojson$",
  full.names = FALSE
)

# ファイル名順に並べる
files <- sort(files)

for (f in files) {
  
  # 例: 3-4-7.geojson
  fname <- tools::file_path_sans_ext(f)
  
  # c("3","4","7")
  parts <- strsplit(fname, "-")[[1]]
  
  animal_num <- parts[1]
  region_num <- parts[2]
  roi_num <- parts[3]
  
  region_name <- region_map[region_num]
  
  # Polygon object 名
  polygon_name <- paste0(
    "Polygon_2DGM_N",
    animal_num,
    "_",
    region_name,
    "_",
    roi_num
  )
  
  # sections 名
  section_name <- paste0("2DGM_",
    animal_num, "_",
    region_num, "_",
    roi_num
  )
  
  # polygon object を取得
  polygon_obj <- get(polygon_name)
  
  # N1 / N2 / N3 に応じて fov 名変更
  fov_name <- switch(
    animal_num,
    "1" = "fov.2DGM",
    "2" = "fov.2DGM",
    "3" = "fov.2DGM"
  )
  
  # list(fov=...) を動的生成
  polygon_list <- setNames(
    list(polygon_obj),
    fov_name
  )
  
  # PolygonAddSections
  merged.obj_SPN_woQIHLPS_afterQC2 <- PolygonAddSections(
    merged.obj_SPN_woQIHLPS_afterQC2,
    polygons = polygon_list,
    pixel_size = 0.2125,
    sections = section_name
  )
  
  cat("Added:", section_name, "\n")
}




#segmentを更新
merged.obj_SPN_woQIHLPS_afterQC2$segments <- NA


# セクションリスト
target_sections_T1_4 <- c("2DGF_1_1_2","2DGF_1_1_3","2DGF_1_1_4","2DGF_1_1_5","2DGF_1_1_6","2DGF_1_1_7","2DGF_1_1_8","2DGF_1_1_9","2DGF_2_1_10","2DGF_2_1_3","2DGF_2_1_4","2DGF_2_1_5","2DGF_2_1_6","2DGF_2_1_7","2DGF_2_1_8","2DGF_2_1_9","2DGF_3_1_2","2DGF_3_1_3","2DGF_3_1_4","2DGF_3_1_5","2DGF_3_1_6","2DGF_3_1_7","2DGF_3_1_8",
                          "2DGM_1_1_10","2DGM_1_1_4","2DGM_1_1_5","2DGM_1_1_6","2DGM_1_1_7","2DGM_1_1_8","2DGM_1_1_9","2DGM_2_1_10","2DGM_2_1_2","2DGM_2_1_3","2DGM_2_1_4","2DGM_2_1_5","2DGM_2_1_6","2DGM_2_1_7","2DGM_2_1_8","2DGM_2_1_9","2DGM_3_1_10","2DGM_3_1_3","2DGM_3_1_4","2DGM_3_1_5","2DGM_3_1_6","2DGM_3_1_7","2DGM_3_1_8","2DGM_3_1_9",
                          "ColdF_1_1_3","ColdF_1_1_4","ColdF_1_1_5","ColdF_1_1_6","ColdF_1_1_7","ColdF_1_1_8","ColdF_1_1_9","ColdF_2_1_3","ColdF_2_1_4","ColdF_2_1_5","ColdF_2_1_6","ColdF_2_1_7","ColdF_2_1_8","ColdF_2_1_9","ColdF_3_1_1","ColdF_3_1_2","ColdF_3_1_3","ColdF_3_1_4","ColdF_3_1_5","ColdF_3_1_6","ColdF_3_1_7","ColdF_3_1_8",
                          "ColdM_1_1_1","ColdM_1_1_2","ColdM_1_1_3","ColdM_1_1_4","ColdM_1_1_5","ColdM_1_1_6","ColdM_2_1_2","ColdM_2_1_3","ColdM_2_1_4","ColdM_2_1_5","ColdM_2_1_6","ColdM_3_1_1","ColdM_3_1_2","ColdM_3_1_3",
                          "CtrlF_1_1_3","CtrlF_1_1_4","CtrlF_1_1_5","CtrlF_1_1_6","CtrlF_1_1_7","CtrlF_1_1_8","CtrlF_1_1_9","CtrlF_1_1_10","CtrlF_2_1_4","CtrlF_2_1_5","CtrlF_2_1_6","CtrlF_2_1_7","CtrlF_2_1_8","CtrlF_2_1_9","CtrlF_2_1_10","CtrlF_2_1_11","CtrlF_3_1_1","CtrlF_3_1_4","CtrlF_3_1_5","CtrlF_3_1_6","CtrlF_3_1_7","CtrlF_3_1_8","CtrlF_3_1_9","CtrlF_3_1_10","CtrlF_3_1_11",
                          "CtrlM_1_1_1", "CtrlM_1_1_2", "CtrlM_1_1_3", "CtrlM_2_1_3", "CtrlM_2_1_4", "CtrlM_2_1_5", "CtrlM_2_1_6", "CtrlM_3_1_3", "CtrlM_3_1_4", "CtrlM_3_1_5", "CtrlM_3_1_6")

target_sections_T5_7 <- c("2DGF_1_2_1","2DGF_1_2_2","2DGF_1_2_3","2DGF_1_2_4","2DGF_1_2_5","2DGF_1_2_6","2DGF_1_2_7","2DGF_1_2_8","2DGF_1_2_9","2DGF_2_2_1","2DGF_2_2_2","2DGF_2_2_3","2DGF_2_2_4","2DGF_2_2_5","2DGF_2_2_6","2DGF_2_2_7","2DGF_2_2_8","2DGF_3_2_1","2DGF_3_2_2","2DGF_3_2_3","2DGF_3_2_4","2DGF_3_2_5","2DGF_3_2_6","2DGF_3_2_7","2DGF_3_2_8","2DGF_3_2_9",
                          "2DGM_1_2_1","2DGM_1_2_2","2DGM_1_2_3","2DGM_1_2_4","2DGM_1_2_5","2DGM_1_2_6","2DGM_1_2_7","2DGM_1_2_8","2DGM_1_2_9","2DGM_2_2_1","2DGM_2_2_2","2DGM_2_2_3","2DGM_2_2_4","2DGM_2_2_5","2DGM_2_2_6","2DGM_2_2_7","2DGM_2_2_8","2DGM_3_2_1","2DGM_3_2_10","2DGM_3_2_2","2DGM_3_2_3","2DGM_3_2_4","2DGM_3_2_5","2DGM_3_2_6","2DGM_3_2_7","2DGM_3_2_8","2DGM_3_2_9",
                          "ColdF_1_2_1","ColdF_1_2_2","ColdF_1_2_3","ColdF_1_2_4","ColdF_1_2_5","ColdF_1_2_6","ColdF_1_2_7","ColdF_1_2_8","ColdF_2_2_1","ColdF_2_2_2","ColdF_2_2_3","ColdF_2_2_4","ColdF_2_2_5","ColdF_2_2_6","ColdF_2_2_7","ColdF_3_2_2","ColdF_3_2_4","ColdF_3_2_5","ColdF_3_2_7","ColdF_3_2_8","ColdF_3_2_1","ColdF_3_2_3","ColdF_3_2_6",
                          "ColdM_1_2_1","ColdM_1_2_2","ColdM_1_2_3","ColdM_1_2_4","ColdM_1_2_5","ColdM_1_2_6","ColdM_3_2_1","ColdM_3_2_2","ColdM_3_2_3","ColdM_3_2_4","ColdM_3_2_5","ColdM_3_2_6",
                          "CtrlF_1_2_1", "CtrlF_1_2_2", "CtrlF_1_2_3", "CtrlF_1_2_4", "CtrlF_1_2_5", "CtrlF_1_2_6", "CtrlF_1_2_7", "CtrlF_1_2_8", "CtrlF_2_2_1", "CtrlF_2_2_2", "CtrlF_2_2_3", "CtrlF_2_2_4", "CtrlF_2_2_5", "CtrlF_2_2_6", "CtrlF_2_2_7", "CtrlF_2_2_8", "CtrlF_2_2_9", "CtrlF_3_2_1", "CtrlF_3_2_2", "CtrlF_3_2_3", "CtrlF_3_2_4", "CtrlF_3_2_5", "CtrlF_3_2_6", "CtrlF_3_2_7", "CtrlF_3_2_8",
                          "CtrlM_1_2_1", "CtrlM_1_2_2", "CtrlM_1_2_3", "CtrlM_1_2_4", "CtrlM_2_2_1", "CtrlM_2_2_2", "CtrlM_2_2_3", "CtrlM_2_2_4", "CtrlM_2_2_5", "CtrlM_2_2_6", "CtrlM_3_2_1", "CtrlM_3_2_2", "CtrlM_3_2_3", "CtrlM_3_2_4", "CtrlM_3_2_5", "CtrlM_3_2_6"
                          )
target_sections_T8_10 <- c("2DGF_1_3_1","2DGF_1_3_2","2DGF_1_3_3","2DGF_1_3_4","2DGF_1_3_5","2DGF_1_3_6","2DGF_1_3_7","2DGF_1_3_8","2DGF_2_3_1","2DGF_2_3_2","2DGF_2_3_3","2DGF_2_3_4","2DGF_2_3_5","2DGF_2_3_6","2DGF_2_3_7","2DGF_3_3_1","2DGF_3_3_2","2DGF_3_3_3","2DGF_3_3_4","2DGF_3_3_5","2DGF_3_3_6","2DGF_3_3_7","2DGF_3_3_8","2DGF_3_3_9",
                           "2DGM_1_3_1","2DGM_1_3_2","2DGM_1_3_3","2DGM_1_3_4","2DGM_1_3_5","2DGM_1_3_6","2DGM_1_3_7","2DGM_1_3_8","2DGM_1_3_9","2DGM_2_3_1","2DGM_2_3_2","2DGM_2_3_3","2DGM_2_3_4","2DGM_2_3_5","2DGM_2_3_6","2DGM_2_3_7","2DGM_2_3_8","2DGM_3_3_1","2DGM_3_3_10","2DGM_3_3_2","2DGM_3_3_3","2DGM_3_3_4","2DGM_3_3_5","2DGM_3_3_6","2DGM_3_3_7","2DGM_3_3_8","2DGM_3_3_9",
                           "ColdF_1_3_1","ColdF_1_3_2","ColdF_1_3_3","ColdF_1_3_4","ColdF_1_3_5","ColdF_1_3_6","ColdF_1_3_7","ColdF_1_3_8","ColdF_2_3_1","ColdF_2_3_2","ColdF_2_3_3","ColdF_2_3_4","ColdF_2_3_5","ColdF_2_3_6","ColdF_2_3_7","ColdF_2_3_8","ColdF_3_3_1","ColdF_3_3_3","ColdF_3_3_5","ColdF_3_3_7","ColdF_3_3_2","ColdF_3_3_4","ColdF_3_3_6",
                           "ColdM_1_3_1","ColdM_1_3_2","ColdM_1_3_3","ColdM_1_3_4","ColdM_1_3_5","ColdM_1_3_6","ColdM_2_3_1","ColdM_2_3_2","ColdM_2_3_3","ColdM_2_3_4","ColdM_2_3_5","ColdM_2_3_6","ColdM_3_3_1","ColdM_3_3_2","ColdM_3_3_3",
                           "CtrlF_1_3_1","CtrlF_1_3_2","CtrlF_1_3_3","CtrlF_1_3_4","CtrlF_1_3_5","CtrlF_1_3_6","CtrlF_1_3_7","CtrlF_1_3_8","CtrlF_1_3_9","CtrlF_2_3_1","CtrlF_2_3_2","CtrlF_2_3_3","CtrlF_2_3_4","CtrlF_2_3_5","CtrlF_2_3_6","CtrlF_2_3_7","CtrlF_2_3_8","CtrlF_2_3_9","CtrlF_3_3_1","CtrlF_3_3_2","CtrlF_3_3_3","CtrlF_3_3_4","CtrlF_3_3_5","CtrlF_3_3_6","CtrlF_3_3_7","CtrlF_3_3_8",
                           "CtrlM_1_3_1", "CtrlM_1_3_2", "CtrlM_1_3_3", "CtrlM_1_3_4", "CtrlM_2_3_1", "CtrlM_2_3_2", "CtrlM_2_3_3", "CtrlM_2_3_4", "CtrlM_2_3_5", "CtrlM_2_3_6", "CtrlM_2_3_7", "CtrlM_2_3_8", "CtrlM_3_3_1", "CtrlM_3_3_2", "CtrlM_3_3_3", "CtrlM_3_3_4", "CtrlM_3_3_5", "CtrlM_3_3_6"
                           
                           
)
target_sections_T11_13 <- c("2DGF_1_4_1","2DGF_1_4_2","2DGF_1_4_3","2DGF_1_4_4","2DGF_1_4_5","2DGF_1_4_6","2DGF_1_4_7","2DGF_1_4_8","2DGF_1_4_9","2DGF_2_4_1","2DGF_2_4_2","2DGF_2_4_3","2DGF_2_4_4","2DGF_2_4_5","2DGF_2_4_6","2DGF_2_4_7","2DGF_2_4_8","2DGF_3_4_1","2DGF_3_4_2","2DGF_3_4_3","2DGF_3_4_4","2DGF_3_4_5","2DGF_3_4_6","2DGF_3_4_7","2DGF_3_4_8",
                            "2DGM_1_4_1","2DGM_1_4_2","2DGM_1_4_3","2DGM_1_4_4","2DGM_1_4_5","2DGM_1_4_6","2DGM_1_4_7","2DGM_1_4_8","2DGM_2_4_1","2DGM_2_4_10","2DGM_2_4_2","2DGM_2_4_3","2DGM_2_4_4","2DGM_2_4_5","2DGM_2_4_6","2DGM_2_4_7","2DGM_2_4_8","2DGM_2_4_9","2DGM_3_4_1","2DGM_3_4_2","2DGM_3_4_3","2DGM_3_4_4","2DGM_3_4_5","2DGM_3_4_6","2DGM_3_4_7","2DGM_3_4_8",
                            "ColdF_1_4_1","ColdF_1_4_3","ColdF_1_4_5","ColdF_1_4_6","ColdF_1_4_7","ColdF_1_4_8","ColdF_2_4_1","ColdF_2_4_2","ColdF_2_4_4","ColdF_2_4_6","ColdF_2_4_8","ColdF_2_4_3","ColdF_2_4_5","ColdF_2_4_7","ColdF_3_4_1","ColdF_3_4_2","ColdF_3_4_3","ColdF_3_4_4","ColdF_3_4_5","ColdF_3_4_6","ColdF_3_4_7","ColdF_1_4_2","ColdF_1_4_4",
                            "ColdM_1_4_1","ColdM_1_4_2","ColdM_1_4_3","ColdM_1_4_4","ColdM_1_4_5","ColdM_1_4_6","ColdM_1_4_7","ColdM_1_4_8","ColdM_2_4_1","ColdM_2_4_2","ColdM_2_4_3","ColdM_2_4_4","ColdM_2_4_5","ColdM_2_4_6","ColdM_2_4_7","ColdM_2_4_8","ColdM_3_4_1","ColdM_3_4_2","ColdM_3_4_3",
                            "CtrlF_1_4_1","CtrlF_1_4_2","CtrlF_1_4_3","CtrlF_1_4_4","CtrlF_1_4_5","CtrlF_1_4_6","CtrlF_1_4_7","CtrlF_1_4_8","CtrlF_2_4_1","CtrlF_2_4_2","CtrlF_2_4_3","CtrlF_2_4_4","CtrlF_2_4_5","CtrlF_2_4_6","CtrlF_2_4_7","CtrlF_2_4_8","CtrlF_3_4_1","CtrlF_3_4_2","CtrlF_3_4_3","CtrlF_3_4_4","CtrlF_3_4_5","CtrlF_3_4_6","CtrlF_3_4_7","CtrlF_3_4_8",
                            "CtrlM_1_4_1", "CtrlM_1_4_2", "CtrlM_1_4_3", "CtrlM_1_4_4", "CtrlM_1_4_5", "CtrlM_2_4_1", "CtrlM_2_4_2", "CtrlM_2_4_3", "CtrlM_2_4_4", "CtrlM_2_4_5", "CtrlM_2_4_6", "CtrlM_2_4_7", "CtrlM_2_4_8", "CtrlM_3_4_1", "CtrlM_3_4_2", "CtrlM_3_4_3", "CtrlM_3_4_4", "CtrlM_3_4_5", "CtrlM_3_4_6", "CtrlM_3_4_7"
)
target_sections_L1_3 <- c("2DGF_1_5_1","2DGF_1_5_10","2DGF_1_5_2","2DGF_1_5_3","2DGF_1_5_4","2DGF_1_5_5","2DGF_1_5_6","2DGF_1_5_7","2DGF_1_5_8","2DGF_1_5_9","2DGF_2_5_1","2DGF_2_5_2","2DGF_2_5_3","2DGF_2_5_4","2DGF_2_5_5","2DGF_2_5_6","2DGF_2_5_7","2DGF_2_5_8","2DGF_3_5_1","2DGF_3_5_2","2DGF_3_5_3","2DGF_3_5_4","2DGF_3_5_5","2DGF_3_5_6","2DGF_3_5_7","2DGF_3_5_8",
                          "2DGM_1_5_1","2DGM_1_5_10","2DGM_1_5_2","2DGM_1_5_3","2DGM_1_5_4","2DGM_1_5_5","2DGM_1_5_6","2DGM_1_5_7","2DGM_1_5_8","2DGM_1_5_9","2DGM_2_5_1","2DGM_2_5_10","2DGM_2_5_2","2DGM_2_5_3","2DGM_2_5_4","2DGM_2_5_5","2DGM_2_5_6","2DGM_2_5_7","2DGM_2_5_8","2DGM_2_5_9","2DGM_3_5_1","2DGM_3_5_2","2DGM_3_5_3","2DGM_3_5_4","2DGM_3_5_5","2DGM_3_5_6","2DGM_3_5_7","2DGM_3_5_8",
                          "ColdF_1_5_1","ColdF_1_5_3","ColdF_1_5_5","ColdF_1_5_7","ColdF_1_5_9","ColdF_1_5_10",
                          "ColdF_1_5_2",
                          "ColdF_1_5_4",
                          "ColdF_1_5_6",
                          "ColdF_1_5_8",
                          "ColdF_2_5_1","ColdF_2_5_3","ColdF_2_5_5","ColdF_2_5_7","ColdF_2_5_8","ColdF_2_5_9",
                          "ColdF_2_5_2",
                          "ColdF_2_5_4",
                          "ColdF_2_5_6",
                          "ColdF_2_5_10",
                          "ColdF_3_5_2","ColdF_3_5_4","ColdF_3_5_5","ColdF_3_5_6","ColdF_3_5_8","ColdF_3_5_1",
                          "ColdF_3_5_3",
                          "ColdF_3_5_7",
                          "ColdM_1_5_1","ColdM_1_5_2","ColdM_1_5_3","ColdM_1_5_4","ColdM_1_5_5","ColdM_1_5_6","ColdM_2_5_1","ColdM_2_5_2","ColdM_2_5_3","ColdM_2_5_4","ColdM_2_5_5","ColdM_2_5_6","ColdM_3_5_1","ColdM_3_5_2","ColdM_3_5_3",
                          "CtrlF_1_5_1","CtrlF_1_5_2","CtrlF_1_5_3","CtrlF_1_5_4","CtrlF_1_5_5","CtrlF_1_5_6","CtrlF_1_5_7","CtrlF_1_5_8","CtrlF_1_5_9","CtrlF_2_5_1","CtrlF_2_5_2","CtrlF_2_5_3","CtrlF_2_5_4","CtrlF_2_5_5","CtrlF_2_5_6","CtrlF_2_5_7","CtrlF_2_5_8","CtrlF_2_5_9","CtrlF_3_5_1","CtrlF_3_5_2","CtrlF_3_5_3","CtrlF_3_5_4","CtrlF_3_5_5","CtrlF_3_5_6","CtrlF_3_5_7","CtrlF_3_5_8",
                          "CtrlM_1_5_1", "CtrlM_1_5_2", "CtrlM_1_5_3", "CtrlM_1_5_4", "CtrlM_1_5_5", "CtrlM_1_5_6", "CtrlM_2_5_1", "CtrlM_2_5_2", "CtrlM_2_5_3", "CtrlM_2_5_4", "CtrlM_2_5_5", "CtrlM_2_5_6", "CtrlM_2_5_7", "CtrlM_2_5_8", "CtrlM_3_5_1", "CtrlM_3_5_10", "CtrlM_3_5_11", "CtrlM_3_5_12", "CtrlM_3_5_13", "CtrlM_3_5_14", "CtrlM_3_5_15", "CtrlM_3_5_2", "CtrlM_3_5_3", "CtrlM_3_5_4", "CtrlM_3_5_5", "CtrlM_3_5_6", "CtrlM_3_5_7", "CtrlM_3_5_8", "CtrlM_3_5_9"
)
target_sections_L45 <- c("ColdF_1_6_1","ColdF_2_6_4","ColdM_2_6_1","ColdM_2_6_2","ColdM_3_6_4","CtrlF_1_6_1","CtrlM_1_6_1"
                         

)
target_sections_L6S1 <- c("2DGF_1_7_6","2DGF_2_7_6","2DGF_2_7_7","2DGF_2_7_8","2DGF_3_7_6","2DGF_3_7_7","2DGM_1_7_5","2DGM_1_7_7","2DGM_2_7_6","ColdF_1_7_5","ColdF_1_7_6","ColdF_3_7_4","ColdF_3_7_5","ColdF_3_7_6","ColdM_3_7_3",
                          "CtrlF_1_7_1","CtrlF_2_7_5","CtrlF_2_7_6","CtrlF_3_7_6","CtrlM_1_7_3", "CtrlM_2_7_5", "CtrlM_2_7_6", "CtrlM_3_7_4", "CtrlM_3_7_5", "CtrlM_3_7_6"
                          

                          
                          )
target_sections_S2_4 <- c("2DGF_1_8_1","2DGF_1_8_2","2DGF_1_8_3","2DGF_1_8_4","2DGF_1_8_7","2DGF_1_8_8","2DGF_2_8_5","2DGF_2_8_6","2DGF_2_8_7","2DGF_2_8_8","2DGF_2_8_9","2DGF_3_8_1","2DGF_3_8_2","2DGF_3_8_3",
                          "2DGM_1_8_1","2DGM_1_8_2","2DGM_1_8_3","2DGM_1_8_4","2DGM_2_8_3","2DGM_2_8_4","2DGM_2_8_5","2DGM_2_8_6","2DGM_2_8_7",
                          "ColdF_1_8_1","ColdF_1_8_2","ColdF_1_8_6","ColdF_2_8_2","ColdF_2_8_3","ColdF_2_8_4","ColdF_2_8_5",
                          "ColdF_3_8_1","ColdF_3_8_2","ColdM_1_8_1","ColdM_1_8_2","ColdM_1_8_3","ColdM_1_8_4","ColdM_1_8_5","ColdM_2_8_4","ColdM_2_8_2","ColdM_2_8_3","ColdM_2_8_5","ColdM_2_8_6",
                          "ColdM_3_8_1","ColdM_3_8_2",
                          "CtrlF_1_8_2","CtrlF_1_8_3","CtrlF_1_8_4","CtrlF_1_8_5","CtrlF_2_8_1","CtrlF_2_8_2","CtrlF_2_8_3","CtrlF_2_8_4","CtrlF_2_8_5","CtrlF_3_8_1","CtrlF_3_8_2","CtrlF_3_8_3","CtrlF_3_8_4","CtrlF_3_8_5","CtrlF_3_8_6",
                          "CtrlM_1_8_1")


# segments 列に "T1-4" を代入
merged.obj_SPN_woQIHLPS_afterQC2$segments[
  merged.obj_SPN_woQIHLPS_afterQC2$sections %in% target_sections_T1_4
] <- "T1-4"
merged.obj_SPN_woQIHLPS_afterQC2$segments[
  merged.obj_SPN_woQIHLPS_afterQC2$sections %in% target_sections_T5_7
] <- "T5-7"
merged.obj_SPN_woQIHLPS_afterQC2$segments[
  merged.obj_SPN_woQIHLPS_afterQC2$sections %in% target_sections_T8_10
] <- "T8-10"
merged.obj_SPN_woQIHLPS_afterQC2$segments[
  merged.obj_SPN_woQIHLPS_afterQC2$sections %in% target_sections_T11_13
] <- "T11-13"
merged.obj_SPN_woQIHLPS_afterQC2$segments[
  merged.obj_SPN_woQIHLPS_afterQC2$sections %in% target_sections_L1_3
] <- "L1-3"
merged.obj_SPN_woQIHLPS_afterQC2$segments[
  merged.obj_SPN_woQIHLPS_afterQC2$sections %in% target_sections_L45
] <- "L45"
merged.obj_SPN_woQIHLPS_afterQC2$segments[
  merged.obj_SPN_woQIHLPS_afterQC2$sections %in% target_sections_L6S1
] <- "L6S1"
merged.obj_SPN_woQIHLPS_afterQC2$segments[
  merged.obj_SPN_woQIHLPS_afterQC2$sections %in% target_sections_S2_4
] <- "S2-4"

# 確認
table(merged.obj_SPN_woQIHLPS_afterQC2@meta.data$segments)


# これまでに、sections列やsegments列がNAになっている細胞IDの確認
merged.obj_SPN_woQIHLPS_afterQC2@meta.data[is.na(merged.obj_SPN_woQIHLPS_afterQC2$segments), c("sections", "segments") ]


# Xenium Explorerを確認して、細胞ごとに直接segments, sectionsの値を入れる
# sections列を更新
merged.obj_SPN_woQIHLPS_afterQC2$sections[
  rownames(merged.obj_SPN_woQIHLPS_afterQC2@meta.data) == "2DGF_bfnbbeni-1"
] <- "2DGF_2_8_9"

# segments列を更新
merged.obj_SPN_woQIHLPS_afterQC2$segments[
  rownames(merged.obj_SPN_woQIHLPS_afterQC2@meta.data) == "2DGF_bfnbbeni-1"
] <- "S2-4"
# 確認
merged.obj_SPN_woQIHLPS_afterQC2@meta.data[
  "2DGF_bfnbbeni-1",
  c("sections", "segments")
]



