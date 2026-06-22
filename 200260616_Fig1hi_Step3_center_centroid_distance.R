# dfにQIHスライド上のSPNの各点の情報を、cellIDとともに格納

# すべてのfov名を取得
fov_list <- c( "fov","fov.Ctrl2", "fov.Ctrl3", "fov.2", "fov.Cold2", "fov.Cold3", "fov.2.3", "fov.coldF", "fov.2DGM", "fov.2DGF")


# 各fovごとにcell情報を抽出
df_list <- lapply(fov_list, function(fov) {
  coords <- merged.obj_SPN_woQIHLPS_afterQC3@images[[fov]]@boundaries[["centroids"]]@coords
  cells  <- merged.obj_SPN_woQIHLPS_afterQC3@images[[fov]]@boundaries[["centroids"]]@cells
  data.frame(
    Cell = cells,
    coords
  )
})

# 全FOVをまとめる
df <- do.call(rbind, df_list)
df$distance <- NA  # まず全行 NA

# 
# df <- data.frame(merged.obj_SPN_woQIHLPS_afterQC3@images[["fov.LPS"]]@boundaries[["centroids"]]@coords)
# df <- data.frame(Cell = merged.obj_SPN_woQIHLPS_afterQC3@images[["fov.LPS"]]@boundaries[["centroids"]]@cells, df)

# geometry カラムから座標を抽出
# N1
# coords_df <- as.data.frame(st_coordinates(Center_QIH_N1_T1_4_1))
# cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
#   !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
#     merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == "1_1_1", ])
# 
# # distance 列を作り、対象行のみ計算
# df$distance[df$Cell %in% cells] <- sqrt((df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 + 
#                                           (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2)


# N1_T1-4
sections <- paste0("1_1_", 1:9)
centers <- paste0("Center_2DGF_N1_T1_4_", 1:9)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_T5-7
sections <- paste0("1_2_", 1:9)
centers <- paste0("Center_2DGF_N1_T5_7_", 1:9)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_T8-10
sections <- paste0("1_3_", 1:8)
centers <- paste0("Center_2DGF_N1_T8_10_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_T11-13
sections <- paste0("1_4_", 1:9)
centers <- paste0("Center_2DGF_N1_T11_13_", 1:9)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_L1-3
sections <- paste0("1_5_", 1:10)
centers <- paste0("Center_2DGF_N1_L1_3_", 1:10)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells]- coords_df$Y*0.2125)^2
  )
}
# N1_L45
sections <- paste0("1_6_", 1:6)
centers <- paste0("Center_2DGF_N1_L45_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N1_L6S1
sections <- paste0("1_7_", 1:8)
centers <- paste0("Center_2DGF_N1_L6S1_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_S2-4
sections <- paste0("1_8_", 1:8)
centers <- paste0("Center_2DGF_N1_S2_4_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_T1-4
sections <- paste0("2_1_", 1:10)
centers <- paste0("Center_2DGF_N2_T1_4_", 1:10)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_T5-7
sections <- paste0("2_2_", 1:8)
centers <- paste0("Center_2DGF_N2_T5_7_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_T8-10
sections <- paste0("2_3_", 1:7)
centers <- paste0("Center_2DGF_N2_T8_10_", 1:7)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_T11-13
sections <- paste0("2_4_", 1:8)
centers <- paste0("Center_2DGF_N2_T11_13_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N2_L1-3
sections <- paste0("2_5_", 1:9)
centers <- paste0("Center_2DGF_N2_L1_3_", 1:9)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N2_L45
sections <- paste0("2_6_", 1:6)
centers <- paste0("Center_2DGF_N2_L45_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N2_L6S1
sections <- paste0("2_7_", 1:8)
centers <- paste0("Center_2DGF_N2_L6S1_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N2_S2-4
sections <- paste0("2_8_", 1:9)
centers <- paste0("Center_2DGF_N2_S2_4_", 1:9)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T1-4
sections <- paste0("3_1_", 1:8)
centers <- paste0("Center_2DGF_N3_T1_4_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T5-7
sections <- paste0("3_2_", 1:9)
centers <- paste0("Center_2DGF_N3_T5_7_", 1:9)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T8_10
sections <- paste0("3_3_", 1:9)
centers <- paste0("Center_2DGF_N3_T8_10_", 1:9)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T11_13
sections <- paste0("3_4_", 1:8)
centers <- paste0("Center_2DGF_N3_T11_13_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_L1_3
sections <- paste0("3_5_", 1:8)
centers <- paste0("Center_2DGF_N3_L1_3_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_L45
sections <- paste0("3_6_", 1:6)
centers <- paste0("Center_2DGF_N3_L45_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_L6S1
sections <- paste0("3_7_", 1:7)
centers <- paste0("Center_2DGF_N3_L6S1_", 1:7)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_S2_4
sections <- paste0("3_8_", 1:7)
centers <- paste0("Center_2DGF_N3_S2_4_", 1:7)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGF", ])
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_T1-4
sections <- paste0("1_1_", 1:10)
centers <- paste0("Center_2DGM_N1_T1_4_", 1:10)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_T5-7
sections <- paste0("1_2_", 1:9)
centers <- paste0("Center_2DGM_N1_T5_7_", 1:9)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_T8-10
sections <- paste0("1_3_", 1:9)
centers <- paste0("Center_2DGM_N1_T8_10_", 1:9)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_T11-13
sections <- paste0("1_4_", 1:8)
centers <- paste0("Center_2DGM_N1_T11_13_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_L1-3
sections <- paste0("1_5_", 1:10)
centers <- paste0("Center_2DGM_N1_L1_3_", 1:10)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells]- coords_df$Y*0.2125)^2
  )
}
# N1_L45
sections <- paste0("1_6_", 1:6)
centers <- paste0("Center_2DGM_N1_L45_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N1_L6S1
sections <- paste0("1_7_", 1:7)
centers <- paste0("Center_2DGM_N1_L6S1_", 1:7)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_S2-4
sections <- paste0("1_8_", 1:5)
centers <- paste0("Center_2DGM_N1_S2_4_", 1:5)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_T1-4
sections <- paste0("2_1_", 1:10)
centers <- paste0("Center_2DGM_N2_T1_4_", 1:10)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_T5-7
sections <- paste0("2_2_", 1:8)
centers <- paste0("Center_2DGM_N2_T5_7_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_T8-10
sections <- paste0("2_3_", 1:8)
centers <- paste0("Center_2DGM_N2_T8_10_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_T11-13
sections <- paste0("2_4_", 1:10)
centers <- paste0("Center_2DGM_N2_T11_13_", 1:10)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N2_L1-3
sections <- paste0("2_5_", 1:10)
centers <- paste0("Center_2DGM_N2_L1_3_", 1:10)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N2_L45
sections <- paste0("2_6_", 1:5)
centers <- paste0("Center_2DGM_N2_L45_", 1:5)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N2_L6S1
sections <- paste0("2_7_", 1:7)
centers <- paste0("Center_2DGM_N2_L6S1_", 1:7)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_S2-4
sections <- paste0("2_8_", 1:7)
centers <- paste0("Center_2DGM_N2_S2_4_", 1:7)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T1-4
sections <- paste0("3_1_", 1:10)
centers <- paste0("Center_2DGM_N3_T1_4_", 1:10)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T5-7
sections <- paste0("3_2_", 1:10)
centers <- paste0("Center_2DGM_N3_T5_7_", 1:10)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T8_10
sections <- paste0("3_3_", 1:10)
centers <- paste0("Center_2DGM_N3_T8_10_", 1:10)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T11_13
sections <- paste0("3_4_", 1:8)
centers <- paste0("Center_2DGM_N3_T11_13_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_L1_3
sections <- paste0("3_5_", 1:9)
centers <- paste0("Center_2DGM_N3_L1_3_", 1:9)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_L45
sections <- paste0("3_6_", 1:2)
centers <- paste0("Center_2DGM_N3_L45_", 1:2)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_L6S1
sections <- paste0("3_7_", 1:4)
centers <- paste0("Center_2DGM_N3_L6S1_", 1:4)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_S2_4
sections <- paste0("3_8_", 1:2)
centers <- paste0("Center_2DGM_N3_S2_4_", 1:2)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "2DGM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}


# CtrlF
# N1_T1-4
sections <- paste0("1_1_", 2:10)
centers <- paste0("Center_CtrlF_N1_T1_4_", 2:10)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_T5-7
sections <- paste0("1_2_", 1:8)
centers <- paste0("Center_CtrlF_N1_T5_7_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_T8-10
sections <- paste0("1_3_", 1:9)
centers <- paste0("Center_CtrlF_N1_T8_10_", 1:9)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_T11-13
sections <- paste0("1_4_", 1:8)
centers <- paste0("Center_CtrlF_N1_T11_13_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_L1-3
sections <- paste0("1_5_", 1:9)
centers <- paste0("Center_CtrlF_N1_L1_3_", 1:9)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells]- coords_df$Y*0.2125)^2
  )
}
# N1_L45
sections <- paste0("1_6_", 1:5)
centers <- paste0("Center_CtrlF_N1_L45_", 1:5)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N1_L6S1
sections <- paste0("1_7_", 1:6)
centers <- paste0("Center_CtrlF_N1_L6S1_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_S2-4
sections <- paste0("1_8_", 1:5)
centers <- paste0("Center_CtrlF_N1_S2_4_", 1:5)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_T1-4
sections <- paste0("2_1_", 1:11)
centers <- paste0("Center_CtrlF_N2_T1_4_", 1:11)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_T5-7
sections <- paste0("2_2_", 1:9)
centers <- paste0("Center_CtrlF_N2_T5_7_", 1:9)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_T8-10
sections <- paste0("2_3_", 1:9)
centers <- paste0("Center_CtrlF_N2_T8_10_", 1:9)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_T11-13
sections <- paste0("2_4_", 1:8)
centers <- paste0("Center_CtrlF_N2_T11_13_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N2_L1-3
sections <- paste0("2_5_", 1:9)
centers <- paste0("Center_CtrlF_N2_L1_3_", 1:9)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N2_L45
sections <- paste0("2_6_", 1:5)
centers <- paste0("Center_CtrlF_N2_L45_", 1:5)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N2_L6S1
sections <- paste0("2_7_", 1:6)
centers <- paste0("Center_CtrlF_N2_L6S1_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_S2-4
sections <- paste0("2_8_", 1:6)
centers <- paste0("Center_CtrlF_N2_S2_4_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T1-4
sections <- paste0("3_1_", 1:11)
centers <- paste0("Center_CtrlF_N3_T1_4_", 1:11)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T5-7
sections <- paste0("3_2_", 1:8)
centers <- paste0("Center_CtrlF_N3_T5_7_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T8_10
sections <- paste0("3_3_", 1:8)
centers <- paste0("Center_CtrlF_N3_T8_10_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T11_13
sections <- paste0("3_4_", 1:8)
centers <- paste0("Center_CtrlF_N3_T11_13_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_L1_3
sections <- paste0("3_5_", 1:9)
centers <- paste0("Center_CtrlF_N3_L1_3_", 1:9)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_L45
sections <- paste0("3_6_", 1:5)
centers <- paste0("Center_CtrlF_N3_L45_", 1:5)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_L6S1
sections <- paste0("3_7_", 1:6)
centers <- paste0("Center_CtrlF_N3_L6S1_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_S2_4
sections <- paste0("3_8_", 1:6)
centers <- paste0("Center_CtrlF_N3_S2_4_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}



# ColdF
# N1_T1-4
sections <- paste0("1_1_", 1:9)
centers <- paste0("Center_ColdF_N1_T1_4_", 1:9)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_T5-7
sections <- paste0("1_2_", 1:8)
centers <- paste0("Center_ColdF_N1_T5_7_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_T8-10
sections <- paste0("1_3_", 1:8)
centers <- paste0("Center_ColdF_N1_T8_10_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_T11-13
sections <- paste0("1_4_", 1:8)
centers <- paste0("Center_ColdF_N1_T11_13_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_L1-3
sections <- paste0("1_5_", 1:10)
centers <- paste0("Center_ColdF_N1_L1_3_", 1:10)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells]- coords_df$Y*0.2125)^2
  )
}
# N1_L45
sections <- paste0("1_6_", 1:6)
centers <- paste0("Center_ColdF_N1_L45_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N1_L6S1
sections <- paste0("1_7_", 1:6)
centers <- paste0("Center_ColdF_N1_L6S1_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_S2-4
sections <- paste0("1_8_", 1:5)
centers <- paste0("Center_ColdF_N1_S2_4_", 1:5)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_T1-4
sections <- paste0("2_1_", 1:9)
centers <- paste0("Center_ColdF_N2_T1_4_", 1:9)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_T5-7
sections <- paste0("2_2_", 1:7)
centers <- paste0("Center_ColdF_N2_T5_7_", 1:7)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_T8-10
sections <- paste0("2_3_", 1:8)
centers <- paste0("Center_ColdF_N2_T8_10_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_T11-13
sections <- paste0("2_4_", 1:8)
centers <- paste0("Center_ColdF_N2_T11_13_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N2_L1-3
sections <- paste0("2_5_", 1:10)
centers <- paste0("Center_ColdF_N2_L1_3_", 1:10)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N2_L45
sections <- paste0("2_6_", 1:5)
centers <- paste0("Center_ColdF_N2_L45_", 1:5)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N2_L6S1
sections <- paste0("2_7_", 1:6)
centers <- paste0("Center_ColdF_N2_L6S1_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_S2-4
sections <- paste0("2_8_", 1:5)
centers <- paste0("Center_ColdF_N2_S2_4_", 1:5)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T1-4
sections <- paste0("3_1_", 1:8)
centers <- paste0("Center_ColdF_N3_T1_4_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T5-7
sections <- paste0("3_2_", 1:8)
centers <- paste0("Center_ColdF_N3_T5_7_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T8_10
sections <- paste0("3_3_", 1:7)
centers <- paste0("Center_ColdF_N3_T8_10_", 1:7)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T11_13
sections <- paste0("3_4_", 1:7)
centers <- paste0("Center_ColdF_N3_T11_13_", 1:7)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_L1_3
sections <- paste0("3_5_", 1:9)
centers <- paste0("Center_ColdF_N3_L1_3_", 1:9)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_L45
sections <- paste0("3_6_", 1:5)
centers <- paste0("Center_ColdF_N3_L45_", 1:5)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_L6S1
sections <- paste0("3_7_", 1:6)
centers <- paste0("Center_ColdF_N3_L6S1_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_S2_4
sections <- paste0("3_8_", 1:5)
centers <- paste0("Center_ColdF_N3_S2_4_", 1:5)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdF", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}












# N1_T1-4
sections <- paste0("1_1_", 1:6)
centers <- paste0("Center_Cold_N1_T1_4_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_T5-7
sections <- paste0("1_2_", 1:6)
centers <- paste0("Center_Cold_N1_T5_7_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_T8-10
sections <- paste0("1_3_", 1:6)
centers <- paste0("Center_Cold_N1_T8_10_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_T11-13
sections <- paste0("1_4_", 1:8)
centers <- paste0("Center_Cold_N1_T11_13_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_L1-3
sections <- paste0("1_5_", 1:6)
centers <- paste0("Center_Cold_N1_L1_3_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells]- coords_df$Y*0.2125)^2
  )
}
# N1_L45
sections <- paste0("1_6_", 1:6)
centers <- paste0("Center_Cold_N1_L45_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N1_L6S1
sections <- paste0("1_7_", 1:6)
centers <- paste0("Center_Cold_N1_L6S1_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_S2-4
sections <- paste0("1_8_", 1:7)
centers <- paste0("Center_Cold_N1_S2_4_", 1:7)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_T1-4
sections <- paste0("2_1_", 1:6)
centers <- paste0("Center_Cold_N2_T1_4_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# # N2_T5-7
# sections <- paste0("2_2_", 1:7)
# centers <- paste0("Center_Cold_N2_T5_7_", 1:7)
# 
# for (i in seq_along(sections)) {
#   # 座標を取得
#   coords_df <- as.data.frame(st_coordinates(get(centers[i])))
#   
#   # 該当する cell ID を取得
#   cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
#     !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
#       merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
#       merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
#   
#   # distance 計算して df に追加
#   df$distance[df$Cell %in% cells] <- sqrt(
#     (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
#       (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
#   )
# }

# N2_T8-10
sections <- paste0("2_3_", 1:6)
centers <- paste0("Center_Cold_N2_T8_10_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_T11-13
sections <- paste0("2_4_", 1:8)
centers <- paste0("Center_Cold_N2_T11_13_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N2_L1-3
sections <- paste0("2_5_", 1:6)
centers <- paste0("Center_Cold_N2_L1_3_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N2_L45
sections <- paste0("2_6_", 1:6)
centers <- paste0("Center_Cold_N2_L45_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N2_L6S1
sections <- paste0("2_7_", 1:6)
centers <- paste0("Center_Cold_N2_L6S1_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N2_S2-4
sections <- paste0("2_8_", 1:7)
centers <- paste0("Center_Cold_N2_S2_4_", 1:7)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T1-4
sections <- paste0("3_1_", 1:3)
centers <- paste0("Center_Cold_N3_T1_4_", 1:3)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T5-7
sections <- paste0("3_2_", 1:6)
centers <- paste0("Center_Cold_N3_T5_7_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T8_10
sections <- paste0("3_3_", 1:3)
centers <- paste0("Center_Cold_N3_T8_10_", 1:3)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T11_13
sections <- paste0("3_4_", 1:3)
centers <- paste0("Center_Cold_N3_T11_13_", 1:3)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_L1_3
sections <- paste0("3_5_", 1:3)
centers <- paste0("Center_Cold_N3_L1_3_", 1:3)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_L45
sections <- paste0("3_6_", 1:5)
centers <- paste0("Center_Cold_N3_L45_", 1:5)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_L6S1
sections <- paste0("3_7_", 1:3)
centers <- paste0("Center_Cold_N3_L6S1_", 1:3)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_S2_4
sections <- paste0("3_8_", 1:3)
centers <- paste0("Center_Cold_N3_S2_4_", 1:3)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "ColdM", ])
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N1_T1-4
sections <- paste0("1_1_", 1:3)
centers <- paste0("Center_Ctrl_N1_T1_4_", 1:3)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_T5-7
sections <- paste0("1_2_", 1:4)
centers <- paste0("Center_Ctrl_N1_T5_7_", 1:4)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_T8-10
sections <- paste0("1_3_", 1:4)
centers <- paste0("Center_Ctrl_N1_T8_10_", 1:4)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_T11-13
sections <- paste0("1_4_", 1:5)
centers <- paste0("Center_Ctrl_N1_T11_13_", 1:5)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_L1-3
sections <- paste0("1_5_", 1:6)
centers <- paste0("Center_Ctrl_N1_L1_3_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells]- coords_df$Y*0.2125)^2
  )
}
# N1_L45
sections <- paste0("1_6_", 1:3)
centers <- paste0("Center_Ctrl_N1_L45_", 1:3)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N1_L6S1
sections <- paste0("1_7_", 1:3)
centers <- paste0("Center_Ctrl_N1_L6S1_", 1:3)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N1_S2-4
sections <- paste0("1_8_", 1:1)
centers <- paste0("Center_Ctrl_N1_S2_4_", 1:1)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_T1-4
sections <- paste0("2_1_", 1:6)
centers <- paste0("Center_Ctrl_N2_T1_4_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_T5-7
sections <- paste0("2_2_", 1:6)
centers <- paste0("Center_Ctrl_N2_T5_7_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_T8-10
sections <- paste0("2_3_", 1:8)
centers <- paste0("Center_Ctrl_N2_T8_10_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N2_T11-13
sections <- paste0("2_4_", 1:8)
centers <- paste0("Center_Ctrl_N2_T11_13_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N2_L1-3
sections <- paste0("2_5_", 1:8)
centers <- paste0("Center_Ctrl_N2_L1_3_", 1:8)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N2_L45
sections <- paste0("2_6_", 1:6)
centers <- paste0("Center_Ctrl_N2_L45_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N2_L6S1
sections <- paste0("2_7_", 1:6)
centers <- paste0("Center_Ctrl_N2_L6S1_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# N2_S2-4
sections <- paste0("2_8_", 1:5)
centers <- paste0("Center_Ctrl_N2_S2_4_", 1:5)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T1-4
sections <- paste0("3_1_", 1:6)
centers <- paste0("Center_Ctrl_N3_T1_4_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T5-7
sections <- paste0("3_2_", 1:6)
centers <- paste0("Center_Ctrl_N3_T5_7_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T8_10
sections <- paste0("3_3_", 1:6)
centers <- paste0("Center_Ctrl_N3_T8_10_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_T11_13
sections <- paste0("3_4_", 1:7)
centers <- paste0("Center_Ctrl_N3_T11_13_", 1:7)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_L1_3
sections <- paste0("3_5_", 1:16)
centers <- paste0("Center_Ctrl_N3_L1_3_", 1:16)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_L45
sections <- paste0("3_6_", 1:6)
centers <- paste0("Center_Ctrl_N3_L45_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}

# N3_L6S1
sections <- paste0("3_7_", 1:6)
centers <- paste0("Center_Ctrl_N3_L6S1_", 1:6)

for (i in seq_along(sections)) {
  # 座標を取得
  coords_df <- as.data.frame(st_coordinates(get(centers[i])))
  
  # 該当する cell ID を取得
  cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
      merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "CtrlM", ])
  
  # distance 計算して df に追加
  df$distance[df$Cell %in% cells] <- sqrt(
    (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
      (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
  )
}
# 
# # N3_S2_4
# sections <- paste0("3_8_", 1:6)
# centers <- paste0("Center_Ctrl_N3_S2_4_", 1:6)
# 
# for (i in seq_along(sections)) {
#   # 座標を取得
#   coords_df <- as.data.frame(st_coordinates(get(centers[i])))
#   
#   # 該当する cell ID を取得
#   cells <- rownames(merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
#     !is.na(merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections) &
#       merged.obj_SPN_woQIHLPS_afterQC3@meta.data$sections == sections[i] &
#       merged.obj_SPN_woQIHLPS_afterQC3@meta.data$ConditionSex == "COld", ])
#   # distance 計算して df に追加
#   df$distance[df$Cell %in% cells] <- sqrt(
#     (df$x[df$Cell %in% cells] - coords_df$X*0.2125)^2 +
#       (df$y[df$Cell %in% cells] - coords_df$Y*0.2125)^2
#   )
# }



ggplot(df, aes(x = "", y = df$distance)) +  # x軸は空文字で1本のバイオリンに
  geom_violin(fill = "skyblue", color = "black") +# バイオリン
  geom_jitter(width = 0.1, size = 1.5, color = "red")+
  labs(x = "", y = "x") +
  theme_minimal()



# Seuratオブジェクトのmeta.dataを取得
meta <- merged.obj_SPN_woQIHLPS_afterQC3@meta.data

# metaに新しい列を追加（NAで初期化）
meta$distance <- NA

# 行名とdf$Cellをマッチングして値を代入
meta$distance <- df$distance[match(rownames(meta), df$Cell)]

# 更新したmetaをSeuratオブジェクトに戻す
merged.obj_SPN_woQIHLPS_afterQC3@meta.data <- meta


# VlnPlotで長い距離が出ている細胞を要チェック
VlnPlot(merged.obj_SPN_woQIHLPS_afterQC3, features = "distance", split.by = "ConditionSex",group.by = "segments")
# ->(2026/6/17現在)あとColdFのT1-4のColdF_lfklaljk-1,ColdF_lfkkmjih-1は場所がおかしいのでSPNではない
# CtrlF_dccepcge-1,CtrlF_deiilmij-1もSPNではない


# NAとなっている細胞IDを確認する
na_cells <- rownames(
  merged.obj_SPN_woQIHLPS_afterQC3@meta.data[
    is.na(merged.obj_SPN_woQIHLPS_afterQC3$distance),
  ]
)
na_cells
# →基本的にsectionsでNAとなっている細胞と同じ

# 論文用Fig描画
VlnPlot(merged.obj_SPN_woQIHLPS_afterQC3, features = c("distance"),group.by="segments", pt.size = 0.1)
VlnPlot(merged.obj_SPN_woQIHLPS_afterQC3, features = c("distance"), pt.size = 0.1)

