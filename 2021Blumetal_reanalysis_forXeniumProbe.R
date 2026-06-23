library(Seurat)
library(SeuratObject)
library(ggplot2)
library(patchwork
        )
library(dplyr)
#以下はボルケーノプロット作成に必要
library(EnhancedVolcano)
library(gridExtra)
library(grid)
library(Seurat)
library(org.Mm.eg.db)
library(rvcheck)
library(clusterProfiler)
library(KEGGREST)
library(DOSE)
library(ReactomePA)
library(stringr)

library(EnhancedVolcano)
library(gridExtra)
library(grid)

library(ggplot2)
library(patchwork)

library(SeuratDisk)
library(scCustomize)

install.packages("anndata")
library(anndata)

install.packages("reticulate")
library(reticulate)
use_virtualenv("anndata")
anndata <- reticulate::import('anndata')

# h5ad形式に変換
# anndata package
tmp <- anndata::AnnData(
  X = Matrix::t(SC_merged_cholinergicneuron_Nos1@assays$SCT@counts),
  var = data.frame(row.names = rownames(SC_merged_cholinergicneuron_Nos1[["SCT"]])),
  obs = SC_merged_cholinergicneuron_Nos1@meta.data,
  layers = list(data = Matrix::t(SC_merged_cholinergicneuron_Nos1@assays$SCT@data)),
  obsm = list(
    rna.umap = SC_merged_cholinergicneuron_Nos1@reductions$umap@cell.embeddings
  )
)
anndata::write_h5ad(tmp, filename = "SC_merged_cholinergicneuron_Nos1.h5ad")

tmp2 <- anndata::AnnData(
  X = Matrix::t(SC_merged_cholinergicneuron_Nos1@assays$RNA@layers$counts),
  var = data.frame(row.names = rownames(SC_merged_cholinergicneuron_Nos1[["RNA"]])),
  obs = SC_merged_cholinergicneuron_Nos1@meta.data,
  obsm = list(
    rna.umap = SC_merged_cholinergicneuron_Nos1@reductions$umap@cell.embeddings
  )
)
anndata::write_h5ad(tmp2, filename = "SC_merged_cholinergicneuron_Nos1_RNA.h5ad")





as.anndata(x = SC_merged_cholinergicneuron_Nos1, file_path = "~/singlecellRNAseq/2021_Blum", file_name = "SC_merged_cholinergicneuron_Nos1.h5ad")
SC_merged_cholinergicneuron_Nos1[["RNA"]] <- as(object = SC_merged_cholinergicneuron_Nos1[["RNA"]], Class = "Assay")
SaveH5Seurat(SC_merged_cholinergicneuron_Nos1,"SC_merged_cholinergicneuron_Nos1.h5ad", overwrite = TRUE)
Convert("SC_merged_cholinergicneuron_Nos1.h5seurat", "SC_merged_cholinergicneuron_Nos1.h5ad", assay="RNA", overwrite = TRUE)
options(memory.limit = 32768 * 1024 * 1024) # 32GB のメモリ制限を設定

options(Seurat.object.assay.version = "v5")
class(Spnal1[["RNA"]])
#Xenium提出用のMexファイルを作る

# Define function
writeCounts <- function(out_dir, data, barcodes = colnames(data), gene.id = rownames(data), gene.symbol = rownames(data), feature.type = "Gene Expression", subset = 1:length(barcodes)) {
  require("R.utils")
  require("Matrix")
  
  if (file.exists(out_dir) || (dir.exists(out_dir) && length(list.files(out_dir)) > 0)) {
    stop("The specified output directory already exists! Not overwriting")
  }
  dir.create(out_dir, recursive = TRUE)
  
  if (require("data.table")) {
    data.table::fwrite(
      data.table::data.table(barcode = barcodes),
      file.path(out_dir, "barcodes.tsv.gz"),
      col.names = FALSE
    )
    
    data.table::fwrite(
      data.table::data.table(
        gene_id = gene.id,
        gene_symbol = gene.symbol,
        feature_type = feature.type
      ),
      file.path(out_dir, "features.tsv.gz"),
      col.names = FALSE
    )
  } else {
    write.table(
      data.frame(barcode = barcodes[subset]),
      gzfile(file.path(out_dir, "barcodes.tsv.gz")),
      sep = "\t", quote = FALSE,
      col.names = FALSE, row.names = FALSE
    )
    
    write.table(
      data.frame(
        gene_id = gene.id,
        gene_symbol = gene.symbol,
        feature_type = feature.type
      ),
      file.path(out_dir, "features.tsv.gz"),
      sep = "\t", quote = FALSE,
      col.names = FALSE, row.names = FALSE
    )
  }
  
  Matrix::writeMM(data[, subset], file.path(out_dir, "matrix.mtx"))
  R.utils::gzip(file.path(out_dir, "matrix.mtx"), remove = TRUE)
}

# Run function
writeCounts(
  "reference_data",
  GetAssayData(SC_merged, assay="RNA", slot ="counts"),
  gene.id = rownames(SC_merged),
  gene.symbol = rownames(SC_merged),
  barcodes = colnames(SC_merged),
  subset = subset
)

# List files
list.files("reference_data")

# Output
#[1] "barcodes.tsv.gz" "features.tsv.gz" "matrix.mtx.gz"



#まずは6バッチを順にロード
#Spnal1
Spnal1.data <- Read10X_h5("./singlecellRNAseq/2021_Blum/GSE161621_RAW/GSM4911289_mixed10_11_5.h5")

Spnal1<- CreateSeuratObject(counts = Spnal1.data, project = "Spinal1", min.cells = 3, min.features = 200)

Spnal1[["percent.mt"]] <- PercentageFeatureSet(Spnal1, pattern = "^mt-")
Spnal1[["percent.rp"]] <- PercentageFeatureSet(Spnal1, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")
VlnPlot(Spnal1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 3)
plot1 <- FeatureScatter(Spnal1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Spnal1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

Spnal1 <- subset(Spnal1, subset = nFeature_RNA > 400 & nFeature_RNA < 8000  & percent.mt < 10 & percent.rp < 10)

VlnPlot(Spnal1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2)
Spnal1　　<- Add_Alt_Feature_ID(
  Spnal1,
  hdf5_file = "./singlecellRNAseq/2021_Blum/GSE161621_RAW/GSM4911289_mixed10_11_5.h5",
  assay = 'RNA'
)

plot1 <- FeatureScatter(Spnal1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Spnal1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2




#Spnal2
Spnal2.data <- Read10X_h5("./singlecellRNAseq/2021_Blum/GSE161621_RAW/GSM4911290_mixed4_5_31_a.h5")

Spnal2　　<- CreateSeuratObject(counts
                              = Spnal2.data, project = "Spinal2", min.cells = 3, min.features = 200)


Spnal2[["percent.mt"]] <- PercentageFeatureSet(Spnal2, pattern = "^mt-")
Spnal2[["percent.rp"]] <- PercentageFeatureSet(Spnal2, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")
VlnPlot(Spnal2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 3)
plot1 <- FeatureScatter(Spnal2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Spnal2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

Spnal2 <- subset(Spnal2, subset = nFeature_RNA > 400 & nFeature_RNA < 8000  & percent.mt < 10 & percent.rp < 10)

VlnPlot(Spnal2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2)
plot1 <- FeatureScatter(Spnal2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Spnal2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
Spnal2 <- Add_Alt_Feature_ID(seurat_object = Spnal2,hdf5_file ="./singlecellRNAseq/2021_Blum/GSE161621_RAW/GSM4911290_mixed4_5_31_a.h5", assay = "RNA")



#Spnal3
Spnal3.data <- Read10X_h5("./singlecellRNAseq/2021_Blum/GSE161621_RAW/GSM4911291_mixed4_5_31_b.h5")

Spnal3<- CreateSeuratObject(counts = Spnal3.data, project = "Spinal3", min.cells = 3, min.features = 200)


Spnal3[["percent.mt"]] <- PercentageFeatureSet(Spnal3, pattern = "^mt-")
Spnal3[["percent.rp"]] <- PercentageFeatureSet(Spnal3, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")
VlnPlot(Spnal3_, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 3)
plot1 <- FeatureScatter(Spnal3_, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Spnal3_, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

Spnal3 <- subset(Spnal3, subset = nFeature_RNA > 400 & nFeature_RNA < 8000  & percent.mt < 10 & percent.rp < 10)

VlnPlot(Spnal3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2)
plot1 <- FeatureScatter(Spnal3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Spnal3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
Spnal3 <- Add_Alt_Feature_ID(seurat_object = Spnal3,hdf5_file ="./singlecellRNAseq/2021_Blum/GSE161621_RAW/GSM4911291_mixed4_5_31_b.h5", assay = "RNA")


#Spnal4
Spnal4.data <- Read10X_h5("./singlecellRNAseq/2021_Blum/GSE161621_RAW/GSM4911292_female5_6_25.h5")

Spnal4<- CreateSeuratObject(counts = Spnal4.data, project = "Spinal4", min.cells = 3, min.features = 200)


Spnal4[["percent.mt"]] <- PercentageFeatureSet(Spnal4, pattern = "^mt-")
Spnal4[["percent.rp"]] <- PercentageFeatureSet(Spnal4, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")
VlnPlot(Spnal4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 3)
plot1 <- FeatureScatter(Spnal4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Spnal4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

Spnal4 <- subset(Spnal4, subset = nFeature_RNA > 400 & nFeature_RNA < 8000  & percent.mt < 10 & percent.rp < 10)

VlnPlot(Spnal4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2)
plot1 <- FeatureScatter(Spnal4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Spnal4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
Spnal4 <- Add_Alt_Feature_ID(seurat_object = Spnal4,hdf5_file ="./singlecellRNAseq/2021_Blum/GSE161621_RAW/GSM4911292_female5_6_25.h5", assay = "RNA")

Spnal4 <- NormalizeData(Spnal4, normalization.method = "LogNormalize", scale.factor = 10000)
WhichCells(object = Spnal4, expression = Xist>0.5)
#Spnal5
Spnal5.data <- Read10X_h5("./singlecellRNAseq/2021_Blum/GSE161621_RAW/GSM4911293_male4_pilot_4_25.h5")

Spnal5<- CreateSeuratObject(counts = Spnal5.data, project = "Spinal5", min.cells = 3, min.features = 200)


Spnal5[["percent.mt"]] <- PercentageFeatureSet(Spnal5, pattern = "^mt-")
Spnal5[["percent.rp"]] <- PercentageFeatureSet(Spnal5, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")
VlnPlot(Spnal5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 3)
plot1 <- FeatureScatter(Spnal5, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Spnal5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

Spnal5 <- subset(Spnal5, subset = nFeature_RNA > 400 & nFeature_RNA < 8000  & percent.mt < 10 & percent.rp < 10)

VlnPlot(Spnal5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2)
plot1 <- FeatureScatter(Spnal5, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Spnal5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
Spnal5 <- Add_Alt_Feature_ID(seurat_object = Spnal5,hdf5_file ="./singlecellRNAseq/2021_Blum/GSE161621_RAW/GSM4911293_male4_pilot_4_25.h5", assay = "RNA")

Spnal5 <- NormalizeData(Spnal5, normalization.method = "LogNormalize", scale.factor = 10000)
WhichCells(object = Spnal5, expression = Xist>1)

#Spnal6
Spnal6.data <- Read10X_h5("./singlecellRNAseq/2021_Blum/GSE161621_RAW/GSM4911294_female5_7_30.h5")

Spnal6<- CreateSeuratObject(counts = Spnal6.data, project = "Spinal6", min.cells = 3, min.features = 200)


Spnal6[["percent.mt"]] <- PercentageFeatureSet(Spnal6, pattern = "^mt-")
Spnal6[["percent.rp"]] <- PercentageFeatureSet(Spnal6, pattern = "^Rps|^Rpl|^Mrpl|^Mrps")
VlnPlot(Spnal6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 3)
plot1 <- FeatureScatter(Spnal6, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Spnal6, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

Spnal6 <- subset(Spnal6, subset = nFeature_RNA > 400 & nFeature_RNA < 8000  & percent.mt < 10 & percent.rp < 10)

VlnPlot(Spnal6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2)
plot1 <- FeatureScatter(Spnal6, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Spnal6, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

Spnal6 <- Add_Alt_Feature_ID(seurat_object = Spnal6,hdf5_file ="./singlecellRNAseq/2021_Blum/GSE161621_RAW/GSM4911294_female5_7_30.h5", assay = "RNA")
Spnal6 <- NormalizeData(Spnal6, normalization.method = "LogNormalize", scale.factor = 10000)
WhichCells(object = Spnal6, expression = Xist>1)
#merge
SC_merged <- merge(Spnal1, y = c(Spnal2,Spnal3,Spnal4,Spnal5,Spnal6), add.cell.ids = c("sc1","sc2","sc3","sc4","sc5","sc6"))


#version5のSeurat Objectにする
SC_merged <- UpdateSeuratObject(SC_merged)

# SCTransform以下、いつもの処理
SC_merged <- SCTransform(SC_merged, vars.to.regress = c("percent.mt","percent.rp"))
SC_merged <- RunPCA(SC_merged,npcs=50)
ElbowPlot(SC_merged, ndims = 50)
SC_merged <- IntegrateLayers(object = SC_merged, method = CCAIntegration, normalization.method = "SCT",
                              verbose = FALSE)

SC_merged <- FindNeighbors(SC_merged, dims = 1:50, reduction = "integrated.dr", verbose = FALSE)
SC_merged <- FindClusters(SC_merged, resolution = 0.6)
SC_merged <- RunUMAP(SC_merged, dims = 1:50, reduction = "integrated.dr", verbose = FALSE)

DimPlot(SC_merged, reduction = "umap",label=TRUE,label.size=3)
DimPlot(SC_merged, reduction = "umap",label=TRUE, split.by = "orig.ident",label.size=3)
DimPlot(SC_merged, group.by='orig.ident', label=TRUE)


SC_merged[["RNA"]] <- JoinLayers(SC_merged[["RNA"]])

VlnPlot(SC_merged, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rp"), ncol = 2)
VlnPlot(SC_merged, features = c("nCount_RNA"), y.max = 25000)
VlnPlot(SC_merged, features = c("nFeature_RNA"), y.max = 5000)


plots <- VlnPlot(SC_merged, features = c("Camk2a", "Rbfox3", "Slc17a6", "Slc17a8","Slc5a7","Chat", "Gad1", "Gad2","Slc6a5"), slot="count",
                 pt.size = 0.3)
wrap_plots(plots = plots, ncol = 1)

plots <- VlnPlot(SC_merged, features = c("Mog","Plp1","Bcas1","Pdgfra","Csf1r","Cx3cr1","P2ry12","Aqp4",'Col3a1',"Pecam1","Rarres2","Atp13a5"),slot="count", 
                 pt.size = 0.3)
wrap_plots(plots = plots, ncol = 1)


# くらすたー情報の保存
SC_merged <- RenameIdents(SC_merged_Norm,'12'='LQ','24'='NFO','32'='LQ','34'='IO')
SC_merged[["cellltype"]] <- Idents(SC_merged)

SC_merged_SpinalX <- subset(SC_merged, subset = orig.ident == "Spinal6")
write.csv(SC_merged_SpinalX[[]], 'Spinal6_total.csv')

#クラスターごとの細胞数生データ
SC_merged_cellnum <- table(SC_merged$orig.ident)
SC_merged_cellnum


# Neuronのみを抽出
SC_merged_neuron <- subset(SC_merged, idents = c(0,1,3,4,6,7,8,9,10,11,13,14,16,17,18,19,20,21,22,23,25,26,27,28,29,30,31,33,35,36,37,38,39,40,41,42))
DimPlot(SC_merged_neuron)

#いつもの処理

SC_merged_neuron[["RNA"]] <- split(SC_merged_neuron[["RNA"]], f = SC_merged_neuron$orig.ident)

# SCTransform以下、いつもの処理
SC_merged_neuron <- SCTransform(SC_merged_neuron, vars.to.regress = c("percent.mt","percent.rp"))
SC_merged_neuron <- RunPCA(SC_merged_neuron,npcs=50)
ElbowPlot(SC_merged_neuron, ndims = 50)
SC_merged_neuron <- IntegrateLayers(object = SC_merged_neuron, method = CCAIntegration, normalization.method = "SCT",
                             verbose = FALSE)

SC_merged_neuron <- FindNeighbors(SC_merged_neuron, dims = 1:50, reduction = "integrated.dr", verbose = FALSE)
SC_merged_neuron <- FindClusters(SC_merged_neuron, resolution = 0.8)
SC_merged_neuron <- RunUMAP(SC_merged_neuron, dims = 1:50, reduction = "integrated.dr", verbose = FALSE)

DimPlot(SC_merged_neuron, reduction = "umap",label=TRUE,label.size=3)
DimPlot(SC_merged_neuron, reduction = "umap",label=TRUE, split.by = "orig.ident",label.size=3)
DimPlot(SC_merged_neuron, group.by='orig.ident', label=TRUE)


SC_merged_neuron@active.assay='SCT'
SC_merged_neuron[["RNA"]] <- JoinLayers(SC_merged_neuron[["RNA"]])

VlnPlot(SC_merged_neuron, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rp"), ncol = 2)
VlnPlot(SC_merged_neuron, features = c("nCount_RNA"), y.max = 25000)
VlnPlot(SC_merged_neuron, features = c("nFeature_RNA"), y.max = 5000)


plots <- VlnPlot(SC_merged_neuron, features = c("Camk2a", "Rbfox3", "Slc17a6", "Slc17a8","Slc5a7","Chat", "Gad1", "Gad2","Slc6a5"), slot="count",
                 pt.size = 0.3)
wrap_plots(plots = plots, ncol = 1)

plots <- VlnPlot(SC_merged_neuron, features = c("Mog","Plp1","Bcas1","Pdgfra","Csf1r","Cx3cr1","P2ry12","Aqp4",'Col3a1',"Pecam1","Rarres2","Atp13a5"),slot="count", 
                 pt.size = 0.3)
wrap_plots(plots = plots, ncol = 1)



# くらすたー情報の保存
SC_merged_neuron <- RenameIdents(SC_merged_neuron,'0'='Ex_interneuron','5'='Ex_interneuron','15'='Ex_interneuron','16'='Ex_interneuron','17'='Ex_interneuron','19'='Ex_interneuron','24'='Ex_interneuron','25'='Ex_interneuron','27'='Ex_interneuron','29'='Ex_interneuron','30'='Ex_interneuron','32'='Ex_interneuron','40'='Ex_interneuron','42'='Ex_interneuron','44'='Ex_interneuron','47'='Ex_interneuron',
                                 '1'="In_interneuron",'2'="In_interneuron",'6'="In_interneuron",'9'="In_interneuron",'11'="In_interneuron",'23'="In_interneuron",'31'="In_interneuron",'33'="In_interneuron",'34'="In_interneuron",'35'="In_interneuron",'37'="In_interneuron",'38'="In_interneuron",'39'="In_interneuron",'43'="In_interneuron",'45'="In_interneuron" )
SC_merged_neuron[["cellltype"]] <- Idents(SC_merged_neuron)





# cholinergic Neuronのみを抽出
SC_merged_cholinergicneuron <- subset(SC_merged_neuron, idents = c(3,4,7,8,10,12,13,14,18,20,21,22,26,28,36,41,46,48))
DimPlot(SC_merged_cholinergicneuron)

#いつもの処理

SC_merged_cholinergicneuron[["RNA"]] <- split(SC_merged_cholinergicneuron[["RNA"]], f = SC_merged_cholinergicneuron$orig.ident)

# SCTransform以下、いつもの処理
SC_merged_cholinergicneuron <- SCTransform(SC_merged_cholinergicneuron, vars.to.regress = c("percent.mt","percent.rp"))
SC_merged_cholinergicneuron <- RunPCA(SC_merged_cholinergicneuron,npcs=50)
ElbowPlot(SC_merged_cholinergicneuron, ndims = 150)
SC_merged_cholinergicneuron <- IntegrateLayers(object = SC_merged_cholinergicneuron, method = CCAIntegration, normalization.method = "SCT",
                                    verbose = FALSE)

SC_merged_cholinergicneuron <- FindNeighbors(SC_merged_cholinergicneuron, dims = 1:50, reduction = "integrated.dr", verbose = FALSE)
SC_merged_cholinergicneuron <- FindClusters(SC_merged_cholinergicneuron, resolution = 0.5)
SC_merged_cholinergicneuron <- RunUMAP(SC_merged_cholinergicneuron, dims = 1:50, reduction = "integrated.dr", verbose = FALSE)

DimPlot(SC_merged_cholinergicneuron, reduction = "umap",label=TRUE,label.size=3)
DimPlot(SC_merged_cholinergicneuron, reduction = "umap",label=TRUE, split.by = "orig.ident",label.size=3)
DimPlot(SC_merged_cholinergicneuron, group.by='orig.ident', label=TRUE)

SC_merged_cholinergicneuron.hierarchical<-BuildClusterTree(object = SC_merged_cholinergicneuron, assay = 'SCT',slot = "data")
PlotClusterTree(SC_merged_cholinergicneuron.hierarchical)

SC_merged_cholinergicneuron@active.assay='RNA'
SC_merged_cholinergicneuron[["RNA"]] <- JoinLayers(SC_merged_cholinergicneuron[["RNA"]])

VlnPlot(SC_merged_cholinergicneuron, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rp"), ncol = 2)
VlnPlot(SC_merged_cholinergicneuron, features = c("nCount_RNA"), y.max = 25000)
VlnPlot(SC_merged_cholinergicneuron, features = c("nFeature_RNA"), y.max = 5000)


plots <- VlnPlot(SC_merged_cholinergicneuron, features = c("Camk2a","Slc17a6", "Slc17a8","Slc5a7","Chat", "Gad1", "Gad2","Slc6a5","Bcl6","Pax2","Nos1","Pitx2"),layer='count',
                 stack = T,flip = T,pt.size = 0.3)
wrap_plots(plots = plots, ncol = 1)

plots <- VlnPlot(SC_merged_cholinergicneuron_Norm, features = c("Mog","Plp1","Bcas1","Pdgfra","Csf1r","Cx3cr1","P2ry12","Aqp4",'Col3a1',"Pecam1","Rarres2","Atp13a5"),slot="count", 
                 pt.size = 0.3)
wrap_plots(plots = plots, ncol = 1)


# くらすたー情報の保存
SC_merged_cholinergicneuron <- RenameIdents(SC_merged_cholinergicneuron,'11'='LQ','24'='LQ','0'='Skel_MN','6'='Skel_MN','12'='Skel_MN','17'='Skel_MN','22'='Skel_MN','23'='Skel_MN')
SC_merged_cholinergicneuron[["cellltype"]] <- Idents(SC_merged_cholinergicneuron)


# Nos1+cholinergic Neuronのみを抽出
SC_merged_cholinergicneuron_Nos1 <- subset(SC_merged_cholinergicneuron, idents = c(1,2,3,5,7,8,10,13,14,15,16,18,19,20,21))
DimPlot(SC_merged_cholinergicneuron_Nos1)

#いつもの処理

SC_merged_cholinergicneuron_Nos1[["RNA"]] <- split(SC_merged_cholinergicneuron_Nos1[["RNA"]], f = SC_merged_cholinergicneuron_Nos1$orig.ident)

# SCTransform以下、いつもの処理
SC_merged_cholinergicneuron_Nos1 <- SCTransform(SC_merged_cholinergicneuron_Nos1, vars.to.regress = c("percent.mt","percent.rp"))
SC_merged_cholinergicneuron_Nos1 <- RunPCA(SC_merged_cholinergicneuron_Nos1,npcs=50)
ElbowPlot(SC_merged_cholinergicneuron_Nos1, ndims = 150)
SC_merged_cholinergicneuron_Nos1 <- IntegrateLayers(object = SC_merged_cholinergicneuron_Nos1, method = CCAIntegration, normalization.method = "SCT",
                                               verbose = FALSE)

SC_merged_cholinergicneuron_Nos1 <- FindNeighbors(SC_merged_cholinergicneuron_Nos1, dims = 1:50, reduction = "integrated.dr", verbose = FALSE)
# MPN_merged <- FindClusters(MPN_merged, resolution = 0.05)
SC_merged_cholinergicneuron_Nos1 <- FindClusters(SC_merged_cholinergicneuron_Nos1, resolution = 0.9)
SC_merged_cholinergicneuron_Nos1 <- RunUMAP(SC_merged_cholinergicneuron_Nos1, dims = 1:50, reduction = "integrated.dr", verbose = FALSE)



SC_merged_cholinergicneuron_Nos1<-PrepSCTFindMarkers(SC_merged_cholinergicneuron_Nos1)




plots <- VlnPlot(SC_merged_cholinergicneuron_Nos1, features = c("Camk2a","Slc17a6", "Slc17a8","Slc5a7","Chat", "Gad1", "Gad2","Slc6a5","Bcl6","Pax2","Nos1","Pitx2"),slot="data", pt.size = 0.3,stack = T,flip = T)
wrap_plots(plots = plots, ncol = 1)

plots <- VlnPlot(SC_merged_cholinergicneuron_Nos1, features = c("Mog","Plp1","Bcas1","Pdgfra","Csf1r","Cx3cr1","P2ry12","Aqp4",'Col3a1',"Pecam1","Rarres2","Atp13a5"),slot="data", 
                 pt.size = 0.3,stack = T,flip = T)
wrap_plots(plots = plots, ncol = 1)


# EDFig1a
DimPlot(SC_merged_cholinergicneuron_Nos1, reduction = "umap",label=TRUE,label.size=3)



# くらすたー情報の保存
SC_merged_cholinergicneuron_Nos1 <- RenameIdents(SC_merged_cholinergicneuron_Nos1,'0'='sn1','1'='sn2','2'='sn3','3'='sn4','4'='sn4','5'='sn5','6'='sn6','7'='ns7','8'='sn8','9'='sn9','10'='sn10','11'='sn11','12'='sn12'
                                                 ,'13'='sn13','14'='sn14','15'='sn15','16'='sn16','17'='sn17','18'='sn18','19'='sn19','20'='sn20','21'='sn21','22'='sn22','23'='sn23','24'='sn24','25'='sn25','26'='sn26')
SC_merged_cholinergicneuron_Nos1[["cellltype"]] <- Idents(SC_merged_cholinergicneuron_Nos1)



SC_merged_cholinergicneuron_Nos1@active.ident <- SC_merged_cholinergicneuron_Nos1$cellltype

# EDFig1b
markers.to.plot =  c("Tmem132c", "Sst", "Slc26a7", "Penk", "Mdga1", "Lmo3", "Frem1", 
                     "Gpc3", "Rorb", "Palld", "Slc16a12", "Foxp2", "Trhr", "Nts", "Nfib",
                     "Fam19a1", "Pde1c", "Nek10", "Sall3", "Postn", "Col12a1",'Pkib' ,"Il1r1",
                     "Vav3", "Ccbe1", "Chodl", "Fam20a", "Ndst4", "Fam198b", "Tnc", 
                     "Cartpt","Fam163a","Spon1", "Lpar1", "Daam2","Ano1","Kif6",
                     "Reln","Col24a1","Tll2","Vegfc","Nms","Galnt14","Plch1","Lama3","Glp2r",
                     "Plcz1", "Lgr5","Agbl1","Syt10","Mpped2","Bcl11a","Dscaml1",
                     "Piezo2","Pappa2","Kcnmb2","Cd40","Fstl4","Ror1","Trpc3","Pde11a","Qrfpr",
                     "Etv1","Col25a1","Mamdc2","Fgf10","Agtr1b","Calca","Fn1","Ccdc3",
                     "Aldh1a1","Grin2c","Sox5","Tspan9","Eln","Grpr","Igf1","Lama2","Mctp2","Lhfpl2",
                     "Cdh23","Pcsk5","Scube2","Npffr2","Rspo2","Sema5a","Bdnf","Phldb2","Pde3a")


DotPlot(SC_merged_cholinergicneuron_Nos1, features = markers.to.plot.short, cols =c("white","tomato2"),dot.scale = 8) +RotatedAxis()

