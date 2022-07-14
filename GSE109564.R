## Load in library packages
suppressPackageStartupMessages({
  library(dplyr)
  library(plyr)
  library(skimr)
  library(ggplot2)
  library(ggbeeswarm)
  library(ggthemes)
  library(rtracklayer)
  library(Seurat)
  library(harmony)
  library(remotes)
  library(scAnnotatR)
  library(SingleR)
  library(celldex)
  library(SingleCellExperiment)
  library(sceasy)
  library(scater)
  library(reticulate)
  library(SeuratDisk)
  library(SeuratWrappers)
  library(monocle3)
  library(anndata)
  library(LEAP)
  library(destiny)
  library(phateR)
  library(readr)
  library(TSCAN)
  library(slingshot)
  library(celda)
  library(Signac)
  library(motifmatchr)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg19)

  # Load knitr to keep code within margin
  library(knitr)
})


## Install Bioconductor packages if need be...

#remotes::install_github("campbio/celda")
#remotes::install_github("milescsmith/ReductionWrappers")
#remotes::install_github("mojaveazure/seurat-disk")
#remotes::install_github('satijalab/seurat-wrappers')
#install_github("cellgeni/sceasy")
#remotes::install_github('cole-trapnell-lab/leidenbase')
#BiocManager::install("rtracklayer")
#BiocManager::install("scater")
#BiocManager::install("scAnnotatR")
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
#BiocManager::install("monocle3")
#remotes::install_github('cole-trapnell-lab/monocle3', dependencies = TRUE)
#BiocManager::install("TSCAN")
#BiocManager::install("slingshot")
#BiocManager::install("SingleCellExperiment")
#BiocManager::install("destiny")
#BiocManager::install("motifmatchr")
#BiocManager::install("JASPAR2020")
#BiocManager::install("TFBSTools")
#BiocManager::install("BSgenome")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
#BiocManager::install("monocle")
#BiocManager::install(c('DelayedArray', 'DelayedMatrixStats', 'org.Hs.eg.db', 'org.Mm.eg.db'))
#install_github("cole-trapnell-lab/garnett")
#library(garnett)
#library(org.Hs.eg.db)



# Load in Data
GSE109564 <- read.table("GSE109564_Kidney.biopsy.dge.txt.gz")

# Convert data into Seurat objects
kidney_biopsy <- CreateSeuratObject(GSE109564)

# Calculate the percentage of mitochondrial genes in kidney_biopsy
mito_biopsy <- grep(pattern = "^MT-", x = rownames(kidney_biopsy@assays[["RNA"]]), value = T)
percent_mito_biopsy <- Matrix::colSums(kidney_biopsy@assays[["RNA"]][mito_biopsy,]) / Matrix::colSums(kidney_biopsy@assays[["RNA"]])
kidney_biopsy <- AddMetaData(object = kidney_biopsy, metadata = percent_mito_biopsy, col.name = "percent.mito")



# Exploratory Data Analysis
VlnPlot(kidney_biopsy, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
par(mfrow = c(1, 2))
FeatureScatter(kidney_biopsy, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(kidney_biopsy, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")



# Subset and normalize kidney_biopsy
kidney_biopsy <- subset(kidney_biopsy, subset = nFeature_RNA < 4000 & percent.mito < 0.30)
kidney_biopsy <- NormalizeData(kidney_biopsy, normalization.method = "LogNormalize", scale.factor = 10000)
kidney_biopsy <- FindVariableFeatures(kidney_biopsy, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 6), dispersion.cutoff = c(1, Inf))
top20_biopsy <- head(VariableFeatures(kidney_biopsy), 20)
feature_plot_biopsy <- VariableFeaturePlot(kidney_biopsy, pt.size = 0.5, log = T, selection.method = "mean.var.plot")
top20_plot_biopsy <- LabelPoints(plot = feature_plot_biopsy, points = top20_biopsy, repel = T)
top20_plot_biopsy
kidney_biopsy <- ScaleData(object = kidney_biopsy, vars.to.regress = c("nCount_RNA", "percent.mito"))



# Principal Component Analysis
kidney_biopsy <- RunPCA(object = kidney_biopsy, npcs = 50, verbose = F)

# Plot PCA of all 4 tubes
DimPlot(object = kidney_biopsy, reduction = "pca")



# Distribution of p-values for each PC
kidney_biopsy <- JackStraw(kidney_biopsy, num.replicate = 100)
kidney_biopsy <- ScoreJackStraw(kidney_biopsy, dims = 1:20)
JackStrawPlot(kidney_biopsy, dims = 1:15)

# Observe how many PCs to choose according to the "elbow"
ElbowPlot(kidney_biopsy)

PCHeatmap(kidney_biopsy, dims = 1:6, cells = 500, balanced = T)



# Compute nearest neighbor
kidney_biopsy <- FindNeighbors(kidney_biopsy, reduction = "pca", dims = 1:20)

# Obtain number of clusters
kidney_biopsy <- FindClusters(kidney_biopsy, algorithm = 4, resolution = 0.6)



# tSNE plot
kidney_biopsy <- RunTSNE(kidney_biopsy, dims = 1:20)
DimPlot(object = kidney_biopsy, reduction = "tsne")



# Look for batch effect (there are none per authors' notes)
kidney_biopsy <- RunUMAP(kidney_biopsy, reduction = "pca", dims = 1:20)
DimPlot(kidney_biopsy, reduction = "umap", group.by = "orig.ident")



# Differentially Expressed Features (Cluster Biomarkers)
kidney_biopsy_markers <- FindAllMarkers(kidney_biopsy, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
kidney_biopsy_markers %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC)
kidney_biopsy_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC) -> top3
DoHeatmap(kidney_biopsy, features = top3$gene) + NoLegend()



# Manual Cell Type Identification
manual_cell_types <- c("LOH (DL)","PT","EC","Mono-1","LOH (AL)","T Cells","Mono-2","Plasma-1","Collecting Duct","Pericyte","Fibroblast","B Cells","Mast Cells","Myofibroblast","Plasma-2","Cycling")
names(manual_cell_types) <- levels(kidney_biopsy)
kidney_biopsy <- RenameIdents(kidney_biopsy, manual_cell_types)
DimPlot(kidney_biopsy, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

TSNEPlot(object = kidney_biopsy, label = T, pt.size = 0.5)
UMAPPlot(object = kidney_biopsy, label = T, pt.size = 0.5)
DoHeatmap(kidney_biopsy, features = top3$gene, size = 3, angle = 90) + NoLegend()


### END OF REPLICATION
