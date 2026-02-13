#-------------------------------------------------------scDblFinder---------------------------------------------------
#Filter out doublets: scDblFinder

#setpath
setwd("D:/AMAN_PhD/Script/scRNA_seq/Seurat/Doublets")

#install package
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("scDblFinder")
BiocManager::install("DoubletFinder")
BiocManager::install("scater")
BiocManager::install("bluster")
BiocManager::install("scDblFinder", force = TRUE)


#load libraries
library(bluster)
library(scDblFinder)
library(DoubletFinder)
library(Seurat)
library(tidyverse)
library(scater)

#set directory
data_dir <- "D:/AMAN_PhD/Script/scRNA_seq/Seurat/Doublets/raw_feature_bc_matrix"

#create counts matrix
cts <- Read10X(data.dir = data_dir)




cts[1:10,1:10]

# create Seurat object
pbmc.seurat <- CreateSeuratObject(counts = cts)
str(pbmc.seurat)

#Mitochondrial reads
pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = '^MT-')

#Filtered reads
pbmc.seurat.filtered <- subset(pbmc.seurat, subset = nCount_RNA > 800 &
                                 nFeature_RNA > 500 &
                                 mitoPercent < 10)

sce <- as.SingleCellExperiment(pbmc.seurat.filtered)

sce <- logNormCounts(sce)

sce <- runPCA(sce, ncomponents = 50)

sce <-scDblFinder(sce)

 

colData(sce)

table(sce$scDblFinder.class)

pbmc.seurat.filtered$scDblFinder_class <- sce$scDblFinder.class
pbmc.seurat.filtered$scDblFinder_score <- sce$scDblFinder.score

pbmc.seurat.filtered <- NormalizeData(pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindVariableFeatures(pbmc.seurat.filtered)
pbmc.seurat.filtered <- ScaleData(pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunPCA(pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunUMAP(pbmc.seurat.filtered, dims = 1:30)

DimPlot(
  pbmc.seurat.filtered,
  group.by = "scDblFinder_class",
  cols = c("singlet" = "grey70", "doublet" = "red")
)

seurat_singlets <- subset(
  pbmc.seurat.filtered,
  subset = scDblFinder_class == "singlet"
)
