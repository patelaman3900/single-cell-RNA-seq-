#----------------------------------Integrating scRNA seq dataset in R using Harmony----------------------------------
#Steps
#1.Soft assign cells to clusters, favoring mixed dataset representation
#2.Get cluster centroids for each dataset
#3.Get dataset correction factors for each cluster
#4.Move cells based on soft cluster membership
#Repeat Again through step 1 until convergence

#GOAL of STUDY:
#To access cell-type-specific changes in gene expression with IFN beta treatment. 
#STUDY DESIGN: 
#PBMC from 8 Lupus patients split into control and stimulated with interferon beta   
#GOAL of ANALYSIS
#To integrate data by condition, overlay cells that are similar in both conditions

#Requirements
#R packages:- 
#Harmony
#Seurat
#SeuratData
#tidyverse
#ggplot2



#SCRIPT
#Script to integrate across conditions using Harmony
setwd("D:/AMAN_PhD/Script/scRNA_seq/Seurat/Integration")

#install packages
install.packages("harmony")
install.packages("Seurat")
install.packages("remotes")
remotes::install_github("satijalab/seurat-data")
install.packages("tidyverse")
install.packages("ggplot2")

# set seed for reproducibility
set.seed(1234)

library(harmony)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(ggplot2)



#get data -------------------------
AvailableData()
# install dataset
InstallData("ifnb")

# load dataset
LoadData("ifnb")
str(ifnb)
#update the seurat object structure
ifnb <- UpdateSeuratObject(ifnb)

# QC and filtering
ifnb$mito.percent <- PercentageFeatureSet(ifnb, pattern = '^MT-')
View(ifnb@meta.data)

# explore QC

# filter

ifnb.filtered <- subset(ifnb, subset = nCount_RNA > 800 &
                          nFeature_RNA > 200 & 
                          mito.percent < 5)

# standard workflow steps
ifnb.filtered <- NormalizeData(ifnb.filtered)
ifnb.filtered <- FindVariableFeatures(ifnb.filtered)
ifnb.filtered <- ScaleData(ifnb.filtered)
ifnb.filtered <- RunPCA(ifnb.filtered)
ElbowPlot(ifnb.filtered)
ifnb.filtered <- RunUMAP(ifnb.filtered, dims = 1:20, reduction = 'pca')

#DimPlot to check the structure of data
before <- DimPlot(ifnb.filtered, reduction = 'umap', group.by = 'stim')


# run Harmony -----------
ifnb.harmony <- ifnb.filtered %>%
  RunHarmony(group.by.vars = 'stim', plot_convergence = FALSE)

ifnb.harmony@reductions

ifnb.harmony.embed <- Embeddings(ifnb.harmony, "harmony")
ifnb.harmony.embed[1:10,1:10]



# Do UMAP and clustering using ** Harmony embeddings instead of PCA **
ifnb.harmony <- ifnb.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# visualize 
after <- DimPlot(ifnb.harmony, reduction = 'umap', group.by = 'stim')

before|after



#1) Harmony Paper:  https://doi.org/10.1038/s41592-019-0619-0
#2) Harmony Vignette: https://portals.broadinstitute.org/harmony/articles/quickstart.html
#3) Link to code: https://github.com/kpatel427/YouTubeTutorials/blob/main/singleCell_integrate_harmony.R
#4) Additional resources: https://www.singlecellcourse.org/scrna-seq-dataset-integration.html#harmony-3-vs-5-10k-pbmc
#https://swaruplab.bio.uci.edu/tutorial/integration/integration_tutorial.html

#----------------------------------Integrating scRNA seq dataset in R using Harmony----------------------------------
