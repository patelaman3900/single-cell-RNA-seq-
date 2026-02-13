#-------------------------------------------Markers and cluster Identification---------------------------------------------
#Finding differentially expressed features and cluster identification in scRNA Seq data in R
#Find markers between specific clusters{ findMarkers() }
#Find markers in a cluster compared to all other clusters{ findAllMarkers() }
#Find markers conserved across conditions { findConservedMArkers }
#CONDITION
 #Control
 #Treated
#CELL TYPE
 #Cell type A
 #Cell type B

#GOAL of STUDY: To assess cell-type-specific changes in gene expression in 
                #eight lupus patient samples treated with interferon (IFN)-B
#STUDY DESIGN: Peripheral blood mononuclear cells (PBMCs) from eight lupus patients were split into a 
                #stimulated and control group and the stimulated group was treated with interferon beta.
#GOAL of ANALYSIS: To identify clusters
                  #To get genes differentially expressed between conditons in a particular celltype

#Required packages
 #Seurat
 #tidyverse

#SCRIPT
# script to identify cluster identity -----------------
# Finding markers in every cluster
# Finding conserved markers 
# Finding markers DE between conditions

setwd("D:/AMAN_PhD/Script/scRNA_seq/Seurat/Markers_and_cluster_identification")

set.seed(1234)
#install packages
install.packages("Seurat")
install.packages("tidyverse")
install.packages('BiocManager')
BiocManager::install('multtest')
install.packages('metap')
chooseCRANmirror()
install.packages("metap", type = "binary")

#load library
library(Seurat)
library(tidyverse)
library(multtest)
library("metap")


# load data
ifnb_harmony <- readRDS('D:/AMAN_PhD/Script/scRNA_seq/Seurat/Markers_and_cluster_identification/ifnb_harmony.rds')
str(ifnb_harmony)
View(ifnb_harmony@meta.data)

# visualize data
clusters <- DimPlot(ifnb_harmony, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
condition <- DimPlot(ifnb_harmony, reduction = 'umap', group.by = 'stim')

condition|clusters

# findAll markers -----------------

FindAllMarkers(ifnb_harmony,            #seurat object
               logfc.threshold = 0.25,  #min avg value of cluster one wrt to all other
               min.pct = 0.1,           #min percentage of gene tested
               only.pos = TRUE,
               test.use = 'DESeq2',
               slot = 'counts')


# findConserved markers -------------

# Notes:
# slot depends on the type of the test used, 
# default is data slot that stores normalized data
# DefaultAssay(ifnb_harmony) <- 'RNA'

DefaultAssay(ifnb_harmony)

markers_cluster3 <- FindConservedMarkers(ifnb_harmony,
                                         ident.1 = 3,         #cluster want to identify wrt to other clusters
                                         grouping.var = 'stim')

head(markers_cluster3)
#gene name= FCGR3A
#control p value = 0 
#control log2FC= 4.099 (upregulated in cluster 3 wrt all other clusters)
#control pct1= 0.977(gene is expressed in 97.7% of cells in cluster 3)
#control pct2= 0.206 (gene is expressed in 20.6% of cells in all other  cluster)
#control p adj val= 0

# let's visualize top features
FeaturePlot(ifnb_harmony, features = c('FCGR3A'), min.cutoff = 'q10') #q10= 10th quantile
#grey= expression is lower than given quantile
#purple= expression is higher  than given quantile

# min-cut off explanation:
seq(1,5)
SetQuantile('q50', seq(1,5))
SetQuantile('q10', seq(1,5))





# rename cluster 3 ident (CD16 Mono marker= 'FCGR3A')
Idents(ifnb_harmony)
ifnb_harmony <- RenameIdents(ifnb_harmony, `3` = 'CD16 Mono')
Idents(ifnb_harmony)

DimPlot(ifnb_harmony, reduction = 'umap', label = T)

# cells already have annotations provided in the metadata
View(ifnb_harmony@meta.data)

# Settings cluster identities is an iterative step(loop)
# multiple approaches could be taken - automatic/manual annotations (sometimes both)
# need to make sure each cell type forms a separate cluster

# setting Idents as Seurat annotations provided (also a sanity check!)
Idents(ifnb_harmony) <- ifnb_harmony@meta.data$seurat_annotations
Idents(ifnb_harmony)

DimPlot(ifnb_harmony, reduction = 'umap', label = TRUE)


# findMarkers between conditions ---------------------
ifnb_harmony$celltype.cnd <- paste0(ifnb_harmony$seurat_annotations,'_', ifnb_harmony$stim)
View(ifnb_harmony@meta.data)
Idents(ifnb_harmony) <- ifnb_harmony$celltype.cnd

DimPlot(ifnb_harmony, reduction = 'umap', label = TRUE)

# find markers (compare betweenthe clusters)
b.interferon.response <- FindMarkers(ifnb_harmony, ident.1 = 'CD16 Mono_STIM', ident.2 = 'CD16 Mono_CTRL')

head(b.interferon.response)

# plotting conserved features vs DE features between conditions
head(markers_cluster3)


FeaturePlot(ifnb_harmony, features = c('FCGR3A', 'AIF1', 'IFIT1'), split.by = 'stim', min.cutoff = 'q10')




#DATA 
#1) Data: 
 # https://drive.google.com/file/d/13I2250SX-vQfV41th0JmF6BIyx0FGjzO/view

#2) Link to code:
 # https://github.com/kpatel427/YouTubeTutorials/blob/main/singleCell_CI_markers.R

#3) Vignettes:
 #▸ https://satijalab.org/seurat/articles...
 #▸ https://satijalab.org/seurat/articles...
 #▸ https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionIV/lessons/SC_marker_identification.html

#4) Marker databases:
 #a. SCSig: https://www.gsea-msigdb.org/gsea/msigdb/supplementary_genesets.jsp
 #b. PangloDB: https://panglaodb.se/
 #c. CellMarker: http://bio-bigdata.hrbmu.edu.cn/CellMarker/help.jsp

#-------------------------------------------Markers and cluster Identification---------------------------------------------