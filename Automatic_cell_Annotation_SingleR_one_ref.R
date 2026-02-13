#---------------------------------------------Automatic cell Annotation (SingleR) one Reference-------------------------------------------------------------------
#Cell Annotation Workflow
#Automatic Annotation
 # BY comparison of data with the annotated reference data
 # By using known marker gene for known cell type
#Manual Annotation
 #Identify and label the unknown / unlabelled cell cluster
#Verification  
 #Experiment 
 #Independent data integration 
 # Statistical support
 # Consultant with experts
 
#Marker-based Annotation
 #Labels cells or cell clusters on the basis of the characteristic expression of known marker genes.
 #Known relationships between marker genes and cell types are obtained from databases, such as
   #MSigDB
   #PanglaoDB
   #CellMarker
   #manually from the literature

#Reference-based Annotation
 #Transfer labels from a reference cell or cluster (from experty annotated scRNA-seq data) to a sufficiently 
 #similar one in the query (data to be annotated)
 #Reference single-cell data are obtained from sources such as -
   #Gene Expression Omnibus (GEO)
   #The Single Cell Expression Atlas
   #Cell atlas projects

#TOOLS FOR MARKER BASED ANNOTATION
#Tool 1
 #AUCell
#Type
 #Marker based
#Language
 #R
#Resolution
 #Single cells
#Approach
 #Area under the curve to estimate marker gene set enrichment
#Allows 'None'
 #Yes
#Notes
#Because of low detection rates at the level of single cells, it requires many markers for every cell type


#Tool 2
 #SCINA
#Type
 #Marker based
#Language
 #R
#Resolution
 #Single cells
#Approach
 #Expectation maximization, Gaussian mixture model
#Allows 'None'
 #(Optional)
#Notes
 #Simultaneously clusters and annotates cells; robust to the inclusion of incorrect marker genes


#Tool 3
 #GSEA/GSVA
#Type
 #Marker based
#Language
 #R/Java
#Resolution
 #Clusters of cells
#Approach
 #Enrichment test
#Allows 'None'
 #Yes
#Notes
 #Marker gene lists must be reformatted in GMT format. Markers must all be differentially expressed in the same direction in the cluster

#Strengths:
  #These methods will assign labels only to cells associated with known markers, and other cells will remain unlabeled.
#Pitfalls:
 #Markers are not easily accessible for all cell types.
 #The marker gene or gene set (a collection of marker genes) should be specifically and consistently expressed in a given cell or cluster
 #Work well once a relevant and sufficiently large set of marker genes are collected
 #These methods work better for annotating major cell types and may not be able to effectively distinguish subtypes. 


#TOOLS FOR REFERENCE BASED ANNOTATION
#Tool 1 
 #singleCell Net
#Type
 #Reference based
#Language
 #R 
#Resolution
 #Single cells
#Approach
 #Relative-expression gene pairs + random forest
#Allows 'None'
 #Yes, but rarely does so even when it should
#Notes
 #10-100x slower than other methods; high accuracy

#Tool 2
 #scmap-cluster
#Type
 #Reference based
#Language
 #R 
#Resolution
 #Single cells
#Approach
 #Consistent correlations
#Allows 'None'
 #Yes  
#Notes
 #Fastest method available; balances false-positives and false-negatives; includes web interface for use with a large pre-built reference 
 #or custom reference set

#Tool 3
 #scmap-cell
#Type
 #Reference based
#Language
 #R 
#Resolution
 #Single cells
#Approach
 #Approximate nearest neighbors
#Allows 'None'
 #Yes  
#Notes
 #Assigns individual cells to nearest neighbor cells in reference; allows mapping of cell trajectories; fast and scalable

#Tool 4
 #singleR
#Type
 #Reference based
#Language
 #R 
#Resolution
 #Single cells
#Approach
 #k-nearest neighbors, support vector machine, random forest, nearest mean classifier and linear discriminant analysis
#Allows 'None'
 #(optional)
#Notes
 #Expertise required for correct design and appropriate training of classifier while avoiding overtraining


#Strengths:
  #Accuracy of assigned labels and avoiding incorrect labeling of novel cell types
#Pitfalls:
  #Approach is feasible only if high-quality and relevant annotated reference single-cell data are available.
  #Some tools have low accuracy if the reference data are incomplete or represent a poor match 


#SingleR - How does it work?
  #Input:Unannotated scRNA-seq data
  #Reference transcriptomes of pure cell types
  #Output:Annotated single cells

#Step 1:
  #Identifying variable genes among cell types in the reference set
#Step 2:
  #Correlating each single-cell transcriptome with each sample in the reference set
#Step3: 
  #Iterative fine-tuning-reducing the reference set to only top cell types


#Study design and Goal of Analysis
#Goal: 
 #Annotate cell types in 20k human peripheral blood mononuclear cells (PBMCs)
#Study Design: 
 #PBMCs of a healthy female donor aged 25-30 were obtained by 10x Genomics.

#Required R packages
 #singleR (v1.6.1)
 #celldex (v1.2.0)  reference data set 
 #Seurat (v4.3.0)
 #tidyverse (v2.0.0)
 #pheatmap (v1.0.12)


#SCRIPT
# script to annotate cell types from 20k Human PBMCs from a healthy female donor
setwd("D:/AMAN_PhD/Script/scRNA_seq/Seurat/Automatic_cell_Annotation_SingleR")


#install packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SingleR")
BiocManager::install("celldex")
install.packages("pheatmap")
install.packages("hdf5r")

#load libraries
library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
library(pheatmap)
library(hdf5r)

# Input Data 10X CellRanger .HDF5 format --------------
hdf5_obj <- Read10X_h5(filename = 'D:/AMAN_PhD/Script/scRNA_seq/Seurat/Automatic_cell_Annotation_SingleR/20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5',
                       use.names = TRUE,
                       unique.features = TRUE)
pbmc.seurat <- CreateSeuratObject(counts = hdf5_obj)
pbmc.seurat

# QC and Filtering -----------
# explore QC
pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = '^MT-')
pbmc.seurat.filtered <- subset(pbmc.seurat, subset = nCount_RNA > 800 &
                                 nFeature_RNA > 500 &
                                 mitoPercent < 10)


# It is a good practice to filter out cells with non-sufficient genes identified and genes with non-sufficient expression across cells.


# pre-process standard workflow ---------------
pbmc.seurat.filtered <- NormalizeData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindVariableFeatures(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- ScaleData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunPCA(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindNeighbors(object = pbmc.seurat.filtered, dims = 1:20)
pbmc.seurat.filtered <- FindClusters(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:20)

# running steps above to get clusters
View(pbmc.seurat.filtered@meta.data)
DimPlot(pbmc.seurat.filtered, reduction = 'umap')

# get reference data -----------
ref <- celldex::HumanPrimaryCellAtlasData()
View(as.data.frame(colData(ref)))

# expression values are log counts (log normalized counts)


# run SingleR (default mode) ---------
# default for SingleR is to perform annotation of each individual cell in the test dataset

pbmc_counts <- GetAssayData(pbmc.seurat.filtered, layer = 'counts')

pred <- SingleR(test = pbmc_counts,
                ref = ref,
                labels = ref$label.main)

pred

pbmc.seurat.filtered$singleR.labels <- pred$labels[match(rownames(pbmc.seurat.filtered@meta.data), rownames(pred))]
DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = 'singleR.labels')


# Annotation diagnostics ----------


# ...Based on the scores within cells -----------
pred
pred$scores

plotScoreHeatmap(pred)


# ...Based on deltas across cells ----------

plotDeltaDistribution(pred)




# ...Comparing to unsupervised clustering ------------

tab <- table(Assigned=pred$labels, Clusters=pbmc.seurat.filtered$seurat_clusters)
pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10))



#DATA RESOURCES
#▸ Link to Data:
 # https://www.10xgenomics.com/datasets/20-k-human-pbm-cs-3-ht-v-3-1-chromium-x-3-1-high-6-1-0

#▸ Link to Code:
 # https://github.com/kpatel427/YouTubeTutorials/blob/main/annotateSingleR.R

#▸ Resources/Vignettes:
 #1. https://www.nature.com/articles/s41596-021-00534-0#sec3
 #2. https://bioconductor.org/books/release/SingleRBook/
 #3. https://bioconductor.org/books/3.14/OSCA.basic/cell-type-annotation.html#assigning-cell-labels-from-reference-data



#---------------------------------------------Automatic cell Annotation (SingleR) one Reference-------------------------------------------------------------------


