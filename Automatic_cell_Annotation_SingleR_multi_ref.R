#---------------------------------------------Automatic cell Annotation (SingleR) multiple Reference-------------------------------------------------------------------

#Strategies for using multiple references
 #Strategy 1: Using reference-specific labels in a combined reference
  #Cons:
   #Batch effects are not taken care of
   #Loss of precision due to noise during the calculation of the score in each reference
   #Risk of technical variation dominating classification results
 #Strategy 2: Combining scores across multiple references (default approach implemented in SingleR)
  #Cons:
   #Lack of consistency in labels across references complicates interpretation
 #Strategy 3: Using harmonized labels in a combined reference
  #Cons:
   #Assumes that harmonized labels are available
   #Mapping process also runs the risk of discarding relevant information about the biological status (e.g., activation status, disease condition) if there is 
   #no obvious counterpart for that state in the ontology.

#Study design and Goal of Analysis
 #Goal: Annotate cell types in 20k human peripheral blood mononuclear cells (PBMCs)
 #Study Design: PBMCs of a healthy female donor aged 25-30 were obtained by 10x Genomics.


#Required R packages
 #singleR (v1.6.1)
 #celldex (v1.2.0)   reference data set
 #Seurat (v4.3.0)
 #tidyverse (v2.0.0)
 #pheatmap (v1.0.12)

#SCRIPT
# script to annotate cell types from 20k Human PBMCs from a healthy female donor
#Set path 
setwd("D:/AMAN_PhD/Script/scRNA_seq/Seurat/Automatic_cell_Annotation_SingleR")

#Install packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SingleR")
BiocManager::install("celldex")
install.packages("pheatmap")
install.packages("hdf5r")


#Load libraries
library(SingleR)
library(celldex)  
library(Seurat)
library(tidyverse)
library(pheatmap)

# 10X CellRanger .HDF5 format ---------
hdf5_obj <- Read10X_h5(filename = "D:/AMAN_PhD/Script/scRNA_seq/Seurat/Automatic_cell_Annotation_SingleR/20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5",
                       use.names = TRUE,
                       unique.features = TRUE)
pbmc.seurat <- CreateSeuratObject(counts = hdf5_obj)
pbmc.seurat

# QC and Filtering ----------
# explore QC
pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = '^MT-')
VlnPlot(pbmc.seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)
pbmc.seurat.filtered <- subset(pbmc.seurat, subset = nCount_RNA > 800 &
                                 nFeature_RNA > 500 &
                                 mitoPercent < 10)

# It is a good practice to filter out cells with non-sufficient genes identified and genes with non-sufficient expression across cells.


# pre-process standard workflow ---------------
pbmc.seurat.filtered <- NormalizeData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindVariableFeatures(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- ScaleData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunPCA(object = pbmc.seurat.filtered)
ElbowPlot(pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindNeighbors(object = pbmc.seurat.filtered, dims = 1:20)
pbmc.seurat.filtered <- FindClusters(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:20)

# running steps above to get clusters
DimPlot(pbmc.seurat.filtered, reduction = "umap")
View(pbmc.seurat.filtered@meta.data)


# run SingleR with multiple reference datasets (default mode) ---------

# for pbmc data, we will use two datasets
hpca <- celldex::HumanPrimaryCellAtlasData()
dice <- celldex::DatabaseImmuneCellExpressionData()

# ...1. Strategy 1: Using reference-specific labels ----------
hpca$label.main
dice$label.main

# adding ref info to labels
hpca$label.main <- paste0('HPCA.', hpca$label.main)
dice$label.main <- paste0('DICE.', dice$label.main)

# create a combined ref based on shared genes
shared <- intersect(rownames(hpca), rownames(dice))
combined <- cbind(hpca[shared,], dice[shared,])
combined
combined$label.main

# run singleR using combined ref
# savings counts into a separate object
pbmc_counts <- GetAssayData(
  pbmc.seurat.filtered,
  assay = "RNA",
  layer = "counts"
)

com.res1 <- SingleR(test = pbmc_counts, ref = combined, labels = combined$label.main)
table(com.res1$labels)

pbmc.seurat.filtered$com.res1.labels <- com.res1[match(rownames(pbmc.seurat.filtered@meta.data), rownames(com.res1)), 'labels']
View(pbmc.seurat.filtered@meta.data)

DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = 'com.res1.labels', label = TRUE)

# ...2. Strategy 2: Comparing scores across references ----------

hpca$label.main
dice$label.main
hpca$label.main <- gsub('HPCA\\.','', hpca$label.main)
dice$label.main <- gsub('DICE\\.','', dice$label.main)

com.res2 <- SingleR(test = pbmc_counts, 
                    ref = list(HPCA = hpca, DICE = dice),
                    labels = list(hpca$label.main, dice$label.main))

# Check the final label from the combined assignment.
table(com.res2$labels)

# which reference scored best for which label?
grouping <- paste0(com.res2$labels,'.', com.res2$reference)
best_ref <- as.data.frame(split(com.res2, grouping))

# get de. genes from each individual references
metadata(com.res2$orig.results$HPCA)$de.genes
metadata(com.res2$orig.results$DICE)$de.genes

# Combined diagnostics
plotScoreHeatmap(com.res2)


# ...3. Strategy 3: Using Harmonized Labels ----------

hpca.ont <- celldex::HumanPrimaryCellAtlasData(cell.ont = 'nonna')
dice.ont <- celldex::DatabaseImmuneCellExpressionData(cell.ont = 'nonna')

# Using the same sets of genes:
shared <- intersect(rownames(hpca.ont), rownames(dice.ont))
hpca.ont <- hpca.ont[shared,]
dice.ont <- dice.ont[shared,]

# Showing the top 10 most frequent terms:
tail(sort(table(hpca.ont$label.ont)),10)
tail(sort(table(dice.ont$label.ont)), 10)

# using label.ont instead on label.main while running SingleR

com.res3 <- SingleR(test = pbmc_counts,
                    ref = list(HPCA = hpca.ont, DICE = dice.ont),
                    labels = list(hpca.ont$label.ont, dice.ont$label.ont))


table(com.res3$labels)



# How to map cell ontology terms? ----------------

colData(hpca.ont)
colData(dice.ont)

hpca.fle <- system.file("mapping","hpca.tsv", package = "celldex")
hpca.mapping <- read.delim(hpca.fle, header = F)





#DATA SOURCE

#▸ Link to Data:
  #https://www.10xgenomics.com/datasets/20-k-human-pbm-cs-3-ht-v-3-1-chromium-x-3-1-high-6-1-0

#▸ Link to Code:
  #https://github.com/kpatel427/YouTubeTutorials/blob/main/annotateSingleR_multipleRefs.R

#▸ Resources/Vignettes:
  #1. https://www.nature.com/articles/s41596-021-00534-0#sec3
  #2. https://bioconductor.org/books/release/SingleRBook/
  #3. https://bioconductor.org/books/3.14/OSCA.basic/cell-type-annotation.html#assigning-cell-labels-from-reference-data



#---------------------------------------------Automatic cell Annotation (SingleR) multiple Reference-------------------------------------------------------------------
