#-------------------------------------------------------Doublet Finder---------------------------------------------------

#Homotypic doublets: doublets derived from transcriptionally similar cells
#Heterotypic doublets: doublets derived from transcriptionally distinct cells
#Doublet finder is more sensitive to Heterotypic doublet than homotypic

#DoubletFinder needs 3 parameters:
#pN = the number of artificial doublets(default 0.25 or 25%)
#pK = the neighborhood size (pK) used to compute the number of artificial nearest neighbors
#Exp = the number of expected real doublets

#For example, in a dataset with 15,000 real cells, a pN of 0.25 would represent the integration 
# of 5,000 artificial doublets, and a pk of 0.01 would represent a pK of 200 cells.

#How does DoubletFinder work?
#1.Simulate artificial doublets from existing scRNA seq data by averaging the gene expression profiles of random pairs of cells 

#2.Merges and pre-processes real and artificial data using the "Seurat" single-cell analysis pipeline

#3.Performs dimensionality reduction on the merged real-artificial data using PCA, producing a low dimensional space that 
# describes the similarity between real and artificial cells.

#4.Detects the nearest neighbors for every real cell in principal component (PC) space, and is used to compute each cell's 
# proportion of artificial nearest neighbors (pANN) (Highly dependent on pK)

#5.Predicts real doublets as cells with the top n pANN values, where in is set to the total number of expected doublets

#Original Data → Simulate Doublets → PCA → Define Neighbours → Threashold pANN → Doublets Removed

#Expected number of doubletts(Exp):-to be given by the reagent provider(eg. 10X genomics)
#Multiplet Rate(%)     No.of Cells Loaded     No. of Cells Recovered
#0.4%	                   800	                     500
#0.8%	                   1,600	             1,000
#1.6%	                   3,200	             2,000
#2.3%                       4,800	             3,000
#3.1%	                   6,400	             4,000
#3.9%	                   8,000	             5,000
#4.6%	                   9,600	             6,000
#5.4%	                   11,200	             7,000
#6.1%	                   12,800	             8,000
#6.9%	                   14,400	             9,000
#7.6%	                   16,000	             10,000

##########Best practices##########
#DoubletFinder not to be applied on aggregated scRNA-seq data
#Not preferable to run on merged data(different samples may have different proportions of cell types and merged objects can be large in size)
#Should be run on distinct samples separately
#Input data should be cleared of low-quality cells
#Remove clusters with low RNA UMIs, high mitochondrial read % and uninformative marker genes 


#R packages:
#DoubletFinder  
#Seurat  
#tidyverse  
#ggplot2


#1)Data: 
# https://www.10xgenomics.com/datasets/10k-human-pbmcs-3-v3-1-chromium-controller-3-1-high
# (Feature / cell matrix (raw)):- 10k_PBMC_3p_nextgem_Chromium_Controller_raw_feature_bc_matrix.tar.gz


#Uncompress the file 
#cd /mnt/d/AMAN_PhD/Script/scRNA_seq/Seurat/Doublets
#tar -xvzf 10k_PBMC_3p_nextgem_Chromium_Controller_raw_feature_bc_matrix.tar.gz




#Filter out doublets: DoubletFinder

#setpath
setwd("D:/AMAN_PhD/Script/scRNA_seq/Seurat/Doublets")

#install package
install.packages("remotes")
install.packages("DoubletFinder")
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")

# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(DoubletFinder)

#create counts matrix
cts <- ReadMtx(mtx = 'D:/AMAN_PhD/Script/scRNA_seq/Seurat/Doublets/raw_feature_bc_matrix/matrix.mtx.gz',
               features = 'D:/AMAN_PhD/Script/scRNA_seq/Seurat/Doublets/raw_feature_bc_matrix/features.tsv.gz',
               cells = 'D:/AMAN_PhD/Script/scRNA_seq/Seurat/Doublets/raw_feature_bc_matrix/barcodes.tsv.gz')

cts[1:10,1:10]

# create Seurat object
pbmc.seurat <- CreateSeuratObject(counts = cts)
str(pbmc.seurat)

# QC and Filtering
# explore QC

#Mitochondrial reads
pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = '^MT-')

#Filtered reads
pbmc.seurat.filtered <- subset(pbmc.seurat, subset = nCount_RNA > 800 &
                                 nFeature_RNA > 500 &
                                 mitoPercent < 10)

pbmc.seurat
pbmc.seurat.filtered

# pre-process standard workflow
pbmc.seurat.filtered <- NormalizeData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindVariableFeatures(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- ScaleData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunPCA(object = pbmc.seurat.filtered)
ElbowPlot(pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindNeighbors(object = pbmc.seurat.filtered, dims = 1:20)
pbmc.seurat.filtered <- FindClusters(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:20)


## pK Identification (no ground-truth)(Require R version 4) -------------------------------------------------
sweep.res.list_pbmc <- paramSweep_v3(pbmc.seurat.filtered, PCs = 1:20, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)

ggplot(bcmvn_pbmc, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_pbmc %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))


## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- pbmc.seurat.filtered@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.076*nrow(pbmc.seurat.filtered@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
pbmc.seurat.filtered <- doubletFinder_v3(pbmc.seurat.filtered, 
                                         PCs = 1:20, 
                                         pN = 0.25, 
                                         pK = pK, 
                                         nExp = nExp_poi.adj,
                                         reuse.pANN = FALSE, sct = FALSE)


# visualize doublets
DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = "DF.classifications_0.25_0.21_691")


# number of singlets and doublets
table(pbmc.seurat.filtered@meta.data$DF.classifications_0.25_0.21_691)




#Data set
#1)Data: 
# https://www.10xgenomics.com/datasets/10k-human-pbmcs-3-v3-1-chromium-controller-3-1-high
# (Feature / cell matrix (raw))

#2)Link to code
# https://github.com/kpatel427/YouTubeTutorials/blob/main/singleCell_doublets.R

#3)DoubletFinder github:
# https://github.com/chris-mcginnis-ucsf/DoubletFinder

#4)DoubletFinder paper:
# https://doi.org/10.1016/j.cels.2019.03.003



#-------------------------------------------------------Doublet Finder---------------------------------------------------
