#----------------------------Pseudo-bulk differential expression analysis using scRNA Seq data---------------------------------------------
#Pseudo bulk analysis
#Aggregating the counts and metadata to the sample/replicate level.
#To leverage existing robust bulk RNA-seq DE frameworks, such as DESeq2, edgeR and limma
#Why perform pseudo-bulk analysis?
  #Single-cell RNAseq data tend to exhibit an abundance of zero counts, a complicated distribution, and huge heterogeneity
  #The heterogeneity within and between cell populations manifests major challenges to the differential gene expression 
    #analysis in scRNAseq data.
  #Single-cell methods to identify highly expressed genes as DE and exhibit low sensitivity for genes having low expression
  #Single cell methods often inflate the p-values as each cell is treated as a sample
  #If cells are treated as samples, then variation across a population is not truly investigated
  #DE testing performed on "pseudo-bulk" expression profiles leverages the resolution offered by single-cell technologies to 
    #define the labels and combines it with the statistical rigor of existing methods for DE analyses
  #Each sample is represented no more than once for each condition, avoiding problems from unmodelled correlations between samples.
  #To infer which genes might be important for a condition at the population level (not the individual level), samples need to be 
   #acquired from different organisms/samples, not different cells.

#Required packages
 #ExperimentHub
 #Seurat
 #SESeq2
 #tidyverse


#SCRIPT

#install packages
install.packages("experiment")
install.packages("BiocManager")
BiocManager::install("SparseArray")
BiocManager::install(c(
  "SummarizedExperiment",
  "S4Vectors",
  "IRanges",
  "GenomeInfoDb",
  "Biobase"
), force=TRUE)
BiocManager::install("DESeq2")

#set path
setwd("D:/AMAN_PhD/Script/scRNA_seq/Seurat/Pseudo_bulk_DEG")

#load library
library(ExperimentHub)
library(Seurat)
library(DESeq2)
library(tidyverse)


# get data
eh <- ExperimentHub()
query(eh, "Kang")

sce <- eh[["EH2259"]]
seu.obj <- as.Seurat(sce, data = NULL)
View(seu.obj@meta.data)


# QC and filtering
# explore QC


# get mito percent
seu.obj$mitoPercent <- PercentageFeatureSet(seu.obj, pattern = '^MT-')
View(seu.obj@meta.data)

# filter
seu.filtered <- subset(seu.obj, subset = nFeature_originalexp > 200 & nFeature_originalexp < 2500 &
                         nCount_originalexp > 800 & 
                         mitoPercent < 5 &
                         multiplets == 'singlet')

seu.obj
seu.filtered

# run Seurat's standard workflow steps
seu.filtered <- NormalizeData(seu.filtered)
seu.filtered <- FindVariableFeatures(seu.filtered)
seu.filtered <- ScaleData(seu.filtered)
seu.filtered <- RunPCA(seu.filtered)
ElbowPlot(seu.filtered)
seu.filtered <- RunUMAP(seu.filtered, dims = 1:20)

# visualize 
cell_plot <- DimPlot(seu.filtered, reduction = 'umap', group.by = 'cell', label = TRUE)
cond_plot <- DimPlot(seu.filtered, reduction = 'umap', group.by = 'stim')

cell_plot|cond_plot

# pseudo-bulk workflow -----------------
# Acquiring necessary metrics for aggregation across cells in a sample
# 1. counts matrix - sample level
# counts aggregate to sample level

View(seu.filtered@meta.data)
seu.filtered$samples <- paste0(seu.filtered$stim, seu.filtered$ind)

DefaultAssay(seu.filtered)

cts <- AggregateExpression(seu.filtered, 
                           group.by = c("cell", "samples"),
                           assays = 'originalexp',
                           slot = "counts",
                           return.seurat = FALSE)

cts <- cts$originalexp


# transpose
cts.t <- t(cts)


# convert to data.frame
cts.t <- as.data.frame(cts.t)

# get values where to split
splitRows <- gsub('_.*', '', rownames(cts.t))


# split data.frame
cts.split <- split.data.frame(cts.t,
                              f = factor(splitRows))

# fix colnames and transpose

cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
})

#gsub('.*_(.*)', '\\1', 'B cells_ctrl101')



# Let's run DE analysis with B cells
# 1. Get counts matrix
counts_bcell <- cts.split.modified$`B cells`


# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_bcell))

colData <- colData %>%
  mutate(condition = ifelse(grepl('stim', samples), 'Stimulated', 'Control')) %>%
  column_to_rownames(var = 'samples')

# get more information from metadata




# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_bcell,
                              colData = colData,
                              design = ~ condition)

# filter
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds)



# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds, name = "condition_Stimulated_vs_Control")
res



#DATA
#1) Link to code:
  #https://github.com/kpatel427/YouTubeTutorials/blob/main/singleCell_pseudoBulk.R


#2) Vignettes:
  #▸ https://hbctraining.github.io/scRNA-s...
  #▸ http://biocworkshops2019.bioconductor...
  #▸ https://bioconductor.org/books/3.14/OSCA.multisample/multi-sample-comparisons.html#motivation-2
  #▸ https://bioconductor.org/books/3.14/OSCA.multisample/multi-sample-comparisons.html#creating-pseudo-bulk-samples

#3) Papers:
  #▸ https://doi.org/10.1038/s41467-021-25960-2
  #▸ https://doi.org/10.1186/s12859-019-2599-6
  #▸ https://doi.org/10.1186/s13059-018-1406-4


#----------------------------Pseudo-bulk differential expression analysis using scRNA Seq data---------------------------------------------

