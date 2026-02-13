#----------------------------Reading Single-cell data into a Seurat object in R----------------------------
#Requirement
#R(v4.1.1)
#RStudio(v1.4.1717)
#R packages
#Seurat(v4.0.6)
#SeuratDisk(v0.0.0.9019)

#Feature -carcode sparse matrix
#Rows= feature(gene)
#column= cell barcode  
#values= counts

#Different Input formats
#Name            #Extension 
#10x hdf5               .hdf5      
#R Data Fomat           .rds
#AnnData Object         .h5ad 
#Loom                   .loom
#text based MArket      .mtx
#Exchange Format(MEX)   .mtx 


#Link to scRNA-sequencing technologies paper: https://www.ncbi.nlm.nih.gov/pmc/arti...
#----------------------------Reading Single-cell data into a Seurat object in R----------------------------


##--------------------------------The Script for reading the file--------------------------------------------
#Reading single cell matrices in various format
#and converting to Seurat object
setwd("D:/AMAN_PhD/Script/scRNA_seq/Seurat")

#instal packages
install.packages("Seurat")
#install.packages("SeuratDisk")
install.packages("remotes")
remotes::install_github("mojaveazure/seurat-disk")
install.packages("BiocManager")
BiocManager::install("rhdf5")
install.packages("hdf5r", type = "binary")
install.packages("tidyverse")
install.packages("umap")


#load library
library("Seurat")
library(SeuratDisk)
library(tidyverse)
library("ggplot2")
library(umap)
#Reading each type of file
## .rds format (R format)
rds_obj = readRDS("ependymal_cells.rds")
str(rds_obj)
head(rds_obj)


## .HDF5 format (10X CellRanger) 
hdf5_obj = Read10X_h5(filename="20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5",
           use.names = TRUE,
           unique.features = TRUE)
#check the data in the count matrix by
hdf5_obj[1:10,1:10]
#convert this to seurat object
seurat_hdf_5 = CreateSeuratObject(counts = hdf5_obj)
#Confirm the Seurat object 
str(seurat_hdf_5)


## .mtx file (CellRanger Output)
mtx_obj = ReadMtx(mtx = "path/file_name.mtx.gz",
        Features="path/features.tsv.gz",
        Cells="path/barcode.tsv.gz")
#convert this count matrix to Seurat object
seurat_mtx = CreateSeuratObject(counts = mtx_obj)
#check the data in the count matrix by
seurat_mtx[1:10,1:10]
#Confirm the Seurat object 
str(seurat_mtx)


## .loom files(large sc dataset)
loom_obj = Connect(filename="file_name.loom",
        mode ='r')
#convert this count matrix to Seurat object
seurat_loom = as.Seurat(loom_obj)
#check the data in the count matrix by
seurat_loom[1:10,1:10]
#Confirm the Seurat object 
str(seurat_loom)


## .h5ad format (scanpy package)
#step 1: convert AnnData object to an h5Seurat file
Convert("adata_SS2_for_download.h5ad",
        dest="h5seurat",overwrite=TRUE)
 #step 2: Load h5Seurat file into a Seurat object
seurat_andata = LoadH5Seurat("adata_SS2_for_download.h5seurat")
#check the data in the count matrix by
seurat_andata[1:10,1:10]
#Confirm the Seurat object 
str(seurat_andata)


###The standard workflow to analyze scRNA seq data
#Count Matrix → QC and filtering → Normalization → Identify highly variable genes → Scale data → 
#Linear dimensionality reduction (PCA) → Clustering → Non-linear dimensionality reduction (UMAP/t-SNE)
#R packages
  #Seurat
  #tidyverse
#(Gene Expression - Feature / cell matrix HDF5 (raw))- dataset used below
#https://www.10xgenomics.com/datasets/20-k-mixture-of-nsclc-dt-cs-from-7-donors-3-v-3-1-3-1-standard-6-1-0
#Batch file can also be used 

#Load the NSCLC datasets
nsclc.sparse.m = Read10X_h5("20k_NSCLC_DTC_3p_nextgem_donor_1_count_sample_feature_bc_matrix.h5")
str(nsclc.sparse.m)
head(nsclc.sparse.m)
cts = nsclc.sparse.m$`Gene Expression`  #extract gene expression counts

#initilize the seurat object with the raw(non-normalized data)
nsclc.seurat.obj = CreateSeuratObject(counts = cts, 
                                      project="NSCLC", 
                                      min.features =200 )
nsclc.seurat.obj
str(nsclc.seurat.obj)

# 1. QC------------
View(nsclc.seurat.obj@meta.data)

# % Mitochondrial Reads
#PercentageFeatureSet()= percent.mt= (UMI count from mt genes/total UMI counts)*100
# ^MT= Annotation used for mitochondrial genes is seurat
#nsclc.seurat.obj[["percent.mt"]] = stores the vector inside the seurat object's metadata

nsclc.seurat.obj[["percent.mt"]] = PercentageFeatureSet(nsclc.seurat.obj,pattern = "^MT")
View(nsclc.seurat.obj@meta.data)
#Violin Plot of 
VlnPlot(nsclc.seurat.obj,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#scatter plot(x= feature1,y=feature2)
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
  geom_smooth(method='lm')

#(Use DoubletFinder package to remove duplicate cells)
# 2. Filtering------------ 
nsclc.seurat.obj = subset(nsclc.seurat.obj,
                          subset= nFeature_RNA>200 & 
                            nFeature_RNA <2500 &
                            percent.mt < 5)

# 3. Normalize data----------
#nsclc.seurat.obj = NormalizeData(nsclc.seurat.obj, normalization.method = "LogNormalize", scale.factor  = 1000)
#OR as above are default value
nsclc.seurat.obj = NormalizeData(nsclc.seurat.obj)
#To see which commands have been used (check in commands slot)
str(nsclc.seurat.obj)

#4. Identify highly variable features--------
#vst= Variance Stabilizing Transformation (mean variance relationship for significantly biological genes)
nsclc.seurat.obj = FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)

#Identify the 10 most highly variable genes
top10 = head(VariableFeatures(nsclc.seurat.obj),10)

 

#Plot variable features with and without labels
plot1 = VariableFeaturePlot(nsclc.seurat.obj)
#repel = avoid text overlapping
plot2 = LabelPoints(plot= plot1, points = top10, repel= TRUE, xnudge = 0, ynudge = 0)

# 5. Scaling the data to remove the sources of variation
all.genes = rownames(nsclc.seurat.obj) 
nsclc.seurat.obj= ScaleData(nsclc.seurat.obj, features = all.genes)
str(nsclc.seurat.obj)

# 6. Perform Linear dimensionality reduction------
nsclc.seurat.obj= RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))

# Visualize PCA results
print(nsclc.seurat.obj[["pca"]],dims=1:5, nfeatures = 5)
DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE)

# Determine dimensionality of the data(remove PC with low variance)
ElbowPlot(nsclc.seurat.obj)

# 7. Clustering----------
nsclc.seurat.obj = FindNeighbors(nsclc.seurat.obj, dims = 1:15)

#Understanding Resolution (low res=few clusters)
nsclc.seurat.obj = FindClusters(nsclc.seurat.obj, resolution = c(0.1,0.3,0.5,0.7,1))
View(nsclc.seurat.obj@meta.data)

DimPlot(nsclc.seurat.obj, group.by =  "RNA_snn_res.0.5", label = TRUE)

#Setting Identity of clusters
Idents(nsclc.seurat.obj)
Idents(nsclc.seurat.obj)= "RNA_snn_res.0.5"
Idents(nsclc.seurat.obj)

# Non-linear Dimensionality reduction
nsclc.seurat.obj= RunUMAP(nsclc.seurat.obj, dims = 1:15)
#Label can be set as TRUE or use LabelCluster function to help label individual clusters
DimPlot(nsclc.seurat.obj, reduction = "umap")
  
##--------------------------------The Script for reading the file--------------------------------------------