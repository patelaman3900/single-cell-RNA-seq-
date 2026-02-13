#Integrating scRNA seq dataset in R using Seurat(CCA)
#Goal of the study: Identification of distinct tumor cell populations and key genetic 
#mechanisms in Hepatoblastoma (HB) through single cell RNA sequencing (scRNA-seq) - GSE180665
#Study design- sample are collected from 3 patients (ID- HB17,HB30,HB53) from different location
#Background liver   #Tumor   #PDX Tumor
#HB17      yes            yes       yes   
#HB30      No             yes       yes
#HB53      yes            yes       No

#Goal of analysis: To integrate data from different patients and correct for batch effects.
#When to integrate?
#multiple scRNA seq data from different sample,condition,Treatment
#cell label transfer from reference to query dataset
#multimodel sc data(scRNA Seq, scATAC seq) to sc multi-omics dataset, signal collected from separate assays
#scRNA Seq and spatial expression data (integrate topological arrangement of cells in tissue with gene expression data)
#Types of Integration 
#Horizontal Integration
#same modality from independent cells
# Eg-scRNA seq from same tissue from different patients/sequencing technology
#Assay are anchored by common gene set
#Vertical Integration
#Multiple modalities profiled simultaneously from same cells
# Eg- scRNA seq and scATAC Seq from same cells
#Assays are anchored by cells
#Diagonal Integration
#Different modalities from the different cells 
# Eg- scRNA and scATAC seq performed on separete group of cells
#Batch correction methods
#MNN, ##Seurat v3, #LIGER, #HArmony, #BBKNN, #scVI, #Conos, #Scmap, #Scanorama, #scAlign

#Required R packages
#Seurat, #tidyverse, #ggplot2, gridExtra

#NCBI- Gene Expression Omnibus (GSE180665)
#choose custom in GSE180665_RAW.tar and download all the filtered_feature_bc_matrix.tar.gz(7 files) 
#Uncompress this GSE180665_RAW.tar file in ubuntu(tar -xvf GSE180665_RAW.tar)
#now uncompress each of the .tar.gz files (for i in *.gz;do tar -xvzf $i; done) to folders (each have 3 .mtz.ga fles)

#setpath
setwd("D:/AMAN_PhD/Script/scRNA_seq/Seurat/Integration")

#Install packages
install.packages("gridExtra")

#load library
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

#get data  location
dirs = list.dirs(path= 'D:/AMAN_PhD/Script/scRNA_seq/Seurat/Integration/',recursive = F,full.names = F)

#create the count matrix
for (x in dirs) {
  #extract clean sample name
  name = gsub('_filtered_feature_bc_matrix','',basename(x))
  #read the 10x files
  cts = ReadMtx(
   mtx = paste0(x,'/matrix.mtx.gz'),
   features = paste0(x,'/features.tsv.gz'),
   cells = paste0(x,'/barcodes.tsv.gz')
   )
#create seurat object
 assign(name,CreateSeuratObject(counts=cts))
 }
 
#merge dataset to perform quality control into one go (no integration)
merged_seurat = merge(HB17_background, y =  c(HB17_PDX,HB17_tumor,HB30_PDX,HB30_tumor,
                              HB53_background,HB53_tumor),
      add.cell.ids=ls()[3:9],
      project= 'HB')

#QC & filtering---------------
View(merged_seurat@meta.data)

#create a sample column
merged_seurat$sample= rownames(merged_seurat@meta.data)

#split the sample column
merged_seurat@meta.data= separate(merged_seurat@meta.data,col = 'sample',
         into = c('Patient','Type','Barcode'),
         sep='_')

#ensure all patients sample merged 
unique(merged_seurat@meta.data$Patient)
#ensure all patients cell source merged
unique(merged_seurat@meta.data$Type)
#calculate mitochondrial percentage

#Calculate the mitochondrial percentage 
merged_seurat$mitoPercent= PercentageFeatureSet(merged_seurat,pattern= '^MT-')

#explore QC

#filtering
merged_seurat_filtered = subset(merged_seurat,subset=nCount_RNA>800 &
        nFeature_RNA>500 &
         mitoPercent<10)

#To check if any batch effects present (standard workflow)
merged_seurat_filtered = NormalizeData(object = merged_seurat_filtered)               #Normalize data
merged_seurat_filtered = FindVariableFeatures(object = merged_seurat_filtered)        #Find variable genes
merged_seurat_filtered = ScaleData(object = merged_seurat_filtered)                   #data scaling 
merged_seurat_filtered = RunPCA(object = merged_seurat_filtered)                      #PCA analysis 
ElbowPlot(merged_seurat_filtered)                                                     #Plot PCA by elbow plot      
merged_seurat_filtered = FindNeighbors(object = merged_seurat_filtered, dims = 1:20)  #Find neighbour features
merged_seurat_filtered = FindClusters(object = merged_seurat_filtered)                #clustering
merged_seurat_filtered =  RunUMAP(object = merged_seurat_filtered, dims = 1:20)       #Run UMAP

#plot on UMAP
p1 = DimPlot(merged_seurat_filtered, reduction ='umap', group.by = 'Patient' )         
p2 = DimPlot(merged_seurat_filtered, reduction ='umap', group.by = 'Type',
             cols = c('red', 'green', 'blue'))
#to see UMAP side by side 
grid.arrange(p1,p2,nrow = 2, ncol = 2)

#plot shows the presence on batch effect
#batch effect correction

obj.list = SplitObject(merged_seurat_filtered , split.by = 'Patient')

#Run Normalization for each object(3) in the list 
#for(i in 1:length(obj.list)){
  #obj.list[[i]] = NormalizeData(object = obj.list[[i]])
  #obj.list[[i]] = FindVariableFeatures(object= obj.list[[i]])
#}

#use lapply for better speed 
obj.list <- lapply(obj.list, NormalizeData)
obj.list <- lapply(obj.list, FindVariableFeatures)

#select integration features 
features = SelectIntegrationFeatures(object.list = obj.list)

#find integration anchors(CCA) to integrate data across different patients
anchors = FindIntegrationAnchors(object.list= obj.list,
                       anchor.features = features)

#integrate data
seurat.integrated = IntegrateData(anchorset = anchors)
        
 #scale data, run PCA and UMAP and visualize integrated data
seurat.integrated = ScaleData(object = seurat.integrated) 
seurat.integrated = RunPCA(object = seurat.integrated)
seurat.integrated = RunUMAP(object = seurat.integrated, dims = 1:50 )

p3 = Dimplot(seurat.integrated, reduction = 'umap', group.by = 'Patient')
p4 = Dimplot(seurat.integrated, reduction = 'umap', group.by = 'Type',
             cols = c("red", "green", "blue"))

#plot on UMAP
grid.arrange(p3, p4, ncol = 2)

#to visualize before and after integration and batch effect correction
grid.arrange(p1, p2, p3, p4, ncol = 2, nrow= 2)



#1) Paper: 
#Modeling Hepatoblastoma: Identification of Distinct Tumor Cell Populations and Key Genetic 
#Mechanisms through Single Cell Sequencing (scRNA-seq): https://pmc.ncbi.nlm.nih.gov/articles/PMC8426487/

#2) Data: 
#GSE180665: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180665

#3) Link to code:
#https://github.com/kpatel427/YouTubeTutorials/blob/main/singleCell_integration.R

#4) Seurat Integration Vignette:
#https://satijalab.org/seurat/articles/integration_introduction.html

#5) Additional resources:
#‣ CCA method: https://doi.org/10.1016/j.cell.2019.05.031
#‣ MNN method: https://doi.org/10.1038/nbt.4091
#‣ https://doi.org/10.1038/s41587-021-00895-7
#‣ https://satijalab.org/seurat/articles...
#‣ https://satijalab.org/seurat/archive/v3.0/integration.html
#‣ https://github.com/hbctraining/scRNA-seq/blob/master/lessons/06_SC_SCT_and_integration.md
