##-------------------------------
## Single cell analysis
## Author: Peng Yang
##-------------------------------
library(dplyr)
library(Matrix)
library(Seurat)
library(patchwork)

## Sample list
Samples <- list.files('FilteredData_txt/')
Samples <- unique(unlist(lapply(Samples, function(x) paste0(strsplit(x, '_')[[1]][c(1,2)], collapse = '_'))))

## Sample list for seurat object
Sample_list <- list()
sam = Samples[1]


## Read all the samples
for(sam in Samples){
  
  ## Load data
  dat <- Matrix::readMM(gzfile(paste0('Raw_Data/', sam,'_matrix.mtx.gz')))
  barcodes <- read.table(gzfile(paste0("Raw_Data/", sam,"_barcodes.tsv.gz")))   
  feathers <- read.table(gzfile(paste0('Raw_Data/', sam,'_features.tsv.gz')))
  
  colnames(dat) <- barcodes$V1
  rownames(dat) <- feathers$V2
  
  ## Create Seurat Object
  dat_seurat <- CreateSeuratObject(counts = dat, min.cells = 3, min.features = 200, project = strsplit(sam, '_')[[1]][2])
  
  ## mitochondrial counts
  dat_seurat[["percent.mt"]] <- PercentageFeatureSet(dat_seurat, pattern = "^MT-")
  
  ## QC
  # VlnPlot(dat_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  dat_seurat <- subset(dat_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 15)
  
  ## Pre-processing
  dat_seurat <- NormalizeData(dat_seurat)
  dat_seurat <- FindVariableFeatures(dat_seurat,  selection.method = "vst", nfeatures = 2000)
  
  Sample_list[[strsplit(sam, '_')[[1]][2]]] <- dat_seurat
  
  
}

merged_combined <- reduce(seurat_object_list, merge)


## Run Harmony Integration
options(repr.plot.height = 2.5, repr.plot.width = 6)
merged_combined <- merged_combined %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
## Run clustering
merged_combined <- RunUMAP(merged_combined, reduction = 'harmony', dims = 1:30,
                           reduction.name = 'umap.harmony')
merged_combined <- FindNeighbors(merged_combined, reduction = 'harmony', dims = 1:30)
merged_combined <- FindClusters(merged_combined, resolution = 0.5)
merged_combined <- FindClusters(merged_combined, resolution = 0.1)

save(merged_combined, file = 'Data/merged_data_Harmony_new.RData')