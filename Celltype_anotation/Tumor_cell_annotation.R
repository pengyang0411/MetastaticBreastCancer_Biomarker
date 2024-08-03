##---------------------------------------
## Run inferCNV sample by sample
## Author: Peng Yang
##---------------------------------------
library(Seurat)
library(SeuratData)
library(infercnv)
rm(list=ls(all=TRUE))
args <- commandArgs(trailingOnly = TRUE)
id = as.numeric(args[1])
# id = 6

# id = 12
## Read the file list
Files <- list.files('FilteredData_txt/')
sampleName <- strsplit(Files[id], split = '_filteredData.txt')[[1]]

if(file.exists(paste0('res_inferCNV/', sampleName, '/obj.RData'))){
  
  load(paste0('res_inferCNV/', sampleName, '/obj.RData'))

}else{
  
  # sampleName <- strsplit(Files[id], split = '_filteredData.txt')[[1]]
  
  ## Read the sample
  dat <- read.table(paste0('FilteredData_txt/', sampleName, '_filteredData.txt'))
  ## Create Seurat Object
  obj <- CreateSeuratObject(dat, project = sampleName, assay = "RNA",
                            min.cells = 0, min.features = 0, names.field = 1,
                            names.delim = "_", meta.data = NULL)
  
  ## Run clustering and label immune cells as reference
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, npcs = 45)
  obj <- RunUMAP(obj, reduction = 'pca', dims = 1:30,
                 reduction.name = 'umap.pca')
  obj <- FindNeighbors(obj, reduction = 'pca', dims = 1:30)
  obj <- FindClusters(obj, resolution = 1)
  
  Idents(obj) <- obj$RNA_snn_res.1
  obj$seurat_clusters <- obj$RNA_snn_res.1
  
  dir.create('res_inferCNV')
  dir.create(path = paste0('res_inferCNV/',sampleName))
  
  ## Save the clustering results
  res_cluster <- obj$seurat_clusters
  names(res_cluster) <- paste(sampleName, names(res_cluster), sep = '_')
  save(res_cluster, file = paste0('res_inferCNV/', sampleName, '/res_cluster.RData'))
  
  ## Identify the reference cluster
  Index     <- which(obj@assays$RNA@counts['PTPRC', ] > 5)
  reference <- rep('Unknown', length(obj$orig.ident))
  reference[Index] <- 'Immune'
  obj$reference <- as.factor(reference)
  
  ## Merge the cell type annoations
  obj$cellTypes <- cellTypes[names(res_cluster), 2]
  Epithelial <- which(names(res_cluster) %in% rownames(epithelial_cluster))
  obj$cellTypes[Epithelial] <- 'Epithelial'
  
  ## Visualization
  pt1 <- DimPlot(obj, reduction = 'umap.pca', group.by = 'seurat_clusters', label = TRUE, repel = TRUE)
  
  pt2 <- DimPlot(obj, reduction = 'umap.pca', group.by = 'cellTypes', label = TRUE, repel = TRUE)
  
  ## Making dotplot for certain genes
  pt3 <- DotPlot(obj, features = c('CD3D', 'CD2','PTPRC', ## Tcells
                                   'CD14', 'AIF1', ## Marcrophage
                                   'CD79A',  ## B cell
                                   'COL1A2', 'COL1A1', ## Fibroblast
                                   'VWF', 'CDH5', ## Endothelial
                                   'KRT19', 'EPCAM')) +  ## Epithelial
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  Idents(obj) <- obj$cellTypes
  pt4 <- DotPlot(obj, features = c('CD3D', 'CD2','PTPRC', ## Tcells
                                   'CD14', 'AIF1', ## Marcrophage
                                   'CD79A',  ## B cell
                                   'COL1A2', 'COL1A1', ## Fibroblast
                                   'VWF', 'CDH5', ## Endothelial
                                   'KRT19', 'EPCAM')) +  ## Epithelial
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ## Top marker genes for epithelial cells
  # pt4 <- DotPlot(obj, features = unique(topMarker$`Cell type`)) +  ## Epithelial
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  
  pdf(paste0('res_inferCNV/', sampleName, '/dim_dot_plot.pdf'),
      width = 14, height = 10)
  print(pt1 + pt2 + pt3 + pt4)
  dev.off()
  
  
  pt5 <- DimPlot(obj, reduction = 'umap.pca', group.by = 'reference', label = TRUE, repel = TRUE)
  pt6 <- FeaturePlot(obj, reduction = 'umap.pca', features = 'EPCAM')
  pt7 <- FeaturePlot(obj, reduction = 'umap.pca', features = 'PTPRC')
  
  
  ##
  pdf(paste0('res_inferCNV/', sampleName, '/dim_feature_plot.pdf'),
      width = 10, height = 8)
  print(pt1 + pt4 + pt5 + pt6)
  dev.off()
  
  # pt1 + pt2 + pt3 
  
  table(obj$seurat_clusters)

  
  save(obj, file = paste0('res_inferCNV/', sampleName, '/obj.RData'))
  
}

load('Data/merged_reference.RData')


## Merge the merged_refernce with obj
obj.combined <- merge(x = obj, y = merged_reference, 
                      add.cell.ids = c(sampleName, 'Reference'))

Indx <- which(is.na(obj.combined$seurat_clusters))
obj.combined$seurat_clusters[Indx] <- 'Reference'


# ## Merge the reference to copykat
# obj$copyKat_cluster[obj$reference == 'Immune'] <- 'Immune'
# obj$copyKat_cluster <- as.factor(obj$copyKat_cluster)


## Run inferCNV
gglist <- read.table(file = "Data/geneOrderingFile.txt")
tmp <- cbind(colnames(obj.combined@assays$RNA@counts),
             as.character(obj.combined$seurat_clusters))
rownames(tmp) <- NULL
head(tmp)
write.table(tmp, file = paste0('res_inferCNV/', sampleName, '/cellAnnotationFile.txt'), row.names = FALSE,
            col.names = FALSE, quote = FALSE, sep = "\t")

# library(infercnv)
# out_dir = "/inferCNV_PCA_epithelial"
out_dir = paste0('/res_inferCNV/', sampleName)
infercnv_obj = infercnv::CreateInfercnvObject(
  raw_counts_matrix=obj.combined@assays$RNA@counts,
  annotations_file=paste0('res_inferCNV/', sampleName, '/cellAnnotationFile.txt'),
  delim="\t",
  gene_order_file="Data/geneOrderingFile.txt",
  ref_group_names=c("Reference")) ## PTPRC positive

if(id == 12){
   if.HMM = FALSE
}else{
   if.HMM = TRUE
}

infercnv_obj_default = infercnv::run(
  infercnv_obj,
  cutoff = 0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
  out_dir=paste0('res_inferCNV/', sampleName),
  cluster_by_groups=TRUE,
  plot_steps=FALSE,
  denoise=TRUE,
  HMM=if.HMM,
  num_threads = 18,
  HMM_report_by = c("subcluster"),
  analysis_mode = "subclusters",
  tumor_subcluster_partition_method = "random_trees",
  no_prelim_plot=TRUE,
  png_res=300,
  hclust_method='ward.D2'
)




