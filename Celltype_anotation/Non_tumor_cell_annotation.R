##---------------------------------------
## Non-tumor cell annotations
## Author: Peng Yang
##---------------------------------------

library(purrr)
library(Seurat)
library(SeuratData)
library(dplyr)
library(patchwork)
library(harmony)

library(SingleR)
library(SingleCellExperiment)
library(scuttle)
library(scran)
load('Data/merged_data_Harmony_new.RData')

Idents(merged_combined) <- merged_combined$Tumor_cells
combined_nontumor <- subset(merged_combined, idents = 'non-tumor cells')

## Differeniate early vs late progressor
Status_new <- as.character(combined_nontumor$Status)
Status_new[Status_new == 'Progressor'] = 'Late Progressor'
Status_new[combined_nontumor$orig.ident %in% c('PA11', 'PA46', 'PA144')] = 'Early Progressor'
combined_nontumor$Status_new <- as.factor(Status_new)

save(combined_nontumor, file = 'Data/combined_nonTumor.RData')


## Supervised annotate cell types from pbmc
if(file.exists('Data/pbmc_combine_nonTumor_new.RData')){
  load('Data/pbmc_combine_nonTumor_new.RData')
}else{
  # rm(merged_combined)
  library(SeuratDisk)
  ## Read the reference panel
  reference <- LoadH5Seurat("Data/pbmc_multimodal.h5seurat")
  
  anchors <- FindTransferAnchors(
    reference = reference,
    query = combined_nontumor,
    normalization.method = "LogNormalize",
    reference.reduction = "spca",
    dims = 1:50
  )
  
  combined_nontumor <- MapQuery(
    anchorset = anchors,
    query = combined_nontumor,
    reference = reference,
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      predicted_ADT = "ADT"
    ),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )
  
  
  save(combined_nontumor, file = 'Data/pbmc_combine_nonTumor_new.RData')
}

## Supervised annotate immune cell types from breast
if(file.exists('Data/pbmc_combine_nonTumor_new_breat_ref.RData')){
  load('Data/pbmc_combine_nonTumor_new_breat_ref.RData')
}else{
  # rm(merged_combined)
  library(SeuratDisk)
  ## Read the reference panel
  
  load('Data/BRCA_SCdata_asRef.RData')
  reference <- Seurat::CreateSeuratObject(counts = BRCA_countmat, 
                                          meta.data = BRCA_metadata)
  reference <- FindVariableFeatures(reference)
  reference <- ScaleData(reference)
  reference <- RunPCA(reference, npcs = 50)
  reference <- RunUMAP(reference, reduction = 'pca', dims = 1:30,
                       reduction.name = 'umap.pca.reference')
  reference <- FindNeighbors(reference, reduction = 'pca', dims = 1:30)
  reference <- FindClusters(reference, resolution = 0.5)
  
  reference$celltype_major <- as.factor(BRCA_metadata$celltype_major)
  reference$celltype_minor <- as.factor(BRCA_metadata$celltype_minor)
  
  anchors <- FindTransferAnchors(
    reference = reference,
    query = combined_nontumor,
    normalization.method = "LogNormalize",
    reference.reduction = "pca",
    dims = 1:50
  )
  
  combined_nontumor <- MapQuery(
    anchorset = anchors,
    query = combined_nontumor,
    reference = reference,
    refdata = list(
      celltype_major = "celltype_major",
      celltype_minor = "celltype_minor"
    ),
    reference.reduction = "pca"#,
    # reduction.model = "umap.pca"
  )
  
  predict_celltypes <- data.frame(celltype_major = combined_nontumor$predicted.celltype_major,
                                  celltype_minor = combined_nontumor$predicted.celltype_minor)
  
  combined_nontumor$celltype_major <- as.factor(predict_celltypes$celltype_major)
  combined_nontumor$celltype_minor <- as.factor(predict_celltypes$celltype_minor)
  
  save(combined_nontumor, file = 'Data/pbmc_combine_nonTumor_new_breat_ref.RData')
  
}


##---------------------------------------
## Final anaotation for non-tumor cells 
##---------------------------------------

if(file.exists('Data/combined_nonTumor_final.RData')){
  
  load('Data/combined_nonTumor_final.RData')
  
}else{
  
  load('Data/pbmc_combine_nonTumor_new.RData')
  load('Data/pbmc_combine_nonTumor_new_breat_ref.RData')
  
  combined_nontumor$celltype_major <- as.factor(predict_celltypes$celltype_major)
  combined_nontumor$celltype_minor <- as.factor(predict_celltypes$celltype_minor)
  
  ## ------------------------------------------------------------------------------------
  # Regarding the annotation I2, change the name “CD14 Mono” to “Macrophage.” 
  # In addition, we do not want the “HSPC” but want the “CAFs,” “PVL, “ and “Endothelial” 
  # annotations (from major referenced by breast). 
  ## ------------------------------------------------------------------------------------
  
  table(combined_nontumor$predicted.celltype.l1)
  table(combined_nontumor$predicted.celltype.l2)
  
  CT_L2 <- as.character(combined_nontumor$predicted.celltype.l2)
  CT_L2[which(CT_L2 == 'CD14 Mono')] = 'Macrophage'
  
  CT_L2[which(combined_nontumor$celltype_major == 'CAFs')] = 'CAFs'
  CT_L2[which(combined_nontumor$celltype_major == 'Endothelial')] = 'Endothelial'
  CT_L2[which(combined_nontumor$celltype_major == 'PVL')] = 'PVL'
  
  table(CT_L2)
  
  CT_L1 <- as.character(combined_nontumor$predicted.celltype.l1)
  
  CT_L1[which(combined_nontumor$celltype_major == 'CAFs')] = 'CAFs'
  CT_L1[which(combined_nontumor$celltype_major == 'Endothelial')] = 'Endothelial'
  CT_L1[which(combined_nontumor$celltype_major == 'PVL')] = 'PVL'
  
  
  table(CT_L1)
  
  combined_nontumor$CT_L1 <- as.factor(CT_L1)
  combined_nontumor$CT_L2 <- as.factor(CT_L2)
  
  save(combined_nontumor, file = 'Data/combined_nonTumor_final.RData')
  
  
}