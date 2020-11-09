library(VISION)
library(data.table)
library(cowplot)
library(Seurat)
library(SeuratDisk)


SEURAT_HDF5_FILE <- here::here(
  '..','..','scrnaseq_tcsl154',
  'rds','scrna.wnn.all.clean.h5Seurat')

SEURAT_RDS_FILE <- here::here(
    '..','..','scrnaseq_tcsl154',
    'rds','scrna.hashed.rds')

source(here::here('r','scrna','scrna_analysis_fxns.R'))
source(here::here('r','fig_colors.R'))

this_t_type <- 'CD4'

# DO ALL CELLS IN T_TYPE
####################################

scrna_subset <- subset(
  LoadH5Seurat(file = SEURAT_HDF5_FILE),
  subset=(CD4v8 == this_t_type))

stopifnot(length(unique(scrna_subset@meta.data$CD4v8)) == 1)

sct_output <- sct_umap(scrna_subset, title=paste(this_t_type,' ALL'))
scrna_subset <- sct_output$obj


SEURAT_HDF5_ALL_FILE <- here::here(
  '..','..','scrnaseq_tcsl154',
  'rds',
  paste0('scrna.wnn.', this_t_type, '.clean.h5Seurat'))

SaveH5Seurat(
  DietSeurat(
    scrna_subset,
    counts=F, data=T, scale.data=T,
    assays=c('SCT','SCT_ADT'),
    dimreducs=c('wnn.umap','adt_pca','rna_pca','umap.hto')),
  filename = SEURAT_HDF5_ALL_FILE, overwrite=T)
sct_output$obj <- NULL

# DO ACTIVATED CLUSTERS FROM MIXED ONLY
##########################################
act_clust_data <- plot_clust_by_act(scrna_subset)
activated_clusters <- act_clust_data$car_clust_data[
  !(car %in% c('KLRG1','Untransduced')),
  list(pct_clust_car, k_type=k_type[which.max(pct_clust_car)]), c('car', 'seurat_clusters')][,
    sum(k_type == 'cd19+')/.N, by='seurat_clusters'][V1 ==1, seurat_clusters]

scrna_subset <- subset(scrna_subset,
  subset=(CD4v8 == this_t_type & seurat_clusters %in% activated_clusters))

sct_output <- sct_umap(scrna_subset, title=paste(this_t_type,' ACT'))
scrna_subset <- sct_output$obj

SEURAT_HDF5_ACT_FILE <- here::here(
  '..','..','scrnaseq_tcsl154',
  'rds',
  paste0('scrna.wnn.', this_t_type, '.act.clean.h5Seurat'))

SaveH5Seurat(
  DietSeurat(
    scrna_subset,
    counts=F, data=T, scale.data=T,
    assays=c('SCT','SCT_ADT'),
    dimreducs=c('wnn.umap','adt_pca','rna_pca','umap.hto')),
  filename = SEURAT_HDF5_ACT_FILE, overwrite=T)

# DO STIM ONLY
####################################

scrna_subset <- subset(
  LoadH5Seurat(file = SEURAT_HDF5_FILE),
  subset=(CD4v8 == this_t_type & k_type == 'cd19+'))

stopifnot(length(unique(scrna_subset@meta.data$CD4v8)) == 1)

sct_output <- sct_umap(scrna_subset, title=paste(this_t_type,' STIM'))
scrna_subset <- sct_output$obj

SEURAT_HDF5_STIM_FILE <- here::here(
  '..','..','scrnaseq_tcsl154',
  'rds',
  paste0('scrna.wnn.', this_t_type, '.stim.clean.h5Seurat'))

SaveH5Seurat(
  DietSeurat(
    scrna_subset,
    counts=F, data=T, scale.data=T,
    assays=c('SCT','SCT_ADT'),
    dimreducs=c('wnn.umap','adt_pca','rna_pca','umap.hto')),
  filename = SEURAT_HDF5_STIM_FILE, overwrite=T)

sct_output$obj <- NULL

