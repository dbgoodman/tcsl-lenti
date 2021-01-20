library(VISION)
library(data.table)
library(cowplot)
library(Seurat)
library(SeuratDisk)


SEURAT_HDF5_FILE <- here::here(
  '..','..','scrnaseq_tcsl154',
  'rds','scrna.wnn.int.all.clean.h5Seurat')

SEURAT_RDS_FILE <- here::here(
     '..','..','scrnaseq_tcsl154',
     'rds',
     paste0('combined.integrated.umap.data.rds'))

source(here::here('r','scrna','scrna_analysis_fxns.R'))
source(here::here('r','fig_colors.R'))

this_t_type <- 'CD8'


# DO ALL CELLS IN T_TYPE
####################################

scrna_subset <- subset(
  LoadH5Seurat((SEURAT_HDF5_FILE)),
  subset=(t_type == this_t_type))

stopifnot(length(unique(scrna_subset@meta.data$t_type)) == 1)

sct_output <- sct_umap(scrna_subset, title=paste(this_t_type,' ALL'))
scrna_subset <- sct_output$obj


SEURAT_HDF5_ALL_FILE <- here::here(
  '..','..','scrnaseq_tcsl154',
  'rds',
  paste0('scrna.wnn.int.', this_t_type, '.clean.h5Seurat'))

SaveH5Seurat(
  DietSeurat(
    scrna_subset,
    counts=F, data=T, scale.data=T,
    assays=c('SCT_INT','SCT_ADT_INT'),
    dimreducs=c('wnn.umap','adt_pca','rna_pca')),
  filename = SEURAT_HDF5_ALL_FILE, overwrite=T)

sct_output$obj <- NULL

# DO ACTIVATED CLUSTERS FROM MIXED ONLY
##########################################

act_clust_data <- plot_clust_by_act(scrna_subset)
activated_clusters <- act_clust_data$car_clust_data[
  !(car == 'KLRG1' | sample %in% c('Untransduced.stim','Untransduced.baseline')),
  list(pct_clust_car, k_type=k_type[which.max(pct_clust_car)]), c('car', 'seurat_clusters')][,
    sum(k_type == 'cd19+')/.N, by='seurat_clusters'][V1 ==1, seurat_clusters]

scrna_subset <- subset(scrna_subset,
  subset=(t_type == this_t_type & seurat_clusters %in% activated_clusters))

sct_output <- sct_umap(scrna_subset, title=paste(this_t_type,' ACT'))
scrna_subset <- sct_output$obj

cd8_umap_plots <- umap_plot(scrna_cd8, assay='SCT_INT', cluster_name='seurat_clusters', reduction='wnn.umap')
cd8_umap_plots_adt <- umap_plot(scrna_cd8, assay='SCT_ADT_INT', cluster_name='seurat_clusters', reduction='wnn.umap')

SEURAT_HDF5_ACT_FILE <- here::here(
  '..','..','scrnaseq_tcsl154',
  'rds',
  paste0('scrna.wnn.int.', this_t_type, '.act.clean.h5Seurat'))

SaveH5Seurat(
  DietSeurat(
    scrna_subset,
    counts=F, data=T, scale.data=T,
    assays=c('SCT','SCT_ADT'),
    dimreducs=c('wnn.umap','adt_pca','rna_pca')),
  filename = SEURAT_HDF5_ACT_FILE, overwrite=T)

# DO STIM ONLY
####################################

scrna_subset <- subset(
  readRDS(SEURAT_RDS_FILE),
  subset=(t_type == this_t_type & (k_type %in% c('cd19+','bead_stim'))))

#stopifnot(length(unique(scrna_subset@meta.data$t_type)) == 1)

sct_output <- sct_umap(scrna_subset, title=paste(this_t_type,' STIM'))
scrna_subset <- sct_output$obj

SEURAT_HDF5_STIM_FILE <- here::here(
  '..','..','scrnaseq_tcsl154',
  'rds',
  paste0('scrna.wnn.int.', this_t_type, '.stim.clean.h5Seurat'))

SaveH5Seurat(
  DietSeurat(
    scrna_subset,
    counts=F, data=T, scale.data=T,
    assays=c('SCT_INT','SCT_ADT_INT'),
    dimreducs=c('wnn.umap','adt_pca','rna_pca')),
  filename = SEURAT_HDF5_STIM_FILE, overwrite=T)

sct_output$obj <- NULL

# DO ACT CARS ONLY
####################################

scrna_subset <- subset(
  LoadH5Seurat((SEURAT_HDF5_FILE)),
  subset=(t_type == this_t_type & (
    (k_type =='cd19+' & car %in% c('CD28','CD40','TACI','BAFF-R','4-1BB','Zeta')) |
    (k_type == 'bead_stim') & car == 'Untransduced')))

sct_output <- sct_umap(scrna_subset, title=paste(this_t_type,' STIM'))
scrna_subset <- sct_output$obj

# put plots and marker data inside seurat obj
sct_output$obj <- NULL
scrna_subset@misc['sct_output'] <- sct_output

SEURAT_HDF5_ACT_CAR_FILE <- here::here(
  '..','..','scrnaseq_tcsl154',
  'rds',
  paste0('scrna.wnn.int.', this_t_type, '.stim.clean.h5Seurat'))

SaveH5Seurat(
  DietSeurat(
    scrna_subset,
    counts=F, data=T, scale.data=T,
    assays=c('SCT_INT','SCT_ADT_INT'),
    dimreducs=c('wnn.umap','adt_pca','rna_pca')),
  filename = SEURAT_HDF5_STIM_FILE, overwrite=T)

sct_output$obj <- NULL

