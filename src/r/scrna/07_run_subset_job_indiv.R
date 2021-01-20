library(VISION)
library(data.table)
library(cowplot)
library(Seurat)
library(SeuratDisk)

source(here::here('r','scrna','scrna_analysis_fxns.R'))
source(here::here('r','fig_colors.R'))

# set from global environment if being run as job
# SEURAT_HDF5_INFILE <- here::here(
#   '..','..','scrnaseq_tcsl154',
#   'rds','scrna.wnn.int.all.clean.h5Seurat')
# 
# SCRNA_SUBSET_FXN <- function(scrna)
#   subset(scrna,
#   subset=(t_type == 'CD8' & (
#     (k_type =='cd19+' & car %in% c('CD28','CD40','TACI','BAFF-R','4-1BB','Zeta')) |
#     (k_type == 'bead_stim') & car == 'Untransduced')))
# 
# SEURAT_HDF5_OUTFILE <- here::here(
#   '..','..','scrnaseq_tcsl154',
#   'rds',
#   paste0('scrna.wnn.int.', 'CD8.act_car', '.clean.h5Seurat'))
# TITLE <- ''

message(paste(c('Outfile:', SEURAT_HDF5_OUTFILE)))

# get rds filename for saving plots and markers
SEURAT_RDS_OUTFILE <- paste0(tools::file_path_sans_ext(SEURAT_HDF5_OUTFILE),'.rds')

scrna_subset <- SCRNA_SUBSET_FXN(SEURAT_OBJ, THIS_T_TYPE)

# Run scTransform UMAP on subset
sct_output <- sct_umap(scrna_subset, title=TITLE)
scrna_subset <- sct_output$obj

# put plots and marker data inside seurat obj
sct_output$obj <- NULL

message(paste(c('RDS Size:', object.size(sct_output))))

saveRDS(sct_output, file=SEURAT_RDS_OUTFILE)

SaveH5Seurat(
  {
    DietSeurat(
      scrna_subset,
      counts=F, data=T, scale.data=T,
      assays=c('SCT_INT','SCT_ADT_INT'),
      dimreducs=c('wnn.umap','adt_pca','rna_pca'))
  },
  filename = SEURAT_HDF5_OUTFILE, overwrite=T)