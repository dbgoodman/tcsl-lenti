library(VISION)
library(data.table)
library(cowplot)
library(Seurat)
library(SeuratDisk)
source(here::here('r','scrna','scrna_analysis_fxns.R'))


SEURAT_RDS_FILE <- here::here(
    '..','..','scrnaseq_tcsl154',
    'rds','scrna.hashed.rds')

SEURAT_RDS_WNN_FILE <- here::here(
    '..','..','scrnaseq_tcsl154',
    'rds','scrna.wnn.all.rds')

SEURAT_HDF5_FILE <- here::here(
  '..','..','scrnaseq_tcsl154',
  'rds','scrna.wnn.all.h5Seurat')

scrna_output <- sct_umap(subset(
  readRDS(SEURAT_RDS_FILE),
  subset=(car != 'K_only')), title='All cells')

scrna <- scrna_output$obj

# use only raw counts, put them in .X in adata obj
scrna@assays$RNA@data <- scrna@assays$RNA@counts
scrna <- DietSeurat(scrna, counts = TRUE, data = TRUE, scale.data = FALSE)

SaveH5Seurat(scrna, filename = SEURAT_HDF5_FILE, overwrite=T)
Convert(SEURAT_HDF5_FILE, dest = "h5ad", assay='RNA', overwrite=T)
