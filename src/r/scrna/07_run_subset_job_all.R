library(VISION)
library(data.table)
library(cowplot)
library(Seurat)
library(SeuratDisk)
library(rstudioapi)

# {CD4, CD8} x {All, Act, Act CARs}
# CLEAR ENVIRONMENT FIRST: rm(list=ls())

t_types <- c('CD4','CD8')

SEURAT_HDF5_INFILE <- here::here(
  '..','..','scrnaseq_tcsl154',
  'rds','scrna.wnn.int.all.clean.h5Seurat')

SEURAT_OBJ <- LoadH5Seurat(SEURAT_HDF5_INFILE)
act_clust_data <- plot_clust_by_act(SEURAT_OBJ)
act_clusters <- unique(act_clust_data$car_clust_data[cd19_clust_act_l2fc > 1, seurat_clusters])
  
scrna_subset_fxns <- list(
  'all' = function(scrna, this_t_type) subset(scrna, subset=(t_type == this_t_type)),
  'act' = function(scrna, this_t_type) subset(scrna, subset=(t_type == this_t_type & 
    (k_type =='cd19+' | k_type == 'bead_stim'))),
  'act_cars'= function(scrna, this_t_type) subset(scrna, subset=(t_type == this_t_type & 
    ((k_type =='cd19+' & car %in% c('CD28','CD40','TACI','BAFF-R','4-1BB','Zeta','CD30')) |
    (k_type == 'bead_stim' & car == 'Untransduced')))),
  'act_clust'= function(scrna, this_t_type) subset(scrna, subset=(t_type == this_t_type & 
    ((k_type =='cd19+' & car %in% c('CD28','CD40','TACI','BAFF-R','4-1BB','Zeta','CD30')) |
    (k_type == 'bead_stim' & car == 'Untransduced')) & seurat_clusters %in% act_clusters))
)



for (this_t_type in t_types) {
  for (this_subset_fxn in names(scrna_subset_fxns)) {
  
    SCRNA_SUBSET_FXN <- scrna_subset_fxns[[this_subset_fxn]]
    
    SEURAT_HDF5_OUTFILE <- here::here(
      '..','..','scrnaseq_tcsl154',
      'rds',
      paste0('scrna.wnn.int.', paste(this_t_type, this_subset_fxn, sep='.'), '.clean.h5Seurat'))
    
    THIS_T_TYPE <- this_t_type
    
    TITLE <- paste(this_t_type, this_subset_fxn)
    
    if (!(file.exists(SEURAT_HDF5_OUTFILE))) {
      jobRunScript(
        here::here('r','scrna','07_run_subset_job_indiv.R'),
        name = NULL,
        encoding = "unknown",
        workingDir = NULL,
        importEnv = TRUE,
        exportEnv = ""
      )
    }
  }
}