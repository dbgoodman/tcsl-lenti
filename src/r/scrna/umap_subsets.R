library(dplyr)
library(tidyr)
library(tidyr)
library(Seurat)
library(patchwork)
library(cluster)
library(parallelDist)
library(ggplot2)
library(data.table)
library(cowplot)
library(doParallel)

# BiocManager::install("clusterProfiler", version = "3.8")
# BiocManager::install("pathview")
# BiocManager::install("enrichplot")
# BiocManager::install(organism, character.only = TRUE)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(ReactomePA)

reload_data=F
source(here::here('r','fig_colors.R'))
source(here::here('r','scrna','scrna_analysis_fxns.R'))


if (!('scrna.hashed' %in% ls()))
  scrna.hashed <- readRDS(here::here(
    '..','..','scrnaseq_tcsl154',
    'rds','scrna.hashed.rds'))


stimulatory_cars <- c('4-1BB','BAFF-R','CD40','TACI','CD28','Zeta')

subset_fxns= list(
  
  # all cars
  function(obj) subset(obj, subset = (car != 'K_only')),
  function(obj) subset(obj, subset = (CD4v8 == 'CD4' & car != 'K_only')),
  function(obj) subset(obj, subset = (CD4v8 == 'CD8' & car != 'K_only')),
  function(obj) subset(obj, subset = (k_type == 'cd19+' & car != 'K_only')),
  
  # stimulated 4s and 8s
  function(obj) subset(obj, subset = (k_type == 'cd19+' & CD4v8 == 'CD4' & car != 'K_only')),
  function(obj) subset(obj, subset = (k_type == 'cd19+' & CD4v8 == 'CD8' & car != 'K_only')),
  
  #unstimulated 4s and 8s
  function(obj) subset(obj, subset = (k_type == 'baseline' & car != 'K_only')),
  function(obj) subset(obj, subset = (k_type == 'baseline' & CD4v8 == 'CD4' & car != 'K_only')),
  function(obj) subset(obj, subset = (k_type == 'baseline' & CD4v8 == 'CD8' & car != 'K_only')),
  
  #stimulatory cars only
  function(obj) subset(obj, subset = (k_type == 'cd19+' & CD4v8 == 'CD4' & car %in% stimulatory_cars)),
  function(obj) subset(obj, subset = (k_type == 'cd19+' & CD4v8 == 'CD8' & car %in% stimulatory_cars)),

  #stimulatory cars, central memory
  function(obj) subset(obj, subset = (k_type == 'cd19+' & CD4v8 == 'CD4' & car %in% stimulatory_cars &
    CD45RO.adt > 0.25 & CD62L.adt > 0.6)),
  function(obj) subset(obj, subset = (k_type == 'cd19+' & CD4v8 == 'CD8' & car %in% stimulatory_cars &
    CD45RO.adt > 0.25 & CD62L.adt > 0.6)),
  #stimulatory cars, effector memory
  function(obj) subset(obj, subset = (k_type == 'cd19+' & CD4v8 == 'CD4' & car %in% stimulatory_cars &
    CD45RO.adt > 0.25 & CD62L.adt < 0.6)),
  function(obj) subset(obj, subset = (k_type == 'cd19+' & CD4v8 == 'CD8' & car %in% stimulatory_cars &
    CD45RO.adt > 0.25 & CD62L.adt < 0.6))
)


subset_names= c('All','CD4','CD8','CD19',
              'CD4.CD19','CD8.CD19',
              'baseline','CD4.baseline', 'CD8.baseline',
              'CD4.CD19.stimulatory','CD8.CD19.stimulatory',
              'CD4.CD19.stimulatory.tcm','CD8.CD19.stimulatory.tcm',
              'CD4.CD19.stimulatory.teff','CD8.CD19.stimulatory.teff')
  
umap_plots <- data.table(subset_name= subset_names, subset_fxn_num=1:length(subset_fxns))
  
  
umap_plots[c(1:11), RNA := list(list(
    umap_subset(subset_fxns[[subset_fxn_num]](scrna.hashed), 
      plot_title_text=paste0(.BY[1],' RNA'), assay='RNA', vision=T, save_plots=T, rerun=F))), 
    by='subset_name']

