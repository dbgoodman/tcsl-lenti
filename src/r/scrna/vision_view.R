library(VISION)
library(data.table)

SEURAT_RDS_FILE <- here::here(
    '..','..','scrnaseq_tcsl154',
    'rds','scrna.hashed.rds')

VISION_RDS_FILE <- here::here(
    '..','..','scrnaseq_tcsl154',
    'rds','scrna.vision.rds')

STIMULATORY_CARS <- c('4-1BB','BAFF-R','CD40','TACI','CD28','Zeta')

vis <- readRDS(here::here(
    '..','..','scrnaseq_tcsl154',
    'rds','scrna.vision.cd8.stim.act.rds'))

viewResults(vis, port=3838)
