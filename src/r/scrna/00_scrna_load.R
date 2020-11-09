source(here::here('r','scrna','HTOdemux.R'))


#####################

PERCENT.MT.CUTOFF <- 13

SEURAT_RDS_FILE <- here::here(
    '..','..','scrnaseq_tcsl154',
    'rds','scrna.hashed.rds')

SEURAT_HDF5_FILE <- here::here(
  '..','..','scrnaseq_tcsl154',
  'rds','scrna.raw.h5Seurat')

## Load the dataset 
######################
# 
# scrna.data <- Read10X(data.dir = here::here(
#     '..','..','scrnaseq_tcsl154',
#     'cellranger','outs','raw_feature_bc_matrix'))
# 
# #rename some adt features
# hto_rows <- paste('Hashtag',1:12,sep='_')
# 
# adt_feature_names <- fread(here::here(
#     '..','..','scrnaseq_tcsl154',
#     'meta','adt_names.tsv'))
# rownames(scrna.data$`Antibody Capture`) <- c(paste0(adt_feature_names$vis_name,'.adt'), hto_rows)
# 
# #rna features
# scrna <- CreateSeuratObject(
#     counts = scrna.data$`Gene Expression`, 
#     project = "tcsl154", 
#     min.cells = 3,
#     min.features = 200)
# 
# adt_rows <- setdiff(rownames(scrna.data$`Antibody Capture`), hto_rows)
# 
# #adt features
# scrna[['ADT']] <- CreateAssayObject(
#     counts = scrna.data$`Antibody Capture`[adt_rows, colnames(scrna)])
# 
# #hto features
# scrna[['HTO']] <- CreateAssayObject(
#   counts = scrna.data$`Antibody Capture`[hto_rows, colnames(scrna)])
# 
# scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^MT-")
# scrna <- subset(
#   scrna, 
#   subset = nFeature_RNA > 800 & nFeature_RNA < 7000 & percent.mt < PERCENT.MT.CUTOFF)
# 
# SaveH5Seurat(scrna, filename = SEURAT_HDF5_FILE, overwrite=T)
# 
# rm(scrna.data)
# gc()

## HTO 
######################


janky_doublet_removal <- function(scrna) {

  DefaultAssay(scrna) <-  'HTO'
  scrna <- SCTransform(scrna, assay='HTO', verbose = FALSE, new.assay.name='SCT_HTO')
  
  DefaultAssay(scrna) <-  'SCT_HTO'
  scrna <- RunUMAP(scrna, assay='SCT_HTO', features = rownames(scrna[["SCT_HTO"]]), reduction.name='umap.hto')
  scrna <- assign_HTO_sample(scrna)
  
  HTO_sample_metadata <- fread(here::here(
      '..','..','scrnaseq_tcsl154','meta','sample_metadata.csv'))
  names(HTO_sample_metadata)[2] <- 'HTO_call'
  scrna@meta.data$HTO_keep <- data.table(scrna@meta.data)[, !HTO_is_doublet & HTO_call %in% HTO_sample_metadata$HTO]
  
  metadata.rows <- data.frame(
    HTO_sample_metadata[data.table(scrna@meta.data), on='HTO_call', nomatch=NA],
    row.names= rownames(scrna@meta.data))
  scrna <- AddMetaData(scrna, metadata=metadata.rows)
  
  scrna@meta.data$sample <- paste(
    scrna@meta.data$car,
    ifelse(scrna@meta.data$k_type=='baseline', 'unstim','stim'),
    sep='.')

}