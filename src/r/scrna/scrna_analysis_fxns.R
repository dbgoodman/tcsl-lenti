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


# get gene ids and entrez ids for gsea analysis
get_gene_id_list <- function(subset_obj) {

  gene_symbols <- rownames(scrna.hashed@assays$RNA@data)
  gene_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
  
  # get missing genes, try to match using synonym table
  unknown_symbols <- gene_symbols[!(gene_symbols %in% gene_ids$SYMBOL)]
  gene_synonyms <- fread(here::here('..','data','Homo_sapiens.gene_info.gz'))
  gene_to_syn <- gene_synonyms[, list(synonym=unlist(strsplit(as.character(Synonyms),"\\|"))), by='Symbol']
  gene_to_syn <- gene_to_syn[gene_synonyms[, list(GeneID, Symbol)], on='Symbol']
  gene_ids <- rbind(
    gene_ids, 
    gene_to_syn[Symbol %in% unknown_symbols, list(SYMBOL=Symbol, ENTREZID=GeneID)],
    gene_to_syn[synonym %in% unknown_symbols, list(SYMBOL=synonym, ENTREZID=GeneID)])
  names(gene_ids) <- c('gene','entrez')
  
  # list of remaining unknown genes, avoid unnamed A* genes
  unknown_symbols <- grep('^[A[CFLP]\\d+', unknown_symbols[!(unknown_symbols %in% gene_ids$gene)], value=T, invert=T)
  
  return(gene_ids)
}

# run gsea analysis on GO, Reactome, Kegg
pathway_enrichment <- function(marker_set,
    logfc.thresh=0.25, min.pct=0.2, minGSSize=3, maxGSSize=100, nperm=1000, ident='car') {
  
  organism = "org.Hs.eg.db"

  gene_data <- data.table(marker_set, keep.rownames = T)
  
  #symbol gene list
  names(gene_data)[1] <- 'gene'
  gene_data <- gene_ids[gene_data, on='gene']
  gene_data <- gene_data[!duplicated(gene_data$gene)]
  gene_data <- gene_data[!duplicated(gene_data$entrez)]
  gene_list <- gene_data$avg_logFC
  names(gene_list) <- gene_data$gene
  gene_list = sort(gene_list, decreasing = TRUE)
  
  #entrez gene list
  entrez_gene_list <- gene_data$avg_logFC
  # Name vector with ENTREZ ids
  names(entrez_gene_list) <- gene_data$entrez
  entrez_gene_list = sort(entrez_gene_list, decreasing = TRUE)
  
  gse_bp <- gseGO(geneList=gene_list, 
    ont ="BP", 
    keyType = "SYMBOL", 
    nPerm = nperm, 
    minGSSize = minGSSize, 
    maxGSSize = maxGSSize, 
    pvalueCutoff = 0.05, 
    verbose = TRUE, 
    OrgDb = organism, 
    pAdjustMethod = "none", by='DOSE')
  
  gse_mf <- gseGO(geneList=gene_list, 
      ont ="MF", 
      keyType = "SYMBOL", 
      nPerm = nperm, 
      minGSSize = minGSSize, 
      maxGSSize = maxGSSize, 
      pvalueCutoff = 0.05, 
      verbose = TRUE, 
      OrgDb = organism, 
      pAdjustMethod = "none", by='DOSE')
  
  kegg_organism = "hsa"
  gse_kegg <- gseKEGG(geneList= entrez_gene_list,
    organism     = kegg_organism,
    nPerm        = nperm,
    minGSSize = minGSSize, 
    maxGSSize = maxGSSize, 
    pvalueCutoff = 0.05,
    pAdjustMethod = "none",
    keyType       = "ncbi-geneid")
  
  gse_pa <- gsePathway(geneList=entrez_gene_list,
    nPerm        = nperm,
    minGSSize = minGSSize, 
    maxGSSize = maxGSSize, 
    pvalueCutoff=0.05)

  gse_results <- list(
    gse_mf= gse_mf,
    gse_bp= gse_bp,
    gse_kegg= gse_kegg,
    gse_pa= gse_pa)
  
  # make table
  pathway_table <- list()
  for (gse in names(gse_results)) {
    
    this_gse_results <- gse_results[[gse]]@result
    
    if (gse %in% c('gse_kegg','gse_pa') && nrow(this_gse_results) > 0) {
      this_gse_results$core_enrichment <- gene_ids[
        data.table(this_gse_results)[, list(
            entrez=unlist(strsplit(core_enrichment, '/'))), by='ID'], 
                on='entrez'][, list(core_enrichment=paste(gene, collapse='/')), by='ID']$core_enrichment
      gse_results[[gse]]@result$core_enrichment <- this_gse_results$core_enrichment
    }
    
    if (nrow(this_gse_results) > 0) {
      pathway_table <- rbind(pathway_table, 
        data.table('gse'=gse, this_gse_results), fill=T)
    }
  }

  return(list(gse_results=gse_results, pathway_table=pathway_table))
}
