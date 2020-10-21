library(dplyr)
library(tidyr)
library(Seurat)
library(patchwork)
library(cluster)
library(parallelDist)
library(ggplot2)
library(data.table)
library(cowplot)
library(doParallel)
library(VISION)
library(spatstat)
library(dendextend)
library(ggdendro)


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

####### GSEA 

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

####### CITESEQ SUBSETS

make_citeseq_dt <- function(obj_subset, assay='ADT') {
  return(data.table(cbind(t(obj_subset@assays[[assay]]@data), obj_subset@meta.data, list(cell.id=colnames(obj_subset))))[car != 'K_only'])
}

make_rna_dt <- function(obj_subset) {
  return(data.table(cbind(t(obj_subset@assays$RNA@data), obj_subset@meta.data))[car != 'K_only'])
}

make_pct_grid <- function(
    obj_subset, x_axis, x_cutoff, y_axis, y_cutoff, split=c('car','k_type','CD4v8')) 
{
  pct_dt <- obj_subset[,
    quadrant := factor((get(x_axis) > x_cutoff)*2 + (get(y_axis) > y_cutoff), 
      levels=c(0:3), label=c('BL','TL','BR','TR'))][, 
        .N, by=c('quadrant', split)][,
          list(pct= N/sum(N), count= N, quadrant), by=split]
 
  pos_dt <- data.table(
    quadrant=c('BL','TL','BR','TR'),
    x=c(-Inf, -Inf, Inf, Inf),
    y=c(-Inf,Inf,-Inf,Inf))
 
  pct_dt <- pos_dt[pct_dt, on='quadrant']
  
  
  labels <- paste0(
    rep(gsub('(.*)\\.\\w+','\\1',x_axis),4),
    c('-','+','-','+'),
    '/',
    rep(gsub('(.*)\\.\\w+','\\1',y_axis),4),
    c('-','-','+','+'))
  
   pct_dt[, lbl := factor(quadrant, labels=labels)]

  
  return(pct_dt)
}

citeseq_density_2d <- function(obj_subset, x_axis, x_cutoff, y_axis, y_cutoff, 
    title=paste0(x_axis,' vs ',y_axis)) 
{
  ggplot(obj_subset, aes_string(x=x_axis, y=y_axis)) + 
      stat_density_2d(aes(fill = stat(ndensity)), geom = "raster", contour = FALSE) + 
      scale_fill_distiller(palette='PuRd', direction=1) + facet_grid(k_type~car) +  
      geom_hline(yintercept=y_cutoff) + geom_vline(xintercept=x_cutoff) + 
      scale_y_continuous(expand = c(0, 0)) + 
      scale_x_continuous(expand=c(0,0)) + 
      labs(title=title) +
      geom_text(
        data=make_pct_grid(obj_subset, x_axis, x_cutoff, y_axis, y_cutoff), 
        size=3, 
        aes(x=x, y=y, hjust=(x > 0), vjust=(y > 0), label=paste0(as.character(round(pct,3)*100),'%\n','(',count,')'))) +
      theme_minimal()
}

citeseq_bar_plot <- function(obj_subset, x_axis, x_cutoff, y_axis, y_cutoff, 
    title=paste0(x_axis,' vs ',y_axis))
{
  ggplot(make_pct_grid(obj_subset, x_axis, x_cutoff, y_axis, y_cutoff)) + 
      scale_fill_brewer(palette='Set1') +
      geom_bar(aes(fill=lbl, x=car, y=pct), stat='identity', position='dodge') +
      facet_wrap(~k_type, ncol=1) + theme_bw()
}

density_bar_plots <- function (obj_subset, x_axis, x_cutoff, y_axis, y_cutoff, 
    title=paste0(x_axis,' vs ',y_axis))
{
  title <- ggdraw() + draw_label(title, fontface='bold')
  plot_grid(
    title,
    citeseq_density_2d(obj_subset, x_axis, x_cutoff, y_axis, y_cutoff, title=''),
    citeseq_bar_plot(obj_subset, x_axis, x_cutoff, y_axis, y_cutoff, title=''),
    ncol=1, 
    rel_heights=c(0.2,1,1)
  )
}

####### PCA / VOLCANO

plot_pca <- function(obj_subset, axes=c('car','k_type'), min_pct_exp = 10, assay='RNA', color_cars=F) {
  
  if ('cluster_dropped' %in% names(obj_subset@meta.data)) {
    obj_subset <- subset(obj_subset, subset=(cluster_dropped == F))
  }
  
  obj_subset@meta.data[['pca_group']] <- do.call(paste, c(obj_subset@meta.data[axes], sep="_"))
  Idents(obj_subset) <- 'pca_group'
  
  if ('RNA' %in% assay) {
    rna_gene_set <- rownames(obj_subset@assays[[assay]]@data)[
      (tabulate(obj_subset@assays$RNA@data@i + 1) / 
        ncol(obj_subset@assays$RNA@data) * 100) > min_pct_exp]
  } else {
    rna_gene_set <- c()
  }
  
  if ('ADT' %in% assay) {
    adt_gene_set <- rownames(obj_subset@assays[[assay]]@data)[
      (tabulate(obj_subset@assays$RNA@data@i + 1) / 
        ncol(obj_subset@assays$ADT@data) * 100) > min_pct_exp]
  } else {
    adt_gene_set <- c()
  }
  
  
  obj.avg <- AverageExpression(obj_subset, assays=assay)
  obj.avg <- t(rbind(obj.avg$ADT, obj.avg$RNA))
  obj.avg <- obj.avg[, apply(obj.avg,2,sd) > 0]
  obj.avg <- obj.avg[, colnames(obj.avg) %in% c(rna_gene_set, adt_gene_set)]
  combined_pca <- prcomp(scale(obj.avg))
  
  # calculate pca stats
  pca.dt <- data.table(pc = data.table(colnames(combined_pca$rotation))[, 
      PC := as.integer(gsub("[A-Z]", "", V1))][, PC], 
      sd = combined_pca$sdev, 
      var = combined_pca$sdev^2, 
      var.norm = combined_pca$sdev^2/sum(combined_pca$sdev^2), 
      var.acc = cumsum(combined_pca$sdev^2/sum(combined_pca$sdev^2)))
  
  melt.pca.dt <- melt(
    pca.dt, measure.vars = c("var.norm", "var.acc"), variable.name = "metric")
  
  pca_metrics_plot <- ggplot(melt.pca.dt) +
    geom_line(aes(x = pc, y = value, color = metric)) +
    geom_point(aes(x = pc, y = value, color = metric)) +
    scale_x_continuous(limits = c(1, NA)) +
    labs(title = "Fraction of Variance Captured by Principal Components, Every Sort Group", 
         x = "Principal Component", y = "Fraction of Variance") +
    theme_bw()
  
  projected.metrics <- data.table(
    scale(obj.avg, combined_pca$center, combined_pca$scale) %*% combined_pca$rotation)
  
  projected.metrics <- projected.metrics[,
      car := rownames(obj.avg)]
  
  if (color_cars) {
    pca_plot <- ggplot(
        projected.metrics, aes(x=PC1, y=PC2, label=car, color='k_type', shape='CD4vCD8')) +
      geom_point() + 
      geom_text_repel()
  } else {
    pca_plot <- ggplot(
        projected.metrics, aes(x=PC1, y=PC2, label=car)) +
      geom_point() + 
      geom_text_repel()    
  }
  
  
  return(
    list(pca_plot= pca_plot, combined_pca= combined_pca, pca_metrics_plot= pca_metrics_plot)
  )
}

plot_gene_assoc_bar <- function(pca_data, PC_name, n_genes=10) {
  
  gene_data <- data.table(pca_data$combined_pca$rotation, keep.rownames=TRUE)[
    c(order(get(PC_name))[1:n_genes], order(-get(PC_name))[1:n_genes]), 
    list(association= get(PC_name), gene= rn)][, 
      gene := factor(gene, levels=gene[order(association)])]

  ggplot(gene_data) + 
    geom_bar(aes(x=association, y=gene, fill=factor(sign(association))), 
      stat='identity', color='black') + 
    theme_minimal() + 
    theme(legend.position='none') +
    labs(x='',y='', title=PC_name)
}

plot_pca_cars <- function(obj_subset, title, axes=c('car','k_type'), min_pct_exp=10, assay='RNA', color_cars=F) {
  
  if ('cluster_dropped' %in% names(obj_subset@meta.data)) {
    obj_subset <- subset(obj_subset, subset=(cluster_dropped == F))
  }
  
  pca_data <- plot_pca(obj_subset, axes, min_pct_exp, assay, color_cars)
  return(plot_grid(
    pca_data$pca_plot + theme_minimal() + labs(title=title),
    plot_gene_assoc_bar(pca_data, 'PC1', n_genes=25),
    plot_gene_assoc_bar(pca_data, 'PC2', n_genes=25),
    ncol=3, nrow=1,
    rel_widths=c(5, 1, 1)))
}

comp_car_volcano <- function(obj_subset, ident.1, ident.2, logfc.thresh=0.2, min.pct=0.1,
    assay='RNA', lbl_pval_cutoff = 10, title='', do_pathway_enrichments=F) {
  
  title <- paste(title,
    paste('[',paste(ident.2,collapse=', '),'] <- vs -> [',paste(ident.1,collapse=', '),']'),
    sep='\n')
  
  Idents(obj_subset) <- 'car'
  markers <- FindMarkers(
    obj_subset, assay=assay,
    ident.1 = ident.1, ident.2=ident.2,
    logfc.threshold = logfc.thresh, min.pct=min.pct)
  
  if (do_pathway_enrichments) {
    pathways <- pathway_enrichment(markers)
    gse_results <- pathways[[1]]
    pathway_table <- pathways[[2]] }
  else {
    pathways <- NULL
    gse_results <- NULL
    pathway_table <- NULL
  }
  
  volcano_ggplot <- ggplot(markers, aes(x=avg_logFC, y=-log10(p_val_adj))) + geom_point() +
    geom_hline(aes(yintercept=-log10(0.05))) + 
    scale_x_continuous(limits=c(-max(abs(markers$avg_logFC)), max(abs(markers$avg_logFC)))) +
    geom_text_repel(
      data=data.table(markers, keep.rownames = T)[-log10(p_val_adj) > lbl_pval_cutoff],
      aes(label=rn)) + 
    theme_minimal() +
    labs(x=title, y='Adjusted -log10(p)')
  
  return(list(markers=markers, plot=volcano_ggplot, gse_results=gse_results, pathway_table=pathway_table))
}

####### UMAP PLOTTING

umap_plot <- function(this_subset, plot_title_text="", assay, cluster_resolution=NULL, cluster_name=NULL, 
    reduction='umap', feat_logfc_thresh=0.25, feat_min_pct=0.5) {
  
  stopifnot(!is.null(cluster_name) || !is.null(cluster_name))
  
  if (!is.null(cluster_resolution))
    cluster_col_name <- paste0(assay,'_snn_res.',as.character(cluster_resolution))
  if (!is.null(cluster_name))
    cluster_col_name <- cluster_name

  if ('cluster_dropped' %in% names(this_subset@meta.data)) {
    this_subset <- subset(this_subset, subset=(cluster_dropped == F))
  }
  
  subset.markers <- FindAllMarkers(this_subset,
  features = VariableFeatures(this_subset),
  min.pct = feat_min_pct, logfc.threshold = feat_logfc_thresh)
  
  top.subset.markers <- subset.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  
  print(names(this_subset@meta.data))
  
  top_genes_plot <- ggplot(top.subset.markers) + 
    geom_bar(aes(y=avg_logFC, x=reorder(gene, avg_logFC)), stat='identity') + 
    facet_wrap(~cluster, scales='free', ncol=3) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  cluster_pct_plot <- ggplot(
      data.table(this_subset@meta.data)[
        CD4v8 != 'intermediate', .N, 
        by=c('car', cluster_col_name, 'CD4v8', 'k_type')][,
          list(cluster=get(cluster_col_name), t_type=CD4v8, frac=N/sum(N)), 
          by=c('car', 'CD4v8', 'k_type')]) + 
    geom_bar(aes_string(x='cluster', y='frac', fill='frac'), stat='identity', color='grey50') +
    scale_fill_distiller(palette='YlGnBu', direction=1) +
    facet_grid(k_type+t_type~car) + 
    theme_bw() +
    coord_flip()
  
  plot_title <- ggdraw() + 
    draw_label(
      plot_title_text,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  
  gene_heatmap <- DoHeatmap(this_subset,
      features = unique(
        (top.subset.markers %>% group_by(cluster) 
          %>% top_n(n=5, wt= avg_log2FC))$gene))
      
  umap_plots <- plot_grid(
    plot_grid(
      plot_title,
      DimPlot(this_subset, reduction = reduction, label=T),
      cluster_pct_plot,
      rel_heights=c(0.2,1.25,1),
      ncol=1),
    gene_heatmap,
    ncol=2)
  
  return(list(umap_plots=umap_plots, marker_genes=subset.markers))
}

####### Run UMAP
calculate_umap <- function(this_subset, 
    plot_title_text="", assay='RNA', cluster_resolution=0.8,
    save=T, save_plots=T, rerun=F, return_all=F, parallel=F) {
  
  this_subset_umap <- subset(this_subset, 
    features=rownames(this_subset@assays[[assay]]))
    
  if (assay == 'RNA') {
    
      this_subset_umap <- NormalizeData(
          this_subset_umap,
          normalization.method = "LogNormalize", 
          scale.factor = 10000)
    } else if (assay == 'ADT') {
      this_subset_umap <- NormalizeData(
        this_subset_umap,
        normalization.method = "CLR")
    }
    
    this_subset_umap <- FindVariableFeatures(
      this_subset_umap, 
      selection.method = "vst",
      nfeatures = 2000)
    
    this_subset_umap <- ScaleData(
      this_subset_umap, assay=assay,
      features= VariableFeatures(this_subset_umap))
  
    this_subset_umap <- RunPCA(
      this_subset_umap, features = VariableFeatures(this_subset_umap))

    this_subset_umap <- FindNeighbors(this_subset_umap, reduction = "pca", dims = 1:8)
    this_subset_umap <- FindClusters(this_subset_umap, resolution = cluster_resolution, algorithm=1)
    this_subset_umap <- RunUMAP(this_subset_umap, dims = 1:8)
  
    this_subset_umap@meta.data[['car.ttype.ktype']] <- paste(
      this_subset_umap@meta.data[['car']],
      this_subset_umap@meta.data[['CD4v8']],
      this_subset_umap@meta.data[['k_type']], sep='_')
    
    return(this_subset_umap)
}


###### UMAP/CORR PLOTS

cluster_car_pct_plot <- function(this_subset, cluster_col_name) {
  return(ggplot(
      data.table(this_subset@meta.data)[
        CD4v8 != 'intermediate', .N, 
        by=c('car', cluster_col_name, 'CD4v8', 'k_type')][,
          list(cluster=get(cluster_col_name), t_type=CD4v8, frac=N/sum(N)), 
          by=c('car', 'CD4v8', 'k_type')]) + 
    geom_bar(aes_string(x='cluster', y='frac', fill='frac'), stat='identity', color='grey50') +
    scale_fill_distiller(palette='YlGnBu', direction=1) +
    facet_grid(k_type+t_type~car) + 
    theme_bw() +
    coord_flip())
}

umap_plot_by_car <- function(this_subset, plot_title_text="", assay, cluster_resolution, cluster_name=NULL) {

  
  if ('cluster_dropped' %in% names(this_subset@meta.data)) {
    this_subset <- subset(this_subset, subset=(cluster_dropped == F))
  }
  Idents(this_subset) <- 'car'
  subset.markers <- FindAllMarkers(this_subset,
  features = VariableFeatures(this_subset),
  min.pct = 0.5, logfc.threshold = 0.25)
  
  top.subset.markers <- subset.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  
  top_genes_plot <- ggplot(top.subset.markers) + 
    geom_bar(aes(y=avg_logFC, x=reorder(gene, avg_logFC)), stat='identity') + 
    facet_wrap(~cluster, scales='free', ncol=3) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  cluster_col_name <- paste0(assay,'_snn_res.',as.character(cluster_resolution))
  
  cluster_pct_plot <- ggplot(
      data.table(this_subset@meta.data)[
        CD4v8 != 'intermediate', .N, 
        by=c('car', cluster_col_name, 'CD4v8', 'k_type')][,
          list(cluster=get(cluster_col_name), t_type=CD4v8, frac=N/sum(N)), 
          by=c('car', 'CD4v8', 'k_type')]) + 
    geom_bar(aes_string(x='cluster', y='frac', fill='frac'), stat='identity', color='grey50') +
    scale_fill_distiller(palette='YlGnBu', direction=1) +
    facet_grid(k_type+t_type~car) + 
    theme_bw() +
    coord_flip()
  
  plot_title <- ggdraw() + 
    draw_label(
      plot_title_text,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  
  gene_heatmap <- DoHeatmap(this_subset,
        features = unique(
          (top.subset.markers %>% group_by(cluster) 
            %>% top_n(n=5, wt= avg_logFC))$gene))
  
  car_density_map <- ggplot(
      cbind(this_subset@meta.data, Embeddings(object = this_subset[["umap"]]))) + 
    geom_density2d(aes(x=UMAP_1, y=UMAP_2)) + facet_wrap(~car) +
    theme_minimal()
      
  umap_plots <- plot_grid(
    plot_grid(
      plot_title,
      car_density_map,
      cluster_pct_plot,
      rel_heights=c(0.2,1.25,1),
      ncol=1),
        plot_grid(
      gene_heatmap,
      top_genes_plot,
      ncol=1,
      rel_heights=c(3,2)),
    ncol=2)
  
  return(umap_plots)
}

corr_plots <- function(subset_obj, assay='ADT') {
  
  DefaultAssay(subset_obj) <- assay
  
  citeseq_cor <- make_citeseq_dt(subset_obj, assay=assay)[, 
      data.table(cor(.SD[, c(rownames(subset_obj@assays[[assay]]@data)), with=F]), keep.rownames=T)]
  
  #remove NA rows
  citeseq_cor <- na.omit(
    citeseq_cor[, which(colSums(is.na(citeseq_cor)) != (nrow(citeseq_cor) - 1)), with=F])
   
  # generate heirarchical clustering
  dendro_m <- as.matrix(citeseq_cor[, -1])
  rownames(dendro_m) <- unlist(citeseq_cor[,1])
  dendro <- as.dendrogram(hclust(d = dist(x = dendro_m)))
  dendro <- rotate(dendro, names(sort(rowMeans(dendro_m))))
  
  # row order from clustering
  citeseq_cor_melt <- melt(citeseq_cor, id.vars='rn')
  
  # melt and reorder rows
  citeseq_cor_melt$rn <- factor(
      citeseq_cor_melt$rn,
      levels = labels(dendro))
  citeseq_cor_melt$variable <- factor(
      citeseq_cor_melt$variable,
      levels = labels(dendro))
  
  # remove markers that are not strongly correlated with anything
  var_has_corr <- citeseq_cor_melt[variable != rn, max(abs(value)) > 0.075, by='variable'][V1 == T]$variable
  citeseq_cor_melt_subset <- citeseq_cor_melt[variable %in% var_has_corr & rn %in% var_has_corr]
  
  average_corr <- ggplot(copy(citeseq_cor_melt_subset)[value > 0.5, value := 0.5]) + 
    geom_tile(aes(x=rn, y=variable, fill=value)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    scale_fill_distiller(palette='PRGn', limits=c(-0.5, 0.5))
  
  return(average_corr)
}

####### RUN UMAP ON SCRNA SUBSET

umap_subset <- function(this_subset, 
    plot_title_text="", assay='RNA', cluster_resolution=0.8,
    save=T, save_plots=T, rerun=F, return_all=F, parallel=F, vision=F) {
  
  DefaultAssay(this_subset) <- assay
  
  subset_name <- gsub(' ','.',plot_title_text)
  
  print(subset_name)
  
  rds_file <- here::here(
    '..','..','scrnaseq_tcsl154',
    'rds',
    paste0(
        'scrna.hashed.',
        subset_name,
        '.rds'))

  vis_rds_file <- here::here(
    '..','..','scrnaseq_tcsl154',
    'rds',
    paste0(
        'vision.',
        subset_name,
        '.rds'))
  
  plot_dir <- here::here(
    '..','figs','scrna',
    'subset_plots',
    subset_name)
  
  if (save_plots) dir.create(plot_dir, recursive=T)
  
  ##### RUN UMAP
  
  if (rerun || !file.exists(rds_file)) {
    
    this_subset_umap <- calculate_umap(this_subset, 
      plot_title_text="", assay='RNA', cluster_resolution=0.8,
      save=T, save_plots=save_plots, rerun=F, return_all=F, parallel=F)
  } else {
    this_subset_umap <- readRDS(rds_file)
  }
  
  cluster_col_name <- paste0(assay,'_snn_res.',as.character(cluster_resolution))
  print(names(this_subset_umap@meta.data))
  
  ##### CLUSTER PRUNING
  
  # prune UMAP clusters:
  cluster_stats <- data.table(
    this_subset_umap@meta.data)[, list(get(cluster_col_name), n_clust_car=.N), 
      by=c(cluster_col_name, 'car', 'CD4v8', 'k_type')][, 
        list(cluster=get(cluster_col_name), n_clust_car, total=sum(n_clust_car), pct=n_clust_car/sum(n_clust_car)), 
          by=c( 'car', 'CD4v8', 'k_type')]
  
  #criteria to drop clusters, all must be true to keep:
  # at least 100 cells
  cluster_stats[, min_cells := sum(n_clust_car) > 100, by='cluster']
  
  # at least 10% of the populaton for one CAR OR at least 5% of total population
  cluster_stats[, max_car_pct := max(pct) > 0.05, by='cluster']
  cluster_stats[, max_tot_pct := sum(n_clust_car)/cluster_stats[, sum(n_clust_car)] > .1, by='cluster']
  
  clusters_kept <- as.numeric(cluster_stats[, 
    all(min_cells, max_car_pct & (max_car_pct | max_tot_pct)), by='cluster'][V1 == T, cluster])
  
  this_subset_umap@meta.data$cluster_dropped <- !(this_subset_umap@meta.data[[cluster_col_name]] %in% clusters_kept)
  
  ##### SAVE UMAP RDS
  
  if ((rerun | (!file.exists(rds_file))) & save)
    saveRDS(this_subset_umap, file = rds_file)
  
  ##### DO VISION ANALYSIS
  if (vision & (rerun || !file.exists(vis_rds_file))) {
    load_vision_data(
      seurat_rds_file=NULL,
      subset_fxn=NULL,
      seurat_obj=this_subset_umap,
      vision_rds_file=vis_rds_file, 
      cluster_name=cluster_col_name,
      cores=8)
  }
  
  ##### MAKE PLOTS
  
  fn_pref = paste0(plot_dir,'/',subset_name,'.')

  # plots
  umap_plot_output <- umap_plot(this_subset_umap, plot_title_text, assay, cluster_resolution)
  subset_umap_plot <- umap_plot_output$umap_plots
  car_umap_plot <- umap_plot_by_car(this_subset_umap, plot_title_text, assay, cluster_resolution)
  car_pca_plot <- plot_pca_cars(this_subset, title=plot_title_text, axes=c('car','k_type','CD4v8'), color_cars=T)
  umap_pca_plot <- plot_pca_cars(this_subset_umap, title=plot_title_text, axes=c(cluster_col_name))


  if (save_plots) {
    ggsave(paste0(fn_pref,'subset_umap.pdf'), subset_umap_plot, width=15, height=9)
    ggsave(paste0(fn_pref,'car_umap_plot.pdf'), car_umap_plot, width=15, height=9)
    ggsave(paste0(fn_pref,'umap_pca_plot.pdf'), umap_pca_plot, width=15, height=9)
    ggsave(paste0(fn_pref,'car_pca_plot.pdf'), car_pca_plot, width=10, height=10)
  }
  
  if (save) {
    cluster_marker_genes <- umap_plot_output$marker_genes
    fwrite(cluster_marker_genes, paste0(fn_pref, 'cluster_genes.csv'))
  }

  # volcano plots
  if (assay == 'RNA') {
    
    Idents(this_subset) <- 'car'
    
    a_sets = list(c('BAFF-R'), c('TACI'), c('CD40'), c('BAFF-R','TACI'), c('BAFF-R','TACI','CD40'),
      c('BAFF-R','TACI','CD40','4-1BB'), c('BAFF-R','TACI','CD40','4-1BB','CD28'), c('CD28'))
    
    b_sets = list(c('CD28'), c('CD28','4-1BB','Zeta'), c('Zeta'), c('4-1BB'))
    
    volcano_plots <- list()
    pathway_enrichments <- data.table()
    marker_table <- data.table()
    
    for (ident_a in a_sets) {
      ident_a_name <- paste(ident_a, collapse='_')
      volcano_plots[[ident_a_name]] = list()
      
      for (ident_b in b_sets) {
        
        #skip if sets overlap
        if (length(intersect(ident_a, ident_b)) > 0) next

        ident_b_name <- paste(ident_b, collapse='_' )
        
        print(paste("Running CAR comparison:", subset_name," - ", ident_a_name,' VS ',ident_b_name))
        
        volcano_plots[[ident_a_name]][[ident_b_name]] <- comp_car_volcano(
          this_subset, ident.1=ident_a, ident.2=ident_b, lbl_pval_cutoff = 3, min.pct=0.25, assay=assay,
          do_pathway_enrichments=F)
        
        # pathway_table <- data.table(
        #   ident_a = ident_a_name,
        #   ident_b = ident_b_name,
        #   volcano_plots[[ident_a_name]][[ident_b_name]]$pathway_table)

        marker_genes <- data.table(
          ident_a = ident_a_name,
          ident_b = ident_b_name,
          volcano_plots[[ident_a_name]][[ident_b_name]]$markers, use.rownames=T)

        # pathway_enrichments <- rbind(pathway_enrichments, pathway_table, fill=T)
        
        marker_table <- rbind(marker_table, marker_genes)
        
        if (save_plots) {
          fn_prefix = paste0(plot_dir,'/',subset_name,'.',
            ident_a_name,'_VS_',ident_b_name)
          ggsave(paste0(fn_prefix, '.volcano.pdf'), 
            volcano_plots[[ident_a_name]][[ident_b_name]]$plot, width=8, height=8, units='in')
          # fwrite(pathway_table, paste0(fn_prefix, '.pathways.csv'))
          fwrite(marker_genes, paste0(fn_prefix, '.genes.csv'))
        }
      }
    }
  }
    
  if (return_all) {
    return(list(
      cluster_plot=subset_umap_plot,
      car_plot=car_umap_plot,
      volcano_plots=volcano_plots,
      pca_plot=car_pca_plot,
      rds_file=rds_file))
  } else {
    return(list(
        cluster_plot=subset_umap_plot,
        car_plot=car_umap_plot,
        volcano_plots=volcano_plots,
        pca_plot=car_pca_plot,
        rds_file=rds_file))
  }
}

####### VISION

# Load VISION
load_vision_data <- function(
    seurat_rds_file=SEURAT_RDS_FILE,
    subset_fxn=function (obj) return(obj),
    seurat_obj=NULL,
    vision_rds_file=VISION_RDS_FILE, 
    cluster_name=NULL,
    cores=8) {
  
  if (is.null(seurat_obj)) {
    scrna.hashed <- subset_fxn(readRDS(seurat_rds_file))
  } else {
    scrna.hashed <- seurat_obj
  }
  
  # Read in expression counts (Genes X Cells, gene names as rownames)
  counts <- data.frame(scrna.hashed@assays[["RNA"]]@counts)
  
  # Read in expression counts (Genes X Cells, gene names as rownames)
  n.umi <- colSums(counts)
  scaled_counts <- t(t(counts) / n.umi) * median(n.umi)
  
  rm(n.umi)
  rm(counts)
  gc()
  
  ### Adjust metadata to get signature associations
  
  meta <- scrna.hashed@meta.data[gsub('-','.',rownames(scrna.hashed@meta.data)) %in% colnames(scaled_counts),]
  rownames(meta) <- gsub('-','.',rownames(meta))
  
  # extract useful metadata columns
  meta_columns = c('k_type','car','CD4v8')
  if (!(is.null(cluster_name))) {
    meta_columns <- c(meta_columns, cluster_name)
  }
  
  meta <- meta[, meta_columns]
  # make a combined column with all combinations of the two
  meta$sample <- factor(paste(meta$k_type, meta$CD4v8, meta$car, sep='_'))
  meta$car_k <- factor(paste(meta$k_type, meta$car, sep='_'))
  meta$t_k <- factor(paste(meta$k_type, meta$k_type, sep='_'))
  meta$tnfr_new <- meta$car
  meta$tnfr <- meta$car
  meta$tnfr_new[meta$car %in% c('CD40','TACI','BAFF-R')] <- 'TNFR_new'
  meta$tnfr[meta$car %in% c('CD40','TACI','BAFF-R','4-1BB')] <- 'TNFR'
  
  # Make vision object
  vis <- Vision(scaled_counts,
    signatures = here::here('..','data',paste0(c('h','c2'),'.all.v7.1.symbols.gmt')),
    meta = meta)
  
  gc()
  options(mc.cores = cores)
  vis <- analyze(vis)
  saveRDS(vis, vision_rds_file)
  return(vis)
}
