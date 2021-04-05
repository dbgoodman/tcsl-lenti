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



AddModuleScoreParallel <- function (object, features, pool = NULL, nbin = 24, ctrl = 100, 
    k = FALSE, assay = NULL, name = "Cluster", seed = 1, search = FALSE, 
    ...) 
{
    if (!is.null(x = seed)) {
        set.seed(seed = seed)
    }
    assay.old <- DefaultAssay(object = object)
    assay <- assay %||% assay.old
    DefaultAssay(object = object) <- assay
    assay.data <- GetAssayData(object = object)
    features.old <- features
    if (k) {
        .NotYetUsed(arg = "k")
        features <- list()
        for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
            features[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == 
                i))
        }
        cluster.length <- length(x = features)
    }
    else {
        if (is.null(x = features)) {
            stop("Missing input feature list")
        }
        features <- lapply(X = features, FUN = function(x) {
            missing.features <- setdiff(x = x, y = rownames(x = object))
            if (length(x = missing.features) > 0) {
                warning("The following features are not present in the object: ", 
                  paste(missing.features, collapse = ", "), ifelse(test = search, 
                    yes = ", attempting to find updated synonyms", 
                    no = ", not searching for symbol synonyms"), 
                  call. = FALSE, immediate. = TRUE)
                if (search) {
                  tryCatch(expr = {
                    updated.features <- UpdateSymbolList(symbols = missing.features, 
                      ...)
                    names(x = updated.features) <- missing.features
                    for (miss in names(x = updated.features)) {
                      index <- which(x == miss)
                      x[index] <- updated.features[miss]
                    }
                  }, error = function(...) {
                    warning("Could not reach HGNC's gene names database", 
                      call. = FALSE, immediate. = TRUE)
                  })
                  missing.features <- setdiff(x = x, y = rownames(x = object))
                  if (length(x = missing.features) > 0) {
                    warning("The following features are still not present in the object: ", 
                      paste(missing.features, collapse = ", "), 
                      call. = FALSE, immediate. = TRUE)
                  }
                }
            }
            return(intersect(x = x, y = rownames(x = object)))
        })
        cluster.length <- length(x = features)
    }
    if (!all(Seurat:::LengthCheck(values = features))) {
        warning(paste("Could not find enough features in the object from the following feature lists:", 
            paste(names(x = which(x = !Seurat:::LengthCheck(values = features)))), 
            "Attempting to match case..."))
        features <- lapply(X = features.old, FUN = CaseMatch, 
            match = rownames(x = object))
    }
    if (!all(Seurat:::LengthCheck(values = features))) {
        stop(paste("The following feature lists do not have enough features present in the object:", 
            paste(names(x = which(x = !Seurat:::LengthCheck(values = features)))), 
            "exiting..."))
    }
    pool <- pool %||% rownames(x = object)
    data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
    data.avg <- data.avg[order(data.avg)]
    data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
        n = nbin, labels = FALSE, right = FALSE)
    names(x = data.cut) <- names(x = data.avg)
    ctrl.use <- vector(mode = "list", length = cluster.length)
    for (i in 1:cluster.length) {
        features.use <- features[[i]]
        for (j in 1:length(x = features.use)) {
            ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(x = data.cut == 
                data.cut[features.use[j]])], size = ctrl, replace = FALSE)))
        }
    }
    ctrl.use <- lapply(X = ctrl.use, FUN = unique)
    ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(x = ctrl.use), 
        ncol = ncol(x = object))
    for (i in 1:length(ctrl.use)) {
        features.use <- ctrl.use[[i]]
        ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use, 
            ])
    }
    features.scores <- matrix(data = numeric(length = 1L), nrow = cluster.length, 
        ncol = ncol(x = object))
    for (i in 1:cluster.length) {
        features.use <- features[[i]]
        data.use <- assay.data[features.use, , drop = FALSE]
        features.scores[i, ] <- Matrix::colMeans(x = data.use)
    }
    features.scores.use <- features.scores - ctrl.scores
    rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
    features.scores.use <- as.data.frame(x = t(x = features.scores.use))
    rownames(x = features.scores.use) <- colnames(x = object)
    object[[colnames(x = features.scores.use)]] <- features.scores.use
    DefaultAssay(object = object) <- assay.old
    return(features.scores.use)
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

get_cluster_pct_plot <- function(this_subset, cluster_col_name) {
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
  
  return(cluster_pct_plot)
}

get_cluster_enrich_plot <- function(
  this_subset, cluster_col_name, max_log2_enrich=2, hide=c('KLRG1','Untransduced'),
  bar_cutoff=0.03, reorder=T, na.rm=T, facet_formula=formula(k_type+t_type~car),
  bar_position='stack',
  by_vars=c('car', 'CD4v8', 'k_type'),
  color_by_donor=F) {
  
  rhs <- formula.tools::rhs.vars(facet_formula)
  
  if (color_by_donor)
    by_vars <- c(by_vars, donor)
    
  enrich_data <- data.table(this_subset@meta.data)[
    CD4v8 != 'intermediate' & !(car %in% hide) , .N, 
    by=c(by_vars, cluster_col_name)][,
      list(cluster=get(cluster_col_name), t_type=CD4v8, frac=N/sum(N)), 
      by=by_vars][,
        rel_frac := log2(frac/mean(frac)),
        by=c('cluster', by_vars[!(by_vars %in% rhs)])]
  
  if (na.rm==T) {
    enrich_data <- enrich_data[!is.na(cluster)]
  }
  
  enrich_data[, mean_frac := mean(frac), by=c('cluster', by_vars[!(by_vars %in% rhs)])]
  
  enrich_data[, rel_frac_disp := ifelse(abs(rel_frac) > max_log2_enrich, 
    sign(rel_frac)*max_log2_enrich,
    rel_frac)]
  
  if (reorder)
    enrich_data[, cluster := reorder(cluster, mean_frac)]
  
  plot <- ggplot(enrich_data[mean_frac > bar_cutoff]) + 
    geom_bar(aes(x=cluster, y=frac, fill=rel_frac_disp), stat='identity', color='grey50', position=bar_position) +
    scale_fill_distiller('Log2 Enrichment\nvs mean CAR', palette='PuOr', direction=1,
      limits=c(-max_log2_enrich, max_log2_enrich)) +
    facet_grid(facet_formula, space='free', scales='free_y') + labs(x='Log2 enrichment vs mean car') +
    theme_bw() +
    coord_flip()
  
  return(list(plot=plot, data=enrich_data))
}
  
umap_plot <- function(this_subset, plot_title_text="", assay, cluster_resolution=NULL, cluster_name=NULL, 
    reduction='umap', feat_logfc_thresh=0.25, feat_min_pct=0.5, subset.markers=NULL) {
  
  DefaultAssay(this_subset) <- assay
  stopifnot(!is.null(cluster_name) || !is.null(cluster_name))
  
  if (!is.null(cluster_resolution))
    cluster_col_name <- paste0(assay,'_snn_res.',as.character(cluster_resolution))
  if (!is.null(cluster_name))
    cluster_col_name <- cluster_name

  if ('cluster_dropped' %in% names(this_subset@meta.data)) {
    this_subset <- subset(this_subset, subset=(cluster_dropped == F))
  }
  
  if (is.null(subset.markers)) {
    subset.markers <- FindAllMarkers(this_subset,
    features = VariableFeatures(this_subset),
    min.pct = feat_min_pct, logfc.threshold = feat_logfc_thresh)
  }
  
  top.subset.markers <- subset.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  
  print(names(this_subset@meta.data))
  
  top_genes_plot <- ggplot(top.subset.markers) + 
    geom_bar(aes(y=avg_logFC, x=reorder(gene, avg_logFC)), stat='identity') + 
    facet_wrap(~cluster, scales='free', ncol=3) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  cluster_pct_plot <- get_cluster_enrich_plot(this_subset, cluster_col_name, bar_cutoff=0, hide='')
  
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

plot_sig_by_car <- function(seurat_obj, vis, signature, x_lab, assay='RNA') {

  DefaultAssay(seurat_obj) <- assay
  sig_scores <- getSignatureScores(vis)

  sig.data <- data.table(
    cbind(seurat_obj@meta.data, sig_scores))[, 
      car := factor(
        car,
        levels=c('CD28','4-1BB','BAFF-R','TACI','CD40','TNR8','KLRG1','Zeta','Untransduced'),
        labels=c(car_order,'Untransduced'))]
  
  plot <- ggplot(sig.data) +
    geom_boxplot(aes_string(x=signature, y='car', fill='car', color='car'), alpha=0.7, outlier.shape=NA) + 
    scale_fill_manual(values=c(car_colors, Untransduced='grey80')) +
    scale_color_manual(values=c(car_colors, Untransduced='grey80')) +
    theme_bw() + theme(panel.grid=element_blank(), text=element_text(family="Arial Unicode MS")) + labs(x=x_lab)
  
  return(plot)
}

plot_sig_by_cluster <- function(seurat_obj, vis, signature, x_lab, cluster_col='VISION_Clusters') {

  DefaultAssay(seurat_obj) <- 'RNA'
  sig_scores <- getSignatureScores(vis)

  sig.data <- data.table(
    cbind(seurat_obj@meta.data, sig_scores))
  
  plot <- ggplot(sig.data) +
  geom_boxplot(aes_string(x=signature, y=cluster_col, fill=cluster_col, color=cluster_col), alpha=0.7, outlier.shape=NA) + 
  theme_bw() + theme(panel.grid=element_blank(), text=element_text(family="Arial Unicode MS")) + labs(x=x_lab)
  
  return(plot)
}

car_feature_dotplot <- function(seurat_obj, features, assay='RNA') {

  #Missing features?
  #print(setdiff(features, intersect(features, rownames(seurat_obj@assays$RNA@scale.data))))
  
  DefaultAssay(seurat_obj) <- assay
  
  feature_dt <- melt(
    data.table(cbind(t(seurat_obj@assays[[assay]]@scale.data[features,]), seurat_obj@meta.data))[
      !(car %in% c('K_only','KLRG1','Untransduced'))], 
    measure.vars=features, variable.name='gene')[,
      `:=`(pct_exp_global=sum(value > 0)/.N, mean_exp_global=mean(value[value > 0])), by=c('CD4v8','gene')][,
        list(pct_exp=(sum(value > 0)/.N), mean_exp_fldch=mean(value[value > 0])/mean_exp_global), 
        by=c('car','CD4v8','gene')]
  
  return(ggplot(feature_dt[CD4v8 != 'intermediate']) + 
    geom_point(
      aes(x=factor(car, levels=c(car_order,'4-1BB', 'Zeta')), 
          y=factor(gene, levels=features), 
          color=mean_exp_fldch, 
          size=pct_exp)) + 
    facet_grid(~CD4v8) + 
    scale_color_distiller('Relative Expression Change\nfrom Mean',direction=1, palette='PiYG') + 
    theme_minimal() + 
    theme(panel.grid=element_blank()) + 
    scale_size_continuous('% Cells Expressing',range=c(2,8)) +
    labs(x='CAR',y='Gene'))

}

signature_plot <- function(seurat_obj, vis_obj, this_signature) {

  this_sig_geneset <- vis_obj@SigGeneImportance[[this_signature]]
  ordered_filtered_sig_geneset <- sort(
    this_sig_geneset[intersect(VariableFeatures(seurat_obj, assay='RNA'), names(this_sig_geneset))],
    decreasing=T)
  
  gene_score <- ggplot(data.table(gene=names(ordered_filtered_sig_geneset), score=(ordered_filtered_sig_geneset))) + 
    geom_tile(aes(x='score', y=reorder(gene, score), fill=score)) + 
    scale_fill_distiller(
      palette='RdBu', 
      limits=c(-max(abs(ordered_filtered_sig_geneset)), max(abs(ordered_filtered_sig_geneset)))) + 
    labs(x='',y='') + theme_minimal() + 
    theme(panel.grid=element_blank()) + 
    scale_x_discrete(expand=c(0,0))
  
  plot_grid(
    DoHeatmap(seurat_obj,
    features=names(ordered_filtered_sig_geneset), 
    assay='RNA', group.by = 'VISION_Clusters'),
    plot_grid(
        cluster_car_pct_plot(seurat_obj, 'VISION_Clusters'),
        plot_sig_by_car(seurat_obj, vis_obj, this_signature, x_lab=this_signature), 
        plot_sig_by_cluster(seurat_obj, vis_obj, this_signature, x_lab=this_signature), ncol=1),
    ncol=2,
    align = 'h', axis = 'tb')
}

### NEW SCT/CLUSTERING Functions

run_sct <- function(seurat_obj, pct.mt.cutoff=13) {
  
  # remove high MT-RNA cells (dying/dead?)
  DefaultAssay(seurat_obj) <- 'RNA'
  seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")
  plots$mt.cutoff <- FeatureScatter(
    seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
    geom_hline(yintercept=pct.mt.cutoff)
  seurat_obj <- subset(seurat_obj, subset=(percent.mt < pct.mt.cutoff))
  
  # get cell cycle scores
  DefaultAssay(seurat_obj) <- 'RNA'
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  #https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html
  seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  seurat_obj <- RunPCA(seurat_obj, features = c(s.genes, g2m.genes), reduction.name = 'cc_pca')
  
  #SCTransform for RNA
  DefaultAssay(seurat_obj) <- 'RNA'
  seurat_obj <- SCTransform(seurat_obj, 
    vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, reduction.name = 'rna_pca')
  
  #SCTransform for ADT
  DefaultAssay(seurat_obj) <- 'ADT'
  seurat_obj <- SCTransform(
    seurat_obj, assay='ADT', new.assay.name = 'SCT_ADT', 
    vars.to.regress = c("S.Score", "G2M.Score"))
  seurat_obj <- RunPCA(
    seurat_obj, features = rownames(seurat_obj[["SCT_ADT"]]), 
    nfeatures.print = 10, reduction.name='adt_pca')
  
  return(seurat_obj)
}

# SCTransform and UMAP of a subset of cells
sct_umap <- function(seurat_obj, pct.mt.cutoff=13, cluster_res=1.3, redo_sct=F, title="") {
  
  plots = list()
  
  # # perform sctransform on 
  # if (redo_sct | !('SCT' %in% names(seurat_obj@assays) | !('SCT_ADT') %in% names(seurat_obj@assays))) 
  #   seurat_obj <- run_sct(seurat_obj, pct.mt.cutoff)
  
  #DefaultAssay(seurat_obj) <- 'SCT_INT'
  #seurat_obj <- RunHarmony(seurat_obj, 'donor', reduction='rna_pca',reduction.save='harmony_rna_pca')
  
  #WNN graph for both PCAs
  seurat_obj <- FindMultiModalNeighbors(
    seurat_obj, reduction.list = list("rna_pca", "adt_pca"), 
    dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
  )
  
  cluster_resolution <- 1.3
  cluster_name <- paste0('wsnn_res.',as.character(cluster_resolution))
  
  seurat_obj <- RunUMAP(seurat_obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  seurat_obj <- FindClusters(seurat_obj, 
    graph.name = "wsnn", algorithm = 3, 
    resolution = cluster_resolution, verbose = FALSE)
  
  pl_clusters <- DimPlot(seurat_obj, reduction = 'wnn.umap', label = TRUE,
    repel = TRUE, label.size = 2.5) + NoLegend()
  pl_cars <- DimPlot(seurat_obj, reduction = 'wnn.umap', group.by = 'sample', label = TRUE,
    repel = TRUE, label.size = 2.5)
  pl_ktype <- DimPlot(seurat_obj, reduction = 'wnn.umap', group.by = 'k_type', label = TRUE,
    repel = TRUE, label.size = 2.5) + NoLegend()
  pl_ccphase <- DimPlot(seurat_obj, reduction = 'wnn.umap', group.by = 'Phase', label = TRUE,
    repel = TRUE, label.size = 2.5) + NoLegend()
  pl_t_type <- DimPlot(seurat_obj, reduction = 'wnn.umap', group.by = 't_type', label = TRUE,
    repel = TRUE, label.size = 2.5) + NoLegend()
  pl_donors <- DimPlot(seurat_obj, reduction = 'wnn.umap', group.by = 'donor', label = TRUE,
    repel = TRUE, label.size = 2.5) + NoLegend()
  plots$combined_umap <- pl_clusters + pl_ktype + pl_ccphase + pl_cars + pl_t_type + pl_donors
  
  DefaultAssay(seurat_obj) <- 'SCT_ADT_INT'
  adt_umap_plots <- umap_plot(seurat_obj, 
    plot_title_text=paste(title,"ADT WNN"), 
    assay='SCT_ADT_INT', cluster_name=cluster_name, reduction= 'wnn.umap')
  
  DefaultAssay(seurat_obj) <- 'SCT_INT'
  rna_umap_plots <- umap_plot(seurat_obj, 
    plot_title_text=paste(title,"RNA WNN"), assay='SCT_INT', cluster_name=cluster_name, reduction= 'wnn.umap')

  plots$adt_umap_plots <- adt_umap_plots$umap_plots
  plots$rna_umap_plots <- rna_umap_plots$umap_plots
  
  return(list(
    obj=seurat_obj, plots=plots, 
    adt_markers= adt_umap_plots$marker_genes,
    rna_markers= rna_umap_plots$marker_genes))
}

# Map clusters from one object to another
map_clusters <- function(from_obj, to_obj, from_cluster_name, to_cluster_name, overwrite.rm=F) {
  from_clust_data <- data.table(from_obj@meta.data, keep.rownames=T)[,
  list(rn, clust= get(from_cluster_name))][
    data.table(to_obj@meta.data, keep.rownames=T), on=('rn')]

  if (overwrite.rm || !(to_cluster_name %in% names(to_obj@meta.data))) {
    to_obj@meta.data[[to_cluster_name]] <- from_clust_data$clust
  } else {
    to_obj@meta.data[[to_cluster_name]][
      !is.na(from_clust_data$clust)] <- from_clust_data$clust[
        !is.na(from_clust_data$clust)]
  }
  return(to_obj)
}

# Map activation % to each cluster (CD19+ vs baseline), of CAR clusters, (not KLRG1/Untransduced)
map_clust_act <- function(seurat_obj, cluster_name) {
  
  car_clust_data <- data.table(
  seurat_obj@meta.data)[, list(n_clust_car= .N), 
  by=c(cluster_name, 'car', 'k_type')][, 
    pct_clust_car := n_clust_car/sum(n_clust_car), 
    by=c('car', 'k_type')][,
      cd19_clust_car_l2fc := log2(pct_clust_car[k_type == 'cd19+']/pct_clust_car[k_type == 'baseline']), 
      by=c(cluster_name, 'car')][, 
        cd19_clust_act_l2fc := mean(cd19_clust_car_l2fc[!(car %in% c('KLRG1','Untransduced'))], na.rm=T), 
        by=cluster_name]
  
  #if NaN, assume only CD19+, make Inf.
  car_clust_data[is.nan(cd19_clust_act_l2fc), cd19_clust_act_l2fc := Inf]
  
  car_clust_data[, pct_clust_car_max := max(pct_clust_car), by=c(cluster_name, 'car')]
  car_clust_data[, is_act := max(pct_clust_car) > pct_clust_car[k_type == 'baseline'], by=c(cluster_name, 'car')]
  car_clust_data[is.na(is_act), is_act := T]
  car_clust_data[, is_act_avg := sum(is_act) > length(unique(car)), by=cluster_name]

  act_plot <- ggplot(car_clust_data) + 
    geom_point(aes_string(x=cluster_name, y='car', fill='cd19_clust_car_l2fc', size='pct_clust_car_max'), shape=21) + 
    scale_fill_distiller(palette='PRGn', direction=1) + scale_size_area(max_size=12) + theme_minimal() + 
    facet_grid(~is_act_avg, scales='free_x')
  
  activation_colors <- ggplot_build(
    ggplot(unique(car_clust_data[, list(cd19_clust_act_l2fc, get(cluster_name))])) + 
    geom_tile(aes_string(x='cluster_name', y=1, fill='cd19_clust_act_l2fc')) + 
    scale_fill_distiller(palette='PRGn', direction=1))$data[1][[1]][c('fill','x')]
  
  act_umap_plot <- DimPlot(
    seurat_obj, group.by=cluster_name,
    cols=activation_colors$fill[order(activation_colors$x)], label=T)
  
  return(list(
    car_clust_data=car_clust_data,
    act_plot=act_plot,
    act_umap_plot= act_umap_plot,
    activation_colors=activation_colors))
}

# Map two cluster types and make plots
compare_clusters <- function(seurat_obj, clust_a, clust_b) {
  
  clust_co_counts <- data.table(seurat_obj@meta.data)[, 
    list(n_clust_car=.N), by=c(clust_a, clust_b)]
  
  clust_co_counts[!is.na(get(clust_a)) & !is.na(get(clust_b)), c(clust_a, 'pct_a_in_b','best_a_in_b') := list(
      get(clust_a), n_clust_car/sum(n_clust_car), get(clust_a)[which.max(n_clust_car/sum(n_clust_car))]), 
    by=c(clust_b)]

  clust_co_counts[!is.na(get(clust_b)) & !is.na(get(clust_a)), c(clust_b, 'pct_b_in_a', 'best_b_in_a') := list(
      get(clust_b), n_clust_car/sum(n_clust_car), get(clust_b)[which.max(n_clust_car/sum(n_clust_car))]), 
    by=c(clust_a)]
  
  a_v_b <- clust_co_counts[,
    list(get(clust_a), pct_clust=n_clust_car/sum(n_clust_car)), by=c(clust_b)]
  b_v_a <- clust_co_counts[,
    list(get(clust_b), pct_clust=n_clust_car/sum(n_clust_car)), by=c(clust_a)]
  
  names(b_v_a)[c(1,2)] <- c(clust_a, clust_b)
  names(a_v_b)[c(1,2)] <- c(clust_b, clust_a)

  plot_a_in_b <- ggplot(a_v_b) + 
    geom_tile(aes_string(x=clust_a, y=clust_b, fill='pct_clust')) +
    scale_fill_distiller(palette='Greens', direction=1) +
    theme_bw()
  plot_b_in_a <- ggplot(b_v_a) + 
    geom_tile(aes_string(x=clust_b, y=clust_a, fill='pct_clust')) +
    scale_fill_distiller(palette='Purples', direction=1) +
    theme_bw()
  
  return(list(clust_co_counts=clust_co_counts, plot_a_in_b=plot_a_in_b, plot_b_in_a=plot_b_in_a))
}

# cd4/cd8 clustering
reassign_t_type <- function(scrna, n_clust=6, plot_obj) {
  
  t_type_adt <- data.table(t(scrna@assays$SCT_ADT@scale.data[c('CD8A.adt','CD4.adt'),]))
  fit <- Mclust(t_type_adt[, c('CD8A.adt','CD4.adt')], G=n_clust)
  t_type_adt$clust <- fit$classification
  doublet_clust <- which.min(table(t_type_adt$clust))
  cd4_clust <- which(fit$parameters$mean[1,] < fit$parameters$mean[2,] & 1:n_clust != doublet_clust)
  cd8_clust <- which((1:n_clust != doublet_clust) & !(1:n_clust %in% cd4_clust))
  
  t_type_adt$doublet_score <- -log(fit$z[,doublet_clust])
  t_type_adt[, t_type := ifelse(clust %in% cd8_clust, 'CD8', 'CD4')]
  t_type_adt[clust == doublet_clust, t_type := 'DP']
  #t_type_adt[, t_type_old := scrna@meta.data$CD4v8]
  
  plot_obj$t_type_mmodel <- plot_grid(
    ggplot(t_type_adt) + 
      geom_point(aes(x=CD8A.adt, y=CD4.adt, shape=factor(clust), color=doublet_score)) + 
      scale_color_distiller(palette='Spectral') +
      theme_bw(),
    ggplot(t_type_adt) + 
      geom_point(aes(x=CD8A.adt, y=CD4.adt, color=factor(clust))) +
      theme_bw(),
    ggplot(t_type_adt) + 
      geom_point(aes(x=CD8A.adt, y=CD4.adt, color=t_type)) +
      theme_bw())
    # ggplot(t_type_adt) + 
    #   geom_point(aes(x=CD8A.adt, y=CD4.adt, color=t_type_old)) +
    #   theme_bw())
  
  # pick cluster with fewest cells, assign them as double positives
  scrna@meta.data$t_type_dp <- t_type_adt[, clust == doublet_clust]
  scrna@meta.data$t_type <- t_type_adt$t_type
  
  # plot_obj$t_type_reassign_table <- gt(as.data.frame.matrix(
  #   table(scrna@meta.data$t_type, scrna@meta.data$CD4v8)))
  
  #scrna@meta.data$t_type_old <- scrna@meta.data$CD4v8
  scrna@meta.data$CD4v8 <- scrna@meta.data$t_type
  
  return(list(seurat_obj=scrna, plot_obj=plot_obj))
}

get_ident_avg <- function(seurat_obj, ident='seurat_clusters', assays=c('SCT_ADT','SCT')) {
  Idents(seurat_obj) <- ident
  av.exp <- AverageExpression(seurat_obj, assays=assays, use.scale=T)
  return(rbindlist(lapply(av.exp,as.data.table)))
}

make_dendroheatmap <- function(
  cor, row_labs=NULL, col_labs=NULL, title='', clust_y=T, clust_x=T,
  dendro_units=grid::unit(0.2, "null"),
  limits=c(0,1), palette='YlOrBr',
  legend_theme=theme(legend.position = "none"))
{

  if (is.null(row_labs)) {
    row_labs <- rownames(cor)
  }
  if (is.null(col_labs)) {
    col_labs <- colnames(cor)
  }
  
  cor.exp <- as.data.frame(cor)
  names(cor.exp) <- col_labs
  rownames(cor.exp) <- row_labs
  cor.exp.dt <- data.table(cor.exp, keep.rownames=T)
  cor.df <- melt(cor.exp.dt, id.vars='rn', variable.name='x')
  names(cor.df)[1] <- 'y'
  
  png("/dev/null");
  hmap_data <- heatmap(cor, keep.dendro=T)
  dev.off()
  
  cor_dendro_y <- hmap_data$Rowv
  cor_dendro_x <- hmap_data$Colv

  # row order
  cor.df[, y := factor(y, levels=rownames(cor.exp)[hmap_data$rowInd])]
  cor.df[, x := factor(x, levels=names(cor.exp)[hmap_data$colInd])]

  # heatmap
  var_heatmap <- ggplot(cor.df) + 
    geom_tile(aes(x = x, y = y, fill = value),
      #colour = "black", size = .20) + 
    )+
    scale_fill_distiller('Correlation', 
      palette=palette, limits=limits, oob=scales::squish, direction=1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill='white'),
        axis.title.y = element_blank())
  
  # col dendro
  dendro_data_col_x <- dendro_data(cor_dendro_x, type = "rectangle")
  dendro_data_col_y <- dendro_data(cor_dendro_y, type = "rectangle")
  
  plot_dendro_x <- axis_canvas(var_heatmap, axis = "x", coord_flip=F) + 
    geom_segment(data = segment(dendro_data_col_x), 
        aes(x = x, y = y, xend = xend, yend = yend))
  
  plot_dendro_y <- axis_canvas(var_heatmap, axis = "y", coord_flip=T) + 
    geom_segment(data = segment(dendro_data_col_y),
        aes(x = x, y = y, xend = xend, yend = yend)) + coord_flip()

  
  if (clust_x) {
    var_heatmap <- insert_xaxis_grob(var_heatmap, 
      plot_dendro_x, dendro_units, position = "top")
  }
  if (clust_y) {
    var_heatmap <- insert_yaxis_grob(var_heatmap, 
      plot_dendro_y, dendro_units, position = "right")
  }
  return(ggdraw(var_heatmap))
}

plot_clust_by_act <- function(seurat_obj, cluster_name='seurat_clusters') {

    car_clust_data <- data.table(
      seurat_obj@meta.data)[, list(n_clust_car= .N), 
      by=c(cluster_name, "car", "k_type")][, 
        pct_clust_car := n_clust_car/sum(n_clust_car), 
        by=.(car, k_type)][,
          cd19_clust_car_l2fc := log2(
              pct_clust_car[k_type == ifelse(.BY[2] == 'Untransduced', 'bead_stim', 'cd19+')]/
                pct_clust_car[k_type == 'baseline']), 
          by=.(get(cluster_name), car)][, 
            cd19_clust_act_l2fc := mean(
                cd19_clust_car_l2fc[!(car %in% c('KLRG1','Untransduced')) | k_type == 'bead_stim'], na.rm=T), 
            by=get(cluster_name)]
    
    car_clust_data[, pct_clust_car_max := max(pct_clust_car), by=.(get(cluster_name), car)]
    car_clust_data[, is_act := max(pct_clust_car) > pct_clust_car[k_type == 'baseline'], 
        by=.(get(cluster_name), car)]
    car_clust_data[is.na(is_act), is_act := T]
    car_clust_data[, is_act_avg := sum(is_act) > length(unique(car)), by=.(get(cluster_name))]
    
    clust_pt_plot <- ggplot(car_clust_data) + 
        geom_point(
            aes(x=get(cluster_name), y=car, fill=cd19_clust_car_l2fc, size=pct_clust_car_max), shape=21) + 
        scale_fill_distiller(palette='PRGn', direction=1) + scale_size_area(max_size=12) + 
        theme_minimal() + 
        facet_grid(~is_act_avg, scales='free')
    
    activation_colors <- ggplot_build(
      ggplot(unique(car_clust_data[, list(cd19_clust_act_l2fc, seurat_clusters)])) + 
      geom_tile(aes(x=get(cluster_name), y='1', fill=cd19_clust_act_l2fc)) + 
      scale_fill_distiller(palette='PRGn', direction=1))$data[1][[1]][c('fill','x')]
    
    act_dimplot <- DimPlot(seurat_obj,
        cols=activation_colors$fill[order(activation_colors$x)], label=T)
    
    return(list(
      car_clust_data=car_clust_data,
        clust_pt_plot=clust_pt_plot,
        activation_colors=activation_colors,
        act_dimplot=act_dimplot))
}

#dendrogram heatmap of genes to clusters/cars
heatmap_clust_v_genes <- function(
  seurat_obj, split_name, assay, features, use_pct=F, lfc_trunc=1.5, use_slot='scale.data') {
  
  slot_data <- slot(seurat_obj@assays[[assay]], use_slot)
  
  missing_features <- features[which(!(features %in% rownames(slot_data)))]
  kept_features <- features[!(features %in% missing_features)]
  message(paste('These features are missing:', paste(missing_features, collapse=', ')))
  
  feature_dt <- unique(melt(
    data.table(cbind(t(slot_data[kept_features,]), seurat_obj@meta.data)), 
    measure.vars=kept_features, variable.name='gene')[,
      `:=`(pct_exp_global=sum(value > 0)/.N, mean_exp_global=mean(value[value > 0])), by=c('CD4v8','gene')][,
        list(pct_exp=(sum(value > 0)/.N), mean_exp_fldch=log2(mean(value[value > 0])/mean_exp_global)), 
        by=c(split_name,'gene')])
  names(feature_dt)[1] <- 'clusters'
  
  feature_dt[, log2_fc_trunc := ifelse(
    test= abs(mean_exp_fldch) > lfc_trunc, 
    yes= lfc_trunc * sign(mean_exp_fldch), 
    no= mean_exp_fldch)]
  
  feature_dt[, gene := factor(gene, levels=kept_features)]
  
  if (use_pct) {
    limits=c(0,1)
    palette='YlOrBr'
    feat.dt <- dcast(feature_dt[!is.na(clusters)][, list(clusters, gene, pct_exp)], gene ~ clusters)
  } else {
    limits=c(-lfc_trunc,lfc_trunc)
    palette='PiYG'
    feat.dt <- dcast(feature_dt[!is.na(clusters)][, list(clusters, gene, log2_fc_trunc)], gene ~ clusters)
  }
  
  feat.mat <- as.matrix(feat.dt[, -'gene'])
  rownames(feat.mat) <- feat.dt[, gene]
  
  dendro_hm_plot <- make_dendroheatmap(
    feat.mat, row_labs=rownames(feat.mat), col_labs=colnames(feat.mat),
    limits=limits, palette=palette)
  
  return(list(
    plot=dendro_hm_plot,
    data=feature_dt
  ))
}

# compare cars/clusts
diff_genes_in_clusts <- function(
  seurat_obj, ident_a, ident_b, cluster='seurat_clusters', logfc.threshold=0.25, min.pct=0.25,
  method='MAST', label_top_adt=20, label_top_rna=5, adt_only=F) {
  
  old_ident <- Idents(seurat_obj)
  Idents(seurat_obj) <- cluster 
  
  diff_genes_adt <- FindMarkers(
    seurat_obj, assay='SCT_ADT_INT',
    ident.1 = ident_a, ident.2=ident_b, method=method,
    logfc.threshold = logfc.threshold, min.pct=min.pct)
  
  if (!adt_only) {
    diff_genes_rna <- FindMarkers(
      seurat_obj, assay='SCT_INT',
      ident.1 = ident_a, ident.2=ident_b, method=method,
      logfc.threshold = logfc.threshold, min.pct=min.pct)
    
    diff_genes <- rbind(data.table(diff_genes_rna, keep.rownames=T)[, assay := 'RNA'],
          data.table(diff_genes_adt, keep.rownames=T)[, assay := 'ADT'])
  } else {
    diff_genes <- data.table(diff_genes_adt, keep.rownames=T)[, assay := 'ADT']
  }
  
  diff_genes[, diff_dir := sign(avg_log2FC)]
  diff_genes[assay=='RNA' & p_val_adj < 0.05, display := order(log(p_val_adj)) %in% 1:label_top_rna, by=diff_dir]
  diff_genes[assay=='ADT' & p_val_adj < 0.05, display := order(log(p_val_adj)) %in% 1:label_top_adt, by=diff_dir]
  diff_genes[display==T, label := rn, by=diff_dir]
  
  #mean_x_neg <- diff_genes[diff_dir == -1, mean(avg_log2FC)]
  #mean_x_pos <- diff_genes[diff_dir == 1, mean(avg_log2FC)]
  mean_y_max <- diff_genes[, max(-log(p_val_adj))]
  
  
  plot_volcano <- ggplot(diff_genes, aes(x=avg_log2FC, y=-log(p_val_adj), label=label, color=assay)) + 
    geom_point() + 
    geom_text_repel() +
    scale_color_manual(values=c('red','black')) +
    annotate(geom='label', x=-0.5, y=mean_y_max+2, label=paste(as.character(ident_b), collapse=', ')) +
    annotate(geom='label', x=0.5, y=mean_y_max+2, label=paste(as.character(ident_a), collapse=', ')) +
    geom_hline(yintercept=-log(0.05), linetype=2) +
    geom_vline(xintercept=0, linetype=2) +
    theme_minimal()
  
  Idents(seurat_obj) <- old_ident
  
    return(list(diff_genes=diff_genes, plot_volcano=plot_volcano))
}

#function to get GSEA genes
get_geneset <- function(geneset_name) {
  tmp <- tempfile()
  download.file(paste0(
    'http://www.gsea-msigdb.org/gsea/msigdb/download_geneset.jsp?geneSetName=',
    geneset_name,
    '&fileType=txt'), tmp)
  gene_list <- fread(tmp)
  file.remove(tmp)
  return(gene_list[-1][,get(geneset_name)])
}

map_geneset <- function(seurat_obj, geneset_list) {

  for (geneset in geneset_list) {
    if (!(paste0(geneset,'1') %in% names(seurat_obj@meta.data))) {
      message(paste0("Adding Geneset: ",geneset))
      DefaultAssay(seurat_obj) <- 'RNA'
      seurat_obj <- AddModuleScore(seurat_obj, features=list('1'=get_geneset(geneset)), name=geneset)
    }
  }
  return(seurat_obj)
}

map_geneset_parallel <- function(seurat_obj, geneset_list, assay='RNA', nbin=24, cores=8) {

  get_module_score_parallel <- function(geneset_i) {
    geneset <- geneset_list[geneset_i]
    message(paste0("Adding Geneset: ",geneset))
    DefaultAssay(seurat_obj) <- assay
    module_score <- tryCatch({
      AddModuleScoreParallel(
        seurat_obj, features=list('1'=get_geneset(geneset)), name=geneset, nbin=nbin)
    }, error= function(e) {
      print(e)
      warning(paste0(geneset, ' failed.'))
      return(data.frame())
    })
    return(module_score)
  }

  geneset_scores <- bplapply(
    1:length(geneset_list),
    get_module_score_parallel,BPPARAM=MulticoreParam(cores, log=T))
  
  geneset_scores <- do.call('cbind', geneset_scores)
  names(geneset_scores) <- gsub('(.*)1$','\\1',names(geneset_scores))
  return(geneset_scores)
}



#remove bead stim, cd30, censored clusters, DP cells, etc
prune_cells <- function(seurat_obj, cars=NULL, k_types=NULL, subtypes=NULL, rescale=F, t_types=NULL) {
  no_cd30 <- seurat_obj$car != 'CD30' & seurat_obj$car != 'K_only'
  no_beadstim <- seurat_obj$k_type != 'bead_stim'
  no_rmclust <- seurat_obj$Subtype != 'REMOVE' &
    !is.na(seurat_obj$Subtype) &
      seurat_obj$Subtype != 'NA'
  
  if (length(cars)) {
    cars <- seurat_obj$car %in% cars
  } else {
    cars <- rep(T,length(seurat_obj$car))
  }
  
  if (length(subtypes)) {
    subtypes <- seurat_obj$Subtype %in% subtypes
  } else {
    subtypes <- rep(T,length(seurat_obj$Subtype))
  }


  if (length(k_types)) {
    k_types <- seurat_obj$k_type %in% k_types
  } else {
    k_types <- rep(T,length(seurat_obj$k_type))
  }
  
  if (length(t_types)) {
    t_types <- seurat_obj$t_type %in% t_types
  } else {
    t_types <- rep(T,length(seurat_obj$t_type))
  }
  
  new_obj <- seurat_obj[, no_cd30 & no_beadstim & no_rmclust & cars & k_types & subtypes & t_types]
  if (rescale) {
    new_obj <- ScaleData(new_obj, assay='SCT_INT')
    new_obj <- ScaleData(new_obj,assay= 'SCT_ADT_INT')
  }
  return(new_obj)
}

