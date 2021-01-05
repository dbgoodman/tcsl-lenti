
library(mclust)
library(Seurat)

DoubleHTODemux <- function (object, assay = "HTO", positive.quantile = 0.99, 
    init = NULL,  nstarts = 100, kfunc = "clara", nsamples = 100, seed = 42,
    ncenters=NULL, verbose = TRUE) 
{
    
    MaxN <- function(x, N = 2){
      len <- length(x)
      if (N > len) {
        warning('N greater than length(x).  Setting N=length(x)')
        N <- length(x)
      }
      sort(x, partial = len - N + 1)[len - N + 1]
    }
    
    whichpart <- function(x, n=30) {
      nx <- length(x)
      p <- nx-n
      xp <- sort(x, partial=p)[p]
      which(x > xp)
    }
  
    if (!is.null(x = seed)) {
        set.seed(seed = seed)
    }
    
    data <- GetAssayData(object = object, assay = assay)
    
    if (is.null(ncenters)) {
      ncenters = nrow(x = data) + 1
    }

    counts <- GetAssayData(object = object, assay = assay, slot = "counts")[, 
        colnames(x = object)]
    counts <- as.matrix(x = counts)
    init.clusters <- clara(x = t(x = GetAssayData(object = object, 
        assay = assay)), k = ncenters, samples = nsamples)
    Idents(object = object, cells = names(x = init.clusters$clustering), 
        drop = TRUE) <- init.clusters$clustering
    
    average.expression <- AverageExpression(object = object, 
        assays = assay, verbose = FALSE)[[assay]]
    
    if (sum(average.expression == 0) > 0) {
        stop("Cells with zero counts exist as a cluster.")
    }
    
    discrete <- GetAssayData(object = object, assay = assay)
    discrete[discrete > 0] <- 0
    for (iter in rownames(x = data)) {
        values <- counts[iter, colnames(object)]
        values.use <- values[WhichCells(
          object = object, 
          idents = levels(x = Idents(object = object))[[
            which.min(x = average.expression[iter, ])]])]
        fit <- suppressWarnings(expr = fitdistrplus::fitdist(
            data = values.use, 
            distr = "nbinom"))
        cutoff <- as.numeric(x = quantile(
          x = fit, probs = positive.quantile)$quantiles[1])
        discrete[iter, names(x = which(x = values > cutoff))] <- 1
        if (verbose) {
            message(paste0("Cutoff for ", iter, " : ", cutoff, 
                " reads"))
        }
    }
    
    npositive <- colSums(x = discrete)
    classification.global <- npositive
    classification.global[npositive == 0] <- "Negative"
    classification.global[npositive == 1] <- "Singlet"
    classification.global[npositive == 2] <- "Doublet"
    classification.global[npositive == 3] <- "Triplet"
    classification.global[npositive > 3] <- "Quadruplet"
    donor.id = rownames(x = data)
    
    hash.maxIDs <- apply(X = data, MARGIN = 2, 
      FUN = function(row) order(-row)[1:4])
    hash.maxValues <- apply(X = data, MARGIN = 2, 
      FUN = function(row) sort(-row)[1:4])
    
    hash.maxID <- rownames(data)[hash.maxIDs[1,]]
    hash.secondID <- rownames(data)[hash.maxIDs[2,]]
    hash.thirdID <- rownames(data)[hash.maxIDs[3,]]
    hash.fourthID <- rownames(data)[hash.maxIDs[4,]]

    hash.max <- hash.maxValues[1,]
    hash.second <- hash.maxValues[2,]
    hash.third <- hash.maxValues[3,]
    hash.fourth <- hash.maxValues[4,]

    hash.margin <- hash.second - hash.third
    hash.mean.margin <- mean(c(hash.max, hash.second)) - hash.third
    hash.mean.dblmgn <- mean(c(hash.max, hash.second)) - mean(c(hash.third, hash.fourth))
    
    doublet_id <- sapply(X = 1:length(x = hash.maxID), FUN = function(x)
      {
        return(paste(
          sort(c(
            as.numeric(gsub('Hashtag-','',hash.maxID[x])),
            as.numeric(gsub('Hashtag-','',hash.secondID[x])))),
          collapse = "_"))
      }
    )

    classification <- classification.global
    
    classification[classification.global == "Negative"] <- "Negative"
    classification[classification.global == "Singlet"] <- gsub(
      'Hashtag-','', hash.maxID[which(x = classification.global == "Singlet")])
    classification[classification.global == "Doublet"] <- doublet_id[
      which(x = classification.global == "Doublet")]
    classification[classification.global == "Triplet"] <- paste(
      doublet_id[which(x = classification.global ==  "Triplet")],
        gsub('Hashtag-','',
          hash.thirdID[which(x = classification.global ==  "Triplet")]),
        sep='+')
        
    classification.metadata <- data.frame(hash.maxID, hash.secondID,
        hash.margin, hash.mean.margin, hash.mean.dblmgn,
        classification, init.clusters$clustering,
        classification.global)
    
    colnames(x = classification.metadata) <- paste(assay, c("maxID", 
        "secondID", "margin", "mean_margin", "mean_dblmargin", "classification","clusters",
        "classification.global"), 
        sep = "_")
    
    object <- AddMetaData(object = object, metadata = classification.metadata)
    Idents(object) <- paste0(assay, "_classification")
    
    triplets <- rownames(x = object[[]])[which(object[[paste0(assay, 
        "_classification.global")]] == "Triplet")]
    singlets <- rownames(x = object[[]])[which(object[[paste0(assay, 
        "_classification.global")]] == "Singlet")]
    
    Idents(object = object, cells = singlets) <- "Singlets"
    Idents(object = object, cells = triplets) <- "Triplets"
    object$hash.ID <- Idents(object = object)
    return(object)
}

assign_HTO_sample <- function(seurat_obj, metadata=FALSE, assay='HTO', slot='counts',
    min_cells=200, doublet_cutoff=6, save_plots=F) {
  
  plots <- list()
  
  counts <- melt(
    data.table(
      t(as.matrix(GetAssayData(seurat_obj, assay = assay, slot = slot))),
      keep.rownames = T), id.vars='rn')
  names(counts) <- c('cell','hto','count')
  
  counts[count > 0, logcount := log10(count)]
  counts[, hto := factor(gsub('Hashtag-(\\d+)$','\\1', hto))]
  counts[, pct := count/sum(count), by='cell']
  
  fit_hto_count <- function(x, measure='count') {
    fit <- Mclust(
      data=counts[hto == x & !is.na(get(measure))][[measure]],
      G=2,
      na.rm=T)
    counts[hto == x & !is.na(get(measure)), score := -log(fit$uncertainty) * (fit$classification * 2 - 3)]
  }
  lapply(1:12, fit_hto_count)
  counts[abs(score) == Inf, score := sign(score) * counts[abs(score) < Inf, max(abs(score))]]
  
  plots$hto_countdist <- ggplot(counts) + 
    geom_density(aes(x=count)) + facet_grid(~hto)

  plots$hto_countdist_clust <- ggplot(counts) + 
    geom_density(aes(x=count, color=score > 0, group=(score > 0))) + facet_grid(~hto)
  
  best_counts <- counts[, list(
    rank=1:4,
    n=sum(count),
    hto=hto[order(-score, -pct)][1:4],
    pct=pct[order(-score, -pct)][1:4],
    score=score[order(-score, pct)][1:4]), by='cell']
  
  cell_summary <- best_counts[, 
    list(
      HTO_call=paste(sort(as.numeric(as.character((hto[1:2])))),collapse='_'),
      HTO_ab=paste(sort(as.numeric(as.character((hto[c(1,2)])))),collapse='_'),
      HTO_ac=paste(sort(as.numeric(as.character((hto[c(1,3)])))),collapse='_'),
      HTO_bc=paste(sort(as.numeric(as.character((hto[c(2,3)])))),collapse='_'),
      HTO_score_margin_BvC= score[2] - score[3],
      HTO_pct_margin_BvC= pct[2] / pct[3],
      HTO_a.score= score[1], HTO_b.score= score[2], HTO_c.score= score[3], HTO_d.score= score[4],
      HTO_a= hto[1], HTO_b= hto[2], HTO_c= hto[3], HTO_d= hto[4],
      HTO_a.pct= pct[1], HTO_b.pct= pct[2], HTO_c.pct= pct[3], HTO_d.pct= pct[4]),
    by='cell']
  
  cell_summary[, HTO_cd_pct := HTO_c.pct + HTO_d.pct]
  metadata.rows <- data.frame(cell_summary[rownames(seurat_obj@meta.data), on='cell'][, -'cell'])
  rownames(metadata.rows) <- rownames(seurat_obj@meta.data)
  seurat_obj <- AddMetaData(seurat_obj, metadata=metadata.rows)
  
  # plot doublets by ratio of C score to B-C
  plots$doublet_ratio <- ggplot(cell_summary) +
      geom_density_2d_filled(aes(x= HTO_score_margin_BvC, y=HTO_c.score), n=100, contour_var = "ndensity") +
      theme_bw() + scale_color_brewer(palette='Spectral') +
      geom_abline(slope=1, intercept=-doublet_cutoff) +
      facet_wrap(~HTO_call)
  
  # margin way with SCT:
  cell_summary[, HTO_is_doublet := HTO_score_margin_BvC - doublet_cutoff < HTO_c.score]
  cell_summary[, HTO_doublet_score := HTO_score_margin_BvC - HTO_c.score]
  cell_summary[, HTO_c_best := factor(ifelse(HTO_is_doublet, as.character(HTO_c), '0'))]
  cell_summary[, HTO_assign := 'ab']
  cell_summary[, HTO_rescue := F]
  cell_summary[, HTO_runnerup := HTO_c]
  
  #for tag combos that do not exist in the metadata, see if we can recover
  # ambiguous means multiple of ab, ac, bc, are incorrect.
  cell_summary[,
    HTO_ambiguous := (HTO_is_doublet | !(HTO_ab %in% metadata$HTO)) & (
      (HTO_ab %in% metadata$HTO) + (HTO_ac %in% metadata$HTO) + (HTO_bc %in% metadata$HTO) > 1)]
  
  plots$cell_type_count <- ggplot(cell_summary) +
      geom_bar(aes(x=as.character(HTO_c_best), fill=HTO_ambiguous)) +
      facet_wrap(~HTO_call, scales='free_x')
  
  cell_summary[!HTO_ambiguous & !(HTO_call %in% metadata$HTO) & (HTO_ac %in% metadata$HTO),
    `:=`(HTO_assign= 'ac', HTO_call= HTO_ac, HTO_runnerup= HTO_b, HTO_rescue= T)]
  cell_summary[!HTO_ambiguous & !(HTO_call %in% metadata$HTO) & (HTO_bc %in% metadata$HTO),
    `:=`(HTO_assign= 'bc', HTO_call= HTO_ac, HTO_runnerup= HTO_b, HTO_rescue= T)]
  cell_summary[HTO_is_doublet & !HTO_ambiguous,
    `:=`(HTO_rescue= T, HTO_is_doublet= F, HTO_runnerup= HTO_c)]
  
  cell_summary <- metadata[cell_summary, on=c(HTO='HTO_call'), nomatch=NA]
  names(cell_summary)[which(names(cell_summary) == 'HTO')] <- 'HTO_call'
  
  cell_summary[, HTO_status := '']
  cell_summary[HTO_rescue==T, HTO_status := 'Rescue']
  cell_summary[!HTO_is_doublet & HTO_doublet_score > doublet_cutoff, HTO_status := 'MidQual']
  cell_summary[!HTO_is_doublet & HTO_doublet_score < 0, HTO_status := 'LowQual']
  cell_summary[!HTO_is_doublet & HTO_doublet_score > doublet_cutoff+2, HTO_status := 'HighQual']
  cell_summary[!(HTO_call %in% metadata$HTO), HTO_status := 'NoMatch']
  cell_summary[HTO_is_doublet==T, HTO_status := 'Doublet']
  
  
  metadata.rows <- data.frame(cell_summary)
  rownames(metadata.rows) <- rownames(seurat_obj@meta.data)
  seurat_obj <- AddMetaData(seurat_obj, metadata=metadata.rows)
  
  if (save_plots) {
    seurat_obj@misc[['HTO_plots']] <- plots
  }
  return(seurat_obj)
}

# after original doublet removal, second pass to get rid of stubborn doublets. this didn't work by itself. :(
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
  
  return(scrna)
}