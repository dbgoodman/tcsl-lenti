
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

assign_HTO_sample <- function(seurat_obj, assay='HTO', slot='counts', min_cells=200) {
  counts <- melt(
    data.table(
      t(as.matrix(GetAssayData(seurat_obj, assay = assay, slot = slot))),
      keep.rownames = T), id.vars='rn')
  names(counts) <- c('cell','hto','count')
  
  counts[count > 0, logcount := log10(count)]
  counts[, hto := factor(gsub('Hashtag-(\\d+)$','\\1', hto))]
  counts[, pct := count/sum(count), by='cell']
  
  #ggplot(counts) + geom_density(aes(x=count)) + facet_grid(~hto)
  
  fit_hto_count <- function(x, measure='pct') {
    fit <- Mclust(
      data=counts[hto == x & !is.na(get(measure))][[measure]],
      G=2,
      na.rm=T)
    counts[hto == x & !is.na(get(measure)), score := -log(fit$uncertainty) * (fit$classification * 2 - 3)]
  }
  lapply(1:12, fit_hto_count)
  counts[abs(score) == Inf, score := sign(score) * counts[abs(score) < Inf, max(abs(score))]]
  
  best_counts <- counts[, list(
    rank=1:4,
    n=sum(count),
    hto=hto[order(-score, -pct)][1:4],
    pct=pct[order(-score, -pct)][1:4],
    score=score[order(-score, pct)][1:4]), by='cell']
  
  cell_summary <- best_counts[, 
    list(
      HTO_call=paste(sort(as.numeric(as.character((hto[1:2])))),collapse='_'),
      HTO_score_margin_ABvC= mean(score[1], score[2]) - score[3],
      HTO_score_margin_BvC= score[2] - score[3],
      HTO_score_margin_ABvC= mean(score[1], score[2]) - mean(score[3], score[4]),
      HTO_pct_margin_ABvC= mean(pct[1], pct[2]) / pct[3],
      HTO_pct_margin_BvC= pct[2] / pct[3],
      HTO_pct_margin_ABvCD= mean(pct[1], pct[2]) / mean(pct[3], pct[4]),
      HTO_a.score= score[1], HTO_b.score= score[2], HTO_c.score= score[3], HTO_d.score= score[4],
      HTO_a= hto[1], HTO_b= hto[2], HTO_c= hto[3], HTO_d= hto[4],
      HTO_a.pct= pct[1], HTO_b.pct= pct[2], HTO_c.pct= pct[3], HTO_d.pct= pct[4]),
    by='cell']
  
  cell_summary[, HTO_cd_pct := HTO_c.pct + HTO_d.pct]
  metadata.rows <- data.frame(cell_summary[rownames(seurat_obj@meta.data), on='cell'][, -'cell'])
  rownames(metadata.rows) <- rownames(seurat_obj@meta.data)
  seurat_obj <- AddMetaData(seurat_obj, metadata=metadata.rows)
  
  
  # find doublets based on log ratio of second two hashtags to total read counts
  ncount_v_bad_tags <- data.table(seurat_obj@meta.data)[, `:=`(
    HTO_logfc_nCount= log10(nCount_HTO)-mean(log10(nCount_HTO)),
    HTO_logfc_pct_cd= log10(HTO_cd_pct)-mean(log10(HTO_cd_pct))),
    by='HTO_call']
  
  ncount_v_bad_tags[, HTO_is_doublet := -HTO_logfc_nCount + 0.5 < HTO_logfc_pct_cd]
  
  # ggplot(ncount_v_bad_tags) +
  #   geom_point(aes(y=logfc_HTO_pct_cd, x=logfc_nCount, color=HTO_call), alpha=0.2) +
  #   theme_bw() +
  #   geom_abline(slope=-1, intercept=0.5) +
  #   facet_wrap(~HTO_call)
  
  metadata.rows <- data.frame(ncount_v_bad_tags)
  rownames(metadata.rows) <- rownames(seurat_obj@meta.data)
  seurat_obj <- AddMetaData(seurat_obj, metadata=metadata.rows)
  
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