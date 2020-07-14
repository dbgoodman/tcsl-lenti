
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
    classification.global[npositive >= 3] <- "Triplet"
    donor.id = rownames(x = data)
    
    hash.maxIDs <- apply(X = data, MARGIN = 2, 
      FUN = function(row) order(-row)[1:3])
    hash.maxValues <- apply(X = data, MARGIN = 2, 
      FUN = function(row) sort(-row)[1:3])
    
    hash.maxID <- rownames(data)[hash.maxIDs[1,]]
    hash.secondID <- rownames(data)[hash.maxIDs[2,]]
    hash.thirdID <- rownames(data)[hash.maxIDs[3,]]

    hash.max <- hash.maxValues[1,]
    hash.second <- hash.maxValues[2,]
    hash.third <- hash.maxValues[3,]
    
    hash.margin <- hash.second - hash.third
    hash.mean.margin <- mean(c(hash.max, hash.second)) - hash.third
    
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
        hash.margin, hash.mean.margin, classification, init.clusters$clustering,
        classification.global)
    
    colnames(x = classification.metadata) <- paste(assay, c("maxID", 
        "secondID", "margin", "mean_margin", "classification","clusters",
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