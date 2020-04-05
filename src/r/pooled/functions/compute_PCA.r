compute_PCA <- function(data.dt, group) {
  # Using car.bin.pct
  cast.car.bin.pct <- dcast(data.dt[k.type != "NA"], sort.group + batch + 
                              donor + timepoint + assay + t.type + CAR.align + 
                              k.type ~ bin, value.var = 'car.bin.pct')
  
  cast.car.bin.pct <- cast.car.bin.pct[, .(sort.group, batch, donor, 
                                           timepoint, assay, t.type, k.type, 
                                           CAR.align, bin.A = A, bin.B = B, 
                                           bin.C = C, bin.D = D)]
  
  cast.car.bin.pct[is.na(cast.car.bin.pct)] <- 0
  
  # compute principal components
  pca <- prcomp(cast.car.bin.pct[sort.group == group, 
                                 .(bin.A, bin.B, bin.C, bin.D)],
                center = T, scale. = T)
  
  # project data onto principal components
  projected.car.bin.pct <- scale(cast.car.bin.pct[sort.group == group,
                                                  .(bin.A, bin.B, bin.C, 
                                                    bin.D)],
                                 pca$center, pca$scale) %*% pca$rotation
  
  projected.car.bin.pct <- cbind(cast.car.bin.pct[sort.group == group,
                                                  .(sort.group, batch,
                                                    donor, timepoint, 
                                                    assay, t.type, k.type,
                                                    CAR.align)], 
                                 projected.car.bin.pct)
  
  # merge CAR.scores with projected.rel.bin.ratio
  ranks <- data.dt[sort.group == group &
                         bin == "A", 
                       .(CAR.align, CAR.score.rank = rank(CAR.score))]
  
  projected.car.bin.pct <- merge(projected.car.bin.pct, ranks, by = "CAR.align")
  
  return(projected.car.bin.pct)
}