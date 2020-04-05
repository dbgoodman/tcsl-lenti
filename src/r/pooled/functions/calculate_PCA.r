calculate_PCA <- function(data.dt, group, metric) {
  
  print(group)
  
  cast.car.bin.pct <- dcast(data.dt[k.type != "NA"], sort.group + batch + 
                              donor + timepoint + assay + t.type + CAR.align + 
                              k.type ~ bin, value.var = metric)
  
  cast.car.bin.pct <- cast.car.bin.pct[, .(sort.group, batch, donor, 
                                           timepoint, assay, t.type, k.type, 
                                           CAR.align, bin.A = A, bin.B = B, 
                                           bin.C = C, bin.D = D)]
  
  cast.car.bin.pct[is.na(cast.car.bin.pct)] <- 0
  
  # compute principal components
  pca <- prcomp(cast.car.bin.pct[sort.group == group, 
                                 .(bin.A, bin.B, bin.C, bin.D)],
                center = T, scale. = T)
  
  # calculate pca stats
  pca.dt <- data.table(pc = data.table(colnames(pca$rotation))[, 
                                                               PC := as.integer(gsub("[A-Z]", "",
                                                                                     V1))][, PC], 
                       sd = pca$sdev, 
                       var = pca$sdev^2, 
                       var.norm = pca$sdev^2/sum(pca$sdev^2), 
                       var.acc = cumsum(pca$sdev^2/sum(pca$sdev^2)))
  
  pca.dt <- cbind(unique(cast.car.bin.pct[sort.group == group,
                                          .(sort.group, batch,
                                            donor, timepoint, 
                                            assay, t.type, k.type)]), 
                  pca.dt)
  
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
  
  # merge CAR.scores and lengths with projected.rel.bin.ratio
  ranks <- data.dt[sort.group == group &
                     bin == "A", 
                   .(CAR.align, CAR.score.rank = rank(CAR.score))]
  
  projected.car.bin.pct <- merge(projected.car.bin.pct, ranks, by = "CAR.align")
  
  projected.car.bin.pct <- merge(projected.car.bin.pct,
                                 unique(data.dt[, .(CAR.align, len)]), 
                                 by = "CAR.align")
  
  return(list(projected.car.bin.pct, pca.dt))
  
}