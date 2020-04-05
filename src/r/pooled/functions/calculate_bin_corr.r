calculate_bin_corr <- function(data.dt, group){
  
  bin.corr <- cor(dcast(data.dt[sort.group == group], 
                        CAR.align + sort.group ~ bin,
                        value.var='rel.bin.ratio')[
                          is.na(A), A := 0][is.na(B), B := 0][
                            is.na(C), C := 0][is.na(D), D := 0][, -c(1,2)])
  
  melt.bin.corr <- cbind(data.table(melt(bin.corr)), 
                         unique(data.dt[sort.group == group, 
                                        .(sort.group, batch, assay, donor, 
                                          t.type, k.type)]))
  
  return(melt.bin.corr)
  
}