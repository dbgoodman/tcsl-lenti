run_deseq <- function(data.dt, ref_bin, test_bin){
  
  # identify inputs
  assay_input <- as.character(unique(data.dt[, assay]))
  k_type <- as.character(unique(data.dt[, k.type]))
  t_type <- as.character(unique(data.dt[, t.type]))
  
  print(paste(assay_input, k_type, t_type, sep = '_'))
  
  # prepare cts and coldata dataframes
  ref.bin.dt <- dcast(unique(data.dt[bin == ref_bin & assay == assay_input & 
                                       k.type == k_type & t.type == t_type,
                                     .(CAR.align, bin.sort.group = 
                                         paste(sort.group, bin, sep = '_'), 
                                       k.type, t.type, batch, assay, donor, 
                                       sort.group, bin, counts)]), 
                      CAR.align ~ bin.sort.group, value.var='counts')
  
  test.bin.dt <- dcast(unique(data.dt[bin == test_bin & assay == assay_input & 
                                        k.type == k_type & t.type == t_type,
                                      .(CAR.align, bin.sort.group = 
                                          paste(sort.group, bin, sep = '_'), 
                                        k.type, t.type, batch, assay, donor, 
                                        sort.group, bin, counts)]), 
                       CAR.align ~ bin.sort.group, value.var='counts')
  
  cts <- merge(ref.bin.dt, test.bin.dt, by = 'CAR.align')
  cts <- data.frame(cts[, -1], row.names = cts[, CAR.align])
  cts[is.na(cts)] <- 0
  
  coldata <- data.frame(condition = c(rep('reference', ncol(ref.bin.dt) - 1), 
                                      rep('test', ncol(test.bin.dt) - 1)), 
                        row.names = c(names(ref.bin.dt)[-1], names(test.bin.dt)[-1]))
  
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ condition)
  
  # pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  # set reference
  dds$condition <- relevel(dds$condition, ref = "reference")
  
  # run DESeq
  dds <- DESeq(dds)
  res <- results(dds)
  
  # shrink log fold change
  resLFC <- lfcShrink(dds, coef="condition_test_vs_reference", type="apeglm")
  
  # convert to data.table
  results.dt <- as.data.table(resLFC)[, CAR.align := row.names(resLFC)]
  results.dt <- cbind(results.dt[, 6], results.dt[, -6])
  results.dt[, assay := assay_input][, k.type := k_type][, t.type := t_type]
  
  return(results.dt)
  
}