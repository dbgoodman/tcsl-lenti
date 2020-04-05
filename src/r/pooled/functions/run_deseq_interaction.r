run_deseq_interaction <- function(data.dt, ref_bin, test_bin){
  
  # identify inputs
  assay_input <- as.character(unique(data.dt[, assay]))
  k_type <- as.character(unique(data.dt[, k.type]))
  t_type <- as.character(unique(data.dt[, t.type]))
  
  # prepare cts and coldata dataframes
  if(ref_bin == 'baseline'){
    ref.bin.dt <- dcast(unique(data.dt[bin == 'D' & assay == assay_input & 
                                         k.type == k_type & t.type == t_type,
                                       .(CAR.align, bin.sort.group = 
                                           paste(batch, donor, timepoint, assay, 
                                                 t.type, 'base', sep = '_'), 
                                         k.type, t.type, batch, assay, donor, 
                                         sort.group, bin = 'base', counts = baseline.counts)]), 
                        CAR.align ~ bin.sort.group, value.var='counts')
    
    num.ref.reps <- ncol(ref.bin.dt) - 1
    
    ref.bin.dt <- cbind(ref.bin.dt[, 1], 
                        do.call("cbind", replicate(length(test_bin), 
                                                   ref.bin.dt[, -1], simplify = FALSE)))
    
    names(ref.bin.dt) <- c(names(ref.bin.dt[, 1]), paste(names(ref.bin.dt[, -1]), 
                                                         rep(test_bin, each=num.ref.reps), 
                                                         sep = '_'))
    
  }else{
    ref.bin.dt <- dcast(unique(data.dt[bin %in% ref_bin & assay == assay_input & 
                                         k.type == k_type & t.type == t_type,
                                       .(CAR.align, bin.sort.group = 
                                           paste(sort.group, bin, sep = '_'), 
                                         k.type, t.type, batch, assay, donor, 
                                         sort.group, bin, counts)]), 
                        CAR.align ~ bin.sort.group, value.var='counts')
  }
  
  test.bin.dt <- dcast(unique(data.dt[assay == assay_input & bin %in% test_bin &
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
                        bin = sapply(strsplit(c(names(ref.bin.dt)[-1], 
                                                names(test.bin.dt)[-1]),"_"), `[`, 7),
                        row.names = c(names(ref.bin.dt)[-1], names(test.bin.dt)[-1]))
  
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ bin+condition)
  
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