# testing deseq2 from pooled_03_deseq2

#groupname = 'CTV2_CD8_pos'
groupname <- c('CTV1_CD8_pos', 'CTV1_CD4_pos')
set_i = 1

#outer loop setup and params
#--------

test_ref_sets <- c(ref=list(), test=list())

# A/B/AB vs C/D/CD
test_ref_sets$ref <- c(as.list(rep('D',3)), as.list(rep('C',3)), rep(list(c('C','D')),3))
test_ref_sets$test <- rep(list('A','B',c('A','B')),3)
test_ref_sets$interaction <- as.list(rep(F, 9))

# A/B/AB/ABCD vs baseline
test_ref_sets$ref <- c(test_ref_sets$ref, as.list(rep('baseline',3)))
test_ref_sets$test <- c(test_ref_sets$test, list('A',c('A','B'),c('A','B','C','D')))
test_ref_sets$interaction <- c(test_ref_sets$interaction, as.list(rep(T, 3)))

all.deseq.results.dt <- data.table()
    
ref_set <- test_ref_sets$ref[[set_i]]
test_set <- test_ref_sets$test[[set_i]]

ref_str <- paste0(ref_set, collapse='')
test_str <- paste0(test_set, collapse='')

inter <- test_ref_sets$interaction[[set_i]]

#run_deseq args
#--------

data.dt <- read.counts[group %in% groupname]
control_replicates = T
interaction = inter
group.control = T
weight.bins = F
ref_bin <- ref_set
test_bin <- test_set


#run_deseq fxn
#--------
run_deseq <- function(data.dt, ref_bin, test_bin, 
    control_replicates = T,
    interaction = T, group.control = F, weight.bins = F) 
{
  # identify inputs
  assay_input <- as.character(unique(data.dt[, assay]))
  k_type <- as.character(unique(data.dt[, k.type]))
  t_type <- as.character(unique(data.dt[, t.type]))
  
  ## 1. Do bin normalization weights =======
  # bin normalization weights
  data.weights <- unique(data.dt[, list(
    batch, donor, timepoint, assay, t.type, k.type, 
    sort.group, bin, bin.pct, bin.reads)])[, 
      list(bin, bin.reads, bin.pct, sort.group,
        read.weight=bin.pct * bin.reads / sum(bin.pct * bin.reads)),
      by=.(batch, donor, timepoint, assay, t.type, k.type)][, 
        read.weight.norm := read.weight/exp(mean(log(read.weight))), 
        by=.(batch, donor, timepoint, assay, t.type, k.type)]
  
  #stopifnot(nrow(interaction(assay_input, k_type, t_type)) == 1)
  
  # prepare cts and coldata dataframes
  
  ## 2. Prepare Ref Bins ============
  if(length(ref_bin) == 1 & ref_bin[1] == 'baseline') {
    # reference is baseline
    # get baseline counts per donor/assay replicate
    ref.bin.dt <- dcast(
      data.dt[
        bin == 'D', 
        .(
          CAR.align, 
          bin.sort.group = paste(
            batch, donor, timepoint, assay, t.type, 'base', sep = '_'), 
          k.type, t.type, batch, assay, donor, 
          sort.group, 
          bin = 'base',
          counts = baseline.counts)], 
      CAR.align ~ bin.sort.group, value.var='counts')
    
    ref.weights <- rep(1, ncol(ref.bin.dt))
    
  } else {
    # reference is a specified bin
    ref.bin.dt <- dcast(
      data.dt[
        bin %in% ref_bin, 
        .(
          CAR.align, 
          bin.sort.group = paste(sort.group, bin, sep = '_'), 
          k.type, t.type, batch, assay, donor, 
          sort.group, 
          bin,
          counts)], 
      CAR.align ~ bin.sort.group, value.var='counts')
    
    ref.weights <- dcast(data.weights[
        bin %in% ref_bin,
        .(
          bin.sort.group = paste(sort.group, bin, sep = '_'), 
          k.type, t.type, batch, assay, donor, 
          sort.group,
          bin,
          read.weight.norm)],
        . ~ bin.sort.group, 
        value.var = 'read.weight.norm')
    
    stopifnot(nrow(ref.bin.dt) == nrow(unique(ref.bin.dt)))
  }
  
  # copy the ref bin columns for each of the test bin columns
  if (interaction == T) {
    
    num.ref.reps <- ncol(ref.bin.dt) - 1
    
    ref.bin.dt <- cbind(ref.bin.dt[, 1], 
      do.call("cbind", replicate(length(test_bin), 
      ref.bin.dt[, -1], simplify = FALSE)))
    
    names(ref.bin.dt) <- c(names(ref.bin.dt[, 1]), 
      paste(names(ref.bin.dt[, -1]), rep(test_bin, each=num.ref.reps), 
        sep = '_'))
  }
  
  ## 3. Prepare Test Bins ============
  test.bin.dt <- dcast(
      data.dt[
        bin %in% test_bin,
        .(
          CAR.align, 
          bin.sort.group = paste(sort.group, bin, sep = '_'), 
          k.type, t.type, batch, assay, donor, 
          sort.group, 
          bin,
          counts)], 
      CAR.align ~ bin.sort.group, value.var='counts')
  
  # check that replicate counts match
  stopifnot(nrow(ref.bin.dt) == nrow(unique(ref.bin.dt)))
  stopifnot(nrow(test.bin.dt) == nrow(unique(test.bin.dt)))
  
  ## 4. Merge and create design matrix ============
  cts <- merge(ref.bin.dt, test.bin.dt, by = 'CAR.align')
  cts <- data.frame(cts[, -1], row.names = cts[, CAR.align])
  cts[is.na(cts)] <- 0
  
  coldata <- data.frame(
    condition = c(
      rep('reference', ncol(ref.bin.dt) - 1), 
      rep('test', ncol(test.bin.dt) - 1)),
    rep = data.table(t(sapply(strsplit(c(
      names(ref.bin.dt)[-1], 
      names(test.bin.dt)[-1]),"_"), `[`, c(1,2))))[,
        paste(V1, V2, sep='_')],
    bin = sapply(strsplit(c(
      names(ref.bin.dt)[-1], 
      names(test.bin.dt)[-1]),"_"), `[`, 7),
    t.type = gsub('.*(CD[48]).*','\\1',names(ref.bin.dt)[-1]),
    row.names = c(names(ref.bin.dt)[-1], names(test.bin.dt)[-1]))
  
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design =  ~ condition + rep + t.type + bin + condition:t.type)
  
  # set reference
  dds$condition <- relevel(dds$condition, ref = "reference")
  
  print(coldata)
  
  # pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  dds <- DESeq2::DESeq(dds)
  
  # shrink log fold change & convert to data.table, merge contrasts
  results.dt <- rbindlist(
      lapply(resultsNames(dds)[-1], function(result_i) 
          as.data.table(lfcShrink(dds, coef=result_i, type="apeglm"), 
        keep.rownames='CAR.align')[, coef := result_i]))
  
  results.dt[, assay := assay_input][, k.type := k_type]
  
  return(results.dt)
}

# run all sets

for (set_i in seq_along(test_ref_sets$ref)) {
    
  ref_set <- test_ref_sets$ref[[set_i]]
  test_set <- test_ref_sets$test[[set_i]]
  
  ref_str <- paste0(ref_set, collapse='')
  test_str <- paste0(test_set, collapse='')
  
  inter <- test_ref_sets$interaction[[set_i]]
  
  deseq.results.dt <- read.counts[batch != 'post-cytof' & k.type == 'pos',
    {
      message(paste(c(ref_str, test_str, inter, .BY[1],"\n"), collapse= ' - '));
      tryCatch(
        run_deseq(
          data.dt = data.table(.SD)[, assay := .BY[1]], 
          ref_bin = ref_set,
          test_bin = test_set,
          interaction = inter,
          group.control = T),
        error= function(e) {message(e); return(data.table())}
      )
    },
    by = .(assay)]
  
  if (nrow(deseq.results.dt) > 0) {
    deseq.results.dt[,
      `:=`(
        ref_set = ref_str,
        test_set = test_str,
        inter = inter)]
  }
  
  all.deseq.results.dt <- rbind(
    all.deseq.results.dt,
    deseq.results.dt, fill=T)
}

save(list=c('all.deseq.results.dt'), 
   file=file.path(data.output.dir, 'pooled_deseq2_data_new.Rdata'))
