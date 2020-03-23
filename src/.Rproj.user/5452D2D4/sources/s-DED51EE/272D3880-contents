require(data.table)
require(scales)
require(here)
require(flowCore)
require(flowWorkspace)

parent_folder = here::here('..','..')

load_platemaps <- function(ctv_dir) {

  # load platemap filenames
  platemap_paths <- data.table(
    filepath= normalizePath(file.path(ctv_dir, list.files(
      ctv_dir,
      pattern='.*platemap.csv$',
      recursive=TRUE))))
  
  #parse day from file name
  platemap_paths[, day := gsub('.*D(\\d+).*', '\\1', filepath)]
  
  # load and merge platemaps
  platemaps <- platemap_paths[, fread(file=.BY[[1]]), by=filepath]
  
  # put 0 back into well number (ugh)
  platemaps[, well := gsub(' ', '',
                           gsub('([A-H])(\\d{1})$', '\\1 0\\2', well))]
  
  platemaps[car == 'KRLG1', car := 'KLRG1'][]
  
  return(platemaps)
}

load_fcs <- function(ctv_dir) {
  ctv_fcs_dir <- file.path(ctv_dir,'processed')
  
  # load fcs filenames
  fcs_paths <- data.table(
    filepath= normalizePath(file.path(ctv_fcs_dir, list.files(
      ctv_fcs_dir,
      pattern='.*fcs',
      recursive=TRUE))))
  
  # parse well name
  fcs_paths[, well := gsub('.*_([A-H]\\d+).fcs','\\1', filepath)]
  fcs_paths[, plate := gsub('.*(\\d+)_[A-H]\\d+.fcs','\\1', filepath)]
  fcs_paths[, day := gsub('.*_(\\d+)_\\d+_[A-H]\\d+.fcs','\\1', filepath)]
  
  fcs_dt <- fcs_paths[, 
                      as.data.frame(flowCore::exprs(read.FCS(filepath))),
                      by=c('well','plate','day')]
  
  # give each event an ID for tracking after melt
  fcs_dt[, event_id := .I, by=c('well','plate','day')]
  
  channel_map <- list(
    cd="379_28 UV-A",
    ctv1="450_50 Violet-A",
    ctv2="470_15 Violet-A",
    gfp="515_20 Blue-A",
    myc= "610_20 YG-A",
    draq7= "730_45 Red-A")
  
  # replace column names
  chan_names <- data.table(merge(
    names(fcs_dt),
    data.table(
      marker=names(channel_map),
      chan=unlist(channel_map)), 
    by.x='x', by.y='chan',
    all=T))
  
  # fill in original column names
  chan_names[is.na(marker), marker := x]
  
  # rename columns and melt
  names(fcs_dt) <- chan_names[
    names(fcs_dt), marker, on='x']
  
  # merge unmelted with plate map
  fcs_dt[, day := as.numeric(day)]
  fcs_dt[, plate := as.numeric(plate)]
  return(list(fcs_dt= fcs_dt, channel_map= channel_map))
}

load_ctv_data <- function(experiment_dir, marker_pos_threshold, max_sample=Inf) {
  
  ctv_dir <- file.path(parent_folder, experiment_dir)
  
  platemaps <- load_platemaps(ctv_dir)
  
  fcs_out <- load_fcs(ctv_dir)
  fcs_dt <- fcs_out$fcs_dt
  channel_map <- fcs_out$channel_map
  rm(fcs_out)
  
  fcs_dt <- merge(
    fcs_dt, platemaps,
    on=intersect(names(platemaps), names(fcs_dt)), all=TRUE)
  
  # flowJo transform 
  for (chan_name in c(names(channel_map), 
        'SSC-A','SSC-W','SSC-H','FSC-A','FSC-H','FSC-W')) { 
    fcs_dt[, c(paste0(chan_name,'_t')) := flowjo_biexp()(fcs_dt[[chan_name]])] 
  }
  
  melt_cols <- c(names(channel_map), 'SSC-A','SSC-W','SSC-H','FSC-A','FSC-H','FSC-W')
  transformed_vars <- paste0(melt_cols,'_t')
  melt_cols <- c(melt_cols, transformed_vars)
  
  # downsample per well if specified
  fcs_dt <- fcs_dt[, .SD[sample(.N, min(c(.N, max_sample)))], by=c('well','plate','day')]
  
  fcs_melt_dt <- melt(
    fcs_dt,
    measure.vars=melt_cols)
  
  marker_thresh_dt <- data.table(
    variable= paste0(names(marker_pos_threshold)),
    threshold= unlist(marker_pos_threshold))
  
  fcs_melt_dt <- fcs_melt_dt[marker_thresh_dt, on='variable'][, gt_thresh := value > threshold]
  fcs_melt_dt[, threshold := NULL]
  
  return(list(fcs_melt_dt= fcs_melt_dt, channel_map=channel_map, melt_cols=melt_cols))
}

# CD4 ======

cd4_marker_pos_threshold <- list(
  cd_t= 2200,
  ctv1_t=1700,
  ctv2_t=800,
  gfp_t=1250,
  myc_t= 400, 
  draq7_t= 1100)

cd4_dir <- "flow/TCSL105 ALL FILES/CTV"

cd4_out <- load_ctv_data(cd4_dir, cd4_marker_pos_threshold)
  
# CD8 ======

cd8_marker_pos_threshold <- list(
  cd_t= 2200,
  ctv1_t=1700,
  ctv2_t=800,
  gfp_t=1250,
  myc_t= 400, 
  draq7_t= 1100)

cd8_dir <- "flow/2019.06.07 TCSL091 ALL FILES/CTV"

cd8_out <- load_ctv_data(cd8_dir, cd8_marker_pos_threshold)

ctv_dt <- rbind(cd8_out$fcs_melt_dt[, t_type := 'cd8'], cd4_out$fcs_melt_dt[, t_type := 'cd4'])
cd4_out$fcs_melt_dt <- NULL
cd8_out$fcs_melt_dt <- NULL

ctv_opt_cd4 <- cd4_out
ctv_opt_cd8 <- cd8_out

save('ctv_opt_cd4', 'ctv_opt_cd8', 'ctv_dt', file=file.path(here::here('..','data','ctv.Rdata')))

# ------
# make a downsampled version for faster plotting and analysis

n_events_per_well <- 15000
cd8_out <- load_ctv_data(cd8_dir, cd8_marker_pos_threshold, max_sample=n_events_per_well)
cd4_out <- load_ctv_data(cd4_dir, cd4_marker_pos_threshold, max_sample=n_events_per_well)
ctv_dt <- rbind(cd8_out$fcs_melt_dt[, t_type := 'cd8'], cd4_out$fcs_melt_dt[, t_type := 'cd4'])
save('ctv_opt_cd4', 'ctv_opt_cd8', 'ctv_dt', file=file.path(here::here('..','data','ctv.sampled.Rdata')))

