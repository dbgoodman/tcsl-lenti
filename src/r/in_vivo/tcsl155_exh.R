library(data.table)
library(scales)
library(here)
library(flowCore)
library(flowWorkspace)

parent_folder = here::here('..','..')

load_platemaps <- function(exh_dir) {
  
  # load platemap filenames
  platemap_paths <- data.table(
    filepath= normalizePath(file.path(exh_dir, '..', '..', list.files(
      file.path(exh_dir, '..', '..'),
      pattern='.*platemap.csv$',
      recursive=TRUE))))
  

  # load and merge platemaps
  platemaps <- platemap_paths[, fread(file=.BY[[1]]), by=filepath]
  
  # put 0 back into well number (ugh)
  platemaps[, well := gsub(' ', '',
                           gsub('([A-H])(\\d{1})$', '\\1 0\\2', well))]
  
  return(platemaps)
}

load_fcs <- function(exh_dir) {
  exh_fcs_dir <- file.path(exh_dir)
  
  # load fcs filenames
  fcs_paths <- data.table(
    filepath= normalizePath(file.path(exh_fcs_dir, list.files(
      exh_fcs_dir,
      pattern='.*fcs',
      recursive=TRUE))))
  
  # parse well name
  fcs_paths[, well := gsub('.*_([A-H]\\d+).fcs','\\1', filepath)]
 
  fcs_dt <- fcs_paths[, 
                      as.data.frame(flowCore::exprs(read.FCS(filepath))),
                      by=c('well')]
  
  # give each event an ID for tracking after melt
  fcs_dt[, event_id := .I, by=c('well')]
  
  channel_map <- list(
    cd4="FJComp-660_20 Violet-A",
    cd8="FJComp-379_28 UV-A",
    tim3="FJComp-586_15 YG-A",
    gfp="FJComp-515_20 Blue-A",
    lag3= "FJComp-800_30 Violet-A",
    pd1= "FJComp-660_20 Red-A",
    cd39="FJComp-710_50 Violet-A")
  
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
  
  return(list(fcs_dt= fcs_dt, channel_map= channel_map))
}


load_exh_data <- function(experiment_dir, marker_pos_threshold, max_sample=Inf) {
  
  exh_dir <- file.path(parent_folder, experiment_dir)
  
  platemaps <- load_platemaps(exh_dir)
  
  fcs_out <- load_fcs(exh_dir)
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
  
  # downsample per well if specified
  fcs_dt <- fcs_dt[, .SD[sample(.N, min(c(.N, max_sample)))], by=c('well')]
  
  melt_cols <- c(names(channel_map), 'SSC-A','SSC-W','SSC-H','FSC-A','FSC-H','FSC-W')
  transformed_vars <- paste0(melt_cols,'_t')
  melt_cols <- c(melt_cols, transformed_vars)
  
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
  tim3_t=1000,
  gfp_t=1300,
  lag3_t=700,
  pd1_t=1400,
  cd39_t=1000)

#cd4_dir <- "flow/TCSL105 ALL FILES/Exhaust"
#AWS Version:
cd4_dir <- "s3-roybal-tcsl/lenti_screen_compiled_data/data/fcs/tcsl155/exh/cd4"

#cd4_out <- load_exh_data(cd4_dir, cd4_marker_pos_threshold)

# CD8 ======

cd8_marker_pos_threshold <- list(
  tim3_t=1000,
  gfp_t=1300,
  lag3_t=700,
  pd1_t=1000,
  cd39_t=1000)

#cd8_dir <- "flow/2019.06.07 TCSL091 ALL FILES/Exhaust"
#AWS Version:
cd8_dir <- "s3-roybal-tcsl/lenti_screen_compiled_data/data/fcs/tcsl155/exh/cd8"

#cd8_out <- load_exh_data(cd8_dir, cd8_marker_pos_threshold)

# exh_dt <- rbind(
#   cd8_out$fcs_melt_dt[, t_type := 'cd8'],
#   cd4_out$fcs_melt_dt[, t_type := 'cd4'])

# cd4_out$fcs_melt_dt <- NULL
# cd8_out$fcs_melt_dt <- NULL
# 
# exh_opt_cd4 <- cd4_out
# exh_opt_cd8 <- cd8_out
# 
# exh_opt_cd4$marker_thresholds <- cd4_marker_pos_threshold
# exh_opt_cd8$marker_thresholds <- cd8_marker_pos_threshold


#save('exh_opt_cd4', 'exh_opt_cd8', 'exh_dt', file=file.path(here::here('..','data','exh.Rdata')))


# ------
# make a downsampled version for faster plotting and analysis

n_events_per_well <- 2500
cd8_out <- load_exh_data(cd8_dir, cd8_marker_pos_threshold, max_sample=n_events_per_well)
cd4_out <- load_exh_data(cd4_dir, cd4_marker_pos_threshold, max_sample=n_events_per_well)
exh_dt <- rbind(cd8_out$fcs_melt_dt[, t_type := 'cd8'], cd4_out$fcs_melt_dt[, t_type := 'cd4'])

cd4_out$fcs_melt_dt <- NULL
cd8_out$fcs_melt_dt <- NULL

exh_opt_cd4 <- cd4_out
exh_opt_cd8 <- cd8_out

exh_opt_cd4$marker_thresholds <- cd4_marker_pos_threshold
exh_opt_cd8$marker_thresholds <- cd8_marker_pos_threshold

fwrite(exh_dt, 
       compress='gzip', 
       file=file.path(here::here('..','data','tcsl155_exh.sampled.csv.gz')))

save(
  'exh_opt_cd4',
  'exh_opt_cd8',
  file=file.path(here::here('..','data','tcsl155_exh.sampled.Rdata')))


