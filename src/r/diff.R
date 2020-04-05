require(data.table)
require(scales)
require(here)
require(flowCore)
require(flowWorkspace)

parent_folder = here::here('..','..')

load_platemaps <- function(diff_dir) {
  
  # load platemap filenames
  platemap_paths <- data.table(
    filepath= normalizePath(file.path(diff_dir, list.files(
      diff_dir,
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

load_fcs <- function(diff_dir, channel_map) {
  diff_fcs_dir <- file.path(diff_dir)
  
  # load fcs filenames
  fcs_paths <- data.table(
    filepath= normalizePath(file.path(diff_fcs_dir, list.files(
      diff_fcs_dir,
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
    cd="FJComp-450_50 Violet-A",
    ccr7="FJComp-780_60 YG-A",
    gfp="FJComp-515_20 Blue-A",
    zombie="FJComp-586_15 Violet-A",
    myc="FJComp-610_20 YG-A",
    cd45ra= "FJComp-660_20 Red-A",
    cd62l= "FJComp-800_30 Violet-A",
    cd95="FJComp-710_50 Violet-A",
    cd27="FJComp-780_60 Red-A",
    cd45ro="FJComp-379_28 UV-A")
  
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

load_diff_data <- function(
  experiment_dir, marker_pos_threshold, max_sample=Inf) 
{
  
  diff_dir <- file.path(parent_folder, experiment_dir)
  
  platemaps <- load_platemaps(diff_dir)
  
  fcs_out <- load_fcs(diff_dir)
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
  
  melt_cols <- c(names(channel_map), 
                 'SSC-A','SSC-W','SSC-H','FSC-A','FSC-H','FSC-W')
  transformed_vars <- paste0(melt_cols,'_t')
  melt_cols <- c(melt_cols, transformed_vars)
  
  # downsample per well if specified
  fcs_dt <- fcs_dt[, .SD[sample(.N, min(c(.N, max_sample)))], 
                   by=c('well','plate','day')]
  
  fcs_melt_dt <- melt(
    fcs_dt,
    measure.vars=melt_cols)
  
  marker_thresh_dt <- data.table(
    variable= paste0(names(marker_pos_threshold)),
    threshold= unlist(marker_pos_threshold))
  
  fcs_melt_dt <- fcs_melt_dt[marker_thresh_dt, on='variable'][, 
    gt_thresh := value > threshold]
  
  fcs_melt_dt[, threshold := NULL]
  
  return(list(
    fcs_melt_dt= fcs_melt_dt, channel_map=channel_map, melt_cols=melt_cols))
}


# CD4 ======

n_events_per_well <- 10000

cd4_marker_pos_threshold <- list(
  cd62l_t=1300,
  cd27_t=1800,
  ccr7_t=1100,
  cd45ro_t=800,
  cd45ra_t=900,
  cd95_t=2100,
  gfp_t=NA,
  cd_t=NA,
  zombie_t=NA,
  myc_t=NA,
  FSC_t=NA,
  SSC_t=NA)

#Box Version:
#cd4_dir <- "flow/TCSL105 ALL FILES/Diff/processed"
#AWS Version:
cd4_dir <- "s3-roybal-tcsl/lenti_screen_compiled_data/data/fcs/tcsl105/diff"

# cd4_out <- load_diff_data(cd4_dir, 
#   cd4_marker_pos_threshold, 
#   max_sample=n_events_per_well)

# CD8 ======

cd8_marker_pos_threshold <- list(
  cd62l_t=1750,
  cd27_t=1800,
  ccr7_t=1050,
  cd45ro_t=800,
  cd45ra_t=900,
  cd95_t=2100,
  gfp_t=NA,
  cd_t=NA,
  zombie_t=NA,
  myc_t=NA,
  FSC_t=NA,
  SSC_t=NA)

#Box Version:
#cd8_dir <- "flow/2019.06.07 TCSL091 ALL FILES/Diff/processed"
#AWS Version:
cd8_dir <- "s3-roybal-tcsl/lenti_screen_compiled_data/data/fcs/tcsl091/diff"

# cd8_out <- load_diff_data(cd8_dir,
#   cd8_marker_pos_threshold, 
#   max_sample=n_events_per_well)


diff_dt <- rbind(
  cd8_out$fcs_melt_dt[, t_type := 'cd8'], 
  cd4_out$fcs_melt_dt[, t_type := 'cd4'])

cd4_out$fcs_melt_dt <- NULL
cd8_out$fcs_melt_dt <- NULL

diff_opt_cd4 <- cd4_out
diff_opt_cd8 <- cd8_out

# save(
#   'diff_opt_cd4', 'diff_opt_cd8',
#   file=file.path(here::here('..','data','diff.Rdata')))

# fwrite(diff_dt, 
#   compress='gzip', 
#   file=file.path(here::here('..','data','diff.csv.gz')))

# ------
# make a downsampled version for faster plotting and analysis

n_events_per_well <- 2500
cd8_out <- load_diff_data(cd8_dir, 
    cd8_marker_pos_threshold, max_sample=n_events_per_well)
cd4_out <- load_diff_data(cd4_dir, 
    cd4_marker_pos_threshold, max_sample=n_events_per_well)

diff_dt <- rbind(
  cd8_out$fcs_melt_dt[, t_type := 'cd8'], 
  cd4_out$fcs_melt_dt[, t_type := 'cd4'])

save('diff_opt_cd4', 'diff_opt_cd8',
     file=file.path(here::here('..','data','diff.sampled.Rdata')))

fwrite(diff_dt, 
       compress='gzip', 
       file=file.path(here::here('..','data','diff.sampled.csv.gz')))

# ------

# sync output back to s3
system('aws s3 sync \\
  --exclude "*" \\
  --include "*.Rdata" \\
  --include "*.csv*" \\
  ~/local-tcsl/data s3://roybal-tcsl//lenti_screen_compiled_data/data')

