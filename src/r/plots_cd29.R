
require(data.table)
require(scales)
require(here)
require(flowCore)
require(flowWorkspace)

source(here::here('r','fig_colors.R'))


platemap <- fread(
  text='well	car	k_type
D02	41BB	pos
D03	BAFFR	pos
D04	CD28	pos
D05	CD40	pos
D06	KLRG1	pos
D07	TACI	pos
D08	Zeta	pos
D09	Untr	pos
B02	41BB	neg
B03	BAFFR	neg
B04	CD28	neg
B05	CD40	neg
B06	KLRG1	neg
B07	TACI	neg
B08	Zeta	neg
B09	Untr	neg
B10	Untr	bead')

cd29_data_dir <- "/home/ec2-user/s3-roybal-tcsl/lenti_screen_compiled_data/data/fcs/tcsl170"

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
  fcs_paths[, plate := gsub('.*_(\\d+)_[A-H]\\d+.fcs','\\1', filepath)]
  fcs_paths[, t_type := gsub('.*/(cd[48])_.*.fcs','\\1', filepath)]
  
  fcs_dt <- fcs_paths[, 
                      as.data.frame(flowCore::exprs(read.FCS(filepath))),
                      by=c('well','plate','t_type')]
  
  # give each event an ID for tracking after melt
  fcs_dt[, event_id := .I, by=c('well','plate')]
  
  channel_map <- list(
    ifng= "FJComp-800_30 Violet-A",
    cd4= "FJComp-450_50 Violet-A",
    cd8= "FJComp-586_15 YG-A",
    gfp= "FJComp-515_20 Blue-A",
    zombie= "FJComp-586_15 Violet-A",
    il2="FJComp-740_35 UV-A",
    cd29="FJComp-660_20 Red-A",
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
  fcs_dt[, plate := as.numeric(plate)]
  return(list(fcs_dt= fcs_dt, channel_map= channel_map))
}

load_diff_data <- function(
  experiment_dir, marker_pos_threshold, platemaps, max_sample=Inf) 
{
  
  diff_dir <- file.path(experiment_dir)
  
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
                   by=c('well','plate')]
  
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

marker_pos_threshold <- list(
  cd45ro_t=800,
  zombie_t=NA,
  gfp_t=NA,
  il2_t=750,
  cd29_t=2800,
  ifng_t=1200)

fcs_out <- load_diff_data(
  cd29_data_dir, 
  marker_pos_threshold,
  platemap,
  max_sample=Inf)

ggplot(fcs_out$fcs_melt_dt[plate == 0 & variable %in% c('cd29_t','ifng_t','il2_t','cd45ro_t')]) + geom_density(aes(x=value, color=factor(k_type, levels=c('neg','pos','bead')))) + facet_grid(car ~ t_type + variable)
ggplot(fcs_out$fcs_melt_dt[k_type != 'bead' & plate == 0 & variable %in% c('cd29_t','ifng_t','il2_t','cd45ro_t')]) + geom_boxplot(aes(x=car, y=value, color=k_type), position=position_dodge(width=1), outlier.shape=NA) + facet_grid(variable~t_type+k_type, scales='free')
ggplot(fcs_out$fcs_melt_dt[gt_thresh==T & plate == 0 & variable %in% c('ifng_t','il2_t')]) + geom_bar(aes(y = (..count..)/sum(..count..), x=car, fill=factor(k_type, levels=c('neg','pos','bead'))), position='dodge') + facet_grid( ~ t_type + variable) + scale_fill_discrete('stim')
ggplot(fcs_out$fcs_melt_dt[k_type != 'bead' & plate == 0 & variable %in% c('cd29_t')]) + geom_boxplot(aes(x=car, y=value, color=k_type), position=position_dodge(width=1), outlier.shape=NA) + facet_grid(variable~k_type+t_type, scales='free') + theme_bw() + labs(y='CD29 MFI', title='CD29 levels')
ggplot(fcs_out$fcs_melt_dt[k_type != 'bead' & plate == 0 & variable %in% c('ifng_t')]) + geom_boxplot(aes(x=car, y=value, color=k_type), position=position_dodge(width=1), outlier.shape=NA) + facet_grid(variable~k_type+t_type, scales='free') + theme_bw() + labs(y='IFNg MFI', title='IFNg levels')


# car name swaps:
fcs_out$fcs_melt_dt[car == 'BAFFR', car := 'BAFF-R']
fcs_out$fcs_melt_dt[car == 'Zeta', car := 'zeta']
fcs_out$fcs_melt_dt[,
    car := factor(car, levels=c(car_order, 'Untr'))]

fcs_pct <- fcs_out$fcs_melt_dt[, list(N=sum(gt_thresh), pct=sum(gt_thresh)/.N), by=c('t_type','car','k_type','variable','plate','well')]

# good plot, IFNG+ % across cars
ggplot(
  fcs_pct[variable %in% c('ifng_t') & k_type == 'pos' & t_type == 'cd8' & plate == 0]) +
  geom_bar(aes(y = pct, x=car, fill=car), stat='identity') +
  scale_y_continuous('% of IFNG positive cells', label=label_percent()) +
  scale_fill_manual(values=c(car_colors, 'Untr'='grey50')) +
  theme_minimal() + facet_grid(k_type ~ t_type + plate)

ggplot(
  fcs_pct[variable %in% c('cd29_t') & plate == 0 & k_type %in% c('neg','pos')]) +
  geom_bar(aes(y = pct, x=car, fill=car), stat='identity') +
  scale_y_continuous('% of CD29 positive cells', label=label_percent()) +
  scale_fill_manual(values=c(car_colors, 'Untr'='grey50')) +
  theme_minimal() + facet_grid(k_type ~ t_type + plate)

paired_pct <- fcs_out$fcs_melt_dt[, 
  list(
    pos_gamma = gt_thresh[variable=='ifng_t'],
    pos_29 = gt_thresh[variable=='cd29_t']), 
  by=c('event_id','t_type','plate','car','well','k_type')][,
    data.table(unlist(table(pos_gamma, pos_29)), keep.rownames=T),
    by=c('t_type','plate','car','well','k_type')][,
      pct := N/sum(N), by=c('t_type','plate','car','well','k_type')]

ggplot(paired_pct[pos_gamma == T]) +
  geom_point(aes(y = pct, x=car, color=pos_29, group=pos_29), position=position_dodge(1)) +
  scale_y_continuous('% of IFNG+ CD29+ cells', label=label_percent()) +
  scale_fill_manual(values=c(car_colors, 'Untr'='grey50')) +
  facet_wrap(k_type ~ t_type) +
  theme_minimal()

cd29_pct <- ggplot(
  fcs_pct[
    variable %in% c('cd29_t') & 
    plate == 0 &
    k_type %in% c('neg','pos') &
    car %in% c('CD28','41BB','BAFF-R','zeta','Untr')]) +
  geom_bar(aes(y = pct, x=car, fill=car), stat='identity') +
  scale_y_continuous(
    '% of CD29 positive cells', label=label_percent(), expand=c(0,0.01),
    limits=c(0,0.65)) +
  scale_fill_manual(values=c(car_colors, 'Untr'='grey50'), guide=NULL) +
  theme_minimal() +
  theme(
      panel.grid=element_blank(),
      panel.background=element_rect(fill=NA))

ifng_29_pct <- ggplot(paired_pct[
    pos_gamma == T & plate == 0 &
    t_type == 'cd8' & k_type == 'pos' &
    car %in% c('CD28','41BB','BAFF-R','zeta','Untr')]) +
  geom_bar(
    aes(y = pct, x=car, fill=car, color=car, alpha=pos_29, group=pos_29),
    stat='identity', position=position_dodge(0.8), width=0.8) +
  scale_y_continuous(
    '% of IFNG+ cells', label=label_percent(), expand=c(0,0.01),
    limits=c(0,0.24)) +
  scale_alpha_manual('CD29+', values=c(0,1)) +
  scale_fill_manual('', values=c(car_colors, 'Untr'='grey70'), guide=NULL) +
  scale_color_manual('', values=c(car_colors, 'Untr'='grey70'), guide=NULL) +
  guides(alpha=guide_legend(override.aes =list(color='black'))) +
  theme_minimal() +
  theme(
      panel.grid=element_blank(),
      panel.background=element_rect(fill=NA))

(cd29_pct + labs(title='% CD29+', x='')) | (ifng_29_pct  + labs(title='% IFNG+', x='')) 
