library("data.table")
library("gdata")
library("ggplot2")
library('scales')

incucyte_dir <- here::here('..','..','incucyte')

load_single_measurement <- function(parent_dir, data_path, platemap_path) {
  
  data_dt <- fread(file=file.path(parent_dir, data_path))
  data_dt <- data_dt[, 2:ncol(data_dt)]
  data_melt_dt <- suppressWarnings(melt(data_dt, id.vars='Elapsed', variable.name='test', value.name='value'))
  data_melt_dt <- data_melt_dt[, c('well','image') := tstrsplit(test, ", Image")]
  platemap_dt <- fread(file=file.path(parent_dir, platemap_path))
  data_melt_dt <- merge(data_melt_dt, platemap_dt, by=c('well'))
  data_melt_dt[, value := as.numeric(value)]
  data_melt_dt[, donor_map := donor]
  data_melt_dt[, day_map := day]
  
  #unused columns
  if('V10' %in% names(data_melt_dt)) data_melt_dt[, V10 := NULL]
  if('test' %in% names(data_melt_dt)) data_melt_dt[, test := NULL]
  if('vol_Ts' %in% names(data_melt_dt)) data_melt_dt[, vol_Ts := NULL]
  
  data_melt_dt[, day := NULL]
  data_melt_dt[, donor := NULL]

    return(data_melt_dt)
  
  # b/c some platemaps have extra columns, only take the first 10
  #return(data_melt_dt[, names(data_melt_dt)[1:10], with=F ])
}

load_incucyte_data <- function(parent_dir, experiment_name) {
  
  experiment_dir = file.path(parent_dir, experiment_name)
  
  filepaths_dt <- data.table(
    data_path=list.files(path=parent_dir, pattern=paste0('.+.', experiment_name, '.+txt'), recursive=TRUE)
  )
  
  filepaths_dt[,
               `:=`(
                 date = gsub('.+/([\\d\\.]+)_[^/]+$', '\\1', data_path),
                 donor = gsub('.+_D\\d{1,2}_donor(\\d).+', '\\1', data_path),
                 day = gsub(paste0('.+', experiment_name, '_D(\\d+).+'), '\\1', data_path),
                 measurement = gsub('.+(_donor\\d)?_(\\w+).txt', '\\2', data_path)
               )]
  
  #load incucyte plate maps
  platemaps_dt <- data.table(
    platemap_path=list.files(path=parent_dir, pattern=paste0('.*',experiment_name,'.*csv'), recursive=TRUE)
  )
  
  platemaps_dt[,
               `:=`(
                 donor = gsub('.+_D\\d{1,2}_inc_donor(\\d).+', '\\1', platemap_path),
                 day = gsub(paste0('.+_', experiment_name, '_D(\\d+).+'), '\\1', platemap_path)
               )]
  
  #merge filepaths_dt and platemaps_dt
  if (platemaps_dt[, all(donor == platemap_path)]) {
    combined_dt <- merge(filepaths_dt, platemaps_dt, by=c('day'))
    combined_data_dt <- combined_dt[, 
      load_single_measurement(parent_dir, data_path, platemap_path), 
      by=c('date', 'day', 'measurement')]
  } else {
    combined_dt <- merge(filepaths_dt, platemaps_dt, by=c('donor','day'))  
    combined_data_dt <- combined_dt[, 
      load_single_measurement(parent_dir, data_path, platemap_path), 
      by=c('date', 'donor', 'day', 'measurement')]
  }
  
  # replace plate donor with explicit donor on map, since some donors are switched
  combined_data_dt[, donor := donor_map]
  combined_data_dt[, donor_map := NULL]
  
  # get number of replicates per timepoint and condition
  combined_data_dt[, n_rep := .N, by=c('measurement','donor','day','well','Elapsed')]
  
  #if individual images available, throw away mean value (where image == well)
  combined_data_dt <- combined_data_dt[n_rep == 1 | well != image]
  
  # if car is '', then this data is extra k-only wells, remove them
  combined_data_dt <- combined_data_dt[car !='']
  
  #d7 for d3 and d4 is actually d8
  combined_data_dt[day==7, day := 8]
  
  # clean up days
  combined_data_dt$day_f <- factor(combined_data_dt$day, levels=as.character(unique(combined_data_dt$day)))
  combined_data_dt$car_f <- factor(combined_data_dt$car, levels=c('41BB', 'BAFF-R', 'CD28', 'CD40', 'KLRG1', 'TACI', 'TNR8', 'zeta', 'none'))
  
  # remove noise:
  #1. identify large differences between successive measurements
  combined_data_dt[, diff := c(0, abs(value[1:.N-1] - value[2:(.N)])), 
              by=c('measurement','donor','day','k562','well','image','car')]
  #2. set difference cutoff to 5e5 in 1 measure
  combined_data_dt[measurement == 'redintensity', noise := diff > 5e5]
  combined_data_dt[measurement == 'redarea', noise := diff > 4e4]
  combined_data_dt[measurement == 'redcount', noise := diff > 250]
  
  #3. flag adjacent measurements as noisy also (this is slow!)
  remove_around_noise <- function(noise_vector) {
    return(as.numeric(any(as.logical(noise_vector))))
  }
  combined_data_dt[, noise_adj := frollapply(noise, 3, remove_around_noise, align='center', fill=0), 
              by=c('measurement','donor','day','k562','well','image','car')]
  
  #4. noisy measurements map across wells for the same day/donor/t_type (same plate)
  combined_data_dt[, noise := any(noise_adj > 0), by=c('measurement','day','k562','Elapsed','donor')]
  combined_data_dt[, noise_adj := NULL]
  
  # normalize
  combined_data_dt[, value_norm := value/(value[Elapsed==0]+0.001), by=c('well','day','donor','car','measurement','k562','image')]
  combined_data_dt[, mean_value_norm := mean(value_norm), by=c('day','donor','car','Elapsed','measurement','k562')]
  combined_data_dt[, std_dev := sd(value_norm), by=c('day','donor','car','Elapsed','measurement','k562')]
  combined_data_dt[, std_error := std_dev/sqrt(.N), by=c('day','donor','car','Elapsed','measurement','k562')]

  #re-order days
  combined_data_dt[, day_f := factor(day_f, levels=sort(as.numeric(levels(day_f))))]
  
  
  return(combined_data_dt)
}

incucyte_dt <- rbind(
  load_incucyte_data(file.path(incucyte_dir,'TCSL105'), 'TCSL105')[, t_type := 'cd4'],
  load_incucyte_data(file.path(incucyte_dir,'TCSL114'), 'TCSL114')[, t_type := 'cd4'],
  load_incucyte_data(file.path(incucyte_dir,'TCSL091'), 'TCSL091')[, t_type := 'cd8']
)

# reorder days
incucyte_dt[, day_f := factor(day_f, levels=sort(as.numeric(levels(day_f))))]

# for cd4 samples, copy none to d3 and d4 samples
incucyte_dt <- rbind(
  incucyte_dt[!is.na(donor)], 
  copy(incucyte_dt[t_type == 'cd4' & is.na(donor)])[, donor := 3], 
  copy(incucyte_dt[t_type == 'cd4' & is.na(donor)])[, donor := 4])

# add normalization to no T cells
incucyte_dt[, mean_value_none := mean_value_norm/mean(mean_value_norm[car == 'none']), by=c('measurement','donor','day','k562','Elapsed', 't_type')]
incucyte_dt[, value_none := value_norm/mean(mean_value_norm[car == 'none']), by=c('measurement','donor','day','k562','Elapsed', 't_type')]
incucyte_dt[, std_dev_value_none := sd(value_none), by=c('measurement','donor','day','k562','Elapsed', 't_type','car')]
incucyte_dt[, std_err_value_none := std_dev_value_none/sqrt(.N), by=c('measurement','donor','day','k562','Elapsed', 't_type','car')]
incucyte_dt[, std_dev_value_none_log := sd(log(value_none)), by=c('measurement','donor','day','k562','Elapsed', 't_type','car')]
incucyte_dt[, std_err_value_none_log := std_dev_value_none_log/sqrt(.N), by=c('measurement','donor','day','k562','Elapsed', 't_type','car')]


save(incucyte_dt, file=file.path(here::here('..','data','incucyte.Rdata')))