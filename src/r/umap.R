# map day 0 to both plus and minus
map_day_0 <- function(df) {
  return(rbind(
    df[k562!='none'], #gets all cd19+-
    #df[rep(keep_none, .N) & k562=='none' & day > 0], 
    df[k562=='none' & day == 0][, k562 := 'cd19+'], #none and day 0, assign as 19+
    df[k562=='none' & day == 0][, k562 := 'cd19-'], #none and day 0, assign as 19-
    df[car =='untrans' & day != 0][, k562 := 'cd19+'], #add back untrans
    df[car =='untrans' & day != 0][, k562 := 'cd19-']
  ))
}


sample_for_umap <- function(df, sample_size) {
  umap_dt <- exh_dt[!is.na(event_id), 
    .SD[event_id %in% sample(unique(event_id), 
      min(sample_size, length(unique(event_id))))],
    by=c('well', 'plate', 'day')]
}

cast_for_umap <- function(df) {
  #cast it so variables are columns and
  #subset sample umap data on variables
  umap_cast_dt <- dcast(df, 
    event_id + well + plate + day + t_type ~ variable, 
    value.var='value')
  return(umap_cast_dt)
}

scale_for_umap <- function(df, umap_vars, censor_negative_min) {
  umap_dt_in <- df[, ..umap_vars]
  
  # for points that are very negative, trim the to below the cutoff
  umap_dt_in[umap_dt_in < censor_negative_min] <- censor_negative_min
  
  # scale each input column
  umap_dt_in[ , 
    (names(umap_dt_in)) := lapply(.SD, scale), .SDcols = names(umap_dt_in)]
  
  return(umap_dt_in)
}

run_umap <- function(dt, chosen_dist, chosen_n_neighbor, sample_n, umap_vars,
    chosen_learning_rate=0.1) {
  
  umap_dt <- dt[!is.na(event_id), 
    .SD[event_id %in% sample(unique(event_id), 
      min(sample_n, length(unique(event_id))))],
    by=c('well', 'plate', 'day')]
  
  # for now use all variables in the channel map
  #umap_vars <- c(paste0(names(channel_map),'_t'))
  
  #cast it so variables are columns and
  #subset sample umap data on variables
  umap_cast_dt <- data.table(dcast(umap_dt, 
    event_id + well + plate + day + t_type ~ variable, 
    value.var='value'))
  
  umap_dt_in <- umap_cast_dt[, ..umap_vars]
  
  # scale each input column
  umap_dt_in[ , 
    (names(umap_dt_in)) := lapply(.SD, scale), .SDcols = names(umap_dt_in)]
  
  # run umap via uwot library
  umap_out <- umap(umap_dt_in, min_dist=chosen_dist,
      learning_rate=chosen_learning_rate,
      n_sgd_threads=1,
      n_neighbors=chosen_n_neighbor, 
      verbose=T, n_threads=8, n_trees=50,
      ret_nn=T)
  
  # nearest neighbors in original space
  nearest_neighbors <- umap_out$nn$euclidean$idx
  neighbor_dist <- umap_out$nn$euclidean$dist
  edge_weights <- scales::rescale(-neighbor_dist, to=c(0,1))
  edge_weights <- edge_weights - min(edge_weights)
  umap_out <- data.table(umap_out$embedding)
  
  #nearest neighbors in embedding
  # umap_nn <- cbind(1:nrow(umap_fcs_dt), 
  #     nnwhich(umap_fcs_dt[, list(umap_1, umap_2)], k=c(1:3)))
  # umap_nd <- cbind(rep(0, nrow(umap_fcs_dt)), 
  #     nndist(umap_fcs_dt[, list(umap_1, umap_2)], k=c(1:3)))
  # umap_ew <- scales::rescale(-umap_nd, to=c(0,1))
  # umap_ew <- umap_ew - min(umap_ew)
      
  # rename the columns
  names(umap_out) <- c('umap_1','umap_2')
  
  # add the umap output to the input dt
  umap_fcs_dt <- cbind(umap_cast_dt[, 1:length(names(umap_cast_dt))], umap_out)
  umap_fcs_dt[, id := 1:.N]
  
  # add back the annotations
  umap_fcs_dt <- unique(umap_dt[, 
    .(donor, car, k562, t_type, day, well, event_id)])[
      umap_fcs_dt, on=.(t_type,well,day,event_id)]
  
  source_python(here::here('py/leiden_clust.py'))
  clust_membership <- leiden_clust(nearest_neighbors,
    edge_weights=scales::rescale(-neighbor_dist, to=c(0,1)))
  umap_fcs_dt[, cluster := factor(clust_membership+1)]
  
  #renumber the clusters
  cluster_ranks <- umap_fcs_dt[, list(mean_day=mean(day)), by=cluster][, rank_day := rank(mean_day)]
  
  umap_fcs_dt <- umap_fcs_dt[cluster_ranks, on='cluster'][, c('cluster', 'rank_day') := list(rank_day, NULL)]
  
  return(umap_fcs_dt)
}

all_plots <- function(umap_fcs_dt, channel_map, umap_vars, cell_type, dendro=F) {
  
  umap_fcs_dt[, cluster := factor(cluster)]
  
  umap_cluster_dt <- umap_fcs_dt[, list(
    mean_umap_1= mean(umap_1), 
    mean_umap_2= mean(umap_2), 
    size= .N), by=cluster]
  
  #UMAP
  umap_plot <- ggplot() + 
  geom_point(data=umap_fcs_dt, 
    aes(x=umap_1, y=umap_2, color=cluster), size=0.2, alpha=0.2) +
  geom_label(data=umap_cluster_dt, aes(
    x=mean_umap_1, y=mean_umap_2,
    label=paste(cluster,size,sep='\n'), 
    color=cluster), alpha=0.3) +
  theme_minimal()
  
  #UMAP biomarker Heatmap
  umap_heatmap <- visualize_params_umap(umap_fcs_dt)
  
  #Dendrogram / cluster heatmap
  cluster_pct_vars <- c(paste0(
    names(channel_map),'_t'), 
    'FSC-A','SSC-A','donor','cd19')
  
    max_cluster <- max(as.numeric(umap_fcs_dt$cluster), na.rm=T)
  
  umap_fcs_dt[, cluster:= factor(cluster, levels = as.character(1:max_cluster))]
  
  
  umap_var_melt_dt <-melt(
    umap_fcs_dt[!is.na(donor)], measure.vars= c(umap_vars)) #add back in previously removed vars
  
  # Use z scale to plot params
  umap_var_melt_dt[, value_scaled := scale(value), by=variable]
  param_cluster_pct_dt <- umap_var_melt_dt[
    !is.na(cluster) & variable %in% cluster_pct_vars,
    list(
      clust_mean = mean(value_scaled),
      clust_sd = sd(value_scaled)),
        by=c('variable','cluster')]
  
  if (dendro) {
    # Make dendrogram
    cluster_mean_cast <- dcast(
      param_cluster_pct_dt, 
      variable ~ cluster, 
      value.var='clust_mean')
    
    cluster_dendro_m <- as.matrix(cluster_mean_cast[, -c(1)])
    rownames(cluster_dendro_m) <- unlist(cluster_mean_cast[,1])
    dendro_clusters <- as.dendrogram(hclust(d = dist(x = t(cluster_dendro_m))))
    dendro_vars <- as.dendrogram(hclust(d = dist(x = cluster_dendro_m)))
    
    # Create dendrogram plot
    dendro_vars_plot <- ggdendrogram(data = dendro_vars, rotate = TRUE) + 
      theme(axis.text.y = element_text(size = 6))
    dendro_clusters_plot <- ggdendrogram(data = dendro_clusters, rotate = TRUE) + 
      theme(axis.text.y = element_text(size = 6))
    
    # column order
    cluster_order <- order.dendrogram(dendro_clusters)
    cluster_names <- colnames(cluster_dendro_m[, cluster_order])
    param_cluster_pct_dt[, cluster:= factor(cluster, levels = cluster_names)]
    
    # row order
    var_order <- order.dendrogram(dendro_vars)
    var_names <- row.names(cluster_dendro_m[var_order, ])
    param_cluster_pct_dt$variable <- factor(
        param_cluster_pct_dt$variable,
        levels = var_names)
  } else {
    param_cluster_pct_dt[, cluster:= factor(cluster, levels = as.character(1:max_cluster))]
  }
  
   #annotate clusters as per variable values (removed gfp, myc, cd, donor)
  '%notin%' <- Negate('%in%')
  annotations <- param_cluster_pct_dt %>%
  mutate(variable=gsub("\\_t$", "", variable)) %>%
  filter(variable %notin% c('gfp','myc','cd','donor', 'zombie')) %>%
  mutate(cluster=as.character(cluster)) %>%
  group_by(cluster) %>% 
  mutate(rk = rank(-abs(clust_mean))) %>%
  filter(rk <= 3) %>% 
  arrange(cluster, rk) %>%
  mutate(variable = ifelse(clust_mean > 0, paste(variable,'+',sep = ""), paste(variable,'-',sep = ""))) %>%
  summarise(variable = paste(variable, collapse=", ")) %>% 
  data.table::transpose()
  colnames(annotations) <- as.character(unlist(annotations[1,]))
  annotations = annotations[-1, ]
  if (dendro) {
    annotations <- annotations[cluster_names]
  } else {
    annotations <- annotations[as.character(1:max_cluster)]
  }
  row.names(annotations) <- NULL 
  
  # heatmap
  var_heatmap <- ggplot(param_cluster_pct_dt) + 
      geom_tile(aes(x=cluster, y=variable, fill=clust_mean), color='black') + 
      geom_text(aes(x=cluster, y=variable, label=round(clust_mean, 2)), size=3, color='black') +
      scale_fill_distiller(palette='PiYG', limits=c(-3,3), oob=scales::squish) +
      theme_minimal() # + theme(axis.text.x = element_text(angle=90, hjust=1))
  # col dendro
  if (dendro) {
    dendro_data_col <- dendro_data(dendro_clusters, type = "rectangle")
    dendro_col <- axis_canvas(var_heatmap, axis = "x") + 
        geom_segment(data = segment(dendro_data_col), 
            aes(x = x, y = y, xend = xend, yend = yend))
    plot_dendroheat <- insert_xaxis_grob(var_heatmap, 
                dendro_col, grid::unit(0.2, "null"), position = "top")
    dendro_plot <- ggdraw(plot_dendroheat)
  }
  
  # annotation table
  tt <- ttheme_default(core=list(fg_params = list(fontsize=10), colhead=list(fg_params = list(fontsize=10, parse=TRUE))))
  annotations[1,] = str_wrap(annotations[1,], 0)
  tbl <- tableGrob(annotations, rows=NULL, theme=tt)
  # tbl <- strwrap(tbl, width = 10, simplify = FALSE)
  tbl$widths <- unit(rep(1/ncol(tbl), ncol(tbl)), "npc")
  
  # CLUSTERS BY DAY
  day_cluster_total_pct_dt <- umap_fcs_dt[,
    data.table(table(day,cluster, t_type, k562))][,
      list(cluster, pct=N/sum(N)), by=c('k562','day','t_type')]
  
  if (dendro) {
    day_cluster_total_pct_dt[, cluster := factor(cluster, levels = cluster_names)]
  } else {
    day_cluster_total_pct_dt[, cluster := factor(cluster, levels = as.character(1:max_cluster))]
  }
  
  clusters_by_day <- ggplot(day_cluster_total_pct_dt[k562=='cd19+'], aes(
      x=cluster, 
      y=factor(day, levels=c(0,6,15,24,33)), 
      fill=pct)) + 
    geom_tile(color='black') +
    geom_text(aes(label=round(pct*100)), size=3, color='white') +
    facet_grid(k562+t_type~., scales='free') +
    scale_fill_viridis_c('',direction=1, limits=c(0, day_cluster_total_pct_dt[, max(pct)])) +
    labs(title='Percent of cells in each day by cluster', y='Day',
      subtitle='(rows sum to 100%)') +
    theme_minimal() + ggtitle(paste('Percent of ',as.character(cell_type),' cells in each day by cluster'))
  # CARs BY CLUSTER/DAY
  car_cluster_day_indiv_k562_total_pct_dt <- umap_fcs_dt[,
    data.table(table(car,cluster,k562,t_type,day, donor))][,
      list(cluster, pct=N/sum(N)), by=c('car','k562','t_type','day', 'donor')][,
        pct_delta := scale(pct), by=c('car','k562','t_type','day', 'donor')]
  
  if (dendro) {
    car_cluster_day_indiv_k562_total_pct_dt[, cluster := factor(cluster, levels = cluster_names)] 
  } else {
    car_cluster_day_indiv_k562_total_pct_dt[, cluster := factor(cluster, levels = as.character(1:max_cluster))]
  }
  
  car_cluster_day_indiv_k562_total_pct_dt$day <- factor(
    car_cluster_day_indiv_k562_total_pct_dt$day, levels=c('0','6','15','24','33'))
  
  min_col <- brewer.pal(name='BrBG', n=11)[2]
  
  cars_by_cluster <- ggplot(car_cluster_day_indiv_k562_total_pct_dt[k562=='cd19+']) + 
      geom_tile(aes(x=cluster, y=car, fill=pct_delta), color='black') +
      facet_grid(day~., scales='free') +
      scale_fill_distiller('',
        palette='BrBG', 
        direction=1,
        #limits=c(-2.8, 2.8), 
        na.value=min_col) +
      labs(title='CAR enrichment\nper day per cluster', y='CAR',
        subtitle='CAR enrichment in cluster\n(z-score) (summed across rows)') +
      theme_minimal() + ggtitle(paste(as.character(cell_type),' CAR enrichment\nper day per cluster')) 
  
  # print all plots and table
  plot_grid(
    plot_grid(umap_plot, umap_heatmap, nrow=2, labels=c('A','B'), rel_heights=1), 
    plot_grid(var_heatmap, clusters_by_day, cars_by_cluster, tbl, 
      labels = c('C','D','E','F'), align = "h", axis='lr', nrow=4, 
      rel_heights = c(2, 2, 5, 1)), 
    nrow=1)
} 

visualize_clusters_umap <- function(df) {
  label_clusters_umap <- df %>% 
    group_by(cluster, cell_type) %>% 
    select(umap_1, umap_2) %>% 
    summarize_all(mean)
  
  ggplot(df, aes(x=umap_1, y=umap_2, color=as.factor(cluster))) + 
    geom_point(size=0.1)+facet_grid(car ~ day) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) + 
    geom_label_repel(aes(label=cell_type, size=0.5), data=label_clusters_umap) + 
    guides(colour=FALSE) +
    ggtitle('Leiden clusters')
}

ggumap <- function(df, df_cluster) {
  
  df[, cluster := factor(cluster)]
  
  df_cluster <- df[, list(
    mean_umap_1= mean(umap_1), 
    mean_umap_2= mean(umap_2), 
    size= .N), by=cluster]
  
  # whole plot, by cluster
  umap_whole <- ggplot() + 
    geom_point(data=df[sample(.N)],
      aes(x=umap_1, y=umap_2, color=cluster), size=0.2, alpha=0.2) +
    geom_label(data=df_cluster, aes(
      x=mean_umap_1, y=mean_umap_2,
      label=paste(cluster,size,sep='\n'), 
      fill=cluster), color='black',alpha=0.3) +
    theme_minimal() +
    ggtitle('UMAP by Cluster')

  # individual cluster plots
  umap_individual <- ggplot() + 
  geom_point(data=df[sample(.N)], 
      aes(x=umap_1, y=umap_2, color=cluster), size=0.2, alpha=0.2) +
    geom_label(data=df_cluster, aes(
      x=mean_umap_1, y=mean_umap_2,
      label=paste(cluster,size,sep='\n'), 
      color=cluster), alpha=0.3) +
    facet_wrap(~cluster, ncol=8) +
    theme_bw()

  # UMAP by day
  label_days_umap <- df %>% group_by(day) %>% select(umap_1, umap_2) %>% summarize_all(mean)
  umap_by_day <- ggplot(df[sample(.N)], aes(x=umap_1, y=umap_2, color=day)) + 
    geom_point(size=0.1,alpha = 0.5) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) + 
    geom_label_repel(aes(label=day, color=day), fill='white', data=label_days_umap) + 
    guides(colour=FALSE)+ggtitle('UMAP by Day') +
    scale_color_viridis(discrete=F)+
    scale_fill_viridis(discrete=F)+
    theme_minimal()

  plot_grid(
    plot_grid(umap_by_day, umap_whole, nrow=1),
    umap_individual, nrow=2, rel_heights=c(2.6,1))
}

color_plots <- function(umap_df, umap_vars) {
  cd4_colors = brewer.pal('Greens', n=9)[c(2,4,6,8,9)]
  cd8_colors = brewer.pal('Oranges', n=9)[c(2,4,6,8,9)]
  color_time <- ggplot(umap_df[][, cd19 := (k562=='cd19+')], 
      aes(x=umap_1, y=umap_2, color=interaction(day, t_type))) +
    geom_point(size=0.1, alpha=0.5) +
    facet_grid(car~cd19) +
    scale_color_manual(values=c(cd4_colors, cd8_colors)) +
    guides(colour = guide_legend(override.aes = list(alpha=1, size=3)))
  
  color_density <- ggplot(umap_df, aes(x=umap_1, y=umap_2)) +
    geom_hex(bins = 70) +
    facet_grid(car~cd19) +
    scale_fill_continuous(
      type = "viridis", limits=c(0,30), oob=scales::squish) +
    theme_bw()
  
  color_markers <- ggplot(
      melt(umap_df, measure.vars=umap_vars)[, 
        value.scaled := scale(value), by='variable'][,
          cd19 := (k562=='cd19+')], 
      aes(x=umap_1, y=umap_2, z=value.scaled)) +
    stat_summary_hex(bins = 70) +
    facet_grid(variable~cd19) +
    scale_fill_distiller(palette='RdYlBu', limits=c(-3, 3), oob=scales::squish) +
    theme_bw()
  
  plot_grid(color_time, color_density, color_markers, ncol=3)
}

viz.umap <- function(dat,dr,param.name,limits=NULL){
  ColVal <- pull(dat %>% select(param.name))
  if(is.null(limits)){
    Lim <- quantile(ColVal,probs=seq(0,1,0.01))[c(2,100)]
    p <- ggplot(dat, aes(x = umap_1, y =umap_2)) +
      geom_point(aes(color = ColVal), size=0.1)+
      theme_classic()+
      scale_color_distiller(name=param.name, palette = "RdYlBu", limits=Lim, oob=squish)+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+ggtitle(param.name)
  } else {
    p <- ggplot(dat, aes(x = umap_1, y = umap_2)) +
      geom_point(aes(color = ColVal), size=0.1)+
      theme_classic()+
      scale_color_distiller(name=param.name, palette = "RdYlBu", limits=c(limits[1],limits[2]), oob=squish)+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
      ggtitle(param.name)
  }
  p
}

visualize_params_umap <- function(df) {
  parameters <- colnames(df)[c(9, 11:15)] # colnames is how to select the right biomarkers
  p <- list()
  for (i in 1:length(parameters)) {
    p[[i]] <- viz.umap(dat=df, param.name=parameters[i])
  }
  do.call(plot_grid,p)
}

umap_contour_plots <- function(umap_df){
  #car x day (tall, narrow)
  car_by_day <- umap_df[, ] %>% ggplot(., aes(x=umap_1, y=umap_2)) +
    geom_point(size=0.1, color= "#DCDCDC") + 
    geom_density2d(aes(colour = stat(nlevel)), contour_var = "ndensity") + facet_grid(day ~ car) + 
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), legend.title = element_blank()) + 
    ggtitle("CAR Densities by Day")

  # clusters x car
  car_total <- umap_df %>% ggplot(., aes(x=umap_1, y=umap_2)) +
    geom_point(size=0.1, color= "#DCDCDC") +
    geom_density2d(aes(colour = stat(nlevel)), contour_var = "ndensity") + 
    facet_wrap(~ car) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), legend.title = element_blank()) + 
    ggtitle("CAR Densities by Cluster")
  
  # CLUSTERS BY DAY
  max_cluster <- max(as.numeric(umap_df$cluster), na.rm=T)
  day_cluster_total_pct_dt <- umap_df[,
    data.table(table(day,cluster, t_type, k562))][,
      list(cluster, pct=N/sum(N)), by=c('k562','day','t_type')]
  
  day_cluster_total_pct_dt[, cluster := factor(cluster, levels = as.character(1:max_cluster))]
  
  clusters_by_day <- ggplot(day_cluster_total_pct_dt, aes(
      x=cluster, 
      y=factor(day, levels=c(0,6,15,24,33)), 
      fill=pct)) + 
    geom_tile(color='black') +
    geom_text(aes(label=round(pct*100)), size=3, color='white') +
    facet_grid(k562+t_type~., scales='free') +
    scale_fill_viridis_c('',direction=1, limits=c(0, day_cluster_total_pct_dt[, max(pct)])) +
    labs(title='Percent of cells in each day by cluster', y='Day',
      subtitle='(rows sum to 100%)') +
    theme_minimal() + ggtitle(paste('Percent of cells in each day by cluster'))
  
  # clusters
  clusters_total <- umap_df %>% ggplot(., aes(x=umap_1, y=umap_2))+
    geom_point(size=0.1, color= "#DCDCDC") + 
    geom_density2d(aes(colour = stat(nlevel)), contour_var = "ndensity") + 
    facet_wrap(~ cluster, nrow=2) + 
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), legend.title = element_blank()) + 
    ggtitle("Cluster Densities")
  
  plot_grid(
    plot_grid(clusters_total, car_total, ncol=1),
    plot_grid(clusters_by_day, car_by_day, ncol=1, rel_heights=c(1,4)),
    ncol=2)
  
}

all_by_day <- function(umap_df) {
  umap_df %>% ggplot(., aes(x=umap_1, y=umap_2))+geom_point(size=0.1, color= "#DCDCDC") + geom_density2d(aes(colour = stat(nlevel)), contour_var = "ndensity") + facet_grid(. ~ day) + theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), legend.title = element_blank()) + ggtitle("All Peaks by Day")
}

visualize_clusters_umap <- function(df) {
  label_clusters_umap <- df %>% group_by(cluster, cell_type) %>% select(umap_1, umap_2) %>% summarize_all(mean)
  ggplot(df, aes(x=umap_1, y=umap_2, color=as.factor(cluster)))+geom_point(size=0.1) +theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+geom_label_repel(aes(label=cell_type, size=10), data=label_clusters_umap)+guides(colour=FALSE)+ggtitle('Leiden clusters')
}

biomarker_corr <- function(df) {
  biomarkers <- df[, c(9, 11:15)]
  corr <- round(cor(biomarkers),2)
  p.mat <- cor_pmat(biomarkers)
  ggcorrplot(corr, method = "circle", type = "upper", outline.col = "white", p.mat = p.mat, ggtheme= ggplot2::theme_gray, colors = c("#6D9EC1", "white", "#E46726"))
}

biomarker_corr_by_day <- function(df) {
  plot_grid(df %>% filter(day==0) %>% biomarker_corr()
            , df %>% filter(day==6) %>% biomarker_corr(), df %>% filter(day==15) %>% biomarker_corr(), df %>% filter(day==24) %>% biomarker_corr(), labels = c('Day 0', 'Day 6', 'Day 15', 'Day 24'), nrow=3)
}

biomarker_corr_by_car <- function(df) {
  plot_grid(df %>% filter(car=='41BB') %>% biomarker_corr()
            , df %>% filter(car=='BAFF-R') %>% biomarker_corr(), df %>% filter(car=='CD28') %>% biomarker_corr(), df %>% filter(car=='CD40') %>% biomarker_corr(), df %>% filter(car=='TACI') %>% biomarker_corr(), df %>% filter(car=='TNR8') %>% biomarker_corr(), df %>% filter(car=='zeta') %>% biomarker_corr(), labels = c('41BB','BAFF-R', 'CD28', 'CD40', 'TACI', 'TNR8', 'zeta'), nrow=3)
}

