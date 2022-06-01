# ------------------------------------------------------------------------------

# TCSL 166 figs

growth_dt <- fread(paste0(getwd(), "/r/in_vivo/tcsl166_growth.csv"))
growth_dt[, length_1 := as.numeric(length_1)]
growth_dt[, length_2 := as.numeric(length_2)]
growth_dt <- growth_dt[car != '']

growth_dt <- rbindlist(
  list(
    growth_dt,
    growth_dt[, list(day=0, width_1=0, length_1=0), by=c('cage','mouse','car')]),
  fill=T)

growth_dt[, car := mapvalues(
  factor(car),
  from=levels(factor(car)), 
  to=c('41BB','BAFF-R','CD28','CD40','KLRG1','None','TACI','Unt','zeta'))]

growth_dt[, vol_lww := length_1 * width_1^2 * 0.5]
growth_dt[, vol_lwl := length_1^2 * width_1 * 0.5]
growth_dt[, vol_lwlw := length_1 * width_1 * (length_1+width_1)/2 * 0.5]

growth_dt_means <- growth_dt[, list(
  vol_lww = mean(vol_lww, na.rm=T),
  vol_lww_se = sd(vol_lww, na.rm=T)/sqrt(.N),
  vol_lwl = mean(vol_lwl, na.rm=T),
  vol_lwl_se = sd(vol_lwl, na.rm=T)/sqrt(.N),
  vol_lwlw = mean(vol_lwlw, na.rm=T),
  vol_lwlw_se = sd(vol_lwlw, na.rm=T)/sqrt(.N)),
  by=c('car','day')]

left_group <- c('None','Unt','41BB','BAFF-R','CD28','TACI')
right_group <- c('None','Unt','KLRG1','zeta')
new_group <- c('41BB','BAFF-R','CD40','Unt','None')
car_linetypes <- setNames(rep(1, length(car_colors)), names(car_colors))
car_linetypes['None'] <- 2

tumor_new <- ggplot(growth_dt_means[day <50 & car %in% new_group][car == 'Untransduced', car := 'Unt'],
                  aes(x=day, y=vol_lwl, color=car, linetype=car, label=car)) +
  geom_line(size=line_size) +
  point_with_family(
    geom_point(size=point_size, shape="\u25CF", show.legend=F),
    "Arial Unicode MS") + 
  geom_errorbar(aes(ymin=vol_lwl-vol_lwl_se, ymax=vol_lwl+vol_lwl_se), show.legend=F, size=errorbar_size, width=0.3) + 
  scale_linetype_manual("", values=car_linetypes[new_group]) +
  scale_color_manual("", values=car_colors[new_group]) +
  
  #labels
  coord_cartesian(xlim= c(2,48), clip = "off") +
  geom_text_repel(data=unique(growth_dt_means[day < 50 & car %in% new_group][day == max(day)]),
                  force        = 3,
                  nudge_x      = 2,
                  direction    = "y",
                  vjust        = 0.5,
                  hjust        = 0,
                  segment.size = 0,
                  xlim = c(52, Inf),
                  fontface     = 'bold'
  ) +
  
  theme_minimal() + 
  plot_theme + 
  theme(legend.position = "none", plot.margin = unit(c(0.1, 2.5, 0.1, 0.1), "cm")) +
  labs(y=bquote('Tumor Volume'~(mm^3)), x='Day')

# --------------------------------------------------------------------------------------------------

growth_dt <- fread(paste0(getwd(), "/r/in_vivo/tcsl179_growth.csv"))
growth_dt <- growth_dt[!(is.na(day))]
growth_dt[, car := gsub('NT', 'None', car)]
growth_dt[, c('car1','car2') := tstrsplit(car, '-')]

#remove baffr-41bb #3
growth_dt <- growth_dt[mouse_number != '809-3']

#make numeric mouse id
growth_dt <- growth_dt[growth_dt[day==6, list(mouse_number, mouse_id= seq((.N))), by='car'], on=c('mouse_number','car')]

stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

ggplot(growth_dt) + geom_line(aes(x=day, y=Area1, color=car, group=mouse_number)) + facet_grid(~car)

growth_mean_dt <- growth_dt[,
  list(
    car1= car1[1], car2= car2[1], 
    area_mean = mean(Area1), 
    area_sd = sd(Area1), area_se = stderr(Area1)),
  by=c('day','car')]


plot_car_group <- function(growth_mean_dt, group1) {
  
  ggplot(growth_mean_dt[car %in% group1 & day < 43],
         aes(x=day, y=area_mean, group=car, color=car1, 
             ymin=area_mean-area_se, ymax=area_mean+area_se,
             label=car)) +
    
    # errorbars
    geom_errorbar(show.legend=F, width=0.3, size=errorbar_size) +
    
    geom_errorbar(data=growth_mean_dt[car %in% group1 & !is.na(car2) & day < 43], 
      size=errorbar_size, linetype=2, show.legend=F, width=0.3, aes(color=car2)) + 
    
    # solid line for singles
    geom_line(data=growth_mean_dt[car %in% group1 & car != 'None'], size=line_size) + 
    
    geom_line(data=growth_mean_dt[car %in% group1 & car == 'None'], size=line_size, linetype=2) + 
    
    # hashed line for doubles
    geom_line(data=growth_mean_dt[car %in% group1 & !is.na(car2) & day < 43], 
      aes(color=car2), linetype=2, size=line_size,  show.legend=F) + 
    
    # full circle for single cars
    point_with_family(
      geom_point(data=growth_mean_dt[car %in% group1 & is.na(car2) & day < 43],
      size=point_size, shape="\u25CF", aes(color=car1), show.legend=F),
      "Arial Unicode MS") + 
  
    # left half circle for double cars
    point_with_family(
      geom_point(data=growth_mean_dt[car %in% group1 & !is.na(car2) & day < 43],
      size=point_size, shape="\u25D7", aes(color=car1),  show.legend=F),
      "Arial Unicode MS") + 
    
    # right half circle for double cars
    point_with_family(
      geom_point(data=growth_mean_dt[car %in% group1 & !is.na(car2) & day < 43], 
      size=point_size, shape="\u25D6", aes(color=car2), show.legend=F),
      "Arial Unicode MS") + 
  
    #black border 
    # point_with_family(
    #   geom_point(data=growth_mean_dt[car %in% group1 & day < 43], size=point_size,
    #   shape="\u25CB", color='black',  show.legend=F),
    #   "Arial Unicode MS") + 
    
    #labels
    coord_cartesian(xlim= c(2,33), clip = "off") +
    geom_text_repel(data=unique(growth_mean_dt[day < 43 & car %in% group1][day == max(day)]),
      force        = 3,
      nudge_x      = 2,
      direction    = "y",
      vjust        = 0.5,
      hjust        = 0,
      segment.size = 0,
      xlim = c(35, Inf),
      fontface     = 'bold'
    ) +
  
    scale_color_manual("", values=rename(car_colors, c('BAFF-R'='BAFFR'))[group1]) +
    theme_minimal() + 
    plot_theme + 
    theme(legend.position = "none", plot.margin = unit(c(0.1, 2, 0.1, 0.1), "cm")) +
    labs(y=bquote('Tumor Volume'~(mm^3)), x='Day')
}

single_bb_baffr <- plot_car_group(growth_mean_dt[day < 43], c('41BB','BAFFR','CD40','Unt','None'))
single_bb_baffr
doubles_bb_baffr <- plot_car_group(
  growth_mean_dt[day < 43],
  c('41BB','BAFFR','CD40','BAFFR-41BB','CD40-41BB','Unt','None')) +
  theme(plot.margin = unit(c(0.1, 3, 0.1, 0.1), "cm"))
  
doubles_28_40 <- plot_car_group(growth_mean_dt[day < 43], c('CD28','CD40','CD40-CD28','Unt','None'))

ggsave(
  here::here('..','figs','editor_rebuttal','combined_in_vivo-rebuttal.png'),
  (tumor_new / plot_spacer()) | (single_bb_baffr / doubles_bb_baffr),
  device='png', height=7.5, width=9.5, units='in')


# STATS
# --------------------------------------------------------------------------------------------------
# cars <- c('BAFFR','41BB','Unt','NT','BAFFR-41BB')
# growth_test_group <- growth_dt[day < 43 & car %in% cars, list(Area1, mouse_id, car, day, mouse_number)]
# growth_test_group[, lArea1 := log(Area1+0.1)]
# 
# with(growth_test_group, {
#   interaction.plot (day, factor(car), Area1, lty=c(1:2),lwd=2,ylab="mean of Area", xlab="Day", trace.label="car")
#   nestinginfo <- groupedData(Area1 ~ car | mouse_number, data= growth_test_group)
#   fit.compsym <- gls(Area1 ~ factor(car)*factor(day), data=nestinginfo, corr=corCompSymm(, form= ~ 1 | mouse_number))
#   fit.nostruct <- gls(Area1 ~ factor(car)*factor(day), data=nestinginfo, corr=corSymm(, form= ~ 1 | mouse_number), weights = varIdent(form = ~ 1 | day))
#   fit.ar1 <- gls(Area1 ~ factor(car)*factor(day), data=nestinginfo, corr=corAR1(, form= ~ 1 | mouse_number))
#   fit.ar1het <- gls(Area1 ~ factor(car)*factor(day), data=nestinginfo, corr=corAR1(, form= ~ 1 | mouse_number), weights=varIdent(form = ~ 1 | day))
#   
#   # compare models for AIC and loglik
#   print(anova(fit.compsym, fit.nostruct, fit.ar1, fit.ar1het)) #compares the models
# 
#   print(anova(fit.nostruct))
#   print(anova(fit.compsym))
#   print(anova(fit.ar1het))
# })

# Unstructured Variance/Covariance Structure gave lowest AIC and loglik, 
# so using that
compute_pairwise_mixed_effects <- function(data) {
  glsControl(maxIter=1000, msMaxIter=1000, tolerance=1e-3)
  with(data, tryCatch({
    
    print(unique(car))
    nestinginfo <- groupedData(Area1 ~ car | mouse_number,data=data)
    
    # unstructured covariance matrix has slightly lower AIC,
    # but does not converge for all/most pairs, so using symmetric
    # matrix instead. This one has fewer dfs and so is likely a better
    # choice anyway.
    
    # fit.nostruct <- gls(Area1 ~ factor(car)*day, 
    #   data=nestinginfo, 
    #   corr=corSymm(, form= ~ 1 | mouse_number),
    #   weights = varIdent(form = ~ 1 | day))
    
    fit.compsym <- gls(Area1 ~ factor(car)*factor(day), 
      data=nestinginfo,
      corr=corCompSymm(, form= ~ 1 | mouse_number))
        
    return(list(cars=unique(car), model=fit.compsym))},
    
    error=function(cond) {
      
      message("Here's the original error message:")
      message(cond)
      return(NA)

    }))
}

models <- apply(growth_dt[day < 43, combn(unique(car),2)], 2, 
    function(x) compute_pairwise_mixed_effects(
      growth_dt[day < 43 & car %in% x]))
p_values <- data.table(t(sapply(models,
    function(model_i) c(model_i$cars, anova(model_i$model)$`p-value`[4]))))

names(p_values) <- c('car_1','car_2','pvalue')

p_values[, pvalue := as.numeric(pvalue)]
p_values[, car_1 := factor(car_1, levels=growth_dt[, unique(car)])]
p_values[, car_2 := factor(car_2, levels=growth_dt[, unique(car)])]

p_values[, signif := ''][
  pvalue < 0.05, signif := '*'][
    pvalue < 0.001, signif := '**'][
      pvalue < 0.00001, signif := '***']

ggplot(p_values) + 
  geom_tile(aes(x=car_1, y=car_2, fill=-log10(pvalue))) + 
  geom_text(data=p_values[pvalue < 0.05], 
    aes(x=car_1, y=car_2, label=signif), color='white', size=10, vjust=0.9) +
  theme_classic() +
  scale_fill_distiller(palette='Oranges',direction=1) +
  theme(legend.position='left')