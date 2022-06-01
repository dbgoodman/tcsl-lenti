library(data.table)
library(ggplot2)
library(RColorBrewer)
library(nlme)
library(plyr)
library(cowplot)
#library(dplyr)


point_size= 9
line_size= 1.5
errorbar_size= 0.8

# CAR COLORS
# --------------------------------------------------------------------------------------------------

source(here::here('r','fig_colors.R'))

#display all colors: display.brewer.all(type="seq")
blues <- brewer.pal(9,'Blues')[c(6,8)] #41bb, cd28
greens <- brewer.pal(9,'Greens')[c(4,6)] #baffr, taci
oranges <- brewer.pal(9,'Oranges')[c(4,6)] #cd40, tnr8
reds <- brewer.pal(9,'PuRd') #klrg1
grey <- 'grey10' #zeta
purples <- brewer.pal(9, 'Purples')[c(4, 6)] #TNR8, CD40

#current set
car_colors <- c(
  '41BB'=blues[1], 'BAFF-R'=greens[1], 
  'CD28'=blues[2], 'CD40'=purples[2], 
  'KLRG1'=reds[5], 'TACI'=greens[2], 
  'TNR8'=purples[1], 'zeta'=grey)

#current set
car_colors <- c(car_colors, 'None'='#999999', 'Unt'='#999999', 'zeta'='black')

car_order=c('CD28','41BB','BAFF-R','TACI','CD40','TNR8','KLRG1','zeta','Unt','None')

car_colors <- car_colors[car_order]

plot_theme <- theme(
  legend.position = c(.05, .95),
  legend.justification = c("left", "top"),
  legend.box.just = "right",
  legend.box.background = element_rect(fill='white'),
  legend.title=element_blank(),
  legend.spacing.y = unit(2, 'pt'),
  legend.margin = margin(2, 3, 3, 3),
  legend.key.height=unit(15,'pt'),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(fill=NA))

point_with_family <- function(layer, family) {
  old_geom <- layer$geom
  new_geom <- ggproto(
    NULL, old_geom,
    draw_panel = function(self, data, panel_params, coord, na.rm = FALSE) {
      pts <- ggproto_parent(GeomPoint, self)$draw_panel(
        data, panel_params, coord, na.rm = na.rm
      )
      pts$gp$fontfamily <- family
      pts
    },
    draw_key = function(self, data, params, size) {
      pts <- ggproto_parent(GeomPoint, self)$draw_key(
        data, params, size
      )
      pts$gp$fontfamily <- family
      pts
    }
  )
  layer$geom <- new_geom
  layer
}

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
                  xlim = c(52, 60),
                  fontface     = 'bold'
  ) +
  
  theme_minimal() + 
  plot_theme + 
  theme(legend.position = "none", plot.margin = unit(c(0.1, 2, 0.1, 0.1), "cm")) +
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
      xlim = c(35, 38),
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

ggsave(here::here('..','figs','editor_rebuttal','single_bb_baffr.pdf'), (tumor_new | single_bb_baffr), device='pdf', height=5, width=10, units='in') 

###### STATS
#https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r
library('datarium')
set.seed(123)
data("selfesteem2", package = "datarium")
selfesteem2 <- selfesteem2 %>%
  gather(key = "time", value = "score", t1, t2, t3) %>%
  convert_as_factor(id, time)
res.aov <- anova_test(
  data = selfesteem2, dv = score, wid = id,
  within = c(treatment, time)
)
get_anova_table(res.aov)

#normality, qqplot
ggqqplot(growth_dt, "Area1", ggtheme = theme_bw()) +
  facet_grid(day ~ car, labeller = "label_both")

#compute anova
growth_test_group <- growth_dt[day< 43][, list(Area1, mouse_id, car, day, mouse_number)]

res.aov <- anova_test(
  data = growth_test_group,
  dv = Area1, wid = mouse_id,
  within = c(car, day)
)

get_anova_table(res.aov)

one.way <- growth_test_group %>%
  group_by(day) %>%
  anova_test(dv = Area1, wid = mouse_id, within = car) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way

pwc <- growth_test_group %>%
  group_by(day) %>%
  pairwise_t_test(
    Area1 ~ car,
    p.adjust.method = "bonferroni"
  )
pwc


# repeated meaures anova example
# --------------------------------------------------------------------------------------------------
phlebitisdata = read.table("phlebitis.csv", header=T, sep=",")
attach(phlebitisdata)
#phlebitisdata #This isn’t necessary, but you might want to see the data structure
#ggplot(phlebitisdata) + geom_line(aes(x=Time, color=Treatment, group=Animal, y=Y))

aov.p = aov(Y~(factor(Treatment)*factor(Time))+Error(factor(Animal)),phlebitisdata )
summary(aov.p)

library(nlme) #activates the nlme library

interaction.plot (Time, factor(Treatment), Y, lty=c(1:3),lwd=2,ylab="mean of Y", xlab="time", trace.label="Treatment")
nestinginfo <- groupedData(Y ~ Treatment | Animal, data= phlebitisdata)
fit.compsym <- gls(Y ~ factor(Treatment)*factor(Time), data=nestinginfo, corr=corCompSymm(, form= ~ 1 | Animal))
fit.nostruct <- gls(Y ~ factor(Treatment)*factor(Time), data=nestinginfo, corr=corSymm(, form= ~ 1 | Animal), weights = varIdent(form = ~ 1 | Time))
fit.ar1 <- gls(Y ~ factor(Treatment)*factor(Time), data=nestinginfo, corr=corAR1(, form= ~ 1 | Animal))
fit.ar1het <- gls(Y ~ factor(Treatment)*factor(Time), data=nestinginfo, corr=corAR1(, form= ~ 1 | Animal), weights=varIdent(form = ~ 1 | Time))
anova(fit.compsym, fit.nostruct, fit.ar1, fit.ar1het) #compares the models
fit.ar1polytime <- gls(Y ~ factor(Treatment)*poly(Time, degree = 3), data=nestinginfo, corr=corAR1(, form= ~ 1 | Animal))

summary(fit.ar1polytime)

anova(fit.compsym)
anova(fit.ar1)
anova(fit.ar1polytime)
anova(fit.ar1polytime, fit.ar1)

detach(phlebitisdata)

# our data
# --------------------------------------------------------------------------------------------------
cars <- c('BAFFR','41BB','Unt','NT','BAFFR-41BB')
growth_test_group <- growth_dt[day < 43 & car %in% cars, list(Area1, mouse_id, car, day, mouse_number)]
growth_test_group[, lArea1 := log(Area1+0.1)]

with(growth_test_group, {
  interaction.plot (day, factor(car), Area1, lty=c(1:2),lwd=2,ylab="mean of Area", xlab="Day", trace.label="car")
  nestinginfo <- groupedData(Area1 ~ car | mouse_number, data= growth_test_group)
  fit.compsym <- gls(Area1 ~ factor(car)*factor(day), data=nestinginfo, corr=corCompSymm(, form= ~ 1 | mouse_number))
  fit.nostruct <- gls(Area1 ~ factor(car)*factor(day), data=nestinginfo, corr=corSymm(, form= ~ 1 | mouse_number), weights = varIdent(form = ~ 1 | day))
  fit.ar1 <- gls(Area1 ~ factor(car)*factor(day), data=nestinginfo, corr=corAR1(, form= ~ 1 | mouse_number))
  fit.ar1het <- gls(Area1 ~ factor(car)*factor(day), data=nestinginfo, corr=corAR1(, form= ~ 1 | mouse_number), weights=varIdent(form = ~ 1 | day))
  
  # compare models for AIC and loglik
  print(anova(fit.compsym, fit.nostruct, fit.ar1, fit.ar1het)) #compares the models

  print(anova(fit.nostruct))
  print(anova(fit.compsym))
  print(anova(fit.ar1het))
})

# Unstructured Variance/Covariance Structure gave lowest AIC and loglik, so using that
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

models <- apply(growth_dt[day < 43, combn(unique(car),2)], 2, function(x) compute_pairwise_mixed_effects(growth_dt[day < 43 & car %in% x]))
p_values <- data.table(t(sapply(models, function(model_i) c(model_i$cars, anova(model_i$model)$`p-value`[4]))))
names(p_values) <- c('car_1','car_2','pvalue')
p_values[, pvalue := as.numeric(pvalue)] 
p_values[, car_1 := factor(car_1, levels=growth_dt[, unique(car)])]
p_values[, car_2 := factor(car_2, levels=growth_dt[, unique(car)])]
}

pairwise_models <- data.table(
  t(growth_dt[day < 43, combn(unique(car),2)]))[, i := .I][, 
    list(gls_fit= compute_pairwise_mixed_effects(
      growth_dt[day < 43 & car %in% c(V1,V2)])), by=c('V1','V2','i')]

ggplot(p_values) + geom_tile(aes(x=car_1, y=car_2, fill=-log10(pvalue))) + geom_text(data=p_values[pvalue < 0.05], aes(x=car_1, y=car_2), color='white', label='*', size=10, vjust=0.9) + theme_classic() + scale_fill_distiller(palette='Oranges',direction=1) + theme(legend.position='left')

# Auc approach
# -----------------------------------------------

growth_test_group <- growth_dt[, list(Area1, mouse_id, car, day, mouse_number)]


calcauc<-function(data) {
  psum<-function(x) rowSums(embed(x,2))
  sum(psum(data$Area1) * diff(data$day)/ 2)
}

with(growth_test_group[day < 43, list(auc=calcauc(.SD)), by=c('mouse_number','car')], 
     pairwise.wilcox.test(auc, car, p.adjust.method='bonferroni'))

with(growth_test_group[day < 43, list(auc=calcauc(.SD)), by=c('mouse_number','car')], 
     pairwise.t.test(auc, car, p.adjust.method='bonferroni'))

# another attempt with aov
# https://stats.idre.ucla.edu/r/seminars/repeated-measures-analysis-with-r/

# this works:
aov_model <- aov(Area1 ~ car * day + Error(mouse_number), data = copy(growth_dt)[, day := factor(day)])
summary(aov_model)


