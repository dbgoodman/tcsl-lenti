library(data.table)
library(ggplot2)
library(RColorBrewer)
library(nlme)
library(plyr)
library(cowplot)
library(ggrepel)
library(patchwork)
#library(dplyr)

point_size= 5
line_size= 0.8
errorbar_size= 0.4

# CAR COLORS
# --------------------------------------------------------------------------------------------------

#source(here::here('r','fig_colors.R'))

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

car_order=c(
  'CD28','41BB','BAFF-R','TACI',
  'CD40','TNR8','KLRG1','zeta',
  'Unt','None')

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
