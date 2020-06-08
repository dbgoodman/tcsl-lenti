library('RColorBrewer')

#display all colors: display.brewer.all(type="seq")
blues <- brewer.pal(9,'Blues')[c(5,7)] #41bb, cd28
greens <- brewer.pal(9,'Greens')[c(5,7)] #baffr, taci
oranges <- brewer.pal(9,'Oranges')[c(4,6)] #cd40, tnr8
reds <- brewer.pal(9,'PuRd') #klrg1
purple <- brewer.pal(9,'Purples')[c(7)] #zeta
grey <- 'grey30' #unt

# new set
grey <- 'grey40' #zeta
purples <- brewer.pal(9, 'Purples')[c(5, 7)] #TNR8, CD40

# old set
old_car_colors <- c(
    '41BB'=blues[1], 'BAFF-R'=greens[1], 
    'CD28'=blues[2], 'CD40'=oranges[1], 
    'KLRG1'=reds[7], 'TACI'=greens[2], 
    'TNR8'=oranges[2], 'zeta'=purple)

#current set
car_colors <- c(
    '41BB'=blues[1], 'BAFF-R'=greens[1], 
    'CD28'=blues[2], 'CD40'=purples[2], 
    'KLRG1'=reds[5], 'TACI'=greens[2], 
    'TNR8'=purples[1], 'zeta'=grey)

car_order=c('CD28','41BB','BAFF-R','TACI','CD40','TNR8','KLRG1','zeta')

show_colors <- function(car_colors, car_order) {
    ggplot(data.table(
            car_colors, factor(car_order, levels=car_order))) + 
        geom_tile(
            aes(y=factor(car_order, levels=car_order), 
            fill=car_order, x=1)) + 
        scale_fill_manual(values=as.character(car_colors))
}