library('RColorBrewer')

#display all colors: display.brewer.all(type="seq")
blues <- brewer.pal(9,'Blues')[c(5,7)] #41bb, cd28
greens <- brewer.pal(9,'Greens')[c(5,7)] #baffr, taci
oranges <- brewer.pal(9,'Oranges')[c(4,6)] #cd40, tnr8
reds <- brewer.pal(9,'PuRd') #klrg1
purple <- brewer.pal(9,'Purples')[c(7)] #zeta
grey <- 'grey30' #unt

colors_light <- toupper(paste0('#',c("978dc6","75abd8","44d9fa","57f2de",
                             "73f1c3","b9efb2","d7f1aa","f6f3a2","f7de9e","f7c898")))
colors_mid <- toupper(paste0('#',c("54478c","2c699a","048ba8","0db39e",
                           "16db93","83e377","b9e769","efea5a","f1c453","f29e4c")))
colors_dark <- toupper(paste0('#',c("241f3d","132e44","023d4a","064f46",
                            "096040","267f1a","587f15","85810d","85620b","81460a")))

# CAR COLORS
# --------------------------------------------------------------------------------------------------

# new set
grey <- 'grey10' #zeta
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

# SCRNA SUBTYPE COLORS
# --------------------------------------------------------------------------------------------------
all_colors <- list(
    'NAIVE-CD62L'=colors_mid[1],
    #'NAIVE-CD7'=colors_mid[2],
    'NAIVE-CD7'='#A084AF',
    'STAT1-IRF1'=colors_mid[4],
    'STIM-INACT'=colors_mid[4],
    'STIM-TC2'=colors_dark[10],
    'TH2'=colors_dark[10],
    'STIM-UNTR'=colors_mid[3],
    'MEM-CD29'=colors_light[3],
    'STIM-IFNG'=colors_mid[6],
    'STIM-GNLY'=colors_mid[8],
    'STIM-OXPHOS'=colors_mid[9],
    'STIM-GLYC'='#ED5F1C') #brewer.pal(9,'Reds')[6]) #colors_mid[10]

# SCRNA PHASE COLORS
# --------------------------------------------------------------------------------------------------

phase_colors <- list(
    'S'='#B27FBE',
    'G1'='#6CA6DC',
    'G2M'='#EF5F9A')


# OTHER_PALETTES
# --------------------------------------------------------------------------------------------------


outlier_cols <- brewer.pal(11, 'PiYG')[c(2,11)]

outlier_cols_light <- brewer.pal(11, 'PiYG')[c(4,8)]

cd19_val <- brewer.pal(11,'PRGn')[c(3,10)]
