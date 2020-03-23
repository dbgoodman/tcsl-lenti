library(data.table)
library(ggplot2)
library(knitr)
library(ggrepel)
library(RColorBrewer)
library(plotly)
library(DESeq2)
library(rnaseqGene)
library(gtable)
library(plyr)    
library(grid)
library(cowplot)


## SETUP 

ngs.dir <- incucyte_dir <- here::here('..','..','ngs_data')

dataset.name <- '2019.01.31.illumina_06'
data.dir <- paste0(tcsl.dir, dataset.name, '/r/out/')
output.dir <- paste0(tcsl.dir, dataset.name, '/r/out/notebook_06/')
function.dir <- paste0(tcsl.dir, dataset.name, '/r/functions/')
metadata.filename <- paste0(tcsl.dir, dataset.name, '/', dataset.name, '_layout.csv')
car.seq.filename <- paste0(tcsl.dir, dataset.name, '/fa/ref.csv')
min.sort.reads <- 100

# global options
opts_knit$set(fig_dim= c(10,10))
theme_set(theme_bw())

# load data
load(paste0(data.dir, 'read_counts.rdata'))

# load functions
source(paste0(function.dir, 'calculate_bin_corr.r'))
source(paste0(function.dir, 'calculate_lm.r'))
source(paste0(function.dir, 'overlapping_strip_labels.r'))
source(paste0(function.dir, 'run_deseq.r'))

# pvalue threshold
padj.thresh <- .05

# convert to days
read.counts[grep('(\\d+)d', timepoint), day := as.numeric(gsub('(\\d+)d','\\1', timepoint))]
read.counts[grep('(\\d+)h', timepoint), day := as.numeric(gsub('(\\d+)h','\\1', timepoint)) / 24]