# cell type identification
# Author: LeeAnn Ramsay
# Date: 26 October 2022

library(Rcpp)
library(data.table)
source('~/GITHUB/MIMIC_CyTOF/Tools/R/clustering_functions.R')
# ____ load cytofkit code ______
# cytofkit source code available at https://github.com/JinmiaoChenLab/cytofkit
# set paths to cytofkit source code location
source('~/GITHUB/MIMIC_CyTOF/Tools/cytofkit/Rphenograph.R')
source('~/GITHUB/MIMIC_CyTOF/Tools/cytofkit/RcppExports.R')
sourceCpp('~/GITHUB/MIMIC_CyTOF/Tools/cytofkit/jaccard_coeff.cpp')
sourceCpp('~/GITHUB/MIMIC_CyTOF/Tools/cytofkit/RcppExports.cpp')
#______ end cytofkit source loading ______
# Rphenograph required packages
library(ggplot2)
library(RANN,nn2)
# load igraph version 1.2.5
library(devtools)
#install_version("igraph", version = '1.2.5', repos = "http://cran.us.r-project.org")
library(igraph)

#__________________________________
# load data ####
#__________________________________
panel <- as.data.frame(fread('~/WLData/MIMIC_IMC/MIMIC_Panel.csv'))
all_markers <- panel$Target
id_markers <- c('CD14','CD31','CD45','CD4','CD68','CD20','CD8a',
                'Sox10','CD45RA','CD3','S100','CD45RO')

MIM <- as.data.frame(fread('~/WLData/MIMIC_IMC/Processed_image_data/MIMIC_meanItensity_location_cell.tsv.gz'))
head(MIM)

# normalize data
cofactor <- 2
xdata <- asinh(MIM[,all_markers] / cofactor)

xdata2 <- xdata[,id_markers]
dim(xdata2)

#_____________________________________________
# run phenograph clustering k=100 ####
#_____________________________________________
pheno.out1 <- Rphenograph(xdata2, k = 100)
clusters1 <- as.numeric(membership(pheno.out1))
names(clusters1) <- rownames(xdata2)
length(unique(clusters1)) # 29

# heatmap check
plot_clustering_heatmap(expr = xdata2, cell_clustering = clusters1, fontsize = 10)
                       # filename = '~/WLData/MIMIC_IMC/Figures/Cell_type_clustering/heatmap_phenograph_k100.pdf')
