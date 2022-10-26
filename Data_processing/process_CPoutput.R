# process CellProfiler cell data output
# Author: LeeAnn Ramsay
# Date: 26 October 2022
library(data.table)
library(dplyr)

#_______________________________________
# load data ####
#_______________________________________
# cell mean intensities (output by CellProfiler 3_measure_mask_LR)
C.full <- as.data.frame(fread('~/WLData/MIMIC_IMC/analysis/cpout/cell.txt'))
# make mean intensity column names shorter
colnames(C.full) <- sub(pattern = 'Intensity_MeanIntensity_FullStack_', replacement = '', x = colnames(C.full))
head(C.full)
# extract required columns (remove duplicates)

col.keep <- c('ImageNumber','ObjectNumber','Metadata_acid','Metadata_acname',
              'FileName_FullStack','FileName_ProbabStack','FileName_cellmask',
              'AreaShape_Area','AreaShape_Eccentricity','Location_Center_X',
              'Location_Center_Y',paste0('c',seq(1,38,1)))
C1 <- C.full[,col.keep]

# load channel names
channels <- fread('~/WLData/MIMIC_IMC/analysis/cpout/panel.csv')
channels$`Tube Number` <- paste0('c',channels$`Tube Number`)

#_______________________________________
# format data ####
#_______________________________________
# rename mean intensity channels
C2 <- C1 %>%
  rename_with(.cols = channels$`Tube Number`, .fn = function(x) channels$Target[channels$`Tube Number` %in% x])
head(C2)

# scale back up the intensity values
C.scale <- C2[,channels$Target] * (2**16)
C.print <- cbind(C2[,!colnames(C2) %in% channels$Target],C.scale)

# extract and correct patient IDs
# TODO

#_____________________________________
# write mean intensity table ####
#_____________________________________
write.table(x = C.print,file = '~/WLData/MIMIC_IMC/Processed_image_data/Output_files/MIMIC_meanItensity_location_cell.tsv',
            sep='\t',row.names = FALSE, quote = FALSE)

gz <- gzfile("~/WLData/MIMIC_IMC/Processed_image_data/MIMIC_meanItensity_location_cell.tsv.gz", "w")
write.table(C.print, gz, quote = FALSE, row.names = FALSE, sep = '\t')
close(gz)
