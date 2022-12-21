library(magrittr)
library(ggplot2)
library(GenomicRanges)
source('~/NystLib/R/GRanges_methods.R')
source('~/NystLib/R/peakUtils.R')

load('rData/sheets.rda')

### get CUT&RUN peaks
cnr.peaks <- getPeakData(cnr.ss, by = 'id', narrowPeak_colname = 'peak_allFrags')
cnr.byID <- cnr.peaks %>% GRanges() %>% split(., mcols(.)$id)

cnr.union <- cnr.byID %>%
  unlist() %>%
  reduce()

### 
