library(magrittr)
library(ggplot2)
library(GenomicRanges)
source('~/NystLib/R/GRanges_methods.R')
source('~/NystLib/R/peakUtils.R')
source('utils.R')

load('rData/sheets.rda')

### get CUT&RUN peaks
cnr.peaks <- getPeakData(cnr.ss, by = 'id', narrowPeak_colname = 'peak_allFrags')
cnr.byID <- cnr.peaks %>% GRanges() %>% split(., mcols(.)$id)

cnr.union <- cnr.byID %>%
  unlist() %>%
  reduce()

### get FAIRE peaks
faire.peaks <- getPeakData(faire.ss, by = 'id', narrowPeak_colname = 'peaks')
faire.byID <- faire.peaks %>% GRanges() %>% split(., mcols(.)$id) #split by replicates
faire.byGrp <- faire.peaks %>% GRanges() %>% split(., mcols(.)$Grp) #split by pooled replicates
faire.byExp <- faire.peaks %>% GRanges() %>% split(., mcols(.)$experiment) #split by experiment


#filter by replicate specific qValue
faire.byGrp <- lapply(faire.byGrp, function(x) grp_qFilter(x)) %>% GRangesList()

faire.union <- faire.byID %>%
  unlist() %>%
  reduce()


### build annotated FAIRE peak dataframe
#TODO - split out osaDegrad reps
#     - filter FAIRe peaks by qVal - take top X peaks
#     - write function to filter by qVal and keep on reproducible peaks for each 'grp'
#What parameters to select "best" most confident peak calls?
#-- reproducible, high confidence calls, 


faire.union %>%
  data.frame() %>%
  dplyr::mutate(peak = 1:nrow(.),
                assay = 'FAIRE',
                WT.3LW = ifelse(GRanges(.) %over% faire.byGrp$wt.3LW, T, F),
                WT.6h = ifelse(GRanges(.) %over% faire.byGrp$wt.6h, T, F),
                WT.24h = ifelse(GRanges(.) %over% faire.byGrp$wt.24h, T, F),
                WT.44h = ifelse(GRanges(.) %over% faire.byGrp$wt.44h, T, F),
                osaGFP.control = ifelse(GRanges(.) %over% faire.byGrp$osaGFP.deGrad_control_faire, T, F),
                osaGFP.deGrad = ifelse(GRanges(.) %over% faire.byGrp$osaGFP.deGrad_nubG4_faire, T, F),
                )