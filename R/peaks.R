library(magrittr)
library(ggplot2)
library(GenomicRanges)
source('~/NystLib/R/GRanges_methods.R')
source('~/NystLib/R/peakUtils.R')
source('utils.R')

load('rData/sheets.rda')

##### load peaks #####

### get CUT&RUN peaks
cnr.peaks <- getPeakData(cnr.ss, by = 'id', narrowPeak_colname = 'peak_allFrags')
cnr.byID <- cnr.peaks %>% GRanges() %>% split(., mcols(.)$id)
cnr.byGrp <- cnr.peaks %>% GRanges() %>% split(., mcols(.)$Grp)

#filter cnr peaks by q val >= 10 and called in each replicate
cnr.byGrp <- lapply(cnr.byGrp, function(x) grp_qFilter(x, q = 10))

#make union cnr peak list
cnr.union <- cnr.byID %>%
  unlist() %>%
  reduce()

### get FAIRE peaks
faire.peaks <- getPeakData(faire.ss, by = 'id', narrowPeak_colname = 'peaks')
faire.byID <- faire.peaks %>% GRanges() %>% split(., mcols(.)$id) #split by replicates
faire.byGrp <- faire.peaks %>% GRanges() %>% split(., mcols(.)$Grp) #split by pooled replicates
faire.byExp <- faire.peaks %>% GRanges() %>% split(., mcols(.)$experiment) #split by experiment

#filter by replicate specific qValue and only peaks that intersect between replicates
faire.byGrp <- lapply(faire.byGrp, function(x) grp_qFilter(x, quantile = 0.75)) 

#make union FAIRE peak list - includes any/all peak calls for all experiments/reps
faire.union <- faire.byID %>%
  unlist() %>%
  reduce()



##### annotate peak lists #####

### build annotated FAIRE peak dataframe
faire.df <- faire.union %>%
  data.frame() %>%
  dplyr::filter(!seqnames %in% c('chr4','chrY','chrM')) %>% #drop regions from potentially problematic chromosomes
  dplyr::mutate(peak = 1:nrow(.),
                assay = 'FAIRE',
                WT.3LW = ifelse(GRanges(.) %over% faire.byGrp$wt.3LW, T, F),
                WT.6h = ifelse(GRanges(.) %over% faire.byGrp$wt.6h, T, F),
                WT.24h = ifelse(GRanges(.) %over% faire.byGrp$wt.24h, T, F),
                WT.44h = ifelse(GRanges(.) %over% faire.byGrp$wt.44h, T, F),
                osaGFP.control = ifelse(GRanges(.) %over% faire.byGrp$osaGFP.deGrad_control_faire, T, F),
                osaGFP.deGrad = ifelse(GRanges(.) %over% faire.byGrp$osaGFP.deGrad_nubG4_faire, T, F)) %>%
  dplyr::filter(WT.3LW | WT.6h | WT.24h | WT.24h | WT.44h | osaGFP.control | osaGFP.deGrad) %>% #drop regions that don't have a good quality reproducible peak called
#get osaGFP deGrad cnr bigwig scores...
  dplyr::mutate(osaGFP.rep1.spikeNorm = get_bw_score(cnr.ss[cnr.ss$baseName == 'osaGFP-3LW-wing-aGFP-CnR-sup-rep1',]$bigwig_allFrags_sacCer3_spikeNorm[1], GRanges(.)),
                osaGFP.rep2.spikeNorm = get_bw_score(cnr.ss[cnr.ss$baseName == 'osaGFP-3LW-wing-aGFP-CnR-sup-rep2',]$bigwig_allFrags_sacCer3_spikeNorm[1], GRanges(.)),
                yw.rep1.spikeNorm = get_bw_score(cnr.ss[cnr.ss$baseName == 'yw-3LW-wing-aGFP-CnR-sup-rep1',]$bigwig_allFrags_sacCer3_spikeNorm[1], GRanges(.)),
                yw.rep2.spikeNorm = get_bw_score(cnr.ss[cnr.ss$baseName == 'yw-3LW-wing-aGFP-CnR-sup-rep2',]$bigwig_allFrags_sacCer3_spikeNorm[1], GRanges(.)),
#get faire bigwig scores... 
                faire.3LW.zScore = get_bw_score(faire.wt.ssPool[faire.wt.ssPool$baseName == 'wt_3LW',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)),
                faire.6h.zScore = get_bw_score(faire.wt.ssPool[faire.wt.ssPool$baseName == 'wt_6h',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)),
                faire.24h.zScore = get_bw_score(faire.wt.ssPool[faire.wt.ssPool$baseName == 'wt_24h',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)),
                faire.44h.zScore = get_bw_score(faire.wt.ssPool[faire.wt.ssPool$baseName == 'wt_44h',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)))
