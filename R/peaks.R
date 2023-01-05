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
cnr.byGrp <- lapply(cnr.byGrp, function(x) grp_qFilter(x, quantile = 0.75, operation = 'subsetByOverlaps'))

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
faire.byGrp <- lapply(faire.byGrp, function(x) grp_qFilter(x, quantile = 0.75, operation = 'subsetByOverlaps')) 

#make union FAIRE peak list - includes any/all peak calls for all experiments/reps
faire.union <- faire.byID %>%
  unlist() %>%
  reduce()


##### bind and annotate peak list #####

faire.df <- faire.union %>% data.frame() %>% dplyr::mutate(assay = 'faire')
cnr.df <- cnr.union %>% data.frame() %>% dplyr::mutate(assay = 'cnr')

peaks <- dplyr::bind_rows(faire.df, cnr.df) %>%
  dplyr::filter(!seqnames %in% c('chr4','chrY','chrM')) %>% #drop regions from potentially problematic chromosomes
  dplyr::arrange(seqnames, start, end) %>%
  dplyr::mutate(WT.3LW = ifelse(GRanges(.) %over% faire.byGrp$wt.3LW, T, F),
                WT.6h = ifelse(GRanges(.) %over% faire.byGrp$wt.6h, T, F),
                WT.24h = ifelse(GRanges(.) %over% faire.byGrp$wt.24h, T, F),
                WT.44h = ifelse(GRanges(.) %over% faire.byGrp$wt.44h, T, F),
                osaGFP.control = ifelse(GRanges(.) %over% faire.byGrp$osaGFP.deGrad_control_faire, T, F),
                osaGFP.deGrad = ifelse(GRanges(.) %over% faire.byGrp$osaGFP.deGrad_nubG4_faire, T, F),
                osa.cnr = ifelse(GRanges(.) %over% cnr.byGrp$osaGFP.sup, T, F),
                yw.cnr = ifelse(GRanges(.) %over% cnr.byGrp$yw.sup, T, F),
                cnr.peakCat = dplyr::case_when(osa.cnr & yw.cnr ~ 'shared',
                                               osa.cnr & !yw.cnr ~ 'osa',
                                               !osa.cnr & yw.cnr ~ 'control',
                                               T ~ 'other')) %>%
  dplyr::filter(WT.3LW | WT.6h | WT.24h | WT.24h | WT.44h | osaGFP.control | osaGFP.deGrad | osa.cnr | yw.cnr) %>% #drop regions that don't have a good quality reproducible peak called
  #get osaGFP deGrad cnr bigwig scores...
  dplyr::mutate(peak = 1:nrow(.),
                osaGFP.rep1.spikeNorm = get_bw_score(cnr.ss[cnr.ss$baseName == 'osaGFP-3LW-wing-aGFP-CnR-sup-rep1',]$bigwig_allFrags_sacCer3_spikeNorm[1], GRanges(.)),
                osaGFP.rep2.spikeNorm = get_bw_score(cnr.ss[cnr.ss$baseName == 'osaGFP-3LW-wing-aGFP-CnR-sup-rep2',]$bigwig_allFrags_sacCer3_spikeNorm[1], GRanges(.)),
                yw.rep1.spikeNorm = get_bw_score(cnr.ss[cnr.ss$baseName == 'yw-3LW-wing-aGFP-CnR-sup-rep1',]$bigwig_allFrags_sacCer3_spikeNorm[1], GRanges(.)),
                yw.rep2.spikeNorm = get_bw_score(cnr.ss[cnr.ss$baseName == 'yw-3LW-wing-aGFP-CnR-sup-rep2',]$bigwig_allFrags_sacCer3_spikeNorm[1], GRanges(.)),
                #get faire bigwig scores... 
                faire.3LW.zScore = get_bw_score(faire.wt.ssPool[faire.wt.ssPool$baseName == 'wt_3LW',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)),
                faire.6h.zScore = get_bw_score(faire.wt.ssPool[faire.wt.ssPool$baseName == 'wt_6h',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)),
                faire.24h.zScore = get_bw_score(faire.wt.ssPool[faire.wt.ssPool$baseName == 'wt_24h',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)),
                faire.44h.zScore = get_bw_score(faire.wt.ssPool[faire.wt.ssPool$baseName == 'wt_44h',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)))

#TODO - what's missing from these annotations?
# - differential enrichment data
# - closing vs opening


##### deseq #####

# CnR
cnr.Counts <- Rsubread::featureCounts(cnr.ss$bam,
                                      annot.ext = makeSAF(peaks), #trying all peaks - FAIRE and CnR - so it will be easy to combine with main peaks df
                                      allowMultiOverlap = T,
                                      nthreads = 4,
                                      isPairedEnd = T)

colnames(cnr.Counts$counts) <- cnr.ss$id

cnr.dds <- DESeq2::DESeqDataSetFromMatrix(countData = cnr.Counts$counts,  colData = cnr.ss,  design = ~grp) %>% DESeq2::DESeq(.)

cnr.dds.df <- DESeq2::results(cnr.dds, contrast = c('grp','osaGFP.sup','yw.sup')) %>% data.frame() 

names(cnr.dds.df) <- paste0('cnr.',names(cnr.dds.df))


# FAIRE 
faire.ss.WT <- faire.ss[faire.ss$experiment == 'WT FAIRE Wing Timecourse',]

faire.Counts <- Rsubread::featureCounts(faire.ss.WT$bam, #just use Bam from WT timecourse - single-end data
                                        annot.ext = makeSAF(peaks), #trying all peaks - FAIRE and CnR - so it will be easy to combine with main peaks df
                                        allowMultiOverlap = T,
                                        nthreads = 4,
                                        isPairedEnd = F)

colnames(faire.Counts$counts) <- faire.ss.WT$id

faire.dds <- DESeq2::DESeqDataSetFromMatrix(countData = faire.Counts$counts,  colData = faire.ss.WT,  design = ~grp) %>%
  DESeq2::DESeq(.)

faire.dds.df <- DESeq2::results(faire.dds, contrast = c('grp','wt.3LW','wt.24h')) %>% data.frame() 

names(faire.dds.df) <- paste0('faire_3LW_24h.',names(faire.dds.df))

#####

peaks <- dplyr::bind_cols(peaks, cnr.dds.df, faire.dds.df)
# 