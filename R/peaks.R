
library(magrittr)
library(ggplot2)
library(GenomicRanges)
library(org.Dm.eg.db)
library(AnnotationDbi)
source('~/NystLib/R/GRanges_methods.R')
source('~/NystLib/R/peakUtils.R') #TODO move/rewrite getPeakData from NystLib to utils.R
source('utils.R')

load('rData/sheets.rda')



dm6 <- BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6
dm6.TxDb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene



#make 1kb bin whole genomic granges
dm6.1kb <- purrr::map(seqnames(dm6), function(x) {
  if(x %in% c('chr2L', 'chr2R', 'chr3L', 'chr3R', 'chrX')){
    lng <- length(dm6[[x]])
    seqnames <- rep_len(x, lng/1000)    
    start <- seq(0, lng-1000, by = 1000)
    end <- seq(999, lng, by = 1000)
    df <- data.frame(seqnames, start, end)
  }
}) %>% dplyr::bind_rows() %>% dplyr::mutate(assay = '1kb background')


##### load peaks #####

### get CUT&RUN peaks
#TODO - error here with getPeakData when running from command line?
cnr.peaks <- getPeakData(cnr.ss, by = 'id', narrowPeak_colname = 'peak_allFrags')
cnr.byID <- cnr.peaks %>% GRanges() %>% split(., mcols(.)$id)
cnr.byGrp <- cnr.peaks %>% GRanges() %>% split(., mcols(.)$grp)

#filter cnr peaks by q val >= 10 and called in each replicate
#TODO - decide on filter cutoff -- 75% was pretty strict, only ~1000 osa specific peaks vs previous approach with no filtering that resulting in ~2700 specific peaks
cnr.byGrp <- lapply(cnr.byGrp, function(x) grp_qFilter(x, quantile = 0.5, operation = 'subsetByOverlaps'))

#make union cnr peak list
cnr.union <- cnr.byID %>%
  unlist() %>%
  reduce()

### get FAIRE peaks
faire.peaks <- getPeakData(faire.ss, by = 'id', narrowPeak_colname = 'peaks')
faire.byID <- faire.peaks %>% GRanges() %>% split(., mcols(.)$id) #split by replicates
faire.byGrp <- faire.peaks %>% GRanges() %>% split(., mcols(.)$grp) #split by pooled replicates
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

#DONE - what's missing from these annotations?
# - deseq annotation
# - closing vs opening


##### Run deSeq2 #####

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

##### bind deseq results to main peaks dataframe and annotate faire behavior between 3LW and 24hAPF #####

# closing --> must be from FAIRE data, must be called a reproducible peak at 3LW, must have a log2(3LW/24h) >= 1 
# static --> must be from FAIRE data, no peak call requirement, must have -1 > log2(3LW/24h) < 1
# opening --> must be from FAIRE data, must be a reproducible peak call at 24h, must have a log2(3LW/24h) <= -1

#TODO - currently gives ~500 closing, ~1100 opening, ~5400 static -- notably ~50% fewer closing regions than I or SLN previously annotated

peaks <- dplyr::bind_cols(peaks, cnr.dds.df, faire.dds.df) %>%
  dplyr::mutate(faireCat.3LW_24h = dplyr::case_when(faire_3LW_24h.log2FoldChange >= 1 ~ 'closing',
                                                    faire_3LW_24h.log2FoldChange < 1 & faire_3LW_24h.log2FoldChange > -1 ~ 'static',
                                                    faire_3LW_24h.log2FoldChange <= -1 ~ 'opening',
                                                    T ~ 'NA')) %>% 
  dplyr::bind_rows(., dm6.1kb) # testing binding in the 1kb windowed background intervals at this step -- don't have any annotation
                               # that way these regions will also get annotations

#### 
# annotate peaks with dm6.TxDb 
####

peaks.anno <- peaks %>%
  GRanges() %>%
  ChIPseeker::annotatePeak(., TxDb=dm6.TxDb, annoDb = "org.Dm.eg.db")

peaks <- ChIPseeker::as.GRanges(peaks.anno) %>% 
  data.frame() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(anno.new = dplyr::case_when(grepl("3' UTR", annotation) ~ "3' UTR",
                                            grepl("5' UTR", annotation) ~ "5' UTR",
                                            grepl("Distal Intergenic", annotation) ~ "Distal Intergenic",
                                            grepl("Downstream", annotation) ~ "Downstream", 
                                            grepl("Exon", annotation) ~ "Exon",
                                            grepl("Intron", annotation) ~ "Intron",
                                            grepl("Promoter \\(2-3kb\\)", annotation) ~ "Promoter (2-3kb)", #why is this not working???
                                            grepl("Promoter", annotation) ~ "Promoter"))

save(peaks, faire.byGrp, cnr.byGrp, file = 'rData/peaks.rda')
  
# 