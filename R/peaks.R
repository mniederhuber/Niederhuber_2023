
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
dm6.500bp <- purrr::map(seqnames(dm6), function(x) {
  if(x %in% c('chr2L', 'chr2R', 'chr3L', 'chr3R', 'chrX')){
    lng <- length(dm6[[x]])
    seqnames <- rep_len(x, lng/500)    
    start <- seq(0, lng-500, by = 500)
    end <- seq(499, lng, by = 500)
    df <- data.frame(seqnames, start, end)
  }
}) %>% dplyr::bind_rows() %>% dplyr::mutate(assay = '500bp background')

##### load peaks #####

### get CUT&RUN peaks
#TODO - error here with getPeakData when running from command line?
cnr.peaks <- getPeakData(cnr.ss, by = 'id', narrowPeak_colname = 'peak_allFrags')
cnr.byID <- cnr.peaks %>% GRanges() %>% split(., mcols(.)$id)
cnr.byGrp <- cnr.peaks %>% GRanges() %>% split(., mcols(.)$grp)

# just using rep 1 coordinates for now
cnr.summits.1 <- read.table('Peaks/osaGFP-3LW-wing-aGFP-CnR-sup-rep1_dm6_trim_q5_dupsRemoved_allFrags_summits.bed', sep = '\t') %>%
  .[,1:3]
colnames(cnr.summits.1) <- c('seqnames', 'start','end')
cnr.summits <- GRanges(cnr.summits.1)

#cnr.summits.2 <- read.table('Peaks/osaGFP-3LW-wing-aGFP-CnR-sup-rep2_dm6_trim_q5_dupsRemoved_allFrags_summits.bed', sep = '\t') %>%
#  .[,1:3]
#colnames(cnr.summits.2) <- c('seqnames', 'start','end')

# find the summits that overlap - keeping rep1 coordinates
#cnr.summits <- subsetByOverlaps(GRanges(cnr.summits.1), GRanges(cnr.summits.2), maxgap = 100) 


#filter cnr peaks by q val >= 10 and called in each replicate
#TODO - decide on filter cutoff -- 75% was pretty strict, only ~1000 osa specific peaks vs previous approach with no filtering that resulting in ~2700 specific peaks
cnr.byGrp <- lapply(cnr.byGrp, function(x) grp_qFilter(x, quantile = 0.5, operation = 'subsetByOverlaps'))


#make union cnr peak list
cnr.union <- cnr.byID %>%
  unlist() %>%
  reduce()

### get FAIRE peaks
faire.ss.WT <- faire.ss[faire.ss$experiment == 'WT FAIRE Wing Timecourse',]
faire.ss.osaDeg <- faire.ss[faire.ss$experiment == 'osaGFP deGrad FAIRE',]

faire.wt.peaks <- getPeakData(faire.ss.WT, by = 'id', narrowPeak_colname = 'peaks')
faire.wt.byID <- faire.wt.peaks %>% GRanges() %>% split(., mcols(.)$id) #split by replicates
faire.wt.byGrp <- faire.wt.peaks %>% GRanges() %>% split(., mcols(.)$grp) #split by pooled replicates



faire.osaDeg.peaks <- getPeakData(faire.ss.osaDeg, by = 'id', narrowPeak_colname = 'peaks')
faire.osaDeg.byID <- faire.osaDeg.peaks %>% GRanges() %>% split(., mcols(.)$id) #split by replicates
faire.osaDeg.byGrp <- faire.osaDeg.peaks %>% GRanges() %>% split(., mcols(.)$grp) #split by pooled replicates

#filter by replicate specific qValue and only peaks that intersect between replicates
faire.wt.byGrp <- lapply(faire.wt.byGrp, function(x) grp_qFilter(x, quantile = 0.75, operation = 'subsetByOverlaps', with_reduce = T))
faire.osaDeg.byGrp <- lapply(faire.osaDeg.byGrp, function(x) grp_qFilter(x, quantile = 0.5, operation = 'subsetByOverlaps', with_reduce = T))

faire.wt.3LW.summits <- read.table('Peaks/wt_3LW_dm6_trim_q5_sorted_dupsRemoved_POOLED_summits.bed', sep = '\t') %>%
  .[,1:3]
colnames(faire.wt.3LW.summits) <- c('seqnames', 'start', 'end')
faire.wt.3LW.summits <- subsetByOverlaps(GRanges(faire.wt.3LW.summits), faire.wt.byGrp$wt.3LW)
                                         
#make union FAIRE peak list - includes any/all peak calls for all wt reps
faire.wt.union <- faire.wt.byID %>%
  unlist() %>%
  reduce()

faire.osaDeg.union <- faire.osaDeg.byID %>%
  unlist() %>%
  reduce()

### get Rotund peaks
rn.ss <- read.csv('~/McKay/Rotund-ChIP/sampleSheet.tsv', sep = '\t') 

rn.peak <- getPeakData(rn.ss, narrowPeak_colname = 'peakFile', by = 'grp') 

q <- quantile(rn.peak$qValue, probs = 0.25) #find the 75% qVal

rn.peak.qf <- rn.peak %>%
  dplyr::mutate(seqnames = paste0('chr',seqnames)) %>%
  dplyr::filter(qValue >= q & seqnames %in% c('chr2L','chr2R','chr3L','chr3R','chrX')) 

#grp.repShared.1 <- IRanges::subsetByOverlaps(grp.list[[1]], grp.list[[2]]) # get coords in 1st rep that overlap rep 2
rn.peak.1 <- IRanges::subsetByOverlaps(GRanges(rn.peak.qf[rn.peak.qf$rep == 'rep1',]), 
                                       GRanges(rn.peak.qf[rn.peak.qf$rep == 'rep2',]))

rn.peak.2 <- IRanges::subsetByOverlaps(GRanges(rn.peak.qf[rn.peak.qf$rep == 'rep2',]), 
                                       GRanges(rn.peak.qf[rn.peak.qf$rep == 'rep1',]))                                      
rn.union <- union(rn.peak.1, rn.peak.2)

##### bind and annotate peak list #####

faire.wt.df <- faire.wt.union %>% data.frame() %>% dplyr::mutate(assay = 'faire', experiment = 'WT FAIRE Wing Timecourse')
faire.osaDeg.df <- faire.osaDeg.union %>% data.frame() %>% dplyr::mutate(assay = 'faire', experiment = 'osaGFP deGrad Pupal Wing FAIRE')
cnr.df <- cnr.union %>% data.frame() %>% dplyr::mutate(assay = 'cnr', experiment = 'osaGFP 3LW Wing CUT&RUN')
rn.df <- rn.union %>% data.frame() %>% dplyr::mutate(assay = 'chip', experiment = 'Rotund 3LW Wing ChIP-seq')

#TODO rename peaks at each step?
peaks <- purrr::map(list(faire.wt.df, faire.osaDeg.df, cnr.df, rn.df), function(x) {
  x %>% #dplyr::bind_rows(faire.wt.df, faire.osaDeg.df, cnr.df) %>%
    dplyr::filter(!seqnames %in% c('chr4','chrY','chrM')) %>% #drop regions from potentially problematic chromosomes
    dplyr::arrange(seqnames, start, end) %>%
    dplyr::mutate(WT.3LW = ifelse(GRanges(.) %over% faire.wt.byGrp$wt.3LW, T, F),
                  WT.6h = ifelse(GRanges(.) %over% faire.wt.byGrp$wt.6h, T, F),
                  WT.18h = ifelse(GRanges(.) %over% faire.wt.byGrp$wt.18h, T, F),
                  WT.24h = ifelse(GRanges(.) %over% faire.wt.byGrp$wt.24h, T, F),
                  WT.36h = ifelse(GRanges(.) %over% faire.wt.byGrp$wt.36h, T, F),
                  WT.44h = ifelse(GRanges(.) %over% faire.wt.byGrp$wt.44h, T, F),
                  osaDeg.control = ifelse(GRanges(.) %over% faire.osaDeg.byGrp$osaGFP.deGrad_control_faire, T, F),
                  osaDeg = ifelse(GRanges(.) %over% faire.osaDeg.byGrp$osaGFP.deGrad_nubG4_faire, T, F),
                  osa.cnr = ifelse(GRanges(.) %over% cnr.byGrp$osaGFP.sup, T, F),
                  yw.cnr = ifelse(GRanges(.) %over% cnr.byGrp$yw.sup, T, F),
                  cnr.peakCat = dplyr::case_when(osa.cnr & yw.cnr ~ 'shared',
                                                 osa.cnr & !yw.cnr ~ 'osa',
                                                 !osa.cnr & yw.cnr ~ 'control',
                                                 T ~ 'other')) %>%
  #  dplyr::filter(WT.3LW | WT.6h | WT.18h | WT.24h | WT.36h | WT.44h | osaDeg.control | osaDeg | osa.cnr | yw.cnr) %>% #drop regions that don't have a good quality reproducible peak called
    #get osaGFP deGrad cnr bigwig scores...
    ## NOTE you still have to filter for each individual experiment, ie. osaDeg.control | osaDeg to just get good peaks from the degrad experiment
    dplyr::mutate(peak = 1:nrow(.),
                  osaGFP.rep1.spikeNorm = get_bw_score(cnr.ss[cnr.ss$baseName == 'osaGFP-3LW-wing-aGFP-CnR-sup-rep1',]$bigwig_allFrags_sacCer3_spikeNorm[1], GRanges(.)),
                  osaGFP.rep2.spikeNorm = get_bw_score(cnr.ss[cnr.ss$baseName == 'osaGFP-3LW-wing-aGFP-CnR-sup-rep2',]$bigwig_allFrags_sacCer3_spikeNorm[1], GRanges(.)),
                  yw.rep1.spikeNorm = get_bw_score(cnr.ss[cnr.ss$baseName == 'yw-3LW-wing-aGFP-CnR-sup-rep1',]$bigwig_allFrags_sacCer3_spikeNorm[1], GRanges(.)),
                  yw.rep2.spikeNorm = get_bw_score(cnr.ss[cnr.ss$baseName == 'yw-3LW-wing-aGFP-CnR-sup-rep2',]$bigwig_allFrags_sacCer3_spikeNorm[1], GRanges(.)),
                  #cnr zNorm...
                  osaGFP.rep1.zScore = get_bw_score(cnr.ss[cnr.ss$baseName == 'osaGFP-3LW-wing-aGFP-CnR-sup-rep1',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)),
                  osaGFP.rep2.zScore = get_bw_score(cnr.ss[cnr.ss$baseName == 'osaGFP-3LW-wing-aGFP-CnR-sup-rep2',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)),
                  yw.rep1.zScore = get_bw_score(cnr.ss[cnr.ss$baseName == 'yw-3LW-wing-aGFP-CnR-sup-rep1',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)),
                  yw.rep2.zScore = get_bw_score(cnr.ss[cnr.ss$baseName == 'yw-3LW-wing-aGFP-CnR-sup-rep2',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)),
                  #get faire bigwig scores... 
                  faire.3LW.zScore = get_bw_score(faire.wt.ssPool[faire.wt.ssPool$baseName == 'wt_3LW',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)),
                  faire.6h.zScore = get_bw_score(faire.wt.ssPool[faire.wt.ssPool$baseName == 'wt_6h',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)),
                  faire.18h.zScore = get_bw_score(faire.wt.ssPool[faire.wt.ssPool$baseName == 'wt_18h',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)),
                  faire.24h.zScore = get_bw_score(faire.wt.ssPool[faire.wt.ssPool$baseName == 'wt_24h',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)),
                  faire.36h.zScore = get_bw_score(faire.wt.ssPool[faire.wt.ssPool$baseName == 'wt_36h',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)),
                  faire.44h.zScore = get_bw_score(faire.wt.ssPool[faire.wt.ssPool$baseName == 'wt_44h',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)))
}) %>%
  append(GRanges(dm6.500bp)) %>%
  setNames(c('faire.wt.df', 'faire.osaDeg.df', 'cnr.df','rn.df','dm6.500bp'))



#DONE - what's missing from these annotations?
# - deseq annotation
# - closing vs opening


##### Run deSeq2 #####

# CnR

cnr.Counts <- Rsubread::featureCounts(cnr.ss$bam,
                                      annot.ext = makeSAF(peaks$cnr.df), #trying all peaks - FAIRE and CnR - so it will be easy to combine with main peaks df
                                      allowMultiOverlap = T,
                                      nthreads = 4,
                                      isPairedEnd = T)

colnames(cnr.Counts$counts) <- cnr.ss$id

cnr.dds <- DESeq2::DESeqDataSetFromMatrix(countData = cnr.Counts$counts,  colData = cnr.ss,  design = ~grp) %>% DESeq2::DESeq(.)

cnr.dds.df <- DESeq2::results(cnr.dds, contrast = c('grp','osaGFP.sup','yw.sup')) %>% data.frame() 

#names(cnr.dds.df) <- paste0('cnr.',names(cnr.dds.df))

peaks$cnr.df %<>% dplyr::bind_cols(., cnr.dds.df)

# FAIRE 

faire.Counts <- Rsubread::featureCounts(faire.ss.WT$bam, #just use Bam from WT timecourse - single-end data
                                        annot.ext = makeSAF(peaks$faire.wt.df), #trying all peaks - FAIRE and CnR - so it will be easy to combine with main peaks df
                                        allowMultiOverlap = T,
                                        nthreads = 4,
                                        isPairedEnd = F)

colnames(faire.Counts$counts) <- faire.ss.WT$id

faire.dds <- DESeq2::DESeqDataSetFromMatrix(countData = faire.Counts$counts,  colData = faire.ss.WT,  design = ~grp) %>%
  DESeq2::DESeq(.)

faire.dds.df <- DESeq2::results(faire.dds, contrast = c('grp','wt.3LW','wt.24h')) %>% data.frame() 

names(faire.dds.df) <- paste0('faire_3LW_24h.',names(faire.dds.df))

peaks$faire.wt.df %<>% dplyr::bind_cols(., faire.dds.df)
# FAIRE - OsaGFP deGrad

#TODO annotate all the diffbind steps - cite source
#### Diff bind
osaDeg.db_sheet  <- faire.ss.osaDeg %>%
  dplyr::rename(SampleID = id,
                Peaks = peaks,
                bamReads = bam) %>%
  dplyr::mutate(PeakCaller = "narrow", 
                Replicate = ifelse(rep == 'Rep1', 1, 2),
                Condition = ifelse(grepl('control', SampleID), 'control','nubG4'))

osaDeg.db <- DiffBind::dba(sampleSheet = osaDeg.db_sheet)
osaDeg.db <- DiffBind::dba.count(osaDeg.db)

info <- DiffBind::dba.show(osaDeg.db)
libsizes <- cbind(LibReads=info$Reads,
                  FRiP=info$FRiP,
                  PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID 

osaDeg.db <- DiffBind::dba.normalize(osaDeg.db, normalize = 'lib') #library size normalization may be more appropriate for non-rnaseq data
norm <- DiffBind::dba.normalize(osaDeg.db, bRetrieve=T)
normlibs <- cbind(FullLibSize=norm$lib.sizes,
                  NormFacs=norm$norm.factors,
                  NormLibSize=round(norm$lib.sizes/norm$norm.factors))
rownames(normlibs) <- info$ID


osaDeg.db <- DiffBind::dba.contrast(osaDeg.db, contrast = c('Condition', 'nubG4','control'))

osaDeg.db <- DiffBind::dba.analyze(osaDeg.db)

dds <- DiffBind::dba.analyze(osaDeg.db, bRetrieveAnalysis =TRUE) #retreive the deseq2 object

dds.results <- DESeq2::results(dds, contrast = c('Condition', 'nubG4', 'control'))

#### 
osaDeg.db.peaks <- osaDeg.db$peaks %>%
  data.frame() %>%
  dplyr::select(seqnames, start, end) %>%
  dplyr::mutate(peak = 1:nrow(.)) 

dds.peaks <- dds.results %>% data.frame() %>%
  dplyr::mutate(peak = 1:nrow(.)) %>%
  dplyr::left_join(osaDeg.db.peaks, ., by = 'peak') 

# final output
#osaDeg.diffBind <- IRanges::subsetByOverlaps(GRanges(dds.peaks), ) %>%
#  data.frame() %>%
#  GRanges() 

## behaves differently than filtering the annotated DF on osaDeg | osaDeg.control - slight difference in peak number
#osaDeg.qF <- IRanges::subsetByOverlaps(GRanges(dds.peaks), GRanges(peaks$faire.osaDeg.df %>% dplyr::filter(osaDeg | osaDeg.control))) 
osaDeg.qF <- dds.peaks %>%
  dplyr::mutate(osaDeg.control = ifelse(GRanges(.) %over% faire.osaDeg.byGrp$osaGFP.deGrad_control_faire, T, F),
                osaDeg = ifelse(GRanges(.) %over% faire.osaDeg.byGrp$osaGFP.deGrad_nubG4_faire, T, F)) %>%
  dplyr::filter(osaDeg | osaDeg.control)

#TODO move to plots?
osaDeg.osaDep <- osaDeg.qF %>%
  data.frame() %>%
  dplyr::filter(log2FoldChange <= -1) %>%
  GRanges()

osaDep.anno <- osaDeg.osaDep %>%
  data.frame() %>%
  dplyr::mutate(#group = ifelse(log2FoldChange <= -1, "Osa Dependent", "Independent"),
                faire.3LW.zScore = get_bw_score(faire.wt.ssPool[faire.wt.ssPool$baseName == 'wt_3LW',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)),
                faire.6h.zScore = get_bw_score(faire.wt.ssPool[faire.wt.ssPool$baseName == 'wt_6h',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)),
                  faire.18h.zScore = get_bw_score(faire.wt.ssPool[faire.wt.ssPool$baseName == 'wt_18h',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)),
                  faire.24h.zScore = get_bw_score(faire.wt.ssPool[faire.wt.ssPool$baseName == 'wt_24h',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)),
                  faire.36h.zScore = get_bw_score(faire.wt.ssPool[faire.wt.ssPool$baseName == 'wt_36h',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)),
                  faire.44h.zScore = get_bw_score(faire.wt.ssPool[faire.wt.ssPool$baseName == 'wt_44h',]$bigwig_rpgcNorm_zNorm[1], GRanges(.)))

# add an osa.dependent annotation to the main osa deg peak set, that way the full FAIRE coordinates can be used for gviz tracks
peaks$faire.osaDeg.df %<>%
  dplyr::mutate(osa.dependent = ifelse(GRanges(.) %over% osaDeg.osaDep, T, F))
  
######

##### bind deseq results to main peaks dataframe and annotate faire behavior between 3LW and 24hAPF #####

# closing --> must be from FAIRE data, must be called a reproducible peak at 3LW, must have a log2(3LW/24h) >= 1 
# static --> must be from FAIRE data, no peak call requirement, must have -1 > log2(3LW/24h) < 1
# opening --> must be from FAIRE data, must be a reproducible peak call at 24h, must have a log2(3LW/24h) <= -1

#TODO - currently gives ~500 closing, ~1100 opening, ~5400 static -- notably ~50% fewer closing regions than I or SLN previously annotated

     
 # TODO - breaking change! withou this annotation here it will need to be done at the peak annotation step in plot notebook for fig 3F
  #dplyr::bind_rows(., dm6.500bp) # testing binding in the 1kb windowed background intervals at this step -- don't have any annotation
                               # that way these regions will also get annotations

#### 
# annotate peaks with dm6.TxDb 
####
peaks <- purrr::map(peaks, function(x) {
  anno <- x %>%
    GRanges() %>%
    ChIPseeker::annotatePeak(., tssRegion = c(-100,500), TxDb=dm6.TxDb, level = 'gene',  annoDb = "org.Dm.eg.db")

  peaks.anno <- ChIPseeker::as.GRanges(anno) %>% 
    data.frame() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(anno.new = dplyr::case_when(grepl("3' UTR", annotation) ~ "3' UTR",
                                              grepl("5' UTR", annotation) ~ "5' UTR",
                                              grepl("Distal Intergenic", annotation) ~ "Distal Intergenic",
                                              grepl("Downstream", annotation) ~ "Downstream", 
                                              grepl("Exon", annotation) ~ "Exon",
                                              grepl("Intron", annotation) ~ "Intron",
                                              grepl("Promoter \\(2-3kb\\)", annotation) ~ "Promoter (2-3kb)", #why is this not working??? 
                                                ## It looks like it's because nothing gets annoated as "Promoter (2-3kb)" maybe because of how promoter region is defined?
                                              grepl("Promoter", annotation) ~ "Promoter"))

#                  flank_gene = stringr::str_split(flank_geneIds, ';'),
#                  flank_distance = stringr::str_split(flank_gene_distances, ';')) 
})

message('Saving output...')
save(peaks, 
     faire.wt.byGrp, 
     faire.osaDeg.byGrp, 
     faire.osaDeg.byID, 
     osaDeg.db.peaks,
     dds.peaks,
     osaDeg.qF, 
     osaDeg.osaDep, 
     osaDep.anno,
     cnr.byGrp, 
     cnr.byID, 
     cnr.summits,
     faire.wt.3LW.summits,
     file = 'rData/peaks.rda')
  
# 