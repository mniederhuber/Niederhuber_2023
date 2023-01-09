### 
# load packages 
### 

library(magrittr)
library(ggplot2)

library(org.Dm.eg.db)
library(AnnotationDbi)

library(EnrichedHeatmap)
library(ComplexHeatmap)

load('rData/peaks.rda')
load('rData/sheets.rda')

###
# assign variables
###

dm6 <- BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6
dm6.TxDb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene

brD <- data.frame('seqnames' = 'chrX', 'start' = 1565708, 'end' = 1567401) %>%
  GenomicRanges::GRanges()

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
viridis.hex <- c("#440154","#3b528b","#21918c","#5ec962","#5ec962","#fde725")

###
# Fig 3A - BigWig signal in osa specific peak calls 
###

bws <- purrr::map(cnr.ss$id, function(x) {
  fp <- cnr.ss[cnr.ss['id'] == x,]$bigwig_allFrags_sacCer3_spikeNorm
  rtracklayer::import.bw(fp)
})

names(bws) <- cnr.ss$id

osaPeaks <- peaks %>% dplyr::filter(osa.cnr & !yw.cnr & assay == 'cnr') %>% GRanges() %>% resize(width = 1, fix = "center")

mats <- purrr::map(bws, function(x) {
  normalizeToMatrix(x, osaPeaks, value_column = 'score', keep = c(0.01, 0.99), extend = 2000)
}) 
names(mats) <- names(bws)

common_min = min(unlist(mats))
common_max = max(unlist(mats))

col_fun = circlize::colorRamp2(seq(common_min, common_max, length.out = 6), viridis.hex)

osa.count <- osaPeaks %>%
  data.frame() %>%
  nrow()

osa.count[1] <- as.character(osa.count[1])

names(osa.count) <- 'Osa Peak'

groups <- ifelse(osaPeaks$osa.cnr, 'Osa Peak', NA)

#TODO - I don't think there's a good way to functionalize EnrichedHeatmap() because of the unique features of the first map and annotation...
hm <- EnrichedHeatmap(mats$osaGFP.sup.rep1,
                      axis_name = c(-2000, 0, 2000),
                      pos_line = F,
                      column_title = 'Osa Rep1',
                      col = col_fun,
                      name = 'spikeNorm Score',
                      top_annotation = HeatmapAnnotation(enriched = anno_enriched(axis_param = list(side = 'left', facing = 'outside'), ylim = c(0,5))),
                      left_annotation = c(rowAnnotation(textbox = anno_textbox(align_to = groups,
                                                                               as.list(osa.count),
                                                                               side = 'left',
                                                                               background_gp = gpar(fill = NA, col = NA),
                                                                               gp = gpar(col = 'black'), 
                                                                               padding = unit(0, 'mm'))))) +
  EnrichedHeatmap(mats$osaGFP.sup.rep2,
                  axis_name = c(-2000, 0, 2000),
                  pos_line = F,
                  column_title = 'Osa Rep2',
                  col = col_fun, 
                  show_heatmap_legend = FALSE,
                  top_annotation = HeatmapAnnotation(enriched = anno_enriched(axis = F, ylim = c(0,5)))) +
  EnrichedHeatmap(mats$yw.sup.rep1,
                  axis_name = c(-2000, 0, 2000),
                  pos_line = F,
                  column_title = 'Control Rep1',
                  col = col_fun, 
                  show_heatmap_legend = FALSE,
                  top_annotation = HeatmapAnnotation(enriched = anno_enriched(axis = F, ylim = c(0,5)))) +  
  EnrichedHeatmap(mats$yw.sup.rep2,
                  axis_name = c(-2000, 0, 2000),
                  pos_line = F,
                  column_title = 'Control Rep2',
                  col = col_fun, 
                  show_heatmap_legend = FALSE,
                  top_annotation = HeatmapAnnotation(enriched = anno_enriched(axis = F, ylim = c(0,5))))

png('rPlots/bwSignal_OsaPeaks.png', width = 6, height = 6, units = 'in', res = 300)
draw(hm, row_title = NULL, merge_legend = T)
dev.off()


###
# Fig 3B - where are osa peaks distributed?
###

peaks %>%
  dplyr::filter(assay == 'cnr' &  osa.cnr & !yw.cnr | assay == '1kb background') %>%
  dplyr::group_by(anno.new, assay) %>%
  dplyr::count() %>%
  ggplot(aes(x = assay, y = n, fill = anno.new)) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_fill_manual(values = cbPalette) +
  ylab('Fraction') +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),)
ggsave('rPlots/osaPeak_annotation.png', width = 4, height = 4, units = 'in', dpi = 300)

###
# Fig 3C - browser shots
## 




