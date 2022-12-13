### load packages
library(magrittr)
#library(GenomicRanges)

#source('scripts/utils.R')

### prep sample sheets
print('cnr.ss')
cnr.ss <- read.csv('sheets/osaGFP-CnR-sampleSheet.tsv', sep = '\t') %>%
  dplyr::mutate(genotype = stringr::str_split_fixed(sample, '-', n = 6)[,1],
                fraction = stringr::str_split_fixed(sample, '-', n = 6)[,6], 
                id = paste(genotype, fraction, rep, sep = '.'),
                grp = paste(genotype, fraction, sep = '.'), .before = 1)

print('idr.ss')
#TODO -- decide if IDR peak calling is useful. Currently not using the IDR scores at all, just shared peaks.
#ideally have only a single cnr sample sheet, so at least merge these into single sheet
cnr.idr.ss <- read.csv('sheets/osaGFP-CnR-idrSheet.tsv', sep = '\t') %>%
  dplyr::mutate(grp = paste(genotype, fraction, fragType, sep = '.'))
print('faire')
faire.wt.ss <- read.csv('sheets/wt-FAIRE-sampleSheet.tsv', header = T, sep = '\t') %>%
  dplyr::mutate(grp = paste(stringr::str_split_fixed(sample, '_', n = 2)[,1],
                            stringr::str_split_fixed(sample, '_', n = 2)[,2],
                            sep = '.'),
                sample = paste(grp,rep, sep = '.'))
print('faire.pool')
faire.wt.ssPool <-read.csv('sheets/wt-FAIRE-sampleSheetPooled.tsv', header = T, sep = '\t') %>%
  dplyr::mutate(grp = paste(stringr::str_split_fixed(sample, '_', n = 2)[,1],
                            stringr::str_split_fixed(sample, '_', n = 2)[,2],
                            sep = '.'),
                time = factor(stringr::str_split_fixed(sample, '_', n = 2)[,2], levels = c('3LW','6h','18h','24h','44h')))


### save rData
save(cnr.ss, cnr.idr.ss, faire.wt.ss, faire.wt.ssPool, file = 'rData/sheets.rda')
