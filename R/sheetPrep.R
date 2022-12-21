### load packages
library(magrittr)
#library(GenomicRanges)

#source('scripts/utils.R')

### prep sample sheets
cnr.ss <- read.csv('sheets/osaGFP-CnR-sampleSheet.tsv', sep = '\t') %>%
  dplyr::mutate(genotype = stringr::str_split_fixed(sample, '-', n = 6)[,1],
                fraction = stringr::str_split_fixed(sample, '-', n = 6)[,6], 
                id = paste(genotype, fraction, rep, sep = '.'),
                grp = paste(genotype, fraction, sep = '.'), .before = 1)

### WT FAIRE
faire.wt.ss <- read.csv('sheets/wt-FAIRE-sampleSheet.tsv', header = T, sep = '\t') %>%
  dplyr::mutate(grp = paste(stringr::str_split_fixed(sample, '_', n = 2)[,1],
                            stringr::str_split_fixed(sample, '_', n = 2)[,2],
                            sep = '.'),
                sample = paste(grp,rep, sep = '.'), .before = 1)

faire.wt.ssPool <-read.csv('sheets/wt-FAIRE-sampleSheetPooled.tsv', header = T, sep = '\t') %>%
  dplyr::mutate(grp = paste(stringr::str_split_fixed(sample, '_', n = 2)[,1],
                            stringr::str_split_fixed(sample, '_', n = 2)[,2],
                            sep = '.'),
                time = factor(stringr::str_split_fixed(sample, '_', n = 2)[,2], levels = c('3LW','6h','18h','24h','44h')), .before = 1)

### osaGFP deGrad FAIRE
faire.osaDeGrad.ss <- read.csv('sheets/osa-deGrad-sampleSheet.tsv', header = T, sep = '\t') %>%
  dplyr::mutate(grp = paste(stringr::str_split_fixed(sample, '_', n = 2)[,1],
                            stringr::str_split_fixed(sample, '_', n = 2)[,2],
                            sep = '.'),
                sample = paste(grp,rep, sep = '.'), .before = 1)

faire.osaDeGrad.ssPool <-read.csv('sheets/osa-deGrad-sampleSheetPooled.tsv', header = T, sep = '\t') %>%
  dplyr::mutate(grp = paste(stringr::str_split_fixed(sample, '_', n = 2)[,1],
                            stringr::str_split_fixed(sample, '_', n = 2)[,2],
                            sep = '.'), before = 1)
### save rData
save(cnr.ss, faire.wt.ss, faire.wt.ssPool, faire.osaDeGrad.ss, faire.osaDeGrad.ssPool, file = 'rData/sheets.rda')
