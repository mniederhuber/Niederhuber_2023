library(magrittr)

brm_complex <- c('Bap111','brm','mor','Snr1','Act5C','Bap55')
pbap_complex <- c('polybromo', 'bap170')
bap_complex <- c('osa')
ino80 <- c('Arp5','Ino80','pho','Act5C')
dom <- c('dom','Tip60','Nipped-A','pont','Reptin')
atac <- c('Atac2','Atac3','D12','Wds','NC2beta')
NURF <- c('E(bx)', 'Nurf-38','ISWI','Caf1-55')
ACF <- c('Acf','Chrac-16')
NURD <- c('Mi-2','MTA1-like','Caf1-55')
NORC <- c('CtBP','Tou')
demethyl <- c('Jarid2','lid','su(var)3-3')
SNF2_like = c('kis','okr','Chd1','Etl1')
control <- c('lexA')

nucleosome_remod <- c(brm_complex, pbap_complex, bap_complex, ino80, dom, atac, NURF, ACF, NURD, SNF2_like)

outdata <- read.csv('../RNAi_Screen-Image_Quantification/output/resultsDF.csv') %>%
  dplyr::mutate(id = paste(symbol, line, sep = "_"),
                nucleosome_remodeler = ifelse(symbol %in% nucleosome_remod, T, F),
                Negcontrol = ifelse(symbol %in% control, 'Negative-Control', 'experimental'),
                complex_general = dplyr::case_when(symbol %in% brm_complex | symbol %in% pbap_complex | symbol %in% bap_complex ~ 'PBAP/BAP',
                                                   symbol %in% control ~ 'Negative-Control',
                                                   T ~ 'other')) %>%
  dplyr::filter(stack >= 10 & gfp >= 8 & size == "1024x1024") %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(area_ratio = gfp_mask_area / gfpNEG_mask_area,
                area_sum = gfp_mask_area + gfpNEG_mask_area,
                log2 = log2(KDtoWT),
                mean_ratio = median(KDtoWT), 
                meanlog2 = mean(log2), # changing from filter on mean_log2...
                bap_hits = dplyr::case_when(mean_ratio >= 1.2 & complex_general == 'PBAP/BAP' ~ 'Brahma Complex', #increase of atleast 150% 
                                            complex_general == 'Negative-Control' ~ 'Negative-Control',
                                            T ~ 'other'))

rnaiScreen <- outdata %>% dplyr::filter(nucleosome_remodeler | Negcontrol == 'Negative-Control')

save(rnaiScreen, file = '../R/rData/rnaiScreen.rda')
