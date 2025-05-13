
#DONE - add functionality to return only qfiltered peaks that are reproduced in each replicate - subsetByOverlaps() ? 
#wrapper to pass faire peak lists (grangeslist) by experiment with replicate ids to qFilter()
grp_qFilter <- function(x, quantile = NULL, q = NULL, operation = c('subsetByOverlaps', 'intersect'), with_reduce = F)  {
  grp <- x$grp
  ids <- unique(x$id) #unique replicate ids
  
  #for each unique replicate id, pass experiment peak granges and replicate name to qFilter, then bind the output list
  grp.filtered <- lapply(ids, function(y) qFilter(x, y, quantile, q)) 
  
  grp.list <- grp.filtered %>% GenomicRanges::GRangesList()
 
  #TODO - cleaner way to write this?
  # TODO -- potential BUG -- return of subsetByOverlaps is coordinates of query, ie rep 1 of the 2 reps. is this even a problem?
  if(operation == 'subsetByOverlaps'){
      # by using Reduce() subsetByOverlaps can be applied to > 2 replicates - which is true for some WT FAIRE experiments
      # subsetByOverlaps will run A vs B => A[B] overlap; A[B] vs C => A[BC] overlap; etc. -- note that subsetByOverlaps returns ranges of query not a merge, in this example A
    if(with_reduce){
      grp.repShared <- Reduce(IRanges::subsetByOverlaps, grp.list) #takes grangeslist of qval filtered peaks and returns a granges object with only regions shared between replicates
    }else{
      grp.repShared.1 <- IRanges::subsetByOverlaps(grp.list[[1]], grp.list[[2]]) # get coords in 1st rep that overlap rep 2
      grp.repShared.2 <- IRanges::subsetByOverlaps(grp.list[[2]], grp.list[[1]]) # recipricol, get coords in 2nd rep that overlap 1st
      grp.repShared <- GenomicRanges::union(grp.repShared.1, grp.repShared.2) #union the two granges -- this avoids coordinates from only one rep being used in final granges
    }
  }
  if(operation == 'intersect'){
    grp.repShared <- Reduce(GenomicRanges::intersect, grp.list) #takes grangeslist of qval filtered peaks and returns a granges object with only regions shared between replicates
  }
  return(grp.repShared)
}

#takes grouped FAIRE granges object - ie. all replicates for same sample - and filters out bottom 25% peaks by qval per replicate
qFilter <- function(df, id, quantile = NULL, q = NULL) {
  
  df %<>% data.frame() %>% dplyr::filter(id == !!id) #filter peak list to specific replicate, "!!" evaluates variable
  
  if(!is.null(quantile)) {
    q <- quantile(df$qValue, probs = quantile) #find the 75% qVal
  } 
    
  df.filtered <- df %>% 
    dplyr::filter(qValue >= q) %>% 
    dplyr::mutate(qCutoff = q) 
  
  return(df.filtered) #return the filtered peak list
}


#The SAF format includes five required columns for each feature: feature identifier, chromosome name, start position, end position and strand
makeSAF <- function(df) {
  # df.filtered <- dplyr::filter(df, grp_rep == x)  
  # df.renamed <- dplyr::mutate(df, Chr = seqnames, .before = 2)  
  df.saf <- df %>%
    dplyr::select(peak, seqnames, start, end, strand) 
  
  colnames(df.saf) <- c('GeneID', 'Chr', 'Start', 'End','Strand')
  
  return(df.saf)
}

## varient function for single bw
get_cnr_track <- function(bw, ylim, fill, col, ...) {
    track <- Gviz::DataTrack(range = bw,
                    genome = 'dm6',
                    type = 'hist',
                    #windowSize = ,
                    fill.histogram = fill,
                    col.histogram = col, 
                    background.title = 'white',
                    fontcolor.title = 'black',
                    col.axis = 'black',
                    ylim = ylim)  
  return(track)
}

get_cnr_tracks <- function(sheet, ylim, ...) {
  track <- purrr::map(split(sheet, sheet$id), ~{ #split will convert group variable to factor, so levels are default alphabetical which works out to be order I want
    #TODO - explicitly denote factor levels? likely necessary for plotting with faire data
    Gviz::DataTrack(range = .$bigwig_rpgcNorm_zNorm,
                    genome = 'dm6',
                    type = 'hist',
                    #windowSize = ,
                    name = .$id,
                    group = .$grp,
                    fill.histogram = .$color,
                    col.histogram = .$color, 
                    background.title = 'white',
                    fontcolor.title = 'black',
                    col.axis = 'black',
                    ylim = ylim)  
  })
  return(track)
}


#get bws for use in heatmaps
get_bws <- function(sheet, by) {
  bw_list <- purrr::map(sheet[[by]], ~{
    fp <- sheet[sheet[by] == .,]$bigwig_rpgcNorm_zNorm
    rtracklayer::import.bw(fp)
    })  
  names(bw_list) <- sheet[[by]]
  return(bw_list)
}

#TODO - uses fbgn.validated which is set as a global variable in notebook, fix that
#TODO - cleanup filtering of problematic fbgns 
#TODO - comment functions
tomUnnest <- function(x) {
  unnest <- x %>%
    tidyr::unnest(tomtom) %>%
    dplyr::filter(!is.na(match_name)) %>%
    dplyr::select(consensus, match_name, match_altname, match_eval) %>%
#    .[c(6,50,51,56)] %>% #DONE change to dplyr::select specific column names
    unique() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(neg_log_eval = -log10(1 + match_eval),
                  FF_FBgn= unlist(stringr::str_split_fixed(match_name, pattern = '_', n = 2)[1])) %>%
    dplyr::filter(FF_FBgn != 'FBgn0051782' & FF_FBgn != 'FBgn0265182') %>% #filter out these fbgn -- have problematic validated matches
    dplyr::mutate(Geneid = purrr::map_chr(FF_FBgn, function(x) {
      fbgn.validated[fbgn.validated$FlyFactor == x,]$Validated
    })) 
  
  return(unnest)
}

## generate a new dataframe for each consensus sequence and its associated tomtom results
tomResults <- function(x) {
  
  unnest <- tomUnnest(x) 
  
  ## make a list of dfs, one for each consensus in the input results 
  tom <- lapply(unique(unnest$consensus), function(x) {
    tom.seq <- unnest %>% 
      dplyr::filter(consensus == x) %>%
      dplyr::left_join(., wing.rnaseq, id = 'Geneid') %>%
      dplyr::ungroup() %>%
#      dplyr::slice_min(evalue, n = 3, with_ties = F)
      dplyr::top_n(8, neg_log_eval) 
    
    tom.seq$match_altname = forcats::fct_reorder(tom.seq$match_altname, tom.seq$neg_log_eval, .desc = F)
    
    return(tom.seq)
  })
  
  return(tom)
}


format_rna <- function(x) {
  x %>%
    .[c(1,3,8:19)] %>%
    reshape2::melt() %>% 
    dplyr::rowwise() %>%
    dplyr::mutate(variable = stringr::str_replace(variable, 'X', ''),
                  rna_grp = stringr::str_split_fixed(variable, pattern = '_', n = 2)[1],
                  rna_grp = factor(rna_grp, levels = c('L3','6h','18h','24h','36h','44h'))) %>%
    dplyr::group_by(consensus, match_altname, rna_grp) %>%
    dplyr::summarise(mean = mean(value)) %>%
    dplyr::mutate(log_mean = log10(1 + mean)) 
}



