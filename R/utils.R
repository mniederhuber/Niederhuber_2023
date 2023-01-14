
#DONE - add functionality to return only qfiltered peaks that are reproduced in each replicate - subsetByOverlaps() ? 
#wrapper to pass faire peak lists (grangeslist) by experiment with replicate ids to qFilter()
grp_qFilter <- function(x, quantile = NULL, q = NULL, operation = c('subsetByOverlaps', 'intersect'))  {
  grp <- x$grp
  ids <- unique(x$id) #unique replicate ids
  
  #for each unique replicate id, pass experiment peak granges and replicate name to qFilter, then bind the output list
  grp.filtered <- lapply(ids, function(y) qFilter(x, y, quantile, q)) 
  
  grp.list <- grp.filtered %>% GenomicRanges::GRangesList()
  
  #TODO - cleaner way to write this?
  if(operation == 'subsetByOverlaps'){
    grp.repShared <- Reduce(IRanges::subsetByOverlaps, grp.list) #takes grangeslist of qval filtered peaks and returns a granges object with only regions shared between replicates
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






