
#DONE - add functionality to return only qfiltered peaks that are reproduced in each replicate - subsetByOverlaps() ? 
#wrapper to pass faire peak lists (grangeslist) by experiment with replicate ids to qFilter()
grp_qFilter <- function(x, quantile = NULL, q = NULL) {
  grp <- x$grp
  ids <- unique(x$id) #unique replicate ids
  
  #for each unique replicate id, pass experiment peak granges and replicate name to qFilter, then bind the output list
  grp.filtered <- lapply(ids, function(y) qFilter(x, y, quantile, q)) 
  
  grp.list <- grp.filtered %>% GenomicRanges::GRangesList()
  
  grp.repShared <- Reduce(IRanges::subsetByOverlaps, grp.list) #takes grangeslist of qval filtered peaks and returns a granges object with only regions shared between replicates
  
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