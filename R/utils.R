#function to filter qValues - set at 75% cutoff 
#takes grouped FAIRE granges object - ie. all replicates for same sample - and filters out bottom 25% peaks by qval per replicate
qFilter <- function(df, id) {
  df %<>% data.frame() %>% dplyr::filter(id == !!id) #filter peak list to specific replicate, "!!" evaluates variable
  q <- quantile(df$qValue, probs = 0.75) #find the 75% qVal
  df.filtered <- df %>% 
    dplyr::filter(qValue >= q) %>% 
    dplyr::mutate(qCutoff = q) 
  
  return(df.filtered) #return the filtered peak list
} 

#TODO - add functionality to return only qfiltered peaks that are reproduced in each replicate - subsetByOverlaps() ? 
#wrapper to pass faire peak lists by experiment with replicate ids to qFilter()
grp_qFilter <- function(x) {
  grp <- x$grp
  ids <- unique(x$id) #unique replicate ids
  
  #for each unique replicate id, pass experiment peak granges and replicate name to qFilter, then bind the output list
  grp.filtered <- purrr::map(ids, function(y) qFilter(x, y)) #%>% dplyr::bind_rows()
   
  return(grp.filtered)
}
