#' Calculate level of determination of for each taxonomic level
#'
#' Counts the number of NA assigned to sequence variants at each level and devides by the total number of ASVs.
#' @param taxtab Taxonomy table - either a data.frame or phyloseq object created with tax_table()
#' @export
#' @return tld
#'

taxLevelDetermination <- function(taxtab){

  tld <- (1-colSums(is.na(taxtab))/nrow(taxtab))*100

  return(tld)

}
