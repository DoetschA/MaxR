#' Retrieve k top taxa from a phyloseq object
#'
#' Sorts all taxa (at a specified level) by overall abundance and returns the names of the k most abundant taxa.
#' @param physeq Phyloseq object containing the microbiome data
#' @param k number of most abundant taxa to retrieve
#' @param taxlevel taxonomic level at which the taxa should be counted and sorted (default: Genus)
#' @param relative use relative abundances for counting taxa (default: TRUE)
#' @param agglomerate sum up all lower level members of each taxon (otherwise only the most abundant sequence variant will be counted (default: TRUE)))
#' @importFrom phyloseq tax_table taxa_sums tax_glom transform_sample_counts
#' @export
#' @return List of names of the most abundant taxa
#' @examples
#' top20genera <- getTopTaxa(physeq_feces, 20)
#' top5phyla  <- getTopTaxa(physeq_feces, 5, taxlevel = "Phylum")

getTopTaxa <- function(physeq, k, taxlevel = "Genus", relative = TRUE, agglomerate = TRUE){
  # use relative values
  if(relative)    physeq <- transform_sample_counts(physeq, function(x) x/sum(x))

  # agglomerate at specified tax level
  if(agglomerate) physeq <- tax_glom(physeq, taxrank = taxlevel)

  # get top taxa
  sortedTaxTable <- tax_table(physeq)[names(sort(taxa_sums(physeq), decreasing = T)), taxlevel]
  topKtaxa <- as.character(sortedTaxTable)[1:k]

  return(topKtaxa)
}
