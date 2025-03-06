#' Augment a taxonomy table
#'
#' Uncharacterized taxa (name "NA") will be replaced with the names of the last defined taxonomic level.
#' The taxa can additionally be numbered
#'
#' @param taxtab taxonomy table, e.g. as produced by the dada2 pipeline
#' @param useBinomSpecies replaces a species name by Genus + species name, if present (default = TRUE)
#' @param numberUnclassified use numbers to distinguish uncharacterized taxa that belong to the same higher-level taxon (default = TRUE)
#' @param replaceIncertaeSedis if any name is "Incertae Sedis", this name is replaced by the next higher tax level with the addition "_i.s." (default = TRUE)
#' @export
#' @return augmented taxonomy table
#' @examples
#' taxtab <- augmentTaxonomy(taxtab)
#' taxtab <- augmentTaxonomy(taxtab, numberUnclassified = F)


augmentTaxonomy <- function(taxtab, useBinomSpecies = TRUE, numberUnclassified = TRUE, replaceIncertaeSedis = TRUE){

  if(is.null(taxtab)){ stop() }

  # define short versions of taxonomic levels
  taxlevels <- colnames(taxtab)
  levelshort <- substr(tolower(taxlevels), 0, 1) # uses first letter as new label

  # variables
  ntaxa <- nrow(taxtab)
  nlevels <- length(taxlevels)

  # scan tax levels one after the other
  augTaxTab <- taxtab

  # apply binomial species names
  if(useBinomSpecies){
    genuslevel <- tolower(taxlevels) == "genus"
    specieslevel <- tolower(taxlevels) == "species"
    if(any(genuslevel) & any(specieslevel)){
      augTaxTab[!is.na(augTaxTab[,specieslevel]), specieslevel] <- paste(augTaxTab[!is.na(augTaxTab[,specieslevel]), genuslevel], augTaxTab[!is.na(augTaxTab[,specieslevel]), specieslevel])
    }
  }

  # replace first level NAs
  NAcount <- sum(is.na(augTaxTab[,1]))
  if(NAcount){ # if at least one NA is already present in first level, use a number
    labelFormats <- paste("%0",ceiling(log10(NAcount)), "d", sep = "")
    augTaxTab[is.na(augTaxTab[,1]),1] <- paste("unclassified", sprintf(labelFormats, 1:NAcount), sep = "_")
  }

  # replace second to last level NAs
  for(i in 2:nlevels){ # start with second tax level
    nlabels <- table(augTaxTab[,i-1])
    labelFormats <- paste("%0",ceiling(log10(nlabels)), "d", sep = "")
    names(labelFormats) <- names(nlabels)
    NAcount <- nlabels*0
    for(j in 1:ntaxa){
      if(is.na(augTaxTab[j,i])){
        # replace NA assigment with the next higher level
        currTaxon <- augTaxTab[j,i-1]
        NAcount[currTaxon] <- NAcount[currTaxon] + 1  # increment NA count for this taxon
        if(numberUnclassified){
          augTaxTab[j,i] <- paste(currTaxon, paste(levelshort[i], sprintf(labelFormats[currTaxon], NAcount[currTaxon]), sep=""), sep="_")
        } else{
          augTaxTab[j,i] <- paste(currTaxon, levelshort[i], sep="_")
        }
      }
      if(augTaxTab[j,i] == "Incertae Sedis"){
        #replace Incertae Sedis with the next higher level
        currTaxon <- augTaxTab[j,i-1]
        augTaxTab[j,i] <- paste(currTaxon, levelshort[i], sep = "_i.s._")
      }
    }
  }



  return(augTaxTab)
}
