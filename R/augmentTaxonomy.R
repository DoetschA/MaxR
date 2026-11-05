#' Augment a taxonomy table (Version 1.3)
#'
#' Uncharacterized taxa (name "NA") will be replaced with the names of the last defined taxonomic level.
#' The taxa can additionally be numbered.
#'
#' Changelog:
#'
#' v1.1:
#' - writing short symbols for taxonomic levels now as prefixes instead of suffixes.
#' - tax symbols are not accumulated
#' - added variable labels for incertae sedis and separator for binomial species names
#'
#' v1.2: (16.10.2025)
#' - fixed bug leading to crashes when counting NAs
#'
#' v1.3: (30.1.2025)
#' - fixed bug leading to crashes with NAs occuring in the first column
#'
#' @param taxtab taxonomy table, e.g. as produced by the dada2 pipeline
#' @param useBinomSpecies replaces a species name by Genus + species name, if present (default = TRUE)
#' @param numberUnclassified use numbers to distinguish uncharacterized taxa that belong to the same higher-level taxon (default = TRUE)
#' @param replaceIncertaeSedis if any name is "Incertae Sedis", this name is replaced by the next higher tax level with the addition "_i.s." (default = TRUE)
#' @param IS_label how incertae sedis is displayed in the taxonomy (Default: 'Incertae Sedis')
#' @param binomsep separator between genus and species term in binomial species name (Default: ' ')
#' @export
#' @return augmented taxonomy table
#' @examples
#' taxtab <- augmentTaxonomy(taxtab)
#' taxtab <- augmentTaxonomy(taxtab, numberUnclassified = F)


augmentTaxonomy <- function(taxtab,
                            useBinomSpecies = TRUE,
                            numberUnclassified = TRUE,
                            replaceIncertaeSedis = TRUE,
                            IS_label = "Incertae Sedis",
                            binomsep = " "
                            ){

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
    # apply only if both levels are actually set (i.e. not NA and not "Incertae Sedis")
    if(any(genuslevel) & any(specieslevel)){
      # select only taxa were species level is set and not incertae sedis
      selected <- !is.na(augTaxTab[,specieslevel]) & (!augTaxTab[,specieslevel] == IS_label)
      augTaxTab[selected, specieslevel] <-
        paste(augTaxTab[selected, genuslevel],
              augTaxTab[selected, specieslevel],
              sep = binomsep)
    }
  }

  # replace first level NAs
  firstColNAs <- is.na(augTaxTab[,1])
  NAcount <- sum(firstColNAs)
  if(NAcount){ # if at least one NA is already present in first level, use a number
    labelFormats <- paste("%0",ceiling(log10(NAcount)), "d", sep = "")
    augTaxTab[firstColNAs, 1] <- paste("unclassified", sprintf(labelFormats, 1:NAcount), sep = "_")
  }

  # replace second to last level NAs
  for(i in 2:nlevels){ # start with second tax level
    # number of times each label occurs
    nlabels <- table(taxtab[,i-1]) # table(augTaxTab[,i-1])
    # set output format for numbers
    labelFormats <- paste("%0",ceiling(log10(nlabels)), "d", sep = "")
    names(labelFormats) <- names(nlabels)

    # cycle taxa in level i
    NAcount <- nlabels*0
    for(j in 1:ntaxa){
      # if NA, replace with lowest available label
      if(is.na(augTaxTab[j,i])){
        # get label of next higher taxon
        currentTaxon <- augTaxTab[j,i-1]
        # write this into original table to propagate changes conveniently
        taxtab[j, i] <- currentTaxon
        # add it to the NA count list
        NAcount <- c(NAcount, setNames(0, currentTaxon))
        # increment NA count for this taxon
        NAcount[currentTaxon] <- NAcount[currentTaxon] + 1
        # replace label
        if(numberUnclassified){
          # use numbers to differentiate between unclassified or incertae sedis taxa
          augTaxTab[j,i] <- paste(levelshort[i],
                                  paste(currentTaxon,
                                        sprintf(labelFormats[currentTaxon],
                                                NAcount[currentTaxon]),
                                        sep = "_"),
                                  sep = "__")
        } else {
          augTaxTab[j,i] <- paste(levelshort[i],
                                  currentTaxon,
                                  sep = "__")
        }

      }

      # in incertae sedis, replace with lowest available level plus IncSed suffix
      if(replaceIncertaeSedis){
        if(taxtab[j,i] == IS_label){
          # get label of next higher taxon
          currentTaxon <- taxtab[j, i-1]
          # write this into original table to propagate changes conveniently
          taxtab[j, i] <- currentTaxon
          # increment NA count for this taxon
          NAcount[currentTaxon] <- NAcount[currentTaxon] + 1
          # replace label
          if(numberUnclassified){
            # i.s. are also numbered
            augTaxTab[j,i] <- paste(levelshort[i],
                                    paste0(currentTaxon,
                                           "_IncSed_",
                                           sprintf(labelFormats[currentTaxon],
                                                   NAcount[currentTaxon])),
                                    sep = "__")
          } else {
            augTaxTab[j,i] <- paste(levelshort[i],
                                    paste0(currentTaxon, "_IncSed"),
                                    sep = "__")
          }
        }
      }
    }
  }


    return(augTaxTab)
}
