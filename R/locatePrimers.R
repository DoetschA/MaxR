#' Locates primers sequences within fastq reads
#'
#' This function to identify primers was adapted from the dada2 ITS pipeline workflow (https://benjjneb.github.io/dada2/ITS_workflow.html).
#' Thanks to benjjneb!
#' @param fwdPrimer forward primer sequence
#' @param revPrimer reverse primer sequence
#' @param fwdFiles fwd FASTQ filenames
#' @param revFiles rev FASTQ filenames
#' @param onlyNfirst only analyse the N first files (for testing purpose)
#' @export
#' @return locatedPrimers
#'
locatePrimers <- function(fwdPrimer, revPrimer, fwdFiles, revFiles, onlyNfirst = NA){
  require(Biostrings)
  require(ShortRead)

  # local functions
  allOrients <- function(primer) {
    # Create all orientations of the input sequence
    dna <- Biostrings::DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = reverse(dna),
                 RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
  }

  primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
  }

  # get possible primer orientations
  fwdPrimer.orients <- allOrients(fwdPrimer)
  revPrimer.orients <- allOrients(revPrimer)

  # get available files


  if(!is.na(onlyNfirst)){
    fwdFiles <- fwdFiles[1:onlyNfirst]
    revFiles <- revFiles[1:onlyNfirst]
  }

  # find locations
  locatedPrimers <- rbind(fwdPrimer.ForwardReads = sapply(fwdPrimer.orients, primerHits, fn = fwdFiles[[1]]),
                          fwdPrimer.ReverseReads = sapply(fwdPrimer.orients, primerHits, fn = revFiles[[1]]),
                          revPrimer.ForwardReads = sapply(revPrimer.orients, primerHits, fn = fwdFiles[[1]]),
                          revPrimer.ReverseReads = sapply(revPrimer.orients, primerHits, fn = revFiles[[1]]))

  # normalize by read count
  readcount  <- countFastq(fwdFiles) # revFiles not needed, since it should be symmetrical at this point
  nreads <- sum(readcount$records)

  locatedPrimers <- locatedPrimers/nreads
  return(locatedPrimers)

}
