#' Displays amount of reads at different steps during dada2 processing
#'
#' Inputs are data objects at various steps in the dada2 workflow and this script summarizes the amounts of reads that remain at each step.
#' Results are displayed as boxplots.
#' @param filtered Output table generate by filterAndTrim()
#' @param denoisedF output of dada() after denoising forward reads
#' @param denoisedR output of dada() after denoising reverse reads
#' @param merged output of mergePairs() after merging denoised forward and reverse reads
#' @param nochim output of removeBimeraDenovo() after removing chimeras
#' @param relative display results relative to initial number of reads (default: T)
#' @export
#' @return track
#'

trackReadsDada <- function(filtered, denoisedF, denoisedR, merged, nochim, relative = T, graphicOutput = T){
  require(reshape2)
  require(dada2)
  require(ggplot2)
  getN <- function(x) sum(getUniques(x))
  if(nrow(filtered)>1){
    track <- cbind(filtered,
                   sapply(denoisedF, getN),
                   sapply(denoisedR, getN),
                   sapply(merged, getN),
                   rowSums(nochim))
  } else {
    track <- cbind(filtout,
                   getN(denoisedF),
                   getN(denoisedR),
                   getN(merged),
                   rowSums(nochim))
  }

  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochim")

  if(relative){
    track <- track/track[,1]
  }

    track <- melt(track)
  colnames(track) <- c("name", "step", "reads")

  if(graphicOutput){
    ggplot(track, aes(x = step,
                      y = reads,
                      group = name,
                      color = name)) +
      geom_line() +
      expand_limits(y = 0) +
      theme(legend.position = "none")
  } else {
    return(track)
  }

}
