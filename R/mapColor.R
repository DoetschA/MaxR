#' Maps values to colors
#'
#' Values are converted into color statements that can e.g. be exported for use in external programs.
#' @param values Vector of values
#' @param map Colormap as defined e.g. by RColorBrewer
#' @param lims two-element vector of lower and upper limit. Defaul: c(min(values), max(values))
#' @export
#' @return color map
#'

mapColor <- function(values, limColors, steps, lims = NULL){
  if(is.null(lims)){
    lims <- c(min(values), max(values))
  }

  map <- colorRampPalette(limColors)

  ii <- findInterval(values, seq(lims[1], lims[2], length.out = steps + 1), all.inside = T)

  return(map(steps)[ii])
}
