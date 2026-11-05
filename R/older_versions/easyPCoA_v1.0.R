#' A wrapper for simplified use of phyloseq's plot_ordination
#'
#' Instead of adding the same old ggplot options over and over, this function sets some options per default.
#' Also the last layer of the plot is removed, preventing double plotting of data points, which is common to plot_ordination.
#' The following ggplot options are used by default:
#'   + coord_fixed() # setting same scale to both axes
#'   + geom_point(size = size) # allowing control of data point sizes per parameter (see below)
#'   + theme(text = element_text(size = textsize), legend.key.size = unit(legend.key.size, "mm")) # allowing control of textsize and legend key size per parameter
#'   + scale_color_manual(values = colorpal) # setting color paletter per parameter
#'   + ggtitle(title) # setting the title per parameter
#' @param physeq a phyloseq object containing your microbiome data set.
#' @param ord a ordination object, output of phyloseq::ordinate()
#' @param color name of the variable to be displayed as colors (optional)
#' @param shape name of the variable to be displayed as colors (optional)
#' @param size data point size (default: 3)
#' @param axes two-element vector selecting the axes to be displayed (default: c(1,2))
#' @param colorpal color palette used for color coding
#' @param title title for the plot (optional)
#' @param textsize base size for text elements (default: 20)
#' @param legend.key.size size of legend key symbols, unit: mm (default: 8)
#' @export
#' @examples
#' myPCoA(physeq_BAL, ord_pcoa_bc, color = "feed", shape = "sampleType", colorpal = dietPalette, title = "PCoA: Bray Curtis")
#'

easyPCoA <- function(physeq, ord, color=NULL, shape=NULL, size=3, axes=c(1,2), colorpal=NULL, title=NULL, textsize = 20, legend.key.size = 8){
  p_pcoa <- plot_ordination(physeq,
                            ord,
                            type="samples",
                            color = color,
                            shape = shape,
                            axes = axes) +
    coord_fixed() +
    geom_point(size = size) +
    theme(text = element_text(size = textsize),
          legend.key.size = unit(legend.key.size, "mm")) +
    scale_color_manual(values = colorpal) +
    ggtitle(title)

  p_pcoa$layers <- p_pcoa$layers[-1]
  p_pcoa
}


