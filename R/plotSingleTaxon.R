#' Display only a single taxon from a phyloseq object
#'
#' @param physeq a phyloseq object containing your microbiome data set.
#' @param taxname name of the taxon to be displayed (should appear in tax_table(physeq))
#' @param param1 grouping variable (will be used for the x axis)
#' @param param2 (optional) secondary grouping variable (will be displayed as different colors)
#' @param taxlevel name of the taxonomic level at which the taxon is characterized (Default: ASV)
#' @param mode Which type of graph should be used? "violin" - violin plot (default), "boxplot" - boxplot, "point" - dot plot
#' @param colorpal Define a color palette to use for the graph.
#' @param relative Display relative abundances if TRUE, absolute counts if FALSE (Default: TRUE).
#' @export
#' @examples
#' plotSingleTaxon(physeq, "ASV0001", "group", colorpal = myPalette21)
#' plotSingleTaxon(physeq, "Pseudomonas", "group", taxlevel = "Genus", colorpal = myPalette21)
#' plotSingleTaxon(physeq, "Firmicutes", "group", "timepoint", taxlevel = "Phylum", colorpal = myPalette21)

plotSingleTaxon <- function(physeq, taxname, param1, param2 = NULL, taxlevel = "ASV", mode = "violin", colorpal, relative = TRUE){
  # reduce phyloseq object (not needed and may cause warning reducing to tree with only 1 tip )
  physeq <- phyloseq(otu_table(physeq), sample_data(physeq), tax_table(physeq))

  # transform to relative values
  if(relative){
    physeq <- transform_sample_counts(physeq, function(x) x/sum(x))
    y_label_abundance <- "rel. abundance"
  } else {
    y_label_abundance <- "absolute counts"
  }
  # gather tax data
  if(taxlevel == "ASV"){
    # select single ASV
    plotdata <- cbind(otu_table(physeq)[,taxname], sample_data(physeq)[, c(param1, param2)])
    taxstring <- paste(taxname, tax_table(physeq)[taxname, "Species"])
  } else {
    # prune unneccessary taxa
    physeq <- prune_taxa(as.logical(tax_table(physeq)[,taxlevel] == taxname), physeq)
    # merge taxa
    options(warn = -1)
    physeq <- tax_glom(physeq, taxrank = taxlevel)
    taxa_names(physeq) <- taxname
    # prepare plot data
    plotdata <- cbind(otu_table(physeq), sample_data(physeq)[, c(param1, param2)])
    taxstring <- taxname
  }

  # assemble plot data
  # check for second parameter
  if(is.null(param2)) param2 <- param1
  # create plot
  if(mode == "violin"){
    ggplot(plotdata, aes_string(x = param1, y = taxname, fill = param2)) +
      geom_violin(scale = "width") +
      stat_summary(fun=median, geom='point', shape = 5) +
      scale_fill_manual(values = colorpal) +
      theme(text = element_text(size=18)) +
      ylab(y_label_abundance) +
      ggtitle(taxstring)
  } else if(mode == "boxplot"){
    ggplot(plotdata, aes_string(x = param1, y = taxname, fill = param2)) +
      geom_boxplot() +
      scale_fill_manual(values = colorpal) +
      theme(text = element_text(size=18)) +
      ylab(y_label_abundance) +
      ggtitle(taxstring)
  } else if(mode == "point"){
    ggplot(plotdata, aes_string(x = param1, y = taxname, color = param2)) +
      geom_point() +
      scale_color_manual(values = colorpal) +
      theme(text = element_text(size=18)) +
      ylab(y_label_abundance) +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                   geom = "crossbar", width = 0.1) +
      ggtitle(taxstring)
  }
}
