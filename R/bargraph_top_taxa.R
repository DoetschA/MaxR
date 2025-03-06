#' Phyloseq-style bargraph of the top k taxa
#'
#' Uses phyloseq bargraphs to display a defined number of most abundant taxa in the data set.
#' @param physeq A phyloseq object containing your microbiome data set.
#' @param k Number of most abundant taxa that should be displayed. Defaults to 17 (corresponding to the default palette.)
#' @param taxlevel Taxonomic level used for coloring the graph. Must use only levels that are defined in tax_table(physeq). Defaults to Genus.
#' @param relative If TRUE, use relative abundances, else use absolute counts. Default: TRUE.
#' @param groupBy (optional) If used, uses this variable to group samples and display sums/average values. Must use variables that are present in sample_data(physeq).
#' @param groupBy2 (optional) Defines a secondary grouping variable.
#' @param mergelevel (optional) taxa are merged on the specified level, displaying only one filled bar per taxon. Always sets colorfill = TRUE.
#' @param colorfill (deprecated, use mergelevel instead) If TRUE, margins of taxa bars are colored as well, displaying the whole taxon on the definedtaxlevel in one color. If FALSE, black margins will allow to distinguish low level taxa (ASVs). Default: FALSE.
#' @param othercol Defines color of "other" bar that sums up all taxa beyond the k most abundant. Set to NULL to disable this bar. Default color is #BBBBBB (medium gray).
#' @param textsize Sets the main text size of the graph. Default is 10.
#' @param colorPalette Define a color palette to use for the graph. Default uses RColorBrewer and combines two palettes: c(brewer.pal(9,"Set1"),brewer.pal(8,"Dark2")).
#' @param maxlegend Maximum number of taxa in legend before legend has two columns (Default: 9)
#' @param xlabel Label to be displayed below X axis (Default: NULL)
#' @seealso phyloseq::plot_bar()
#' @export
#' @examples
#' bargraph_top_taxa(physeq) #simple bargraph of all samples
#'
#' bargraph_top_taxa(physeq, k=10, groupBy = "group") # bargraph of top 10 taxa, grouped by variable "group"
#' bargraph_top_taxa(physeq, k=10, colorfill = T, groupBy = "group") # same with colors filled
#'
#' bargraph_top_taxa(physeq, k=6, colorfill = T, taxlevel = "Phylum", othercol = NULL) # top 6 phylo (no "other")


bargraph_top_taxa <- function(physeq,
                              k = 17,
                              taxlevel = "Genus",
                              relative = T,
                              groupBy = NULL,
                              groupBy2 = NULL,
                              mergelevel = F,
                              colorfill = F,
                              colorPalette = NULL,
                              othercol = "#BBBBBB",
                              textsize = 10,
                              maxlegend = 9,
                              xlabel = NULL){
  require(phyloseq)
  require(RColorBrewer)
  require(ggplot2)

  # define color palette, if none is provided
  if(is.null(colorPalette)){
    colorPalette <- c(brewer.pal(9,"Set1"),brewer.pal(8,"Dark2"))
  }

  # check for empty samples that might cause errors
  if(any(sample_sums(physeq) == 0)){
    physeq <- prune_samples(sample_sums(physeq) > 0, physeq)
  }

  # transform absolute to relative counts (recommended for grouped samples to avoid  means weighted by group size)
  if(relative){
    physeq <- transform_sample_counts(physeq, function(x) 100 * x/sum(x))
  }

  # plot only one bar per taxon by merging taxa before generating the plot
  if(mergelevel){
    colorfill <- TRUE
    physeq <- tax_glom(physeq, taxrank = taxlevel)
  }

  # warning issued if old "colorfill" method is used, since this creates line plots with a lot of objects that get hard to handle
  if(colorfill & !mergelevel){
    warning("The 'colorfill = T' option is deprecated. Use 'mergelevel = T' instead, because it will produce plots with much less objects, which makes them much more convenient to use!")
  }

  # error if groupBy or groupBy2 are not a variable in the data set
  if(!is.null(groupBy) && !(groupBy %in% sample_variables(physeq))){
    stop(paste0("The variable \"", groupBy, "\" is not present in the data set. Please check!"))
  }
  if(!is.null(groupBy2) && !(groupBy2 %in% sample_variables(physeq))){
    stop(paste0("The variable \"", groupBy2, "\" is not present in the data set. Please check!"))
  }

  # group by a variable and determine mean abundances
  if(!is.null(groupBy)){
    if(is.null(groupBy2)){

      groupLevels <- levels(as.factor(sample_data(physeq)@.Data[[which(sample_variables(physeq) == groupBy)]]))

      # check if only one level (or none at all...) is present
      if(length(groupLevels) < 2){warning(paste("Not enough levels for groupBy:", length(groupLevels)))}
      # merge samples and replace the summed abundances (mean is not implemented correctly!)
      physeq <- merge_samples(physeq, groupBy)
      # repeating the transformation gets the values back to percent, effectively calculating the mean
      if(relative){
        physeq <- transform_sample_counts(physeq, function(x) 100 * x/sum(x))
      }
    } else {
      #create new dummy variable containing both grouping factors
      gf1 <- sample_data(physeq)@.Data[[which(names(sample_data(physeq)) == groupBy)]]  # get grouping variable fro sample_data
      gf2 <- sample_data(physeq)@.Data[[which(names(sample_data(physeq)) == groupBy2)]]
      gf1 <- as.factor(gf1)
      gf2 <- as.factor(gf2)

      # check if only one level (or none at all...) is present
      if(length(levels(gf1)) < 2){warning(paste("Not enough levels for groupBy:", length(levels(gf1))))}
      if(length(levels(gf2)) < 2){warning(paste("Not enough levels for groupBy2:", length(levels(gf2))))}

      # create dummy variable for grouping
      gfX <- expand.grid(gf1levels = levels(gf1), gf2levels = levels(gf2))           # get combination of levels
      grouping    <- paste(as.character(gf1), as.character(gf2), sep = ":") # create dummy variable
      groupLevels <- paste(gfX$gf1levels, gfX$gf2levels, sep = ":")
      grouping    <- factor(grouping, levels = groupLevels) # important to keep original order of levels! grouping <- paste(gfX$g1, gfX$g2, sep = ":")       # create dummy variable
      grouping    <- droplevels(grouping) # remove empty levels
      sample_data(physeq)[,"grouping"] <- grouping

      # merge samples and replace the summed abundances (mean is not implemented correctly!)
      physeq <- merge_samples(physeq, "grouping")
      # repeating the transformation gets the values back to percent, effectively calculating the mean
      if(relative){
        physeq <- transform_sample_counts(physeq, function(x) 100 * x/sum(x))
      }
    }
  }


  # find top k unique taxa on requested level, or if k is not numeric use the provided list of taxa names
  list_taxa <- as.character(unique(tax_table(physeq)[names(sort(taxa_sums(physeq), decreasing = T)), taxlevel]))
  if(is.numeric(k)){
    # numeric k => use top k taxa
    if(length(list_taxa) >= k){
      keep_taxa <- list_taxa[1:k]
    } else {
      keep_taxa <- list_taxa
    }
  } else {
    keep_taxa <- list_taxa[list_taxa %in% k]
  }


  physeq_selection <- prune_taxa(tax_table(physeq)[,taxlevel] %in% keep_taxa, physeq)

  # add "others" as a dummy taxon
  if(!(is.na(othercol) || is.null(othercol))){
    other.tax <- tax_table(physeq_selection)[1,]
    row.names(other.tax) <- "_other"
    for(i in 1:length(other.tax)){other.tax[i] <- "_other"}
    other.tax <- merge_phyloseq(tax_table(physeq_selection), tax_table(other.tax))

    other.otu <- as.matrix(sample_sums(physeq) - sample_sums(physeq_selection))
    colnames(other.otu) <- "_other"
    other.otu <- merge_phyloseq(otu_table(physeq_selection), otu_table(other.otu, taxa_are_rows = F))

    physeq_selection <- phyloseq(sample_data(physeq_selection), other.otu, other.tax)
  }

  # barplot of top k taxa
  if(is.numeric(k)){
    colorPalette <- c(othercol, colorPalette[1:k]) # add color for "other"
  } else {
    colorPalette <- c(othercol, colorPalette[1:length(k)]) # add color for "other"
  }
  if(relative) { ytext <- "relative abundance (%)" } else { ytext <- "read counts"}

  if(is.null(groupBy)){
    # no grouping
    p_bar <- plot_bar(physeq_selection, fill=taxlevel) +
      ylab(ytext) +
      theme(text = element_text(size=textsize)) +
      scale_fill_manual(values = colorPalette) +
      scale_color_manual(values = colorPalette) +
      scale_x_discrete(name = xlabel)
  } else {
    # apply grouping
    p_bar <- plot_bar(physeq_selection, fill=taxlevel) +
      ylab(ytext) +
      theme(text = element_text(size=textsize)) +
      scale_fill_manual(values = colorPalette) +
      scale_color_manual(values = colorPalette) +
      scale_x_discrete(limits = groupLevels,
                       name = xlabel)
  }

  # wrap legend if necessary
  if(!is.numeric(k)) k <- length(k)
  if((k + !is.null(othercol)) > maxlegend ) {
    p_bar <- p_bar + guides(fill=guide_legend(ncol=2))
  }

  if(colorfill){
    eval(parse(text=paste("p_bar + geom_bar(aes(color=", taxlevel, ", fill=", taxlevel, "), stat='identity', position='stack')", sep="")))
  } else {
    p_bar #+ facet_grid(cols = vars(phase))
  }

}
