#' Stacked bargraph of data (Version 1.0)
#'
#' Useful for data were you want to show composition of something, e.g.,
#' of a microbiome or of fatty acids. The values are stacked to add up for a sample,
#' but can also be grouped using other factor variables in the dataset.
#' Data should be organized in a long table format with a column for values and a column for the names of the components.
#' @param data Data.frame used for plotting. Use a long format.
#' @param values Column containing the values to be used for the bar graph.
#' @param names Column containing the names for the components that are represented in the values.
#' @param sampleID Column containing the names of samples that should be represented as stacked bar graphs.
#' @param ktop If this is set (i.e. not NULL), only the ktop highest components will be shown and the rest summarized as 'other' (Default: NULL).
#' @param relative If TRUE, colums are normalized to a sum of 1 for each stacked column. Default: FALSE.
#' @param groupBy (optional) If used, uses this factor variable to group samples and display sums/average values.
#' @param groupBy2 (optional) Defines a secondary grouping variable.
#' @param othercol Defines color of "other" bar that sums up all taxa beyond the ktop most abundant. Set to NULL to disable this bar. Default color is #BBBBBB (medium gray).
#' @param textsize Sets the main text size of the graph. Default is 10.
#' @param colorPalette Define a color palette to use for the graph. Default: NULL.
#' @param maxlegend Maximum number of taxa in legend before legend has two columns (Default: 9)
#' @param xlabel Label to be displayed on the X axis (Default: NULL)
#' @param ylabel Label to be displayed on the Y axis (Default: NULL)
#' @seealso MaxR::bargraph_top_taxa
#' @export
#' @examples
#' #simple bargraph of all samples
#' bargraph_stacked(scfa_data,
#'                  values = "concentration",
#'                  names = "SCFA",
#'                  sampleID = "Probenname")
#'
#' #relative values and grouped by variable ("Lauf")
#' bargraph_stacked(scfa_data,
#'                  values = "concentration",
#'                  names = "SCFA",
#'                  sampleID = "Probenname",
#'                  relative = T,
#'                  groupBy = "Lauf")

bargraph_stacked <- function(data,
                             values = NULL,
                             names = NULL,
                             sampleID = NULL,
                             ktop = NULL,
                             relative = F,
                             groupBy = NULL,
                             groupBy2 = NULL,
                             colorPalette = NULL,
                             othercol = "#BBBBBB",
                             textsize = 10,
                             maxlegend = 9,
                             xlabel = NULL,
                             ylabel = NULL){
  require(ggplot2)
  require(dplyr)

  # check for mandatory arguments
  if(is.null(values) | is.null(names) | is.null(sampleID)){
    stop(warning("Mandatory arguments (values, names, sampleID) not set."))
  }

  # apply color for 'other'
  colorPalette <- c(colorPalette[1:ktop], othercol)

  # Group data by sampleID and names, and calculate the sum of values
  data_sum <- data |>
    group_by(!!sym(sampleID), !!sym(names)) |>
    summarise(sum_values = sum(!!sym(values)))

  # calculate ktop components a sum the rest as 'other'
  if(!is.null(ktop)){
    # determine top components across all samples
    include_names <- data_sum |>
      group_by(!!sym(names)) |>
      summarise(sum_values = sum(sum_values)) |>
      arrange(desc(sum_values)) |>
      mutate(rank = row_number()) |>
      filter(rank <= ktop) |>
      select(!!sym(names)) |>
      unlist()
    data_top <- data_sum |>
      filter(!!sym(names) %in% include_names)

    # calculate 'other'
    other_names <- setdiff(levels(include_names), include_names)
    data_other <- data_sum |>
      group_by(!!sym(sampleID)) |>
      filter(!!sym(names) %in% other_names) |>
      summarise(sum_values = sum(sum_values)) |>
      mutate(!!sym(names) := "other")

    # Combine data
    data_plot <- bind_rows(data_top, data_other)
  } else{
    data_plot <- data_sum
  }

  # normalize to relative values
  if(relative) {
    data_plot <- data_plot |>
      group_by(!!sym(sampleID)) |>
      mutate(sum_values = sum_values / sum(sum_values))
  }

  # apply further grouping variables
  if(!is.null(groupBy)){
    if(!is.null(groupBy2)){
      # use two grouping variables
      # re-add grouping variable
      data_plot <- data |>
        select(!!sym(sampleID), !!sym(groupBy), !!sym(groupBy2)) |>
        unique() |>
        right_join(data_plot)

      # group
      data_plot <- data_plot |>
        group_by(!!sym(names), !!sym(groupBy), !!sym(groupBy2)) |>
        summarise(sum_values = sum(sum_values)) |>
        mutate(sampleID = paste0(!!sym(groupBy), ":", !!sym(groupBy2)))
    } else {
      # use only one grouping variable

      # re-add grouping variable
      data_plot <- data |>
        select(!!sym(sampleID), !!sym(groupBy)) |>
        unique() |>
        right_join(data_plot)
      # group
      data_plot <- data_plot |>
        group_by(!!sym(names), !!sym(groupBy)) |>
        summarise(sum_values = sum(sum_values)) |>
        mutate(sampleID = !!sym(groupBy))
    }

    # normalize again within new groups to relative values
    # this makes sure that diffently sized groups are equally weighted
    if(relative) {
      data_plot <- data_plot |>
        group_by(sampleID) |>
        mutate(sum_values = sum_values / sum(sum_values))
    }
  } else {
    # make sure that there is a "sampleID" column
    data_plot <- data_plot |>
      mutate(sampleID = !!sym(sampleID))
  }



  # Plot graph
  p <- ggplot(data_plot,
              aes_string(x = "sampleID", y = "sum_values", fill = names)) +
    geom_col() +
    labs(x = xlabel,
         y = ylabel) +
    theme_bw() +
    theme(text = element_text(size =  textsize))

  # apply color palette
  if(!is.null(colorPalette)){
    p <- p + scale_fill_manual(values = colorPalette)
  }

  # set max number of legend rows if maxlegend is set
  if(n_distinct(data_plot[[names]]) > maxlegend){
    p <- p + guides(fill=guide_legend(ncol=2))
  }

  # return graph
  return(p)
}
