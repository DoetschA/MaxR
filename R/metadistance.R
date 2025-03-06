#' Combines a distance matrix with metadata to a data.frame that can be used for ggplot graphs to e.g. group dissimilarities by certain factors (e.g. belonging to the same group)
#'
#' @param dist distance matrix of class dist()
#' @param meta data.frame containing factors (values will be assigned specfic levels if both samples have the same value and NA otherwise)
#' @param meta_n (optional) data.frame containing numeric metadata (values will be absolute differences)
#' @export
#' @return A data.frame in a format suitable to diplay distance values for certain types of comparisons.
#' Output columns are:
#' - 1st, 2nd: names of the two samples that are compared
#' - 3rd: distance (as defined in 'dist') between the two samples
#' - _comp - TRUE if variable is same in both samples
#' - _labels - displaying compared samples.
#' @examples
#' metadistance(dist_jsd, meta = sample_data(physeq_feces)[,c("proband_ID","kit_purif","lab_purif","MiSeq")], meta_n = ord_pcoa_jsd$vectors[,1:5])
#'
#'

metadistance <- function(dist, meta, meta_n=NULL){

  # remove phyloseq features
  meta <- as.matrix(meta)

  # sample names
  names <- rownames(meta)
  # number of samples, metadata variables and comparisons
  nsamp <- nrow(meta)

  # index comparisons
  comp_i  <- combn(nsamp,2)
  row_i <- comp_i[2,]
  col_i <- comp_i[1,]

  # do comparison
  same <- (meta[row_i,] == meta[col_i,])
  same[is.na(same)] <- FALSE
  same[is.na(meta[row_i,]) & is.na(meta[col_i,])] <- TRUE

    # populate array
  compmat <- meta[row_i,] # fills matrix with values of one of the compared samples

  for(i in 1:ncol(meta)){
    tmplabels <- cbind(meta[row_i[!same[,i]],i], meta[col_i[!same[,i]],i])
    tmplabels <- t(apply(tmplabels, 1, sort))
    compmat[!same[,i],i] <- paste(tmplabels[,1], tmplabels[,2], sep="::")
    #compmat[!same[,i],i] <- paste(meta[row_i[!same[,i]],i], meta[col_i[!same[,i]],i], sep="::")
  }

  # set column names
  colnames(same) <- paste(colnames(meta), "_comp", sep = "")
  colnames(compmat) <- paste(colnames(meta), "_labels", sep = "")

  # numeric vectors
  if(is.null(meta_n)){
    # create data.frame
    output <- data.frame(sample1 = names[row_i], sample2 = names[col_i], dist = as.vector(dist), same, compmat)
  } else {
    meta_d <- abs(meta_n[row_i,] - meta_n[col_i,])
    # create data.frame
    output <- data.frame(sample1 = names[row_i], sample2 = names[col_i], dist = as.vector(dist), meta_d, same, compmat)
  }


  rownames(output) <- paste(names[row_i], names[col_i], sep="::") # rownames now represent both compared samples
  return(output)
}
