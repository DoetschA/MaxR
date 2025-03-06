#' Correct OTU abundances according to rRNA copy numbers
#'
#' THis script uses information from the rrnDB (https://rrndb.umms.med.umich.edu/).
#' Each taxonomic rank present in the input data is screened for taxa present in the data base.
#' If the name is found, the mean or median (depending on parameters) is used to normalize the read counts for this taxon.
#' last changed: A. Dötsch, 21.12.2021
#'
#' @param physeq A phyloseq object containg at least abundance (otu_table) and taxonomy (tax_table) data.
#' @param dbpath Path to the rrndb file.
#' @param median if true, the median copy number is used on the genus level; if FALSE, mean is used for all tax levels (default)
#' @return phyloseq object with copy-number corrected OTU counts
#' @export
#' @examples
#' rrnDB_path <- "/home/bioinf/scripts/rrnDB-5.7_pantaxa_stats_RDP_Silva.tsv"
#' physeq_norm <- rRNA_copy_number_correction(physeq, dbpath = rrnDB_path)
#'


rRNA_copy_number_correction <- function(physeq, dbpath, use_median = F){
  require(phyloseq)

  # load data from database
  rrnDB <- read.table(file = dbpath, sep = "\t", header = T)

  # get taxonomic table from phyloseq object
  taxtab <- tax_table(physeq)

  # prepare read copy number (RCN) matrix
  rcn <- data.frame(asv_id = row.names(taxtab),
                    median_rcn = 1,
                    mean_rcn = 1,
                    std_rcn = 0,
                    lowest_rank = NA,
                    lowest_rank_name = NA)
  row.names(rcn) <- rcn[,1]
  rcn <- rcn[,-1]

  # identify database hits
  for(taxlevel in colnames(taxtab)){
    # check if taxlevel exists in database (e.g. Species is not present!)
    if(!any(rrnDB$rank == tolower(taxlevel))){ next }

    utax <- as.character(unique(taxtab[,taxlevel])) # unique taxa
    for(i in 1:length(utax)){
      dbhit  <- rrnDB$name == utax[i]
      #prevent double hits in database (e.g. for Chloroplast)
      if(sum(dbhit) > 1){dbhit <- which(dbhit)[1]}
      taxhit <- taxtab[,taxlevel] == utax[i]
      if(any(dbhit)){
        rcn[taxhit, c("median_rcn", "mean_rcn", "std_rcn")] <- rrnDB[dbhit, c("median", "mean", "stddev")]
        rcn[taxhit, "lowest_rank"] <- taxlevel
        rcn[taxhit, "lowest_rank_name"] <- as.character(taxtab[i, taxlevel])
      }
    }
  }

  # normalize read counts
  physeq_norm <- physeq
  normfactor <- rcn[,"mean_rcn"]
  if(use_median){
    normfactor[rcn[,"lowest_rank"] == "Genus"] <- rcn[rcn[,"lowest_rank"] == "Genus", "median_rcn"]
  }

  if(taxa_are_rows(physeq)){
    for(i in 1:nsamples(physeq)){
      otu_table(physeq_norm)[,i] <- otu_table(physeq)[,i]/normfactor
    }
  } else {
    for(i in 1:nsamples(physeq)){
      otu_table(physeq_norm)[i,] <- otu_table(physeq)[i,]/normfactor
    }
  }

  isrelative <- !is.integer(otu_table(physeq))
  if(!isrelative){
    otu_table(physeq_norm) <- ceiling(otu_table(physeq_norm)) # ceiling, to prevent zeros
  }

  return(physeq_norm)
}
