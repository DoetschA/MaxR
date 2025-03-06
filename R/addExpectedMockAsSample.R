#' Add theoretical composition of a mock community to phyloseq objects
#'
#' @param physeq a phyloseq object containing your microbiome data set.
#' @param mockName how the sample should be called (e.g., in graphs)
#' @param definition path to the file defining the theoretical composition
#' @param fillFields which variables (metadata) in sample_data(physeq) should be filled with mockName. All the other fields will be filled with NA.
#' @export
#' @return Phyloseq object with added mock community
#' @examples
#' # Add zymo mock composition to objeckt physeq_mock
#' physeq_mock <- addExpectedMockAsSample(physeq_mock, mockName = "zymo_expected", definition = "/data/db/mock_communities/zymo_16s.definition", fillFields = c("sampleType", "sampleName", "templateName"))
#'

addExpectedMockAsSample <- function(physeq, mockName, definition = "/data/db/mock_communities/zymo_16s.definition", fillFields){

  require(phyloseq)
  # read expected mock composition from definition file
  mock.def <- read.csv(file = definition, header = T)
  mock.otu <- data.frame(mock.def$percent * 10)
  colnames(mock.otu) <- mockName
  rownames(mock.otu) <- mock.def$mockID

  #browser()

  # prepare taxonomy
  mock.tax <- as.matrix(mock.def[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")])
  rownames(mock.tax) <- mock.def$mockID
  if(dim(tax_table(physeq))[2] == 6){ # remove Species if absent from physeq object
    mock.tax <- mock.tax[,1:6]
  }

  # prepare sample data
  exp_dat <- sample_data(physeq)[1,]
  rownames(exp_dat) <- mockName
  for(i in 1:length(exp_dat)){ exp_dat[i] <- NA}
  exp_dat[,fillFields] <- mockName

  # merge expected data with phyloseq object (merge_phyloseq does not work well here!)
  physeq_exp <- phyloseq(exp_dat, otu_table(mock.otu, taxa_are_rows = T), tax_table(mock.tax))

  # merge otu tables
  otutab <- rbind(cbind(otu_table(physeq), matrix(0, nrow = nsamples(physeq), ncol = ntaxa(physeq_exp))),
                  cbind(matrix(0, nrow = nsamples(physeq_exp), ncol = ntaxa(physeq)), t(otu_table(physeq_exp))))
  colnames(otutab) <- c(taxa_names(physeq), taxa_names(physeq_exp))

  # merge taxonomy
  taxtab <- rbind(tax_table(physeq), tax_table(physeq_exp))

  # merge sample data
  samdat <- rbind(sample_data(physeq), sample_data(physeq_exp))

  # new phyloseq object
  physeq <- phyloseq(otu_table(otutab, taxa_are_rows = F), tax_table(taxtab), sample_data(samdat))

  return(physeq)
}
