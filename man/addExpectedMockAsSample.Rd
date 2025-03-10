% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addExpectedMockAsSample.R
\name{addExpectedMockAsSample}
\alias{addExpectedMockAsSample}
\title{Add theoretical composition of a mock community to phyloseq objects}
\usage{
addExpectedMockAsSample(
  physeq,
  mockName,
  definition = "/data/db/mock_communities/zymo_16s.definition",
  fillFields
)
}
\arguments{
\item{physeq}{a phyloseq object containing your microbiome data set.}

\item{mockName}{how the sample should be called (e.g., in graphs)}

\item{definition}{path to the file defining the theoretical composition}

\item{fillFields}{which variables (metadata) in sample_data(physeq) should be filled with mockName. All the other fields will be filled with NA.}
}
\value{
Phyloseq object with added mock community
}
\description{
Add theoretical composition of a mock community to phyloseq objects
}
\examples{
# Add zymo mock composition to objeckt physeq_mock
physeq_mock <- addExpectedMockAsSample(physeq_mock, mockName = "zymo_expected", definition = "/data/db/mock_communities/zymo_16s.definition", fillFields = c("sampleType", "sampleName", "templateName"))

}
