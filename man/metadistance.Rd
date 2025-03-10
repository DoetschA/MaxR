% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metadistance.R
\name{metadistance}
\alias{metadistance}
\title{Combines a distance matrix with metadata to a data.frame that can be used for ggplot graphs to e.g. group dissimilarities by certain factors (e.g. belonging to the same group)}
\usage{
metadistance(dist, meta, meta_n = NULL)
}
\arguments{
\item{dist}{distance matrix of class dist()}

\item{meta}{data.frame containing factors (values will be assigned specfic levels if both samples have the same value and NA otherwise)}

\item{meta_n}{(optional) data.frame containing numeric metadata (values will be absolute differences)}
}
\value{
A data.frame in a format suitable to diplay distance values for certain types of comparisons.
Output columns are:
\itemize{
\item 1st, 2nd: names of the two samples that are compared
\item 3rd: distance (as defined in 'dist') between the two samples
\item _comp - TRUE if variable is same in both samples
\item _labels - displaying compared samples.
}
}
\description{
Combines a distance matrix with metadata to a data.frame that can be used for ggplot graphs to e.g. group dissimilarities by certain factors (e.g. belonging to the same group)
}
\examples{
metadistance(dist_jsd, meta = sample_data(physeq_feces)[,c("proband_ID","kit_purif","lab_purif","MiSeq")], meta_n = ord_pcoa_jsd$vectors[,1:5])


}
