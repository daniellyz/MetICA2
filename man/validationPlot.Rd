% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validationPlot.R
\name{validationPlot}
\alias{validationPlot}
\title{Plots to decide number of clusters}
\usage{
validationPlot(M1, cluster_index = T)
}
\arguments{
\item{M1}{The entire list object generated from the function MetICA(X, pcs = 15...)}

\item{cluster_index}{Boolean object. TRUE if geometric indices for clusterinig quality needs to be calculted. The calculation provides an additional criterion for cluster number selection but can be time-consuming.}
}
\value{
Return an overall graph with multiple plots. If cluster_index = TRUE, also return a list object containing clustering quality score.
}
\description{
The function generates multiple plots that help users to decide number of clusters or MetICA components.
}
