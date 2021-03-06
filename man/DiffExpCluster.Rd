% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{DiffExpCluster}
\alias{DiffExpCluster}
\title{Differential expression for different clusters in specific clusters
Selct the differential gene among clusters by t test}
\usage{
DiffExpCluster(Data, label, i)
}
\arguments{
\item{Data}{the scRNA-seq matrix}

\item{label}{the clusters label of the imput data}
}
\value{
P the p_values matrix, the rows correspond to the clusters, the cloumns correspond to the genes

DEgenes, the differential genes for each clusters. the each row is the DE genes for each cluster.
}
\description{
Differential expression for different clusters in specific clusters
Selct the differential gene among clusters by t test
}
