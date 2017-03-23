\name{polyCluster}
\alias{polyCluster}
\title{
Run the polyCluster pipeline.
}
\description{
Runs the entire polyCluster pipeline with the given dataset and clustering algorithms.
}
\usage{
polyCluster(filename, clusterAlg = c("hc", "pm", "km", "nmf"),
    maxK = 7, reps = 100, phenoFile = NULL, ref = "",
    nmfData = NULL, interactive = FALSE)
}
\arguments{
  \item{filename}{
Character string. File name of preprocessed and normalized expression values with genes in rows and samples in columns.
}
  \item{clusterAlg}{
Character vector. Any combination of "hc" (hierarchical clustering), "pm" (partitioning around medoids), "km" (k-means) or "nmf" (nonnegative matrix factorization).
}
  \item{maxK}{
Integer > 2. The maximum number of clusters of samples to evaluate.
}
  \item{reps}{
Integer. The number of resampling (for "hc", "pm" and "km") or starting seed ("nmf") iterations.n.
}
  \item{phenoFile}{
Character string. File name of known phenotypes, with sample name in first column and phenotype in second column.
}
  \item{ref}{
Character string. Reference with which to name the output of the analysis.
}
  \item{nmfData}{
Character string. File path to the result of a previous NMF clustering by this function, in order to save time.
}
  \item{interactive}{
Logical. If FALSE, networks will be plotted with an automatic layout. If TRUE, a users can customise the layout through tkplot.
}
}
\details{
Runs the polyCluster pipeline and outputs files detailing clustering results from each algorithm, as well as reconcilation between the algorithms' results by approaches detailed in Eason et al. (2017).
}
\value{
  \item{NULL}
}
\references{
Eason et al. (2017)
}
\author{
Kate Eason
}
