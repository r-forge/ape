\name{SDM}
\alias{SDM}
\title{construction of consensus distance matrix based on SDM}
\usage{SDM(...)}
\arguments{a list of 2n elements (n>1), the first n elements are the distance matrices, the next n elements are the sequence length from which the matrices have been estimated (can be seen as a degree of confidence in matrices)}
\description{Reconstructs a consensus distance matrix from a set of input distance matrices on overlapping sets of taxa. Potentially missing values in the supermatrix are represented by NA. An error is returned if the input distance matrices can not resolve to a consensus matrix.}
\value{
a 2-element list containing a distance matrix labelled by the union of the set of taxa of the input distance matrices, and a variance matrix associated to the returned distance matrix.}
\references{Criscuolo A, Berry V, Douzery EJP, Gascuel O: SDM: a fast distance-based approach for (super)tree building in phylogenomics.
Syst Biol 2006, 55(5):740-755.}
\author{Andrei Popescu \email{niteloserpopescu@gmail.com}}