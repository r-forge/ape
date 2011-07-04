\name{bionjs}
\alias{bionjs}
\title{Tree Reconstruction from Incomplete Distances Based on bio-NJ*}
\usage{bionjs(X)}
\arguments{\item{X}{a distance matrix}}
\description{Reconstructs a phylogenetic tree from a distance matrix with possibly missing values. Missing values represented by either NA or any negative number.}
\value{
an object of class \code{"phylo"}.
}
\references{\url{http://www.biomedcentral.com/1471-2105/9/166}}
\author{C code and interfacing with \R by Andrei Popescu \email{niteloserpopescu@gmail.com}}