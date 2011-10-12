\name{bionjs}
\alias{bionjs}
\title{Tree Reconstruction from Incomplete Distances Based on bio-NJ*}
\usage{bionjs(X, fS=15)}
\arguments{\item{X}{a distance matrix}\item{fS}{an optional agglomeration criterion parameter}}
\description{Reconstructs a phylogenetic tree from a distance matrix with possibly missing values. Missing values represented by either NA or any negative number. See reference for meaning of fS.}
\value{
an object of class \code{"phylo"}.
}
\references{Criscuolo, A. and Gascuel, O. (2008). Fast NJ-like algorithms to deal with incomplete
distance matrices. BMC Bioinformatics, 9.}
\author{C code and interfacing with \R by Andrei Popescu \email{niteloserpopescu@gmail.com}}