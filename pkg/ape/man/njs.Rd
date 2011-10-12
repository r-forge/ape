\name{njs}
\alias{njs}
\title{Tree Reconstruction from Incomplete Distances Based on NJ*}
\usage{njs(X, fS=15)}
\arguments{\item{X}{a distance matrix}\item{fS}{optional agglomeration criterion parameter}}
\description{Reconstructs a phylogenetic tree from a distance matrix with possibly missing values. Missing values represented by either NA or any negative number. See reference for the meaning of fS.}
\value{
an object of class \code{"phylo"}.
}
\references{Criscuolo, A. and Gascuel, O. (2008). Fast NJ-like algorithms to deal with incomplete
distance matrices. BMC Bioinformatics, 9.}
\author{C code and interfacing with \R by Andrei Popescu \email{niteloserpopescu@gmail.com}}