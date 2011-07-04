\name{triangMtdS}
\alias{triangMtdS}
\title{Tree Reconstruction from Incomplete Distances Based on the Triangles Method}
\usage{triangMtdS(X)}
\arguments{\item{X}{a distance matrix}}
\description{Fast distance-based costruction method. Should only be used when distance measures are fairly reliable or for small number of missing distances.}
\value{
an object of class \code{"phylo"}.
}
\references{\url{http://archive.numdam.org/ARCHIVE/RO/RO_2001__35_2/RO_2001__35_2_283_0/RO_2001__35_2_283_0.pdf}}
\author{C code and interfacing with \R by Andrei Popescu \email{niteloserpopescu@gmail.com}}