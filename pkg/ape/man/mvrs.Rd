\name{mvrs}

\alias{mvrs}

\title{Minimum Variance Reduction*}

\description{MVR version for incomplete distances. See reference for details about the fS parameter}

\usage{
mvr(X, V, fS=15)
}

\arguments{
\item{X}{a distance matrix.}
 \item{V}{a variance matrix.} 

\item{fS}{agglomeration criterion parameter}}

\value{
an object of class \code{"phylo"}.
}

\references{
  Criscuolo, A. and Gascuel, O. (2008). Fast NJ-like algorithms to deal with incomplete
distance matrices. BMC Bioinformatics, 9.
}

\author{Andrei Popescu \email{niteloserpopescu@gmail.com}}

\keyword{models}



