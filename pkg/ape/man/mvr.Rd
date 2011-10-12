\name{mvr}

\alias{mvr}

\title{Minimum Variance Reduction}

\description{
Phylogenetic tree 
construction based on the minimum variance reduction. Can be seen as a version of bio-NJ which is not restricted to the Poison model of variance.}

\usage{
mvr(X, V)
}

\arguments{
\item{X}{a distance matrix.}
\item{V}{a variance matrix.}

}

\value{
an object of class \code{"phylo"}.
}

\references{
Gascuel, O. (2000). Data model and classification by trees: the minimum variance
reduction (mvr) method. Journal of Classification, 17, 67–99.  
}

\author{Andrei Popescu \email{niteloserpopescu@gmail.com}}

\keyword{models}



