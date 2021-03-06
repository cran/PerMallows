% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PerMallows.R
\name{permutations.of}
\alias{permutations.of}
\title{Generate every permutation of perm.length item}
\usage{
permutations.of(perm.length, alert = TRUE)
}
\arguments{
\item{perm.length}{number of items in the permutation}

\item{alert}{optional ask for confirmation when the number of permtuations to show is very large}
}
\value{
A collection of every permutation of the specified number of items
}
\description{
This functions returns a matrix in thich each of rows
is a different permutation of the specified number of items
}
\examples{
permutations.of(3)
permutations.of(10)
}
