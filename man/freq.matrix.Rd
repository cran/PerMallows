\name{freq.matrix}
\alias{freq.matrix}
\title{Compute the frequency matrix}
\usage{
freq.matrix(perm)
}
\arguments{
  \item{perm}{a permutation or a collection of them}
}
\value{
a matrix with n rows and n columns with the proportion of
the permutations in the input that have each item in each
position
}
\description{
Compute the first order marginal probability. In other
words, given at least one permutation, calculate the
proportion of them that have each item in each position
}
\examples{
freq.matrix(c(1,3,2,4,5))
}

