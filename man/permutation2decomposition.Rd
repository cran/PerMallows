\name{permutation2decomposition}
\alias{permutation2decomposition}
\title{Get the decomposition vector}
\usage{
permutation2decomposition(perm, dist.name = "kendall")
}
\arguments{
  \item{perm}{the permutation}

  \item{dist.name}{optional the name of the distance. One
  of: kendall (default), cayley, hamming}
}
\value{
the distance decomposition vector of the given permutation
and distance
}
\description{
Given a permutation and a distance name generate the
decomposition vector
}
\examples{
permutation2decomposition(c(1,2,4,3,5), 'cayley')
}

