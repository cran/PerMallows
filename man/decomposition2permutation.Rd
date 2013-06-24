\name{decomposition2permutation}
\alias{decomposition2permutation}
\title{Get a permutation consistent with a decomposition vector}
\usage{
decomposition2permutation(vec, dist.name = "kendall")
}
\arguments{
  \item{vec}{the permutation}

  \item{dist.name}{optional the name of the distance. One
  of: kendall (default), cayley, hamming}
}
\value{
the distance decomposition vector of the given permutation
and distance
}
\description{
Given a distance decomposition vector and a distance name,
generate uniformly at random a permutation consistent with
the decomposition vector.
}
\examples{
permutation2decomposition(c(1,2,4,3,5), 'cayley')
}

