\name{compose}
\alias{compose}
\title{Compose permutations}
\usage{
compose(perm1, perm2)
}
\arguments{
  \item{perm1}{a permtuation or a collection of
  permutations}

  \item{perm2}{a permtuation or a collection of
  permutations}
}
\value{
the composition of the permutations
}
\description{
This function composes two given permtuations. One of the
arguments can be a collection of permutations, but not both
at a time. In this case, every permutation in the
collection is composed with the other argument
}
\examples{
compose(c(3,1,2,4), c(4,1,3,2))
}

