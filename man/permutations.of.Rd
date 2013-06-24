\name{permutations.of}
\alias{permutations.of}
\title{Generate every permutation}
\usage{
permutations.of(perm.length, alert = TRUE)
}
\arguments{
  \item{perm.length}{number of items in the permutation}

  \item{alert}{optional ask for confirmation when the
  number of permtuations to show is very large}
}
\value{
collection of every permutation of the specified number of
items
}
\description{
Generate every permutation of the specified number of items
}
\examples{
permutations.of(3)
permutations.of(10)
}

