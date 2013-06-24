\name{r.dist.d}
\alias{r.dist.d}
\title{Generate a collection of permutations at a given distance}
\usage{
r.dist.d(n, perm.length, dist.value, dist.name = "kendall")
}
\arguments{
  \item{n}{number of permutations in the sample}

  \item{perm.length}{number of items in the permutations}

  \item{dist.value}{distance value}

  \item{dist.name}{distance name. One of: kendall
  (default), cayley, hamming, ulam}
}
\value{
a sample of permutations at the given distance
}
\description{
Given a number of permutations, the number of items in the
permutations, a distance value and a distance name,
generate a sample of permutations with the specified length
at the given distance
}
\examples{
r.dist.d(1, 4, 2, 'ulam')
}

