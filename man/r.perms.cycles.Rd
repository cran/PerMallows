\name{r.perms.cycles}
\alias{r.perms.cycles}
\title{Generate a collection of permutations with a given number of cycles}
\usage{
r.perms.cycles(n, perm.length, cycles)
}
\arguments{
  \item{n}{number of permutations in the sample}

  \item{perm.length}{number of items in the permutations}

  \item{cycles}{number of cycles}
}
\value{
a sample of permutations with the specified length and
number of cycles
}
\description{
Given a number of permutations, the number of items in the
permutations and number of cycles, generate a sample of
permutations with the specified length and number of cycles
}
\examples{
r.perms.cycles(1, 4, 2)
}

