\name{lmm}
\alias{lmm}
\title{Learn a Mallows Model}
\usage{
lmm(sample, sigma_0_ini = identity.permutation(dim(sample)[2]),
  dist.name = "kendall", estimation = "approx", disk = FALSE)
}
\arguments{
  \item{sample}{the matrix with the permutations to
  estimate}

  \item{sigma_0_ini}{optional the initial guess for the
  consensus permutation}

  \item{dist.name}{optional the name of the distance used
  by the model. One of: kendall (default), cayley, hamming,
  ulam}

  \item{estimation}{optional select the approximated or the
  exact. One of: approx, exact}

  \item{disk}{optional can only be true if estimating a MM
  under the Ulam distance. Insted of generating the whole
  set of SYT and count of permutations per distance, it
  loads the info from a file in the disk}
}
\value{
a list with the parameters of the estimated distribution:
the mode and the dispersion parameter
}
\description{
Learn the parameter of a sample of n permutations comming
from a Mallows Model (MM).
}
\examples{
sample <- matrix(c(1,2,3,4, 1,4,3,2, 1,2,4,3), nrow = 3, ncol = 4, byrow = TRUE)
lmm(sample, dist.name='ulam', estimation='approx')
}

