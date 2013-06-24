\name{lgmm}
\alias{lgmm}
\title{Learn a Generalized Mallows Model}
\usage{
lgmm(sample, sigma_0_ini = identity.permutation(dim(sample)[2]),
  dist.name = "kendall", estimation = "approx")
}
\arguments{
  \item{sample}{the matrix with the permutations to
  estimate}

  \item{sigma_0_ini}{optional the initial guess for the
  consensus permutation}

  \item{dist.name}{optional name of the distance used by
  the GMM. One of: kendall (default), cayley, hamming}

  \item{estimation}{optional select the approximated or the
  exact. One of: approx, exact}
}
\value{
a list with the parameters of the estimated distribution:
the mode and the dispersion parameter vector
}
\description{
Learn the parameter of a sample of n permutations comming
from a Generalized Mallows Model (GMM).
}
\examples{
sample <- matrix(c(1,2,3,4, 1,4,3,2, 1,2,4,3), nrow = 3, ncol = 4, byrow = TRUE)
lgmm(sample, dist.name='cayley', estimation='approx')
}

