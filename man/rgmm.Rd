\name{rgmm}
\alias{rgmm}
\title{Sample a Generalized Mallows Model}
\usage{
rgmm(n, sigma0, theta, dist.name = "kendall", sampling.method = "gibbs")
}
\arguments{
  \item{n}{the number of permutations to be generated}

  \item{sigma0}{central permuation of the GMM}

  \item{theta}{dispersion parameter vector of the GMM}

  \item{dist.name}{optional used name of the distance used
  in the GMM. One of: kendall (default), cayley, hamming}

  \item{sampling.method}{optional name of the sampling
  algorithm. One of: multistage, gibbs (default)}
}
\value{
a sample of permutations from the specified ditribution
}
\description{
Generate a sample of n permutations from a Generalized
Mallows Model (GMM).
}
\examples{
rgmm(2,c(1,2,3,4,5),c(1,1,1,1,1),'hamming')
}

