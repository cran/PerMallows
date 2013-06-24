\name{expectation.gmm}
\alias{expectation.gmm}
\title{Compute the expected distance, GMM under the Hamming distance}
\usage{
expectation.gmm(theta, dist.name = "kendall")
}
\arguments{
  \item{theta}{n dimensional real vector with the
  dispersion parameters}

  \item{dist.name}{optional name of the distance used in
  the GMM. One of: kendall (default), cayley, hamming}
}
\value{
expected distance decomposition vector
}
\description{
Compute the expected distance in the GMM under the Hamming
distance
}
\examples{
expectation.gmm(c(0.38, 0.44, 0.1, 0.2, 1, 0.1))
expectation.gmm(c(2,2,2,2),'hamming')
}

