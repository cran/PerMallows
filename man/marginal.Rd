\name{marginal}
\alias{marginal}
\title{Compute the marginal probability, GMM under the Hamming distance}
\usage{
marginal(h, theta)
}
\arguments{
  \item{h}{n dimensional distance decomposition vector
  where h_j = 0 means that $j$ is a fixed point, h_j = 1
  means that $j$ is an unfixed point and otherwise $j$ is
  not known}

  \item{theta}{n dimensional distance decomposition vector
  with the dispersion parameters}
}
\value{
marginal probability
}
\description{
Compute the marginal probability, GMM under the Hamming
distance, of a distance decomposition vector for which some
positions are known and some are not
}
\examples{
marginal(c(1,0,1,-1,-1), c(0.1, 0.3, 0.7, 0.1, 1))
marginal(c(-1,0,1,-1,-1,0), c(0.1, 0.3, 0.7, 0.1, 0.7, 1))
}

