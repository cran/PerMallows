\name{inversion.at}
\alias{inversion.at}
\title{Inversion operator}
\usage{
inversion.at(perm, i)
}
\arguments{
  \item{perm}{a permutation}

  \item{i}{position of the permutation}
}
\value{
the permutation in the input with an inversion at the
specified position
}
\description{
Given a permutation and a position, swap positions i and
i+1
}
\examples{
inversion.at(c(1,2,3,4,5),2)
}

