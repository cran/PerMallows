\name{insert.at}
\alias{insert.at}
\title{Insert operator}
\usage{
insert.at(perm, i, j)
}
\arguments{
  \item{perm}{a permutation}

  \item{i}{position of the permutation}

  \item{j}{position of the permutation}
}
\value{
the permutation in the input in which th eoperation has
been applied
}
\description{
Given a permutation and two positions i, j, move item in
position i to position j
}
\examples{
insert.at(c(1,2,3,4,5),5,2)
insert.at(c(1,2,3,4,5),2,5)
}

