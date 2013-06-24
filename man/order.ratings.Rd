\name{order.ratings}
\alias{order.ratings}
\title{Convert rating to permutation}
\usage{
order.ratings(ratings)
}
\arguments{
  \item{ratings}{a matrix in which each row is a vector of
  ratings of several items}
}
\value{
a matrix in which each row is the correspondibg permutation
of the items
}
\description{
This function is given a collection of ratings and converts
each row to a permutation
}
\examples{
order.ratings(c(0.1, 4, 0.5, -4))
}

