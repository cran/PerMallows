\name{count.perms.fixed.points}
\alias{count.perms.fixed.points}
\title{Count permutations with a given number of fixed points}
\usage{
count.perms.fixed.points(perm.length, fixed)
}
\arguments{
  \item{perm.length}{number of items in the permutations}

  \item{fixed}{the number of fixed points}
}
\value{
the number of permutations in which the number of fixed
points is equal to the specified argument
}
\description{
Given the number of items in the permutations and the
number of fixed points f, how many permutations are there
with f fixed points?
}
\examples{
count.perms.fixed.points(4,1)
}

