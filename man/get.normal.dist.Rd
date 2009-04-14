\name{get.normal.dist} 
\alias{get.normal.dist}
\title{ Random normal distribution generation from a chain graph }
\description{
  Generates a random normal distribution from a given chain graph structure.
}
\usage{
get.normal.dist(amat)
} 
\arguments{
  \item{amat}{ the adjacency matrix of a chain graph. }
}
\details{
  The function returns a matrix that encodes a mean \eqn{0} block-recursive
  regression model (See Wermuth (1992)).
}
\value{
  A matrix encoding the randomly generated block-recursive regression model.
}
\references{
  Wermuth, N. (1992). Block-recursive regression equations
  (with discussions). \emph{Revista Brasileira de Probabilidade e
    Estatistica,} \bold{6}, 1-56.
}
\author{ Zongming Ma and Xiangrui Meng}
\keyword{ multivariate }
\keyword{ distribution }
