\name{get.multinom.dist} 
\alias{get.multinom.dist}
\title{ Random multinomial distribution generation from a chain graph }
\description{
  Generates a random multinomial distribution from a given chain graph structure. 
}
\usage{
get.multinom.dist(amat, n.state, alpha, beta)
}
\arguments{
  \item{amat}{ the adjacency matrix of a chain graph. }
  \item{n.state}{ a vector of positive integers indicating the desired
    number of states for each variable in the chain graph. }
  \item{alpha, beta}{ shape parameters for Beta distribution in simulating the
    potentials.}
}
\details{
  The function returns a list that encodes a block-recursive conditional
  distribution over the chain graph structure.
}
\value{
  A sequence of block-recursive conditional distributions generated
  according to the graph structure. 
}
\references{
  Cowell, R. G., Dawid, A. P., Lauritzen, S. L. and Spiegelhalter,
  D. J. (1999) \emph{Probabilistic Networks and Expert
    Systems}. pp.77-79. Springer-Verlag, New York.
}
\author{ Zongming Ma and Xiangrui Meng}
\keyword{ multivariate }
\keyword{ distribution }
