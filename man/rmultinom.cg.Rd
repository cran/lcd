\name{rmultinom.cg} 
\alias{rmultinom.cg}
\title{ Random multinomial sample from a chain graph }
\description{
  Generates a desired number of multinomial random samples from a given
  chain graph structure and a sequence of compatible block-recursive
  conditional distributions.
} 
\usage{
rmultinom.cg(n, amat, distn)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ the intended sample size, should be at least \eqn{1}. }
  \item{amat}{ the adjacency matrix of a chain graph. }
  \item{distn}{ a sequence of block-recursive conditional distributions,
  e.g. the one returned by \code{\link{get.multinom.dist}}. }
}
\value{
  An \eqn{n} by \eqn{p} matrix with each row corresponding to
  a random sample.
}
\references{
  Cowell, R. G., Dawid, A. P., Lauritzen, S. L. and Spiegelhalter,
  D. J. (1999) \emph{Probabilistic Networks and Expert
    Systems}. pp.77-79. Springer-Verlag, New York.
}
\author{ Zongming Ma and Xiangrui Meng}
\keyword{ multivariate }
\keyword{ datagen }
\keyword{ distribution }
