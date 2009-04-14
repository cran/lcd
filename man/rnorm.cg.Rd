\name{rnorm.cg} 
\alias{rnorm.cg}
\title{ Random normal sample from a chain graph }
\description{
  Generates a desired number of normal random samples from a given chain graph
  structure and a given block-recursive regression system compatible
  with the chain graph. 
}
\usage{
rnorm.cg(n, amat, Bstar)
}
\arguments{
  \item{n}{ the intended sample size, should be at least \eqn{1}. }
  \item{amat}{ the adjacency matrix of a chain graph. }
  \item{Bstar}{ the matrix representation of the block-recursive
    regression structure, e.g., the one returned by \code{\link{get.normal.dist}}. }
}
\details{
  The function uses a mean \eqn{0} block-recursive
  regression model (See Wermuth (1992)), which is recorded in
  \code{Bstar} and
  normal random samples are generated from the specified block-recursive
  regression model.
}
\value{
  An \eqn{n} by \eqn{p} matrix with each row corresponding to
  a random sample.
}
\references{
   Wermuth, N. (1992). Block-recursive regression equations
  (with discussions). \emph{Revista Brasileira de Probabilidade e
  Estatistica,} \bold{6}, 1-56.
}
\author{ Zongming Ma and Xiangrui Meng}
%\examples{
%data(lcd) 
%tgdata <- rnorm.cg(20000, toy.graph)$data
%}
\keyword{ multivariate }
\keyword{ datagen }
\keyword{ distribution }
