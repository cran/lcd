\name{is.chaingraph}  
\alias{is.chaingraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Chain graph verification }
\description{
  Checks if a given graph is a chain graph.
}
\usage{
is.chaingraph(amat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{amat}{ the adjacency matrix of the graph with dimnames. }
}
\value{
  \item{result }{a logical value, \code{TRUE} if the given graph is a valid
    chain graph and \code{FALSE} otherwise.}
  \item{vert.order }{a topological order of the graph vertices.}
  \item{chain.size }{a vector indicating the size of each chain
    component in the graph, whose order corresponds to \code{vert.order}.}
}
\references{
  Cowell, R. G., Dawid, A. P., Lauritzen, S. L. and Spiegelhalter,
  D. J. (1999) \emph{Probabilistic Networks and Expert
    Systems}. Springer-Verlag, New York.
}
\author{ Zongming Ma and Xiangrui Meng}
\examples{
data(lcd)
is.chaingraph(toy.graph)
}
\keyword{ multivariate }
\keyword{ graphs }
\keyword{ models }
