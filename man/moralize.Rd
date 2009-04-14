\name{moralize} 
\alias{moralize}
\title{ Chain graph moralization }
\description{
  Computes the moral graph of a given chain graph.
} 
\usage{
moralize(amat)
}
\arguments{
  \item{amat}{ the adjacency matrix of the chain graph. }
}
\details{
  Joins the parents of every complex and undirect all the edges.
}
\value{
  The adjacency matrix of the moral graph.
}
\references{
  Lauritzen, S. L. (1996). \emph{Graphical Models.} Clarendon Press, Oxford.
}
\author{ Zongming Ma and Xiangrui Meng}
\examples{
moralize(toy.graph)
}
\keyword{ graphs }
