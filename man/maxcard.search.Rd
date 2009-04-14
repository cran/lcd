\name{maxcard.search} 
\alias{maxcard.search}
\title{ Maximum cardinality search }
\description{
  Performs a maximum cardinality search on an undirected graph to
  determine whether it is triangulated.
}
\usage{
maxcard.search(amat)
}
\arguments{
  \item{amat}{ the adjacency matrix of the undirected graph. }
}
\value{
  \item{is.triangulated }{a logical value indicating whether the input
    graph is triangulated or not.}
  \item{perfect.numbering }{a perfect numbering of the vertices.}
  \item{card }{ number of unlabeled neighbors when labeling each
    variable, with order compatible to the perfect numbering.}
  \item{pi.record }{a record of unlabeled neighbors during the execution
  of the algorithm.}
}
\references{
  Tarjan, R. E. and Yannakakis, M. (1984). Simple linear-time algorithms
  to test chordality of graphs, test acyclicity of hypergraphs, and
  selectively reduce acyclic hypergraphs. \emph{SIAM Journal on
    Computing}, \bold{13}, 566-79.
}
\author{ Zongming Ma and Xiangrui Meng}
\note{
  Only the \code{is.triangulated} and \code{perfect.numbering} are
  supposed to be of interest to the users. The other two output
  components are mainly for track purpose.
}
\keyword{ graphs }
