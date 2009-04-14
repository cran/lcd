\name{tri.ug} 
\alias{tri.ug}
\title{ Triangulation of an undirected graph }
\description{
  Triangulates an undirected graph to a chordal graph.
}
\usage{
tri.ug(amat) 
}
\arguments{
  \item{amat}{ the adjacency matrix of undirected graph. }
}
\details{
  The function implements the `One-step look ahead triangulation'
  algorithm described in Cowell, et al (1999). The criterion \eqn{c(v)}
  is chosen to be the number of edges needed to be added to the graph if
  vertex \eqn{v} is chosen to be labelled.
}
\value{
  The adjacency matrix of the triangulated graph.
}
\references{
  Cowell, R. G., Dawid, A. P., Lauritzen, S. L. and Spiegelhalter,
  D. J. (1999) \emph{Probabilistic Networks and Expert
    Systems}. Springer-Verlag, New York.
}
\author{ Zongming Ma and Xiangrui Meng}
\keyword{ graphs }

