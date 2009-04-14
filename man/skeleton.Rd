\name{skeleton} 
\alias{skeleton}
\title{ Graph skeleton }
\description{
  Returns the skeleton of a graph. The graph can be undirected, directed
  or mixed.
}
\usage{
skeleton(amat)
}
\arguments{ 
  \item{amat}{ the adjacency matrix of the graph. }
}
\value{
  The adjacency matrix of the graph skeleton.
}
\author{ Zongming Ma and Xiangrui Meng}
\examples{
skeleton(toy.graph)
}
\keyword{ graphs }

