\name{draw} 
\alias{draw}
\title{ Draw graph }
\description{
  Draws a graph from its adjacency matrix.
}
\usage{
draw(amat, plain = TRUE)
}
\arguments{
  \item{amat}{ the adjacency matrix of the graph. }
  \item{plain}{ logical value. If \code{TRUE}, then a plain plot is
    drawn for the graph; if \code{FALSE}, a tlkck plot is drawn.}
}
\details{
  The function converts the adjacency matrix to an \code{igraph} object
  and uses the \code{plot.igraph}/\code{tkplot} function in
  \code{igraph} package to draw 
  the graph. Special care is paid to make the directed/undirected edges
  displayed right. Finally, the function returns the vertex list of the
  graph for comparison to the vertex numbering used in the plot.
}
\value{
  The vertex list. The main output is the side effect of the function.
}
\author{ Zongming Ma and Xiangrui Meng}
\keyword{ graphs }
\keyword{ iplot }
