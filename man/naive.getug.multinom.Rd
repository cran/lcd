\name{naive.getug.multinom} 
\alias{naive.getug.multinom}
\title{ A naive function to get an undirected graph for multinomial data }
\description{
  Learns an undirected independence graph from a given data set. The
  data are assumed to be multinomially distributed.
}
\usage{
naive.getug.multinom(freq.tb, p.value, method = "mkb")
}
\arguments{
  \item{freq.tb}{ an object of class \code{freq.tb}, the frequency table
  for the data.}
  \item{p.value}{ the thresholding \eqn{p}-value for each conditional
    independence test. }
  \item{method}{ a string to specify a method, see details. }
}
\details{
  Currently, there are three supported \code{method}'s. \code{"mkb"}
  stands for grow-shrink Markov blanket selection, where greedy Markov
  blanket selection is performed for each vertex and the blankets are the
  neighborhood for vertices of the output graph. The \code{"simple"}
  method simply performs conditional independence test for each vertex
  pair with all the rest vertices as the conditioning set. The
  \code{"fwd"} method uses the forward selection procedure described in
  Edwards (2000).
}
\value{
  The adjacency matrix of an undirected graph.
}
\references{
  Edwards, D. (2000). \emph{Introduction to Graphical Modelling.} 2nd
  Ed. Springer-Verlag, New York.
}
\author{ Zongming Ma and Xiangrui Meng}
\keyword{ multivariate }
\keyword{ graphs }
\keyword{ models }
