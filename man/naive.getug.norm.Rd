\name{naive.getug.norm} 
\alias{naive.getug.norm}
\title{ A naive function to get an undirected graph for normal data }
\description{
  Learns an undirected independence graph from a given data set. The
  data are assumed to be normally distributed.
}
\usage{
naive.getug.norm(data, p.value)
}
\arguments{
  \item{data}{ a data matrix with rows corresponding to observations and
  columns corresponding to random variables.}
  \item{p.value}{ the thresholding \eqn{p}-value for each conditional
    independence test. }
}
\details{
  For each pair of random variables, the function performs a conditional
  independence using the partial correlation coefficient. If the p-value
  of the test is smaller than the given threshold, then there will be an
  edge on the output graph. The function essentially uses the global
  Markov property of the undirected graph.
}
\value{
  The adjacency matrix of the computed undirected graph, with dimnames
  the corresponding random variables.
}
\references{
  Lauritzen, S. L. (1996). \emph{Graphical Models.} Clarendon Press, Oxford.
}
\author{ Zongming Ma and Xiangrui Meng}
\keyword{ multivariate }
\keyword{ graphs }
\keyword{ models }
