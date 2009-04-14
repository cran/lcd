\name{learn.graph} 
\alias{learn.complex.norm}
\alias{learn.complex.multinom}
\alias{learn.mec.norm}
\alias{learn.mec.multinom}
\alias{learn.skeleton.norm}
\alias{learn.skeleton.multinom}
\alias{learn.v}
\title{ Graph learning functions }
\description{
  Learns chain graphs structures via decomposition algorithms with
  normal or multinomial data. 
}
\usage{

learn.complex.norm(skel, cov, n, p.value)

learn.mec.norm(tree, cov, n, p.value, method = "CG")

learn.skeleton.norm(tree, cov, n, p.value, drop = TRUE)

learn.complex.multinom(skel, freq.tb, p.value)

learn.mec.multinom(tree, freq.tb, p.value, method = "CG")

learn.skeleton.multinom(tree, freq.tb, p.value, drop = TRUE)

learn.v(skel, tree)
}
\arguments{

  \item{skel}{ the object returned by \code{learn.skeleton.norm} function. }
  \item{tree}{ an object of class \code{sep.tree}, e.g.,
    the separation tree obtained via \code{ug.to.jtree}. }
  \item{cov}{ the covariance matrix of the data }
  \item{n}{ number of random samples used to obtain the covariance
    matrix. }
  \item{freq.tb}{ the frequency table of the data, should be an object
    from class \code{freq.tb}.}
  \item{p.value}{ thresholding \eqn{p}-value for conditional
    independence test during the learning procedure. }
  \item{method}{ a character string, must be one of "CG", "DAG" and
    "UG", telling the alogrithm whether the underlying graph is a DAG,
    undirected graph or more generally a chain graph.}
  \item{drop}{logical value, with default=\code{TRUE}, whether to drop possibly
    extra edges after recovering and combining local skeletons.}
}
\details{
  \code{learn.mec.norm} is the wrapper function for
  \code{learn.skeleton.norm}, \code{learn.complex.norm} and \code{learn.v}
  for normal data,  which is intended to be mostly called by the
  user. \code{learn.mec.multinom} is the counterpart for discrete
  data. The user specifies
  whether she/he is learning an undirected graph (Markov network), a DAG
  (Bayesian network) or in general a chain graph by supplying the
  specific string to the \code{method} argument. The function itself
  will utilize the correct helper functions to perform the structural
  learning and return the pattern (representing the Markov equivalent
  class) of the graph.

  \code{learn.skeleton.norm} and \code{learn.skeleton.multinom} try to learn the skeleton of the
  graph. \code{learn.v}, \code{learn.complex.norm} (or
  \code{learn.complex.multinom}) are for v-structure/complex discovery
  in DAG/chain graph respectively.

  Some parts of the code of these learning functions are adapted from
  the code for \code{pcAlgo} function in the \code{R} package
  \code{pcalg} written by Markus Kalisch and Martin Maechler. 
}
\value{
  \code{learn.mec.norm}, \code{learn.mec.multinom}, \code{learn.v},
  \code{learn.complex.norm}, \code{learn.complex.multinom} return the
  pattern of the graph, represented in its adjacency matrix.

  \code{learn.skeleton.norm} and \code{learn.skeleton.multinom} return a
    list of the following two items: 
    \item{amat }{the adjacency matrix of the learned skeleton.}
    \item{sep.pairs }{a list of \code{sep.pair} objects, recording the
  separation pairs we obtained during the skeleton learning procedure.}
}
\references{
  Ma, Z., Xie, X. and Geng, Z. (2008). Structural learning of chain
  graphs via decomposition. \emph{J. Mach. Learn. Res.,}
  \bold{9}, 2847-2880.

  Xie, X., Geng, Z. and Zhao, Q. (2006). Decomposition of structural
  learning about directed acyclic graphs. \emph{Artif.
  Intell.,} \bold{170}, 422-439.
}
\author{ Zongming Ma and Xiangrui Meng}
\examples{
set.seed(100)
p.value <- .01
n <- 3000
is.chaingraph(toy.graph)
tgdata <- rnorm.cg(n, toy.graph, get.normal.dist(toy.graph))
tgug <- naive.getug.norm(tgdata, p.value)
tg.jtree <- ug.to.jtree(tgug)
tg.pat <- learn.mec.norm(tg.jtree, cov(tgdata), n, p.value, "CG")
comp.skel(skeleton(toy.graph), skeleton(tg.pat))
comp.pat(pattern(toy.graph), tg.pat)
}

\keyword{ graphs }
\keyword{ multivariate }
\keyword{ models }
