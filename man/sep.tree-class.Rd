\name{sep.tree-class}
\docType{class}
\alias{sep.tree-class}

\title{Class "sep.tree" }
\description{
  Objects representing separation tree as described in Xie, Geng and Zhao
  (2006) and Ma, Xie and
  Geng (2008), which includes junction tree of cliques as a special case.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("sep.tree",
  ...)}. We recommend user to use the functions provided in
this package, such as \code{\link{ug.to.jtree}} to create the separation
tree object.
}
\section{Slots}{
  \describe{
  \item{\code{tree.struct} }{the adjacency matrix of the junction tree.}
  \item{\code{cliques} }{the list of cliques on the junction tree.}
  \item{\code{separators} }{the list of separators on the junction tree, each
    element on the list has two components: the \code{separator}
    compoenent is the
    set of graph vertices in the separator and \code{edge} is the edge
    on the tree that the separator is attahced to.}
  }
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
showClass("sep.tree")
}
\keyword{classes}
