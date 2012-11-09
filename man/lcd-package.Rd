\name{lcd-package}
\alias{lcd-package}
\alias{lcd}
\docType{package}
\title{ Structural learning of chain graphs via the decomposition
  approach }

\description{ This package implements the algorithms for
learning chain graphs (and as a special case, Bayesian networks) via the
decomposition approach as described in Ma, Xie and Geng (2008)and Xie,
Geng and Zhao (2006). The correctness of the algorithms is given in the
above cited paper.

The algorithms are \emph{constraint-based} method: the actual learning
procedure depends on testing statistical significances. See the learning functions for
the particular tests we used in the implementation.

Currently, the package supports learning with continuous and discrete
data, but not the mix of them.

The package also provides some utility functions for graph manipulation, frequency table
construction and compression and random distribution
and random sample generation from chain graphs.
}
\section{List of functions}{
  The package has the following categories of functions.
  \itemize{
    \item Chain graph/Bayesian network learning:
    \code{\link{learn.mec.norm}},
    \code{\link{learn.mec.multinom}},
    \code{\link{learn.skeleton.norm}},
    \code{\link{learn.skeleton.multinom}},
    \code{\link{learn.complex.norm}},
    \code{\link{learn.complex.multinom}},
    \code{\link{learn.v}}.

    \item Graphical model manipulation and graph plotting:
    \code{\link{draw}},
    \code{\link{is.chaingraph}},
    \code{\link{is.separated}},
    \code{\link{moralize}},
    \code{\link{pattern}},
    \code{\link{skeleton}},
    \code{\link{tri.ug}},
    \code{\link{maxcard.search}},
    \code{\link{ug.to.jtree}}.
    
    \item Random distribution and random sample generation from chain graph:
    \code{\link{get.normal.dist}},
    \code{\link{get.multinom.dist}},
    \code{\link{rnorm.cg}},
    \code{\link{rmultinom.cg}}.

    \item Conditional independence tests:
    \code{\link{norm.ci.test}},
    \code{\link{multinom.ci.test}}.

    \item Construction of undirected independence graphs from
    data:
    \code{\link{naive.getug.norm}},
    \code{\link{naive.getug.multinom}}.

    \item Frequency table manipulation:
    \code{\link{as.freq.tb}},
    \code{\link{compress.freq.tb}}.

    \item Graph Comparison:
    \code{\link{comp.pat}},
    \code{\link{comp.skel}}.
  } 
  Moreover, the package include the following two graph structure:
  \itemize{
    \item The ALARM network: \code{\link{alarm.net}};
    \item A toy example: \code{\link{toy.graph}}.
    }
%   The package is intended to be a contribution to the gR-project
%   described in Lauritzen (2002).
}
\author{
  Zongming Ma\cr
  Department of Statistics\cr
  The Wharton School, University of Pennsylvania, USA\cr
  \email{zongming.ma@gmail.com}
  
  
  Xiangrui Meng\cr
  Institute of Computational and Mathematical Engineering\cr
  Stanford University, USA
}
\references{
%   Lauritzen, S. L. (2002). gRaphical Models in R.
%   \emph{R News}, 3(2)39.
  Ma, Z., Xie, X. and Geng, Z. (2008). Structural learning of chain
  graphs via decomposition. \emph{J. Mach. Learn. Res.,}
  \bold{9}, 2847-2880.

  Xie, X., Geng, Z. and Zhao, Q. (2006). Decomposition of structural
  learning about directed acyclic graphs. \emph{Artif.
  Intell.,} \bold{170}, 422-439.
}
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

\keyword{ package }
\keyword{ multivariate }
\keyword{ models }
\keyword{ graphs }

