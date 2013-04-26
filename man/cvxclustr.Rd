\docType{package}
\name{cvxclustr}
\alias{cvxclustr}
\alias{cvxclustr-package}
\title{Convex Clustering via Splitting Methods}
\description{
  This package provides two variable splitting methods
  \itemize{ \item{Alternating Method of Multipliers (ADMM)}
  \item{Alternating Minimization Algorithm (AMA)} } for
  solving a convex formulation of the clustering problem.
  We seek the centroids u_i that minimize \deqn{
  \frac{1}{2} \sum_i || x_i - u_i||_2^2 + \gamma \sum_l
  w_{l} ||u_{l1} - u_{l2} || } Two penalty norms are
  currently supported: 1-norm and 2-norm.
}
\references{
  Eric C. Chi, Kenneth Lange (2013). Splitting Methods for
  convex clustering. arXiv: 1304.0499 [stat.ML].
  \url{http://arxiv.org/abs/1304.0499}.
}

