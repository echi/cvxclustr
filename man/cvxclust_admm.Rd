\name{cvxclust_admm}
\alias{cvxclust_admm}
\title{Convex Clustering via ADMM}
\usage{
  cvxclust_admm(X, Lambda, V, ix, w, gamma, nu = 1,
    type = 2, max_iter = 100, tol = 1e-04,
    accelerate = TRUE)
}
\arguments{
  \item{X}{The q-by-p data matrix whose columns are to be
  clustered.}

  \item{Lambda}{The q-by-k matrix of Lagrange multipliers.}

  \item{V}{The q-by-k matrix of centroid differences.}

  \item{ix}{The k-by-2 matrix of index pairs.}

  \item{w}{The vector of k positive weights.}

  \item{gamma}{The regularization parameter controlling the
  amount of shrinkage.}

  \item{nu}{A positive penalty parameter for quadratic
  deviation term.}

  \item{tol}{The convergence tolerance.}

  \item{max_iter}{The maximum number of iterations.}

  \item{type}{An integer indicating the norm used: 1 =
  1-norm, 2 = 2-norm.}

  \item{accelerate}{If \code{TRUE} (the default),
  acceleration is turned on.}
}
\description{
  \code{cvxclust_admm} performs convex clustering via ADMM.
  This is an R wrapper function around Fortran code.
  Dimensions of various arguments are as follows: \itemize{
  \item{q is the number of data points} \item{p is the
  number of features} \item{k is the number non-zero
  weights.} }
}
\author{
  Eric C. Chi
}

