\name{cvxclust_ama}
\alias{cvxclust_ama}
\title{Convex Clustering via AMA}
\usage{
  cvxclust_ama(X, Lambda, ix, M1, M2, s1, s2, w, gamma, nu,
    eta = 2, type = 2, max_iter = 100, tol = 1e-04,
    accelerate = TRUE, backtracking = TRUE)
}
\arguments{
  \item{X}{The q-by-p data matrix whose columns are to be
  clustered.}

  \item{Lambda}{The q-by-k matrix of Lagrange multipliers}

  \item{ix}{The k-by-2 matrix of index pairs}

  \item{M1}{Index set used to track nonzero weights}

  \item{M2}{Index set used to track nonzero weights}

  \item{s1}{Index set used to track nonzero weights}

  \item{s2}{Index set used to track nonzero weights}

  \item{w}{A vector of k positive weights}

  \item{gamma}{The regularization parameter controlling the
  amount of shrinkage}

  \item{nu}{The initial step size parameter when
  backtracking is applied. Otherwise it is a fixed step
  size in which case there are no guarantees of convergence
  if it exceeds \code{2/ncol(X)}.}

  \item{eta}{The step decrement factor}

  \item{type}{An integer indicating the norm used: 1 =
  1-norm, 2 = 2-norm.}

  \item{max_iter}{The maximum number of iterations.}

  \item{tol}{The convergence tolerance.}

  \item{accelerate}{If \code{TRUE} (the default),
  acceleration is turned on.}

  \item{backtracking}{If \code{TRUE} (the default),
  acceleration is turned on.}
}
\description{
  \code{cvxclust_ama} performs convex clustering via AMA.
  This is an R wrapper function around Fortran code.
  Dimensions of various arguments are as follows: \itemize{
  \item{q is the number of data points} \item{p is the
  number of features} \item{k is the number non-zero
  weights.} }
}

