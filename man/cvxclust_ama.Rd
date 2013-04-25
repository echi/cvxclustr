\name{cvxclust_ama}
\alias{cvxclust_ama}
\title{Convex Clustering via AMA}
\usage{
  cvxclust_ama(X, Lambda, ix, M1, M2, s1, s2, w, gamma, nu,
    eta = 2, type = 2, max_iter = 100, tol = 1e-04,
    accelerate = TRUE)
}
\arguments{
  \item{X}{q-by-p data matrix}

  \item{Lambda}{q-by-k matrix of Lagrange multipliers}

  \item{ix}{k-by-2 matrix of index pairs}

  \item{M1}{}

  \item{M2}{}

  \item{s1}{}

  \item{s2}{}

  \item{w}{vector of k positive weights}

  \item{gamma}{regularization parameter controlling the
  amount of shrinkage}

  \item{eta}{step decrement}

  \item{nu}{positive penalty parameter for quadratic
  deviation term}

  \item{type}{integer indicating the norm used}

  \item{max_iter}{maximum number of iterations}

  \item{tol}{convergence tolerance}

  \item{accelerate}{boolean indicating whether to use
  acceleration}
}
\description{
  \code{cvxclust_ama} performs convex clustering via AMA.
  Let q denote the number of data points and p denote the
  number of covariates. Let k denote the number non-zero
  weights.
}

