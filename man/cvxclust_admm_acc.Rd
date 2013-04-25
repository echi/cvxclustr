\name{cvxclust_admm_acc}
\alias{cvxclust_admm_acc}
\title{Convex Clustering via Accelerated ADMM}
\usage{
  cvxclust_admm_acc(X, Lambda, V, ix, w, gamma, nu = 1,
    type = 2, max_iter = 100, tol = 1e-04)
}
\arguments{
  \item{X}{q-by-p data matrix}

  \item{Lambda}{q-by-k matrix of Lagrange multipliers}

  \item{V}{q-by-k matrix of centroid differences}

  \item{ix}{k-by-2 matrix of index pairs}

  \item{w}{vector of k positive weights}

  \item{gamma}{regularization parameter controlling the
  amount of shrinkage}

  \item{nu}{positive penalty parameter for quadratic
  deviation term}

  \item{type}{integer indicating the norm used}

  \item{max_iter}{maximum number of iterations}

  \item{tol}{convergence tolerance}
}
\description{
  \code{cvxclust_admm_acc} performs convex clustering via
  accelerated ADMM. Let q denote the number of data points
  and p denote the number of covariates. Let k denote the
  number non-zero weights.
}

