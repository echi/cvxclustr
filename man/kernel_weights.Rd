\name{kernel_weights}
\alias{kernel_weights}
\title{Compute Gaussian Kernel Weights}
\usage{
  kernel_weights(X, phi)
}
\arguments{
  \item{X}{The data matrix to be clustered. The rows are
  the features, and the columns are the samples.}

  \item{phi}{The nonnegative parameter that controls the
  scale of kernel weights}
}
\value{
  A vector \cite{w} of weights for convex clustering.
}
\description{
  \code{kernel_weights} computes Gaussian kernel weights
  given a data matrix \code{X} and a scale parameter
  \code{phi}. Namely, the lth weight \code{w[l]} is given
  by \deqn{ w[l] = exp(-phi ||X[,i]-X[,j]||^2) }, where the
  lth pair of nodes is (\code{i},\code{j}).
}
\author{
  Eric C. Chi
}

