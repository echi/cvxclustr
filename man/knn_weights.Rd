\name{knn_weights}
\alias{knn_weights}
\title{"Thin" a weight vector to be positive only for its k-nearest neighbors}
\usage{
  knn_weights(w, k, p)
}
\arguments{
  \item{w}{A vector of nonnegative weights. The ith entry
  \code{w[i]} denotes the weight used between the ith pair
  of centroids. The weights are in dictionary order.}

  \item{k}{The number of nearest neighbors}

  \item{p}{The number of data points.}
}
\value{
  A vector \cite{w} of weights for convex clustering.
}
\description{
  \code{knn_weights} takes a weight vector \code{w} and
  sets the ith component \code{w[i]} to zero if either of
  the two corresponding nodes is not among the other's
  \code{k} nearest neighbors.
}
\author{
  Eric C. Chi
}

