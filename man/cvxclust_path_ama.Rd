\name{cvxclust_path_ama}
\alias{cvxclust_path_ama}
\title{Convex Clustering Path via AMA}
\usage{
  cvxclust_path_ama(X, w, gamma, nu = 10, tol = 1e-04,
    max_iter = 10000, type = 2, accelerate = TRUE,
    backtracking = TRUE)
}
\arguments{
  \item{X}{The data matrix to be clustered. The rows are
  the features, and the columns are the samples.}

  \item{w}{A vector of nonnegative weights. The ith entry
  \code{w[i]} denotes the weight used between the ith pair
  of centroids. The weights are in dictionary order.}

  \item{gamma}{A sequence of regularization parameters.}

  \item{nu}{The initial step size parameter when
  backtracking is applied. Otherwise it is a fixed step
  size in which case there are no guarantees of convergence
  if it exceeds \code{2/ncol(X)}.}

  \item{tol}{The convergence tolerance.}

  \item{max_iter}{The maximum number of iterations.}

  \item{type}{An integer indicating the norm used: 1 =
  1-norm, 2 = 2-norm.}

  \item{accelerate}{If \code{TRUE} (the default),
  acceleration is turned on.}

  \item{backtracking}{If \code{TRUE} (the default),
  backtracking is used.}
}
\description{
  \code{cvxclust_path_ama} estimates the convex clustering
  path via the Alternating Minimization Algorithm. Required
  inputs include a data matrix \code{X} (rows are features;
  columns are samples), a vector of weights \code{w}, and a
  sequence of regularization parameters \code{gamma}. Two
  penalty norms are currently supported: 1-norm and 2-norm.
  AMA is performing proximal gradient ascent on the dual
  function, and therefore can be accelerated with FISTA.
  Additional speed up can be had by employing backtracking.
  Both speed-ups are employed by default.
}
\examples{
## Create a small set of points to cluster.
set.seed(12345)
p = 10
q = 2
X = matrix(rnorm(p*q),q,p)

## Pick some weights and a sequence of regularization parameters.
w = kernel_weights(X,2)
gamma = c()
gamma[1] = 0
gamma[2] = 1.01
i = 3
repeat {
  g = 1.05*gamma[i-1]
  if (g > 50) break
  gamma[i] = g
  i = i + 1
}
nGamma = i-1
gamma = gamma[1:nGamma]

## Perform clustering
sol = cvxclust_path_ama(X,w,gamma)

## Plot the cluster path
library(ggplot2)
df.paths = data.frame(x=c(),y=c(), group=c())
for (j in 1:nGamma) {
  x = sol$U[[j]][1,]
  y = sol$U[[j]][2,]
  df = data.frame(x=x, y=y, group=1:p)
  df.paths = rbind(df.paths,df)
}
data_plot = ggplot(data=df.paths,aes(x=x,y=y))
data_plot = data_plot + geom_path(aes(group=group),colour='grey60',alpha=0.5)
data_plot = data_plot + geom_point(data=data.frame(x=X[1,],y=X[2,]),aes(x=x,y=y))
data_plot + theme_bw()
}
\author{
  Eric C. Chi
}
\seealso{
  \code{\link{cvxclust_path_admm}} for estimating the
  clustering path with ADMM. \code{\link{kernel_weights}}
  and \code{\link{knn_weights}} compute useful weights.
}

