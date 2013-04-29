\name{cvxclust_path_ama}
\alias{cvxclust_path_ama}
\title{Convex Clustering Path via AMA}
\usage{
  cvxclust_path_ama(X, w, gamma, nu = 1, tol = 0.001,
    max_iter = 10000, type = 2, accelerate = TRUE)
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
}
\value{
  \code{U} A list of centroid matrices.

  \code{V} A list of centroid difference matrices.

  \code{Lambda} A list of Lagrange multiplier matrices.
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
  This speed-up is employed by default.
}
\examples{
## Clusterpaths for Mammal Dentition
data(mammals)
X = as.matrix(mammals[,-1])
X = t(scale(X,center=TRUE,scale=FALSE))
p = ncol(X)

## Pick some weights and a sequence of regularization parameters.
k = 5
phi = 0.5
w = kernel_weights(X,phi)
w = knn_weights(w,k,p)
gamma = seq(0.0,43, length.out=100)

## Perform clustering
sol = cvxclust_path_ama(X,w,gamma)

## Plot the cluster path
library(ggplot2)
svdX = svd(X)
pc = svdX$u[,1:2,drop=FALSE]
pc.df = as.data.frame(t(pc)\%*\%X)
nGamma = sol$nGamma
df.paths = data.frame(x=c(),y=c(), group=c())
for (j in 1:nGamma) {
  pcs = t(pc)\%*\%sol$U[[j]]
  x = pcs[1,]
  y = pcs[2,]
  df = data.frame(x=pcs[1,], y=pcs[2,], group=1:p)
  df.paths = rbind(df.paths,df)
}
X_data = as.data.frame(t(X)\%*\%pc)
colnames(X_data) = c("x","y")
X_data$Name = mammals[,1]
data_plot = ggplot(data=df.paths,aes(x=x,y=y))
data_plot = data_plot + geom_path(aes(group=group),colour='grey30',alpha=0.5)
data_plot = data_plot + geom_text(data=X_data,aes(x=x,y=y,label=Name), position=position_jitter(h=0.125,w=0.125))
data_plot = data_plot + geom_point(data=X_data,aes(x=x,y=y),size=1.5)
data_plot = data_plot + xlab('Principal Component 1') + ylab('Principal Component 2')
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

