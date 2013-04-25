\name{cvxclust_path_admm}
\alias{cvxclust_path_admm}
\title{Convex Clustering via ADMM}
\usage{
  cvxclust_path_admm(X, w, gamma, nu = 1, tol = 1e-04,
    max_iter = 10000, type = 2, accelerate = TRUE)
}
\arguments{
  \item{X}{q-by-p data matrix}

  \item{w}{vector of k positive weights. The ith entry w[i]
  denotes the weight used between the ith pair ofcentroids.
  The weights are in dictionary order.}

  \item{gamma}{sequence of regularization parameters}

  \item{nu}{positive penalty parameter for quadratic
  deviation term}

  \item{tol}{convergence tolerance}

  \item{max_iter}{maximum number of iterations}

  \item{type}{integer indicating the norm used: 1 = 1-norm,
  2 = 2-norm}

  \item{accelerate}{boolean indicating whether to use
  acceleration}
}
\description{
  \code{cvxclust_path_admm} estimates the convex clustering
  path via ADMM using warm starts. q denote the number of
  data points and p denote the number of covariates. Let k
  denote the number non-zero weights.
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
sol = cvxclust_path_admm(X,w,gamma)

## Plot the cluster path
library(ggplot2)
df.paths = data.frame(x=c(),y=c(), group=c())
for (j in 1:nGamma) {
  x = sol$UHx[[j]][1,]
  y = sol$UHx[[j]][2,]
  df = data.frame(x=x, y=y, group=1:p)
  df.paths = rbind(df.paths,df)
}
data_plot = ggplot(data=df.paths,aes(x=x,y=y))
data_plot = data_plot + geom_path(aes(group=group),colour='grey60',alpha=0.5)
data_plot = data_plot + geom_point(data=data.frame(x=X[1,],y=X[2,]),aes(x=x,y=y))
data_plot + theme_bw()

## Iris Data
X = as.matrix(iris[,1:4])
X = t(scale(X,center=TRUE,scale=FALSE))
dupix = which(duplicated(t(X)))
X = X[,-dupix]
p = ncol(X)

## Pick some weights and a sequence of regularization parameters.
kgamma = 4
knbs = 5
w = kernel_weights(X,kgamma)
w = knn_weights(w,knbs,p)
gamma = c()
gamma[1] = 0
gamma[2] = 1.01
i = 3
repeat {
  g = 1.05*gamma[i-1]
  if (g > 21) break
  gamma[i] = g
  i = i + 1
}
gamma = gamma[1:(i-1)]

## Perform clustering
sol = cvxclust_path_admm(X,w,gamma,nu=1,tol=1e-1)

## Plot the cluster path
library(ggplot2)
svdX = svd(X)
pc = svdX$u[,1:2,drop=FALSE]
pc.df = as.data.frame(t(pc)\%*\%X)
nGamma = sol$nGamma
Uout = sol$UHx
df.paths = data.frame(x=c(),y=c(), group=c())
for (j in 1:nGamma) {
  pcs = t(pc)\%*\%Uout[[j]]
  x = pcs[1,]
  y = pcs[2,]
  df = data.frame(x=pcs[1,], y=pcs[2,], group=1:p)
  df.paths = rbind(df.paths,df)
}
X_data = as.data.frame(t(X)\%*\%pc)
colnames(X_data) = c("x","y")
X_data$Species = iris[-dupix,5]
data_plot = ggplot(data=df.paths,aes(x=x,y=y))
data_plot = data_plot + geom_path(aes(group=group),colour='grey60',alpha=0.5)
data_plot = data_plot + geom_point(data=X_data,aes(x=x,y=y,colour=Species,shape=Species))
data_plot = data_plot + xlab('Principal Component 1') + ylab('Principal Component 2')
data_plot + theme_bw()
}

