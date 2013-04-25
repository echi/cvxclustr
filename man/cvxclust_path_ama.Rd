\name{cvxclust_path_ama}
\alias{cvxclust_path_ama}
\title{Convex Clustering via AMA}
\usage{
  cvxclust_path_ama(X, w, gamma, nu0 = 10, tol = 1e-04,
    max_iter = 1e+06, type = 2, accelerate = TRUE)
}
\arguments{
  \item{X}{q-by-p data matrix}

  \item{w}{vector of k positive weights}

  \item{gamma}{sequence of regularization parameters}

  \item{nu0}{positive penalty parameter for quadratic
  deviation term}

  \item{tol}{convergence tolerance}

  \item{max_iter}{maximum number of iterations}

  \item{type}{integer indicating the norm used}

  \item{accelerate}{boolean indicating whether to use
  acceleration}
}
\description{
  \code{cvxclust_path_ama} estimates the convex clustering
  path via AMA using warm starts. q denote the number of
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
sol = cvxclust_path_ama(X,w,gamma)

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
}

