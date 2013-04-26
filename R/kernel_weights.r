#' Compute Gaussian Kernel Weights
#' 
#' \code{kernel_weights} computes Gaussian kernel weights given a data matrix \code{X} and a scale parameter \code{phi}. Namely,
#' the lth weight \code{w[l]} is given by
#' \deqn{
#' w[l] = exp(-phi ||X[,i]-X[,j]||^2)
#' }, where the lth pair of nodes is (\code{i},\code{j}).
#' @param X The data matrix to be clustered. The rows are the features, and the columns are the samples.
#' @param phi The nonnegative parameter that controls the scale of kernel weights
#' @author Eric C. Chi
#' @export
#' @return A vector \cite{w} of weights for convex clustering.
kernel_weights = function(X,phi) {
  p = ncol(X)
  k = 1
  w = matrix(0,p*(p-1)/2,1)
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      w[k] = exp(-phi*norm(as.matrix(X[,i,drop=FALSE]-X[,j,drop=FALSE]),'f')^2)
      k = k+1
    }
  }
  return(weights=w)
}
