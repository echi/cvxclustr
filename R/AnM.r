#' Compute Lipschitz bound upper bound on the largest eigenvalue of the Laplacian
#' 
#' \code{AnM} computes the node degrees from the edge information.
#' 
#' @param w vector of weights
#' @param p number of points to cluster
#' @export
#' @examples
#' data(mammals)
#' X <- as.matrix(mammals[,-1])
#' X <- t(scale(X,center=TRUE,scale=FALSE))
#' p <- ncol(X)
#' 
#' ## Pick some weights and a sequence of regularization parameters.
#' k <- 5
#' phi <- 0.5
#' w <- kernel_weights(X,phi)
#' w <- knn_weights(w,k,p)
#' AnM(w,p)
AnM <- cmpfun(function(w,p) {
  ix <- compactify_edges(w,p)$ix
  nEdges <- nrow(ix)
  nVertex <- max(ix)
  deg <- integer(nVertex)
  for (i in 1:nVertex) {
    deg[i] <- length(which(ix[,1] == i)) + length(which(ix[,2] == i))
  }
  bound <- -Inf
  for (j in 1:nEdges) {
    bound <- max(bound,deg[ix[j,1]] + deg[ix[j,2]])
  }
  return(min(bound,nVertex))
})