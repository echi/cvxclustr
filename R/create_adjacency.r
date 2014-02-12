#' Create adjacency matrix from V
#' 
#' \code{create_adjacency} creates an n-by-n sparse adjacency matrix from the matrix of centroid differences.
#' 
#' @param V Matrix of centroid differences
#' @param ix Key mapping centroid difference indices to node pairs
#' @param n Number of points to cluster
#' @export
#' @examples
#' ## Clusterpaths for Mammal Dentition
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
#' gamma <- seq(0.0,43, length.out=100)
#' 
#' ## Perform clustering
#' sol <- cvxclust(X, w, gamma, method = 'ama', accelerate = TRUE)
#' 
#' ## Construct adjacency matrix
#' n <- ncol(X)
#' nK <- n*(n-1)/2
#' ix <- vec2tri(which(w > 0),p)
#' A <- create_adjacency(sol$V[[10]],ix,n)
#' G <- graph.adjacency(A, mode = 'upper')
#' plot(G,vertex.label=as.character(mammals[,1]),vertex.label.cex=0.65,vertex.label.font=2,mark.groups=clist)
create_adjacency <- function(V,ix,n) {
  differences <- apply(V,2,FUN=function(x) {norm(as.matrix(x),'f')})
  connected_ix <- which(differences == 0)
  i <- ix[connected_ix,1]
  j <- ix[connected_ix,2]
  A <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
  A[(j-1)*n + i] <- 1
  return(A)
}