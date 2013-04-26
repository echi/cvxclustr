#' "Thin" a weight vector to be positive only for its k-nearest neighbors
#' 
#' \code{knn_weights} takes a weight vector \code{w} and sets the ith 
#' component \code{w[i]} to zero if either of the two corresponding nodes
#' is not among the other's \code{k} nearest neighbors.
#' 
#' @param w A vector of nonnegative weights. The ith entry \code{w[i]} denotes the weight used between the ith pair of centroids. The weights are in dictionary order.
#' @param k The number of nearest neighbors
#' @param p The number of data points.
#' @author Eric C. Chi
#' @export
#' @return A vector \cite{w} of weights for convex clustering.
knn_weights = function(w,k,p) {
  i = 1
  neighbors = tri2vec(i,(i+1):p,p)
  keep = neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
  for (i in 2:(p-1)) {
    group_A = tri2vec(i,(i+1):p,p)
    group_B = tri2vec(1:(i-1),i,p)
    neighbors = c(group_A,group_B)
    knn = neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
    keep = union(knn,keep)
  }
  i = p
  neighbors = tri2vec(1:(i-1),i,p)
  knn = neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
  keep = union(knn,keep)
  w[-keep] = 0
  return(w)
}
