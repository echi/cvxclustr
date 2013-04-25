#' Compute k-nearest neighbor weights
#' 
#' \code{kernel_weights} computes knn weights
#' 
#' @param X q-by-p data matrix
#' @param k number of nearest neighbors
#' @param p number of data points.
#' @export
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
