#' Find clusters
#' 
#' \code{find_clusters} uses breadth-first search to identify the connected components of the corresponding
#' adjacency graph of the centroid differences vectors.
#' 
#' @param A adjacency matrix
#' @export
#' 
find_clusters <- function(A) {
  G <- graph.adjacency(A, mode = 'upper')
  n <- nrow(A)
  node_seen <- logical(n)
  cluster <- integer(n)
  k <- 1
  for (i in 1:n) {
    if (!node_seen[i]) {
      connected_set <- graph.bfs(G, root=i, unreachable = FALSE)$order
      node_seen[connected_set] <- TRUE
      cluster[connected_set] <- k
      k <- k + 1
    }
  }
  nClusters <- k - 1
  size <- integer(nClusters)
  for (j in 1:nClusters) {
    size[j] <- length(which(cluster == j))
  }
  return(list(cluster=cluster, size=size))
}