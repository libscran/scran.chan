#' Apply graph-based clustering
#'
#' Perform multi-level (i.e., Louvain) clustering on a shared nearest neighbor graph.
#'
#' @param x Numeric matrix where rows are dimensions and columns are cells.
#' @param num.neighbors Integer scalar specifying the number of neighbors to use to construct the graph.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' 
#' @return A list containing:
#' \itemize{
#' \item \code{membership}, a list of integer vectors with cluster assignments for each cell at each level.
#' Assignments are sorted by decreasing resolution (i.e., fewer, larger clusters).
#' \item \code{modularity}, a numeric vector containing the modularity of each level.
#' \item \code{best}, the level with the lowest modularity.
#' }
#' 
#' @author Aaron Lun
#'
#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' clustering <- clusterSNNGraph.chan(x)
#' clustering$modularity
#' lapply(clustering$membership, table)
#' 
#' @export
clusterSNNGraph.chan <- function(x, num.neighbors=10, resolution=1, num.threads=1) {
    neighbors <- build_nn_index(x)
    graph <- build_graph(neighbors, k=num.neighbors, nthreads=num.threads)
    .cluster_snn_graph_internal(graph, resolution=resolution)
}

.cluster_snn_graph_internal <- function(graph, resolution) {
    clustering <- cluster_multilevel(graph, res=resolution)
    clustering$best <- clustering$best + 1L
    clustering$membership <- lapply(clustering$membership, function(x) x + 1L)
    clustering
}
