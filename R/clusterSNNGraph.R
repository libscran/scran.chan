#' Apply graph-based clustering
#'
#' Perform multi-level (i.e., Louvain) clustering on a shared nearest neighbor graph.
#'
#' @param x Numeric matrix where rows are dimensions and columns are cells.
#' @param num.neighbors Integer scalar specifying the number of neighbors to use to construct the graph.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param method String specifying the community detection method to use.
#' @param resolution Numeric scalar specifying the resolution to use for multi-level clustering.
#' Defaults to 1 for \code{method="multilevel"} and 0.05 for \code{method="leiden"}.
#' 
#' @return A list containing \code{membership}, an integer vector with cluster assignments for each cell.
#' Each method may also return additional elements.
#' For \code{method="multilevel"}, we have:
#' \itemize{
#' \item \code{levels}, a list of integer vectors with cluster assignments for each cell at each level.
#' Assignments are sorted by decreasing resolution (i.e., fewer, larger clusters).
#' \item \code{modularity}, a numeric vector containing the modularity of each level.
#' \item \code{best}, the level with the lowest modularity.
#' }
#' For \code{method="leiden"}, we have:
#' \itemize{
#' \item \code{quality}, a numeric scalar containing the quality of the clustering (either the modularity or a related score).
#' }
#' For \code{method="walktrap"}, we have:
#' \itemize{
#' \item \code{merges}, an integer matrix specifying how the clusters were merged to obtain \code{membership}.
#' Each row corresponds to a merge step and contains the IDs of the temporary clusters (not the same as those in \code{membership}). 
#' \item \code{modularity}, a numeric vector containing the modularity before and after each merge step. 
#' }
#' 
#' @author Aaron Lun
#'
#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' clustering <- clusterSNNGraph.chan(x)
#' clustering$modularity
#' table(clustering$membership)
#' 
#' @export
clusterSNNGraph.chan <- function(x,
    num.neighbors=10,
    method=c("multilevel", "leiden", "walktrap"), 
    resolution=NULL, 
    num.threads=1) 
{
    neighbors <- build_nn_index(x)

    method <- match.arg(method)
    resolution <- .default_resolution(method, resolution)
    graph <- build_graph(neighbors, k=num.neighbors, method=method, resolution=resolution, nthreads=num.threads)

    clustering <- cluster_graph(graph)
    .clean_graph_clustering(method, clustering)
}

.default_resolution <- function(method, resolution) {
    if (is.null(resolution)) {
        if (method == "multilevel") {
            resolution <- 1
        } else if (method=="leiden") {
            resolution <- 0.05
        } else {
            resolution <- 1 # dummy value
        }
    } 
    resolution
}

.clean_graph_clustering <- function(method, clustering) { 
    clustering$membership <- factor(clustering$membership + 1L)

    if (method=="multilevel") {
        clustering$best <- clustering$best + 1L
        clustering$levels <- lapply(clustering$levels, function(x) factor(x + 1L))
    } else if (method=="walktrap") {
        clustering$merges <- clustering$merges + 1L
    }

    clustering
}
