#' Run all downstream steps
#'
#' Run all steps in the basic workflow that are downstream of PCA (but before marker detection).
#' 
#' @param x Numeric matrix containing dimensions in the rows and cells in the columns.
#' This is typically a matrix of principal components.
#' @param do.umap Logical scalar, should we perform a UMAP?
#' @param do.tsne Logical scalar, should we perform a t-SNE?
#' @param do.cluster.snn Logical scalar, should we perform graph-based clustering?
#' @param do.cluster.kmeans Logical scalar, should we perform k-means clustering?
#' @param tsne.perplexity Parameters to be used for t-SNE, see \code{\link{runTSNE.chan}} for details.
#' @param umap.num.neighbors,umap.min.dist Parameters to be used for UMAP, see \code{\link{runUMAP.chan}} for details.
#' @param cluster.kmeans.k Parameters to be used for k-means clustering, see \code{\link{clusterKmeans.chan}} for details.
#' @param cluster.snn.num.neighbors,cluster.snn.method,cluster.snn.resolution Parameters to be used for graph-based clustering, see \code{\link{clusterSNNGraph.chan}} for details.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' 
#' @return A list containing \code{"cluster.snn"}, \code{"umap"} and \code{"tsne"},
#' each of which contains the result of their respective function \code{*.chan} functions.
#'
#' @details
#' By running all of these steps together, we can avoid redundant construction of the nearest neighbor index.
#' We can also execute some of the single-threaded steps concurrently for further time savings.
#'
#' It is tempting to re-use the nearest neighbor search results across the different steps,
#' subsetting to the number of neighbors required in each step.
#' However, we do not do so as the approximate nature of the search means that the subsetting may not yield the same results as a direct search for the requested number of neighbors.
#' This can lead to inconsistent behavior where the results of one step depend on whether other steps were run.
#'
#' @author Aaron Lun
#' 
#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' res <- runAllDownstream(x)
#' str(res)
#' @export
runAllDownstream <- function(x,
    do.tsne=TRUE,
    do.umap=TRUE,
    do.cluster.snn=TRUE,
    do.cluster.kmeans=FALSE,
    tsne.perplexity=30, 
    umap.num.neighbors=15, 
    umap.min.dist=0.01,
    cluster.snn.num.neighbors=10, 
    cluster.snn.method=c("multilevel", "leiden", "walktrap"), 
    cluster.snn.resolution=NULL, 
    cluster.kmeans.k=10,
    num.threads=1) 
{
    neighbors <- build_nn_index(x)

    snn.graph <- umap.init <- tsne.init <- kmeans.info <- NULL
    if (do.tsne) {
        tsne.init <- initialize_tsne(neighbors, tsne.perplexity, interpolate=-1, max_depth=7, nthreads=num.threads) # defaults from runTSNE.chan.
    }
    if (do.umap) {
        umap.init <- initialize_umap(neighbors, umap.num.neighbors, umap.min.dist, num.threads)
    }
    if (do.cluster.snn) {
        cluster.snn.method <- match.arg(cluster.snn.method)
        cluster.snn.resolution <- .default_resolution(cluster.snn.method, cluster.snn.resolution)
        snn.graph <- build_graph(neighbors, 
            k=cluster.snn.num.neighbors, 
            method=cluster.snn.method, 
            resolution=cluster.snn.resolution, 
            nthreads=num.threads)
    }
    if (do.cluster.kmeans) {
        kmeans.info <- list(x, cluster.kmeans.k)
    }

    output <- run_all_downstream(snn.graph, umap.init, tsne.init, kmeans.info, num.threads)
    names(output) <- c("cluster.snn", "umap", "tsne", "cluster.kmeans")
    output <- output[!vapply(output, is.null, TRUE)]

    if (do.cluster.snn) {
        output$cluster.snn <- .clean_graph_clustering(cluster.snn.method, output$cluster.snn)
    }
    if (do.umap) {
        output$umap <- t(output$umap)
    }
    if (do.tsne) {
        output$tsne <- t(output$tsne)
    }
    if (do.cluster.kmeans) {
        output$cluster.kmeans <- .clean_kmeans_clustering(output$cluster.kmeans)
    }

    output
}
