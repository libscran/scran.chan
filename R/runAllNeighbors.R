#' Run all neighbor-related steps
#'
#' Run all steps in the basic workflow that require nearest neighbor detection.
#' 
#' @param x Numeric matrix containing dimensions in the rows and cells in the columns.
#' This is typically a matrix of principal components.
#' @param tsne.perplexity Parameters to be used for t-SNE, see \code{\link{runTSNE.chan}} for details.
#' @param umap.num.neighbors Parameters to be used for UMAP, see \code{\link{runUMAP.chan}} for details.
#' @param cluster.snn.num.neighbors,cluster.snn.resolution Parameters to be used for graph-based clustering, see \code{\link{clusterSNNGraph.chan}} for details.
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
#' res <- runAllNeighbors(x)
#' str(res)
#' @export
runAllNeighbors <- function(x,
    tsne.perplexity=30, 
    umap.num.neighbors=15, 
    cluster.snn.num.neighbors=10, 
    cluster.snn.resolution=1, 
    num.threads=1) 
{
    neighbors <- build_nn_index(x)
    tsne.init <- initialize_tsne(neighbors, tsne.perplexity, interpolate=-1, max_depth=7, nthreads=num.threads) # defaults from runTSNE.chan.
    umap.init <- initialize_umap(neighbors, umap.num.neighbors, num.threads)
    snn.graph <- build_graph(neighbors, cluster.snn.num.neighbors, cluster.snn.resolution, num.threads)

    output <- run_all_neighbors(snn.graph, umap.init, tsne.init, num.threads)
    names(output) <- c("cluster.snn", "umap", "tsne")

    # Enforce 1-based indexing, see ?clusterSNNGraph.chan.
    output$cluster.snn$best <- output$cluster.snn$best + 1L
    output$cluster.snn$membership <- lapply(output$cluster.snn$membership, function(x) x + 1L)

    output$umap <- t(output$umap)
    output$tsne <- t(output$tsne)

    output
}
