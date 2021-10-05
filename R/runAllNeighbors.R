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
#' @importFrom parallel clusterCall
runAllNeighbors <- function(x, tsne.perplexity=30, umap.num.neighbors=15, cluster.snn.num.neighbors=10, cluster.snn.resolution=1, num.threads=1) {
    neighbors <- build_nn_index(x)
    tsne.init <- initialize_tsne(neighbors, tsne.perplexity, num.threads)
    umap.init <- initialize_umap(neighbors, umap.num.neighbors, num.threads)
    snn.graph <- build_graph(neighbors, cluster.snn.num.neighbors, num.threads)

    jobs <- list(
        cluster.snn=function() {
            .cluster_snn_graph_internal(snn.graph, resolution=cluster.snn.resolution)
        },
        umap=function() {
            output <- run_umap(umap.init)
            t(output)
        },
        tsne=function() {
            remaining.threads <- max(1L, num.threads - 2L)
            output <- run_tsne(tsne.init, remaining.threads)
            t(output)
        }
    )

    if (num.threads == 1) {
        output <- lapply(jobs, function(f) f())
    } else {
        output <- BiocParallel::bplapply(jobs, function(f) f(), BPPARAM=BiocParallel::MulticoreParam(num.threads))
    }

    output
}
