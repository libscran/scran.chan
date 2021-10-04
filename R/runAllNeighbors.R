#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' res <- runAllNeighbors(x)
#' str(res)
#' @export
#' @importFrom parallel clusterCall
runAllNeighbors <- function(x, tsne.perplexity=30, umap.num.neighbors=15, snn.num.neighbors=10, snn.resolution=1, num.threads=1) {
    neighbors <- build_nn_index(x)
    tsne.init <- initialize_tsne(neighbors, tsne.perplexity, num.threads)
    umap.init <- initialize_umap(neighbors, umap.num.neighbors, num.threads)
    snn.graph <- build_graph(neighbors, snn.num.neighbors, num.threads)

    jobs <- list(
        snn=function() {
            .cluster_snn_graph_internal(snn.graph, resolution=snn.resolution)
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
