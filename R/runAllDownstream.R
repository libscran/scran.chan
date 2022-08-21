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
#' @param tsne.args Named list of t-SNE parameters, see \code{\link{runTSNE.chan}} for details.
#' @param umap.args Named list of UMAP parameters, see \code{\link{runUMAP.chan}} for details.
#' @param cluster.kmeans.args Named list of k-means parameters, see \code{\link{clusterKmeans.chan}} for details.
#' @param cluster.snn.args Named list of SNN parameters, see \code{\link{clusterSNNGraph.chan}} for details.
#' @param tsne.perplexity Deprecated, use \code{tsne.args} instead.
#' @param umap.num.neighbors,umap.min.dist Deprecated, use \code{umap.args} instead.
#' @param cluster.kmeans.k,cluster.kmeans.init Deprecated, use \code{cluster.kmeans.args} instead. 
#' @param cluster.snn.num.neighbors,cluster.snn.method,cluster.snn.resolution Deprecated, use \code{cluster.snn.args} instead.
#' @param drop Logical scalar indicating whether to drop the sweep-based formatting when the parameters are scalars.
#' @param downsample Integer scalar specifying the downsampling level (see \code{k} in \code{\link{downsampleByNeighbors.chan}}.
#' If \code{NULL}, no downsampling is performed.
#' @param approximate Logical scalar specifying whether to perform an approximate neighbor search.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' 
#' @return A list containing \code{"cluster.snn"}, \code{"cluster.kmeans"}, \code{"umap"} and \code{"tsne"},
#' each of which contains the result of their respective function \code{*.chan} functions.
#'
#' @details
#' By running all of these steps together, we can avoid redundant construction of the nearest neighbor index.
#' We can also execute some of the single-threaded steps concurrently for further time savings.
#'
#' @author Aaron Lun
#' 
#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' res <- runAllDownstream(x, num.threads=3)
#' str(res)
#' @export
#' @importFrom parallel stopCluster
runAllDownstream <- function(x,
    do.tsne=TRUE,
    do.umap=TRUE,
    do.cluster.snn=TRUE,
    do.cluster.kmeans=FALSE,
    tsne.args=list(),
    umap.args=list(),
    cluster.snn.args=list(),
    cluster.kmeans.args=list(),
    tsne.perplexity=NULL,
    umap.num.neighbors=NULL, 
    umap.min.dist=NULL,
    cluster.snn.num.neighbors=NULL, 
    cluster.snn.method=NULL,
    cluster.snn.resolution=NULL, 
    cluster.kmeans.k=NULL,
    cluster.kmeans.init=NULL,
    drop=TRUE,
    downsample=NULL,
    approximate=TRUE, 
    num.threads=1) 
{
    # Fixing the parameters.
    legacy_args <- function(name, legacy, formals, args) {
        if (!name %in% names(args)) {
            if (!is.null(legacy)) {
                args[[name]] <- legacy
            } else {
                args[[name]] <- formals[[name]]
            }
        }
        args
    }

    default_args <- function(names, formals, args) {
        for (n in names) {
            if (!n %in% names(args)) {
                args[n] <- formals[n]
            }
        }
        args
    }

    if (do.tsne) {
        tsne.formals <- formals(runTSNE.chan) 
        tsne.args <- legacy_args("perplexity", tsne.perplexity, tsne.formals, tsne.args)
        tsne.args <- default_args(c("interpolate", "max.depth", "seed"), tsne.formals, tsne.args)
    }

    if (do.umap) {
        umap.formals <- formals(runUMAP.chan) 
        umap.args <- legacy_args("num.neighbors", umap.num.neighbors, umap.formals, umap.args)
        umap.args <- legacy_args("min.dist", umap.min.dist, umap.formals, umap.args)
        umap.args <- default_args(c("num.epochs", "seed"), umap.formals, umap.args)
    }

    if (do.cluster.snn) {
        cluster.snn.formals <- formals(clusterSNNGraph.chan) 
        cluster.snn.args <- legacy_args("num.neighbors", cluster.snn.num.neighbors, cluster.snn.formals, cluster.snn.args)
        cluster.snn.args <- legacy_args("method", cluster.snn.method, cluster.snn.formals, cluster.snn.args)
        cluster.snn.args <- legacy_args("resolution", cluster.snn.resolution, cluster.snn.formals, cluster.snn.args)
        cluster.snn.args <- default_args(c("weight.scheme", "steps", "seed"), cluster.snn.formals, cluster.snn.args)
    }

    if (do.cluster.kmeans) {
        cluster.kmeans.formals <- formals(clusterKmeans.chan) 
        cluster.kmeans.args <- legacy_args("k", cluster.kmeans.k, cluster.kmeans.formals, cluster.kmeans.args)
        cluster.kmeans.args <- legacy_args("init.method", cluster.kmeans.init, cluster.kmeans.formals, cluster.kmeans.args)
        cluster.kmeans.args <- default_args("seed", cluster.kmeans.formals, cluster.kmeans.args)
    }

    down <- do_downsample(downsample)
    if (down) {
        chosen <- downsampleByNeighbors.chan(x, downsample, approximate=approximate, num.threads=num.threads)
        x <- x[,chosen$chosen,drop=FALSE]
        m <- match(chosen$assigned, chosen$chosen)
    }

    # Generating the neighbors.
    #
    # It is tempting to re-use the nearest neighbor search results across the
    # different steps, subsetting to the number of neighbors required in each step.
    # However, we do not do so as the approximate nature of the search means that
    # the subsetting may not yield the same results as a direct search for the
    # requested number of neighbors.  This can lead to inconsistent behavior where
    # the results of one step depend on whether other steps were run.
    all.neighbors <- list() 
    if (do.tsne || do.umap || do.cluster.snn) {
        nnbuilt <- build_nn_index(x, approximate=approximate)
        if (do.tsne) {
            all.neighbors <- .find_tsne_neighbors(nnbuilt, perplexity=tsne.args$perplexity, num.threads=num.threads, existing=all.neighbors)
        }
        if (do.umap) {
            all.neighbors <- .find_umap_neighbors(nnbuilt, num.neighbors=umap.args$num.neighbors, num.threads=num.threads, existing=all.neighbors)
        }
        if (do.cluster.snn) {
            # This MUST be last, as distances are not reported for efficiency.
            # We don't want to fill 'all.neighbors' with no-distance results at
            # 'k' that might be used by t-SNE or UMAP. By doing SNN last, we
            # ensure that no-distance results are only cached into
            # 'all.neighbors' for 'k' that isn't needed by t-SNE or UMAP.
            all.neighbors <- .find_snn_neighbors(nnbuilt, num.neighbors=cluster.snn.args$num.neighbors, num.threads=num.threads, existing=all.neighbors)
        }
    }

    # Figuring out how many jobs we have.
    all.params <- list()
    if (do.tsne) {
        all.params$tsne <- do.call(.tsne_sweeper, c(list(all.neighbors), tsne.args, list(.env=NULL)))
    }
    if (do.umap) {
        all.params$umap <- do.call(.umap_sweeper, c(list(all.neighbors), umap.args, list(.env=NULL)))
    }
    if (do.cluster.snn) {
        all.params$cluster.snn <- do.call(.snn_sweeper, c(list(all.neighbors), cluster.snn.args, list(.env=NULL)))
    }
    if (do.cluster.kmeans) {
        all.params$cluster.kmeans <- do.call(.kmeans_sweeper, c(list(x), cluster.kmeans.args, list(.env=NULL)))
    }

    # Setting up the cluster and dispatching the jobs.
    njobs <- sum(vapply(all.params, nrow, 0L))
    if (njobs == 1 || num.threads == 1) {
        env <- new.env()
        env$results <- list()
        env$cluster <- NULL
        env$active <- 0L
        threads.per.job <- num.threads
    } else {
        nnodes <- min(njobs, num.threads)
        env <- spawnCluster(nnodes)
        on.exit(stopCluster(env$cluster))
        threads.per.job <- max(1, floor(num.threads / nnodes))
    }
    common.args <- list(.env=env, num.threads=threads.per.job)

    if (do.tsne) {
        do.call(.tsne_sweeper, c(list(all.neighbors), tsne.args, common.args))
    }
    if (do.umap) {
        do.call(.umap_sweeper, c(list(all.neighbors), umap.args, common.args))
    }
    if (do.cluster.snn) {
        do.call(.snn_sweeper, c(list(all.neighbors), cluster.snn.args, common.args))
    }
    if (do.cluster.kmeans) {
        do.call(.kmeans_sweeper, c(list(x), cluster.kmeans.args, common.args))
    }

    # Cleaning up.
    completed <- finishJobs(env)
    output <- list()
    undown.proj.neighbors <- NULL 

    for (n in names(completed)) {
        new.name <- switch(n,
            clusterKmeans="cluster.kmeans",
            runTSNE="tsne",
            runUMAP="umap",
            clusterSNNGraph="cluster.snn",
            stop("unknown result type '", x, "'")
        )

        store <- list(parameters=all.params[[new.name]], results=completed[[n]])

        if (down) {
            if (new.name == "tsne" || new.name == "umap") {
                store <- .undownsample_embedding(m, store)
            } else if (new.name == "cluster.snn") {
                store <- .undownsample_snn(m, store)
            } else if (new.name == "cluster.kmeans") {
                store <- .undownsample_kmeans(m, store)
            }
        }

        output[[new.name]] <- .drop_sweep(store, drop)
    }

    output[order(names(output))]
}
