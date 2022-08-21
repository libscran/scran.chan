#' Apply graph-based clustering
#'
#' Perform community detection on a shared nearest neighbor graph.
#'
#' @param x Numeric matrix where rows are dimensions and columns are cells.
#' @param num.neighbors Integer scalar specifying the number of neighbors to use to construct the graph.
#' This may also be a vector to perform a parameter sweep.
#' @param weight.scheme String specifying the weighting scheme to use for constructing the SNN graph.
#' This can be \code{"rank"} (default), \code{"jaccard"} or \code{"number"}.
#' @param method String specifying the community detection method to use.
#' Options are multi-level (\code{"multilevel"}), Walktrap (\code{"walktrap"}) or Leiden (\code{"leiden"}).
#' This may also be a vector to perform a parameter sweep.
#' @param resolution Numeric scalar specifying the resolution to use for multi-level or Leiden clustering.
#' This may also be a vector to perform a parameter sweep.
#' @param steps Integer scalar specifying the number of steps to use for Walktrap clustering.
#' This may also be a vector to perform a parameter sweep.
#' @param seed Integer scalar specifying the seed to use for multi-level or Leiden clustering. 
#' This may also be a vector to perform a parameter sweep.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param drop Logical scalar indicating whether to drop the sweep-based formatting when \code{k}, \code{init.method} and \code{seed} are scalars.
#' @param downsample Integer scalar specifying the downsampling level (see \code{k} in \code{\link{downsampleByNeighbors.chan}}.
#' If \code{NULL}, no downsampling is performed.
#' @param approximate Logical scalar specifying whether to perform an approximate neighbor search.
#' 
#' @return By default, a list containing \code{membership}, an integer vector with cluster assignments for each cell.
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
#' If any of the parameters are vectors or \code{drop = FALSE}, we assume that the user is requesting a parameter sweep.
#' A list is returned containing \code{parameters}, a data.frame with each relevant combination of parameters;
#' and \code{results}, a list of length equal to the number of rows of \code{parameters}, where each element is itself a list as described above.
#' 
#' If \code{downsample} is not \code{NULL}, no additional elements are returned for any \code{method}.
#' Only \code{membership} will be present in the list.
#'
#' @author Aaron Lun
#'
#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' clustering <- clusterSNNGraph.chan(x)
#' clustering$modularity
#' table(clustering$membership)
#' 
#' # With a parameter sweep.
#' swept <- clusterSNNGraph.chan(x,
#'    num.neighbors=c(10, 20),
#'    method=c("multilevel", "walktrap"),
#'    weight.scheme=c("rank", "jaccard")
#' )
#' swept$parameters
#' length(swept$results)
#' table(swept$results[[1]]$membership)
#'
#' @export
clusterSNNGraph.chan <- function(x,
    num.neighbors=10,
    weight.scheme="rank",
    method="multilevel",
    resolution=1, 
    steps=4,
    seed=42,
    drop=TRUE,
    approximate=TRUE,
    downsample=NULL,
    num.threads=1)
{
    down <- do_downsample(downsample)
    if (down) {
        chosen <- downsampleByNeighbors.chan(x, downsample, approximate=approximate, num.threads=num.threads)
        x <- x[,chosen$chosen,drop=FALSE]
        m <- match(chosen$assigned, chosen$chosen)
    }

    nnbuilt <- build_nn_index(x, approximate=approximate)
    all.neighbors <- .find_snn_neighbors(nnbuilt, num.neighbors, num.threads=num.threads)

    sweep <- function(...) { 
        .snn_sweeper(all.neighbors, 
            num.neighbors=num.neighbors, 
            weight.scheme=weight.scheme, 
            method=method, 
            resolution=resolution, 
            steps=steps, 
            seed=seed, 
            ...)
    }

    output <- .sweep_wrapper(sweep, "clusterSNNGraph", num.threads=num.threads)

    if (down) {
        output <- .undownsample_snn(m, output)
    }

    .drop_sweep(output, drop)
}

.find_snn_neighbors <- function(nnbuilt, num.neighbors, num.threads, existing = list()) {
    for (k in num.neighbors) {
        kchar <- as.character(k)
        if (!kchar %in% names(existing)) {
            existing[[kchar]] <- find_nearest_neighbor_indices(nnbuilt, k, nthreads=num.threads)
        }
    }
    existing
}

.undownsample_snn <- function(m, output) {
    for (i in seq_along(output$results)) {
        output$results[[i]] <- list(membership = output$results[[i]]$membership[m])
    }
    output
}

.snn.method.choices <- c("multilevel", "leiden", "walktrap")

.snn.weight.choices <- c("rank", "number", "jaccard")

clusterSNNGraph.chan.core <- function(neighbors, weight.scheme, method, resolution, steps, seed, num.threads) {
    clustering <- cluster_snn_graph(
        nnidx=neighbors,
        weight_scheme=weight.scheme,
        method=method,
        resolution=resolution, 
        steps=steps, 
        seed=seed, 
        nthreads=num.threads)

    clustering$membership <- factor(clustering$membership + 1L)
    if (method=="multilevel") {
        clustering$best <- clustering$best + 1L
        clustering$levels <- lapply(clustering$levels, function(x) factor(x + 1L))
    } else if (method=="walktrap") {
        clustering$merges <- clustering$merges + 1L
    }

    clustering
}

.snn_sweeper <- function(neighbors, num.neighbors, weight.scheme, method, resolution, steps, seed, num.threads, .env) {
    counter <- 0L
    preflight <- is.null(.env)
    parameters <- list()

    for (m in seq_along(method)) {
        method[m] <- match.arg(method[m], .snn.method.choices)
    }
    for (w in seq_along(weight.scheme)) {
        weight.scheme[m] <- match.arg(weight.scheme[w], .snn.weight.choices)
    }

    for (k in num.neighbors) {
        curneighbors <- neighbors[[as.character(k)]]$index

        for (w in weight.scheme) {
            for (m in method) {
                if (m == "walktrap") {
                    for (s in steps) {
                        counter <- counter + 1L
                        if (preflight) {
                            parameters[[counter]] <- data.frame(num.neighbors=k, weight.scheme=w, method=m, resolution=NA_real_, seed=NA_integer_, steps=s)
                        } else {
                            submitJob(.env,
                                fun=clusterSNNGraph.chan.core,
                                args=list(curneighbors, weight.scheme=w, method=m, resolution=0, seed=0, steps=s, num.threads=num.threads),
                                type="clusterSNNGraph",
                                index=counter
                            )
                        }
                    }
                } else {
                    for (r in resolution) {
                        for (s in seed) {
                            counter <- counter + 1L
                            if (preflight) {
                                parameters[[counter]] <- data.frame(num.neighbors=k, weight.scheme=w, method=m, resolution=r, seed=s, steps=NA_integer_)
                            } else {
                                submitJob(.env,
                                    fun=clusterSNNGraph.chan.core,
                                    args=list(curneighbors, weight.scheme=w, method=m, resolution=r, seed=s, steps=0, num.threads=num.threads),
                                    type="clusterSNNGraph",
                                    index=counter
                                )
                            }
                        }
                    }
                }
            }
        }
    }

    if (preflight) {
        do.call(rbind, parameters)
    } else {
        NULL
    }
}
