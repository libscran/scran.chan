#' Run UMAP
#'
#' Compute the uniform manifold approximation and projection. 
#'
#' @param x Numeric matrix where rows are dimensions and columns are cells.
#' @param num.neighbors Integer scalar specifying the number of neighbors to use in the UMAP algorithm.
#' This may also be a vector to perform a parameter sweep.
#' @param min.dist Numeric scalar specifying the minimum distance between points.
#' This may also be a vector to perform a parameter sweep.
#' @param seed Integer scalar specifying the seed to use. 
#' This may also be a vector to perform a parameter sweep.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' 
#' @return 
#' By default, a numeric matrix where rows are cells and columns are the two dimensions of the embedding.
#'
#' If any of the parameters are vectors or \code{drop = FALSE}, we assume that the user is requesting a parameter sweep.
#' A list is returned containing \code{parameters}, a data.frame with each relevant combination of parameters;
#' and \code{results}, a list of numeric matrices of length equal to the number of rows of \code{parameters}, containing the embeddings for each combination.
#' 
#' @author Aaron Lun
#'
#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' embedding <- runUMAP.chan(x)
#' plot(embedding[,1], embedding[,2], col=iris[,5])
#' 
#' # Parameter sweep mode.
#' swept <- runUMAP.chan(x, 
#'     num.neighbors=c(10, 20),
#'     seed=c(123, 456))
#' swept$parameters
#' head(swept$results[[1]])
#'
#' @export
runUMAP.chan <- function(x, num.neighbors=15, min.dist=0.01, seed=1234567890, drop=TRUE, num.threads=1) {
    nnbuilt <- build_nn_index(x)
    all.neighbors <- .find_umap_neighbors(nnbuilt, num.neighbors, num.threads=num.threads)
    sweep <- function(...) .umap_sweeper(all.neighbors, num.neighbors=num.neighbors, min.dist=min.dist, seed=seed, ...)
    .sweep_wrapper(sweep, "runUMAP", num.threads=num.threads, drop=drop)
}

.find_umap_neighbors <- function(nnbuilt, num.neighbors, num.threads, existing = list()) {
    for (k in num.neighbors) {
        kchar <- as.character(k)
        if (!kchar %in% names(existing)) {
            existing[[kchar]] <- find_nearest_neighbors(nnbuilt, k, nthreads=num.threads)
        }
    }
    existing
}

runUMAP.chan.core <- function(neighbors, min.dist, seed, num.threads) {
    output <- run_umap(
        nnidx=neighbors$index,
        nndist=neighbors$distance,
        min_dist=min.dist,
        seed=seed,
        nthreads=num.threads
    )
    t(output)
}

.umap_sweeper <- function(neighbors, num.neighbors, min.dist, seed, num.threads, .env) {
    counter <- 0L
    preflight <- is.null(.env)
    parameters <- list()

    for (k in num.neighbors) {
        curneighbors <- neighbors[[as.character(k)]]

        for (d in min.dist) {
            for (s in seed) {
                counter <- counter + 1L
                if (preflight) {
                    parameters[[counter]] <- data.frame(num.neighbors = k, min.dist = d, seed=s)
                } else {
                    submitJob(.env,
                        fun=runUMAP.chan.core,
                        args=list(neighbors=curneighbors, min.dist=d, seed=s, num.threads=num.threads),
                        type="runUMAP",
                        index=counter
                    )
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
