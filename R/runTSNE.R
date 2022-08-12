#' Run an approximate t-SNE
#'
#' Compute the t-stochastic neighbor embedding with assorted approximations for speed.
#'
#' @param x Numeric matrix where rows are dimensions and columns are cells.
#' @param perplexity Numeric scalar specifying the perplexity to use in the t-SNE algorithm.
#' This may also be a vector to perform a parameter sweep.
#' @param interpolate Integer scalar specifying the grid resolution for interpolating repulsive forces (larger is slower but more accurate).
#' A value of zero disables interpolation, while a value of -1 only uses interpolation for large datasets.
#' This may also be a vector to perform a parameter sweep.
#' @param max.depth Integer scalar specifying the maximum depth of the Barnes-Hut quad trees (larger is slower but more accurate).
#' Set to a large integer (e.g., 1000) to eliminate depth restrictions.
#' This may also be a vector to perform a parameter sweep.
#' @param seed Integer scalar specifying the seed to use for generating the initial coordinates.
#' This may also be a vector to perform a parameter sweep.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param drop Logical scalar indicating whether to drop the sweep-based formatting when the parameters are scalars.
#' @param downsample Integer scalar specifying the downsampling level (see \code{k} in \code{\link{downsampleByNeighbors.chan}}.
#' If \code{NULL}, no downsampling is performed.
#' @param approximate Logical scalar specifying whether to perform an approximate neighbor search.
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
#' embedding <- runTSNE.chan(x)
#' plot(embedding[,1], embedding[,2], col=iris[,5])
#'
#' # Parameter sweep mode.
#' swept <- runTSNE.chan(x, 
#'     perplexity=c(10, 20, 50),
#'     seed=c(42, 88))
#' swept$parameters
#' head(swept$results[[1]])
#' 
#' @export
runTSNE.chan <- function(x, perplexity=30, interpolate=-1, max.depth=7, seed=42, drop=TRUE, approximate=TRUE, downsample=NULL, num.threads=1) {
    if (!is.null(downsample)) {
        original <- x
        chosen <- downsampleByNeighbors.chan(x, downsample, approximate=approximate, num.threads=num.threads)
        x <- x[,chosen,drop=FALSE]
    }

    nnbuilt <- build_nn_index(x, approximate=approximate)
    all.neighbors <- .find_tsne_neighbors(nnbuilt, perplexity, num.threads=num.threads)

    sweep <- function(...) {
        .tsne_sweeper(all.neighbors, 
            perplexity=perplexity, 
            interpolate=interpolate, 
            max.depth=max.depth, 
            seed=seed, 
            ...)
    }

    output <- .sweep_wrapper(sweep, "runTSNE", num.threads=num.threads)

    if (!is.null(downsample)) {
        output <- .undownsample_embedding(output, nnbuilt, original, approximate=approximate, num.threads=num.threads)
    }

    .drop_sweep(output, drop)
}

.find_tsne_neighbors <- function(nnbuilt, perplexity, num.threads, existing = list()) {
    for (p in perplexity) {
        k <- perplexity_to_neighbors(p)
        kchar <- as.character(k)
        if (!kchar %in% names(existing)) {
            existing[[kchar]] <- find_nearest_neighbors(nnbuilt, k, nthreads=num.threads)
        }
    }
    existing
}

runTSNE.chan.core <- function(neighbors, perplexity, interpolate, max.depth, seed, num.threads) {
    output <- run_tsne(
        nnidx=neighbors$index, 
        nndist=neighbors$distance, 
        perplexity=perplexity,
        interpolate=interpolate,
        max_depth=max.depth,
        seed=seed,
        nthreads=num.threads
    )
    t(output)
}

.tsne_sweeper <- function(neighbors, perplexity, interpolate, max.depth, seed, num.threads, .env) {
    counter <- 0L
    preflight <- is.null(.env)
    parameters <- list()

    for (p in perplexity) {
        k <- perplexity_to_neighbors(p)
        curneighbors <- neighbors[[as.character(k)]]

        for (i in interpolate) {
            for (d in max.depth) {
                for (s in seed) {
                    counter <- counter + 1L
                    if (preflight) {
                        parameters[[counter]] <- data.frame(perplexity = p, interpolate = i, max.depth = d, seed=s)
                    } else {
                        submitJob(.env, 
                            fun=runTSNE.chan.core, 
                            args=list(neighbors=curneighbors, perplexity=p, interpolate=i, max.depth=d, seed=s, num.threads=num.threads),
                            type="runTSNE", 
                            index=counter
                        )
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
