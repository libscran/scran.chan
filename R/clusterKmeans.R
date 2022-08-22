#' K-means clustering
#'
#' Perform k-means clustering using kmeans++ initialization with the Hartigan-Wong algorithm. 
#'
#' @param x Numeric matrix where rows are dimensions and columns are cells.
#' @param k Integer scalar specifying the number of clusters.
#' This may also be a vector to perform a parameter sweep.
#' @param init.method String specifying the initialization method:
#' PCA partitioning (\code{"pca-part"}), kmeans++ (\code{"kmeans++"}) or random initialization (\code{"random"}).
#' This may also be a vector to perform a parameter sweep.
#' @param seed Integer scalar specifying the seed to use for random or kmeans++ initialization.
#' This may also be a vector to perform a parameter sweep.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param drop Logical scalar indicating whether to drop the sweep-based formatting when \code{k}, \code{init.method} and \code{seed} are scalars.
#' @param downsample Integer scalar specifying the downsampling level (see \code{k} in \code{\link{downsampleByNeighbors.chan}}.
#' If \code{NULL}, no downsampling is performed.
#' 
#' @return 
#' By default, a list is returned containing:
#' \itemize{
#' \item \code{clusters}, a factor containing the cluster assignments for each cell.
#' \item \code{centers}, a numeric matrix with the coordinates of the cluster centroids (dimensions in rows, centers in columns).
#' \item \code{iterations}, an integer scalar specifying the number of Hartigan-Wong iterations used.
#' \item \code{withinss}, a numeric vector containing the within-cluster sum of squares for each cluster.
#' }
#'
#' If any of the parameters are vectors or \code{drop = FALSE}, we assume that the user is requesting a parameter sweep.
#' A list is returned containing \code{parameters}, a data.frame with each relevant combination of parameters;
#' and \code{results}, a list of length equal to the number of rows of \code{parameters}, where each entry contains the result for the corresponding parameter combination.
#'
#' If \code{downsample} is not \code{NULL}, only the \code{clusters} element of the list will be reported.
#' @author Aaron Lun
#'
#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' clustering <- clusterKmeans.chan(x, k=3)
#' table(clustering$clusters, iris[,"Species"])
#'
#' # Parameter sweep mode.
#' swept <- clusterKmeans.chan(x, k=c(2,5,10), 
#'     init.method=c("pca-part", "random"), 
#'     seed=c(1,2,3))
#' swept$parameters
#' length(swept$results)
#' table(swept$results[[1]]$clusters)
#'
#' # Adding some downsampling.
#' down.clust <- clusterKmeans.chan(x, k=3, downsample=5)
#' table(down.clust$clusters, iris[,"Species"])
#' 
#' @export
clusterKmeans.chan <- function(x, k=10, init.method = "pca-part", seed=5489L, drop=TRUE, downsample=NULL, num.threads=1) {
    down <- do_downsample(downsample)
    if (down) {
        chosen <- downsampleByNeighbors.chan(x, downsample, num.threads=num.threads)
        x <- x[,chosen$chosen,drop=FALSE]
        m <- match(chosen$assigned, chosen$chosen)
    }

    sweep <- function(...) .kmeans_sweeper(x, k=k, init.method=init.method, seed=seed, ...)
    output <- .sweep_wrapper(sweep, "clusterKmeans", num.threads=num.threads)

    if (down) {
        output <- .undownsample_kmeans(m, output)
    }

    .drop_sweep(output, drop)
}

.kmeans.init.choices <- c("pca-part", "kmeans++", "random")

clusterKmeans.chan.core <- function(x, k, init.method, seed, num.threads) {
    output <- cluster_kmeans(x, k, init.method, seed=seed, nthreads=num.threads)
    output$clusters <- factor(output$clusters + 1L)
    output 
}

.undownsample_kmeans <- function(m, output) {
    for (i in seq_along(output$results)) {
        output$results[[i]] <- list(clusters=output$results[[i]]$clusters[m])
    }
    output
}

.kmeans_sweeper <- function(x, k, init.method, seed, num.threads, .env) {
    counter <- 0L
    preflight <- is.null(.env)
    parameters <- list()

    for (init in init.method) {
        init <- match.arg(init, .kmeans.init.choices)
        for (k0 in k) {
            if (init == "pca-part") {
                counter <- counter + 1L
                if (preflight) {
                    parameters[[counter]] <- data.frame(k = k0, init.method = init, seed = NA_integer_)
                } else {
                    submitJob(.env, 
                        fun=clusterKmeans.chan.core, 
                        args=list(x=x, init.method=init, k=k0, seed=0, num.threads=num.threads),
                        type="clusterKmeans", 
                        index=counter
                    )
                }
            } else {
                for (s in seed) {
                    counter <- counter + 1L
                    if (preflight) {
                        parameters[[counter]] <- data.frame(k = k0, init.method = init, seed = s)
                    } else {
                        submitJob(.env,
                            fun=clusterKmeans.chan.core, 
                            args=list(x=x, init.method=init, k=k0, seed=s, num.threads=num.threads),
                            type="clusterKmeans", 
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

