#' K-means clustering
#'
#' Perform k-means clustering using kmeans++ initialization with the Hartigan-Wong algorithm. 
#'
#' @param x Numeric matrix where rows are dimensions and columns are cells.
#' @param k Integer scalar specifying the number of clusters.
#' @param init.method String specifying the initialization method:
#' PCA partitioning (\code{"pca-part"}), kmeans++ (\code{"kmeans++"}) or random initialization (\code{"random"}).
#' @param num.threads Integer scalar specifying the number of threads to use.
#' 
#' @return A list containing:
#' \itemize{
#' \item \code{clusters}, an integer vector containing the cluster assignments.
#' \item \code{centers}, a numeric matrix with the coordinates of the cluster centroids (dimensions in rows, centers in columns).
#' \item \code{iterations}, an integer scalar specifying the number of Hartigan-Wong iterations used.
#' \item \code{withinss}, a numeric vector containing the within-cluster sum of squares for each cluster.
#' }
#'
#' @author Aaron Lun
#'
#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' clustering <- clusterKmeans.chan(x)
#' table(clustering$clusters)
#' 
#' @export
clusterKmeans.chan <- function(x, k=10, init.method = "pca-part", num.threads=1) {
    init.choice <- .kmeans_init_choice(init.method)
    output <- cluster_kmeans(x, k, init.choice, nthreads=num.threads)
    .clean_kmeans_clustering(output)
}

.kmeans_init_choice <- function(choice) {
    .choices <- c("pca-part", "kmeans++", "random")
    match(match.arg(choice, .choices), .choices) - 1L
}

.clean_kmeans_clustering <- function(clustering) {
    clustering$clusters <- factor(clustering$clusters + 1L)
    clustering 
}
