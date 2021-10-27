#' K-means clustering
#'
#' Perform k-means clustering using kmeans++ initialization with the Hartigan-Wong algorithm. 
#'
#' @param x Numeric matrix where rows are dimensions and columns are cells.
#' @param k Integer scalar specifying the number of clusters.
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
clusterKmeans.chan <- function(x, k=10, num.threads=1) {
    output <- cluster_kmeans(x, k, nthreads=num.threads)
    .clean_kmeans_clustering(output)
}

.clean_kmeans_clustering <- function(clustering) {
    clustering$clusters <- clustering$clusters + 1L
    clustering 
}
