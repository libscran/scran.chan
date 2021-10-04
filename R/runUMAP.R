#' Run UMAP
#'
#' Compute the uniform manifold approximation and projection. 
#'
#' @param x Numeric matrix where rows are dimensions and columns are cells.
#' @param num.neighbors Integer scalar specifying the number of neighbors to use in the UMAP algorithm.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' 
#' @return A numeric matrix where rows are cells and columns are the two dimensions of the embedding.
#' 
#' @author Aaron Lun
#'
#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' embedding <- runUMAP.chan(x)
#' plot(embedding[,1], embedding[,2], col=iris[,5])
#' 
#' @export
runUMAP.chan <- function(x, num.neighbors=15, num.threads=1) {
    neighbors <- build_nn_index(x)
    init <- initialize_umap(neighbors, num.neighbors, num.threads)
    output <- run_umap(init)
    t(output)
}
