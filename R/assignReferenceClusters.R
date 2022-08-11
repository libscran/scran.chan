#' Assign cells to reference clusters
#'
#' Given a set of clusters in a reference dataset, assign cluster labels to cells in a test dataset using a \pkg{SingleR}-like classification scheme.
#'
#' @param ref Numeric matrix for the reference dataset, where rows are dimensions and columns are cells.
#' This usually contains PC coordinates from \code{\link{runPCA.chan}} or friends.
#' @param clusters Vector or factor of length equal to \code{ncol(ref)}, containing the cluster assignment for each cell in the reference dataset.
#' @param test Numeric matrix for the test dataset, where rows are dimensions and columns are cells.
#' This usually contains PC coordinates similar to \code{ref}.
#' @param quantile Numeric scalar specifying the quantile to use for converting distances to an assignment score for each cluster.
#' Smaller values focus on local similarities to immediate neighbors while larger values focus on similarities to the bulk of the cluster.
#' @param approximate Logical scalar indicating whether an approximate neighbor search should be performed.
#' @param num.threads Integer scalar specifying the number of threads to use.
#'
#' @author Aaron Lun
#'
#' @return
#' Vector or factor of the same type as \code{clusters}, containing cluster assignments for all cells in \code{test}.
#' 
#' @seealso
#' \code{\link{downsampleByNeighbors.chan}}, which is typically used to generate a reference subsample for easier clustering before assignment of the full (test) dataset.
#'
#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' keep <- downsampleByNeighbors.chan(x, 10)
#' sub <- x[,keep]
#'
#' labs <- iris[keep,"Species"]
#' full <- assignReferenceClusters.chan(sub, labs, x)
#' table(full, iris[,"Species"]) # compare to reference
#' @export
assignReferenceClusters.chan <- function(ref, clusters, test, quantile=0.2, approximate=TRUE, num.threads=1) {
    iclusters <- as.integer(factor(clusters)) - 1L
    out <- assign_reference_clusters(ref_data=ref, ref_clusters=iclusters, test_data=test, quantile=quantile, approximate=approximate, nthreads=num.threads)
    clusters[match(out, iclusters)]
}
