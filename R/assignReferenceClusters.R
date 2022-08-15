#' Assign cells to reference clusters
#'
#' Given a set of clusters in a reference dataset, assign cluster labels to cells in a test dataset using a nearest neighbor classification scheme.
#'
#' @param ref Numeric matrix for the reference dataset, where rows are dimensions and columns are cells.
#' This usually contains PC coordinates from \code{\link{runPCA.chan}} or friends.
#' @param clusters Vector or factor of length equal to \code{ncol(ref)}, containing the cluster assignment for each cell in the reference dataset.
#' @param test Numeric matrix for the test dataset, where rows are dimensions and columns are cells.
#' This usually contains PC coordinates similar to \code{ref}.
#' @param num.neighbors Integer scalar specifying the number of neighbors in the reference to use. 
#' Smaller values focus on local similarities to immediate neighbors while larger values focus on similarities to the bulk of the cluster.
#' @param approximate Logical scalar indicating whether an approximate neighbor search should be performed.
#' @param num.threads Integer scalar specifying the number of threads to use.
#'
#' @author Aaron Lun
#'
#' @return A data frame containing:
#' \itemize{
#' \item \code{assigned}, a vector or factor of the same type as \code{clusters}, containing cluster assignments for all cells in \code{test}.
#' \item \code{best.prop}, numeric vector containing the proportion of neighbors supporting the assignment of each cell in \code{assigned}.
#' \item \code{second.prop}, numeric vector containing the second-largest proportion of neighbors. 
#' }
#' 
#' @seealso
#' \code{\link{downsampleByNeighbors.chan}}, which is typically used to generate a reference subsample for easier clustering before assignment of the full (test) dataset.
#'
#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' keep <- downsampleByNeighbors.chan(x, 5)
#' sub <- x[,keep]
#'
#' labs <- iris[keep,"Species"]
#' full <- assignReferenceClusters.chan(sub, labs, x)
#' table(full$assigned, iris[,"Species"]) # compare to reference
#' @export
assignReferenceClusters.chan <- function(ref, clusters, test, num.neighbors=20, approximate=TRUE, num.threads=1) {
    if (is.matrix(ref)) {
        idx <- build_nn_index(ref, approximate=approximate)
    } else {
        idx <- ref
    }

    iclusters <- as.integer(factor(clusters)) - 1L
    out <- assign_reference_clusters(idx, ref_clusters=iclusters, test_data=test, k=num.neighbors, approximate=approximate, nthreads=num.threads)
    out$assigned <- clusters[match(out$assigned, iclusters)]

    df <- data.frame(out)
    rownames(df) <- colnames(test)
    df
}
