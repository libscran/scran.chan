#' Project cells onto an embedding using its neighbors
#'
#' Project cells into a reference embedding (e.g., t-SNE, UMAP) given its neighbors in another space.
#'
#' @param ref Numeric matrix for the reference dataset, where rows are dimensions and columns are cells.
#' This usually contains PC coordinates from \code{\link{runPCA.chan}} or friends.
#' @param embedding Numeric matrix for the embedding of the reference dataset, where rows are cells and columns are dimensions.
#' @param test Numeric matrix for the test dataset, where rows are dimensions and columns are cells.
#' This usually contains PC coordinates similar to \code{ref}.
#' @param k Integer scalar specifying the number of neighbors to use for projection.
#' @param approximate Logical scalar indicating whether an approximate neighbor search should be performed.
#' @param num.threads Integer scalar specifying the number of threads to use.
#'
#' @return Numeric matrix containing the projected embedding coordinates for each cell in the test dataset.
#' Like \code{embedding}, cells are in the rows.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{downsampleByNeighbors.chan}}, which is typically used to generate a reference subsample for easier clustering before assignment of the full (test) dataset.
#'
#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' keep <- downsampleByNeighbors.chan(x, 3) # not too much downsampling
#' sub <- x[,keep]
#'
#' embedding <- runUMAP.chan(sub)
#' full <- projectNeighborEmbedding.chan(sub, embedding, x)
#' plot(full[,1], full[,2], col=iris[,"Species"], 
#'     pch=ifelse(seq_len(ncol(x)) %in% keep, 16, 1)) 
#'
#' @export
projectNeighborEmbedding.chan <- function(ref, embedding, test, k = 20, approximate = TRUE, num.threads = 1) {
    idx <- build_nn_index(ref, approximate=approximate)
    .project_neighbor_embedding(ref, index=idx, embedding=embedding, test=test, k=k, approximate=approximate, num.threads=num.threads)
}

.project_neighbor_embedding <- function(ref, index, embedding, test, k, approximate, num.threads) {
    out <- project_neighbor_embedding(ref_data=ref, ref_index=index, emb_data=t(embedding), test_data=test, k=k, approximate=approximate, nthreads=num.threads)
    out <- t(out)
    rownames(out) <- colnames(test)
    colnames(out) <- colnames(embedding)
    out
}
