#' Marshal results into a SingleCellExperiment
#'
#' Marshal the results of \code{\link{quickBasicAnalysis}} into a SingleCellExperiment object.
#'
#' @param x A SummarizedExperiment or SingleCellExperiment object containing the count matrix used in \code{\link{quickBasicAnalysis}}.
#' @param results A named list containing the output of \code{\link{quickBasicAnalysis}}.
#' @param assay.type String or integer scalar specifying the assay of \code{x} containing the counts.
#' @param include.gene.data Logical scalar indicating whether to include per-gene statistics in the \code{rowData}.
#'
#' @return A SingleCellExperiment decorated with size factors, reduced dimensions and cluster annotations.
#'
#' @author Aaron Lun
#'
#' @details
#' Size factors are stored in the \code{sizeFactors} of the output.
#' We also compute log-transformed normalized expression values and store them in the \code{"logcounts"} assay.
#'
#' All dimensionality reduction results (PCA, t-SNE and UMAP) are stored in the \code{reducedDims}.
#'
#' Clustering results are stored in \code{clusters} columns of the \code{colData}.
#' For multi-level community detection, a clustering result at level \code{N} is stored as \code{clusters.N}.
#' The clustering with the highest modularity is stored in the \code{colLabels}.
#'
#' If \code{include.gene.data=TRUE}, per-gene statistics are stored as nested DataFrames inside the \code{rowData}.
#' Variance modelling statistics are stored in the \code{variances} column.
#'
#' @examples
#' library(scRNAseq)
#' x <- ZeiselBrainData()
#' mito <- grep("^mt-", rownames(x))
#' out <- quickBasicAnalysis(assay(x), qc.subsets=list(Mt=mito), num.threads=1)
#' sce <- marshalToSCE(x, out, assay.type=1)
#' sce
#'
#' @export
marshalToSCE <- function(x, results, assay.type="counts", include.gene.data=FALSE) {
    x <- as(x, "SingleCellExperiment")

    keep <- !results$qc.discard
    x <- x[,keep]
    qc.df <- S4Vectors::DataFrame(results$qc.metrics)
    SummarizedExperiment::colData(x) <- BiocGenerics::cbind(SummarizedExperiment::colData(x), qc.df[keep,,drop=FALSE])

    sf <- results$size.factors
    SingleCellExperiment::sizeFactors(x) <- sf
    ass <- SummarizedExperiment::assay(x, assay.type)
    SummarizedExperiment::assay(x, "logcounts") <- log1p(t(t(DelayedArray::DelayedArray(ass))/sf))/log(2)

    SingleCellExperiment::reducedDim(x, "PCA") <- results$pca$components
    SingleCellExperiment::reducedDim(x, "TSNE") <- results$tsne
    SingleCellExperiment::reducedDim(x, "UMAP") <- results$umap

    clust.info <- results$clusters.snn$membership
    for (i in seq_along(clust.info)) {
        SummarizedExperiment::colData(x, paste0("clusters.", i)) <- results$cluster.snn$membership[[i]]
    }
    m <- results$cluster.snn$best
    SingleCellExperiment::colLabels(x) <- results$cluster.snn$membership[[m]]

    if (include.gene.data) {
        SummarizedExperiment::rowData(x)$variances <- S4Vectors::DataFrame(results$variances)
    }

    x
}
