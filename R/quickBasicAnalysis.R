#' Perform a simple scRNA-seq analysis
#'
#' Analyze a simple single-cell RNA-seq dataset, from quality control to marker gene detection.
#' No frills, no dependencies.
#'
#' @param x A sparse matrix-like object.
#' @param qc.subsets Feature control subsets, see \code{?\link{perCellQCMetrics.chan}} for details.
#' @param qc.nmads Filter threshold stringency, see \code{?\link{perCellQCFilters.chan}} for details.
#' @param trend.span Trend smoothing, see \code{?\link{modelGeneVar.chan}} for details.
#' @param hvg.num Integer scalar specifying the number of HVGs to use for PCA.
#' @param pca.num Number of top PCs, see \code{?\link{runPCA.chan}} for details.
#' @inheritParams runAllNeighbors
#'
#' @return A list containing the elements:
#' \itemize{
#'     \item \code{qc.metrics}, a list containing the output of \code{\link{perCellQCMetrics.chan}}.
#'     \item \code{qc.filters}, a list containing the output of \code{\link{perCellQCFilters.chan}}.
#'     \item \code{qc.discard}, a logical vector specifying which cells to discard during quality control.
#'     This is the same as \code{qc.filters$filters$overall}.
#'     \item \code{size.factors}, a numeric vector containing the size factors for scaling normalization.
#'     \item \code{variances}, a list containing the output of \code{\link{modelGeneVar.chan}}.
#'     This contains a \code{keep} column indicating whether the gene is retained for PCA.
#'     \item \code{pca}, a list containing the output of \code{\link{runPCA.chan}}.
#'     \item \code{cluster.snn}, a list containing the output of \code{\link{clusterSNNGraph.chan}}.
#'     \item \code{umap}, a list containing the output of \code{\link{runUMAP.chan}}.
#'     \item \code{tsne}, a list containing the output of \code{\link{runTSNE.chan}}.
#' }
#' 
#' @author Aaron Lun
#' @examples
#' library(scRNAseq)
#' x <- ZeiselBrainData()
#' mito <- grep("^mt-", rownames(x))
#' out <- quickBasicAnalysis(assay(x), qc.subsets=list(Mt=mito), num.threads=1)
#'
#' @export 
quickBasicAnalysis <- function(x, 
    qc.subsets,
    qc.nmads=3,
    trend.span=0.3,
    hvg.num=2000,
    pca.num=50,
    tsne.perplexity=30, 
    umap.num.neighbors=15, 
    cluster.snn.num.neighbors=10, 
    cluster.snn.resolution=1, 
    num.threads=1
) {
    x <- initializeSparseMatrix(x, num.threads=1)

    results <- list()
    results$qc.metrics <- qc.metrics <- perCellQCMetrics.chan(x, qc.subsets, num.threads=num.threads)
    results$qc.filters <- qc.filters <- perCellQCFilters.chan(qc.metrics$sums, qc.metrics$detected, qc.metrics$subsets, nmads=qc.nmads)
    results$qc.discard <- qc.filters$filters$overall
    x <- filterCells.chan(x, results$qc.discard)

    sf <- qc.metrics$sums[!qc.filters$filters$overall]
    sf <- sf/mean(sf)
    results$size.factors <- sf
    x <- logNormCounts.chan(x, sf)

    results$variances <- modelGeneVar.chan(x, span = trend.span, num.threads=num.threads)
    keep <- rank(-results$variances$residuals, ties.method="first") <= hvg.num
    results$variances$keep <- keep

    results$pca <- runPCA.chan(x, num.comp=pca.num, subset=keep, num.threads=num.threads)
    pcs <- results$pca$components 
    results$pca$components <- t(pcs)

    neighbor.out <- runAllNeighbors(pcs, 
        tsne.perplexity=tsne.perplexity, 
        umap.num.neighbors=umap.num.neighbors, 
        cluster.snn.num.neighbors=cluster.snn.num.neighbors, 
        cluster.snn.resolution=cluster.snn.resolution, 
        num.threads=num.threads)

    marker.out <- scoreMarkers.chan(x, neighbor.out$cluster.snn$membership[[neighbor.out$cluster.snn$best]])

    results <- c(results, neighbor.out)
    results$markers <- marker.out
    results
}
