quickAdtSizeFactors <- function(x, size.factors=NULL, batch=NULL, num.comp=25, num.clusters=20, num.threads=1) {
    norm <- logNormCounts.chan(x, size.factors=size.factors, batch=batch)
    pcs <- runPCA.chan(norm, num.comp=num.comp, num.threads=num.threads)
    clust <- clusterKmeans.chan(pcs$components, num.clusters=num.clusters, num.threads=num.threads)
    groupedSizeFactors(x, clust$clusters)
}

