# A slimmed-down version of scran

⚠️ ⚠️ ⚠️ ⚠️  **This repository is deprecated, use the [scrapper](https://github.com/libscran/scrapper) package instead.** ⚠️ ⚠️ ⚠️ ⚠️  

## Overview

**scran.chan** is a slimmed-down version of the [**scran**](https://bioconductor.org/packages/scran) Bioconductor package,
refactored to use the [**libscran**](https://github.com/LTLA/libscran) C++ library for all of the computation.
This provides methods for an end-to-end analysis of a single-cell RNA-sequencing (scRNA-seq) analysis,
starting from the count matrix and finishes with clusters, markers, and various embeddings (i.e., t-SNE and UMAP).
It's pretty fast and memory-efficient - some casual timings show that 170,000 cells can be analyzed in 5 minutes with 8 threads.

## Installation

Compilation is straightforward provided you have a reasonably up-to-date CMake and a compiler that supports C++17.

```r
devtools::install_github("LTLA/scran.chan")
```

The compilation can take some time, so just be patient.

You'll also want to use a compiler that supports OpenMP to take advantage of parallelization.
Otherwise, the package can still be installed but will not be responsive to any setting of `num.threads`.

## Usage

Once installed, usage is simple:

```r
# Grabbing a test dataset.
library(scRNAseq)
x <- ZeiselBrainData()
mito <- grep("^mt-", rownames(x))

# Running the analysis to get the results.
out <- quickBasicAnalysis(assay(x), qc.subsets=list(Mt=mito)) 

# Optionally packaging into a SingleCellExperiment.
sce <- marshalToSCE(x, out, assay.type=1)
```

Advanced users can use the various `*.chan()` functions to call specific steps.
Note that this requires a matrix to be initialized via `initializeSparseMatrix()`; they will not work with normal R matrices.

## Links

The [**scran.js**](https://www.npmjs.com/package/scran.js) package provides Javascript bindings to the same C++ libraries.

The [**scran-cli**](https://github.com/LTLA/scran-cli) tool provides a command-line interface to a basic scRNA-seq analysis.

## TODO

Add the `quickMergedAnalysis()` pipeline function.
