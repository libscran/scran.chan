#include "config.h"

#include "knncolle.h"
#include "scran/dimensionality_reduction/ProjectNeighborEmbedding.hpp"
#include "Rcpp.h"

//[[Rcpp::export(rng=false)]]
SEXP project_neighbor_embedding(Rcpp::NumericMatrix ref_data, SEXP ref_index, Rcpp::NumericMatrix emb_data, Rcpp::NumericMatrix test_data, int k, bool approximate, int nthreads) {
    KnncollePtr index(ref_index);
    size_t ref_nr = index->ndim(), ref_nc = index->nobs();
    auto ref_ptr = static_cast<const double*>(ref_data.begin());
    if (ref_nr != ref_data.nrow() || ref_nc != ref_data.ncol()) {
        throw std::runtime_error("reference matrix and indices should have the same dimensions");
    }

    size_t test_nr = test_data.nrow(), test_nc = test_data.ncol();
    auto test_ptr = static_cast<const double*>(test_data.begin());
    if (test_nr != ref_nr) {
        throw std::runtime_error("test and reference matrices should have the same number of rows");
    }

    size_t emb_nr = emb_data.nrow(), emb_nc = emb_data.ncol();
    auto emb_ptr = static_cast<const double*>(emb_data.begin());
    if (emb_nc != ref_nc) {
        throw std::runtime_error("reference and embedding matrices should have the same number of columns");
    }

    scran::ProjectNeighborEmbedding runner;
    runner
        .set_num_neighbors(k)
        .set_approximate(approximate)
        .set_num_threads(nthreads);

    Rcpp::NumericMatrix output(emb_nr, test_nc);
    runner.run(ref_nr, ref_nc, ref_ptr, test_nc, test_ptr, emb_nr, emb_ptr, index.get(), static_cast<double*>(output.begin()));

    return output; 
}
