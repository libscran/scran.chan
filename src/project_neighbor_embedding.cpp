#include "config.h"

#include "knncolle.h"
#include "scran/dimensionality_reduction/ProjectNeighborEmbedding.hpp"
#include "Rcpp.h"

//[[Rcpp::export(rng=false)]]
SEXP project_neighbor_embedding(Rcpp::NumericMatrix ref_data, Rcpp::IntegerMatrix nnidx, Rcpp::NumericMatrix nndist, Rcpp::NumericMatrix emb_data, Rcpp::NumericMatrix test_data, int k, bool approximate, int nthreads) {
    size_t ref_nr = ref_data.nrow(), ref_nc = ref_data.ncol();
    auto ref_ptr = static_cast<const double*>(ref_data.begin());

    size_t test_nr = test_data.nrow(), test_nc = test_data.ncol();
    auto test_ptr = static_cast<const double*>(test_data.begin());
    if (test_nr != ref_nr) {
        throw std::runtime_error("test and reference matrices should have the same number of rows");
    }

    auto neighbors = unpack_neighbors<int, double>(nnidx, nndist);
    if (test_nc != neighbors.size()) {
        throw std::runtime_error("number of cells in test matrix and neighbor list should be the same");
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
    runner.run(ref_nr, ref_nc, ref_ptr, test_nc, test_ptr, emb_nr, emb_ptr, neighbors, static_cast<double*>(output.begin()));

    return output; 
}
