#include "config.h"

#include "Rcpp.h"
#include "tatamize.h"

//[[Rcpp::export(rng=false)]]
SEXP combine_matrix(Rcpp::List x, bool byrow) {
    std::vector<std::shared_ptr<tatami::NumericMatrix> > vecs;
    vecs.reserve(x.size());
    for (size_t i = 0, end = x.size(); i < end; ++i) {
        vecs.push_back(extract_NumericMatrix_shared(x[i]));
    }

    std::shared_ptr<tatami::NumericMatrix> output;
    if (byrow) {
        output = tatami::make_DelayedBind<0>(std::move(vecs));
    } else {
        output = tatami::make_DelayedBind<1>(std::move(vecs));
    }

    return new_MatrixChan(std::move(output));
}
