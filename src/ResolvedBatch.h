#ifndef RESOLVED_BATCH_H
#define RESOLVED_BATCH_H

#include "Rcpp.h"
#include <algorithm>

struct ResolvedBatch {
    ResolvedBatch(Rcpp::Nullable<Rcpp::IntegerVector> batch_) {
        if (batch_.isNotNull()) {
            batch = Rcpp::IntegerVector(batch_);
            ptr = batch.begin();
        }
        return;
    }

    size_t number() const {
        size_t nbatches = 1;
        if (ptr != NULL && batch.size()){
            nbatches = *std::max_element(ptr, ptr + batch.size()) + 1;
        }
        return nbatches;
    }

    // Need to carry along the vector to avoid garbage 
    // collection of any allocated memory (e.g., ALTREP'd)
    // from the Rcpp::Integer initialization.
    Rcpp::IntegerVector batch;
    const int* ptr = NULL;
};

#endif
