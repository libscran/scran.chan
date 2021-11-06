#ifndef RESOLVED_FEATURES_H
#define RESOLVED_FEATURES_H

#include "Rcpp.h"

struct ResolvedFeatures {
    ResolvedFeatures(Rcpp::Nullable<Rcpp::LogicalVector> features_) {
        if (features_.isNotNull()) {
            features = Rcpp::LogicalVector(features_);
            ptr = features.begin();
        }
        return;
    }

    Rcpp::LogicalVector features;
    const int* ptr = NULL;
};

#endif
