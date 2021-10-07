#ifndef RUN_UMAP_H
#define RUN_UMAP_H
#include "umappp/umappp.hpp"
#include <vector>
#include "Rcpp.h"

typedef umappp::Umap<float> umap;

struct InitializedUmap {
    InitializedUmap(umap r, umap::Status s, std::vector<float> e, int n) : runner(std::move(r)), status(std::move(s)), embedding(std::move(e)), nobs(n) {}

    void run() {
        runner.run(status, 2, embedding.data());
        return;
    }

    Rcpp::NumericMatrix yield() {
        Rcpp::NumericMatrix output(2, nobs);
        std::copy(embedding.begin(), embedding.end(), output.begin());
        return output;
    }

    umap runner;
    umap::Status status;
    std::vector<float> embedding;
    int nobs;
};

#endif
