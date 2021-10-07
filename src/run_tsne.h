#ifndef RUN_TSNE_H
#define RUN_TSNE_H
#include "qdtsne/qdtsne.hpp"
#include <vector>
#include "Rcpp.h"

typedef qdtsne::Tsne<2, float> MyTsne;

struct InitializedTsne {
    InitializedTsne(MyTsne r, MyTsne::Status<int> s, std::vector<float> e, int n) : runner(std::move(r)), status(std::move(s)), embedding(std::move(e)), nobs(n) {}

    void run() {
        runner.run(status, embedding.data());
        return;
    }

    Rcpp::NumericMatrix yield() {
        Rcpp::NumericMatrix output(2, nobs);
        std::copy(embedding.begin(), embedding.end(), output.begin());
        return output;
    }
    
    MyTsne runner;
    MyTsne::Status<int> status;
    std::vector<float> embedding;
    int nobs;
};

#endif
