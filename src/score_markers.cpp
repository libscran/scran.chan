#include <iostream>
#include "scran/differential_analysis/ScoreMarkers.hpp"
#include "Rcpp.h"
#include "tatamize.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List score_markers(SEXP x, Rcpp::IntegerVector groups) {
    auto mat = extract_NumericMatrix(x);
    const int ngroups = (groups.size() ? *std::max_element(groups.begin(), groups.end()) + 1 : 1);

    std::vector<Rcpp::NumericVector> means, detected;
    means.reserve(ngroups);
    detected.reserve(ngroups);
    std::vector<double*> mptrs(ngroups), dptrs(ngroups);

    for (int g = 0; g < ngroups; ++g) {
        means.emplace_back(mat->nrow());
        detected.emplace_back(mat->nrow());
        mptrs[g] = static_cast<double*>(means.back().begin());
        dptrs[g] = static_cast<double*>(detected.back().begin());
    }

    std::vector<std::vector<Rcpp::NumericVector> > cohens(ngroups), aucs(ngroups);
    std::vector<std::vector<double*> > cptrs(scran::differential_analysis::n_summaries);
    cptrs[scran::differential_analysis::MIN].resize(ngroups);
    cptrs[scran::differential_analysis::MEAN].resize(ngroups);
    cptrs[scran::differential_analysis::MIN_RANK].resize(ngroups);
    auto aptrs = cptrs;

    for (int g = 0; g < ngroups; ++g) {
        cohens[g].reserve(3);
        aucs[g].reserve(3);
        for (int i = 0; i < 3; ++i) {
            cohens[g].emplace_back(mat->nrow());
            aucs[g].emplace_back(mat->nrow());
        }

        cptrs[scran::differential_analysis::MIN][g] = static_cast<double*>(cohens[g][0].begin());
        cptrs[scran::differential_analysis::MEAN][g] = static_cast<double*>(cohens[g][1].begin());
        cptrs[scran::differential_analysis::MIN_RANK][g] = static_cast<double*>(cohens[g][2].begin());
        
        aptrs[scran::differential_analysis::MIN][g] = static_cast<double*>(aucs[g][0].begin());
        aptrs[scran::differential_analysis::MEAN][g] = static_cast<double*>(aucs[g][1].begin());
        aptrs[scran::differential_analysis::MIN_RANK][g] = static_cast<double*>(aucs[g][2].begin());
    }

    scran::ScoreMarkers runner;
    runner.run(mat, static_cast<const int*>(groups.begin()), std::move(mptrs), std::move(dptrs), std::move(cptrs), std::move(aptrs));

    // Organizing the output.
    Rcpp::List cohen_out(ngroups);
    Rcpp::List auc_out(ngroups);
    for (int g = 0; g < ngroups; ++g) {
        cohen_out[g] = Rcpp::List(cohens[g].begin(), cohens[g].end());
        auc_out[g] = Rcpp::List(aucs[g].begin(), aucs[g].end());
    }

    return Rcpp::List::create(
        Rcpp::Named("means") = Rcpp::List(means.begin(), means.end()),
        Rcpp::Named("detected") = Rcpp::List(detected.begin(), detected.end()),
        Rcpp::Named("cohen") = cohen_out,
        Rcpp::Named("auc") = auc_out
    );
}
