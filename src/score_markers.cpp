#include "config.h"

#include "scran/differential_analysis/ScoreMarkers.hpp"
#include "scran/utils/average_vectors.hpp"

#include "Rcpp.h"
#include "tatamize.h"

#include "vector_to_pointers.h"
#include "ResolvedBatch.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List score_markers(SEXP x, Rcpp::IntegerVector groups, Rcpp::Nullable<Rcpp::IntegerVector> batch, bool simple_means_only, double lfc, int nthreads) {
    auto mat = extract_NumericMatrix(x);
    size_t NR = mat->nrow();
    const int ngroups = (groups.size() ? *std::max_element(groups.begin(), groups.end()) + 1 : 1);

    auto batch_info = ResolvedBatch(batch);
    auto bptr = batch_info.ptr;
    size_t nbatches = batch_info.number();

    // Setting up all the damn output vectors (I) - for mean & detected
    std::vector<std::vector<Rcpp::NumericVector> > means(ngroups), detected(ngroups);
    for (int g = 0; g < ngroups; ++g) {
        means[g].reserve(nbatches);
        detected[g].reserve(nbatches);
        for (size_t b = 0; b < nbatches; ++b) {
            means[g].emplace_back(NR);
            detected[g].emplace_back(NR);
        }
    }

    std::vector<std::vector<double*> > mptrs(ngroups), dptrs(ngroups);
    for (int g = 0; g < ngroups; ++g) {
        mptrs[g] = vector_to_pointers<double>(means[g]);
        dptrs[g] = vector_to_pointers<double>(detected[g]);
    }

    // Setting up all the damn output vectors (II) - for log-FC & delta-detected
    std::vector<std::vector<Rcpp::NumericVector> > logfc(ngroups), delta_detected(ngroups);
    std::vector<std::vector<double*> > lptrs(scran::differential_analysis::n_summaries);

    lptrs[scran::differential_analysis::MEAN].resize(ngroups);
    if (!simple_means_only) {
        lptrs[scran::differential_analysis::MIN].resize(ngroups);
        lptrs[scran::differential_analysis::MIN_RANK].resize(ngroups);
    }
    auto ddptrs = lptrs;

    for (int g = 0; g < ngroups; ++g) {
        if (simple_means_only) {
            logfc[g].emplace_back(NR);
            lptrs[scran::differential_analysis::MEAN][g] = static_cast<double*>(logfc[g].back().begin());
            delta_detected[g].emplace_back(NR);
            ddptrs[scran::differential_analysis::MEAN][g] = static_cast<double*>(delta_detected[g].back().begin());
        } else {
            logfc[g].reserve(3);
            delta_detected[g].reserve(3);
            for (int i = 0; i < 3; ++i) {
                logfc[g].emplace_back(NR);
                delta_detected[g].emplace_back(NR);
            }

            lptrs[scran::differential_analysis::MIN][g] = static_cast<double*>(logfc[g][0].begin());
            lptrs[scran::differential_analysis::MEAN][g] = static_cast<double*>(logfc[g][1].begin());
            lptrs[scran::differential_analysis::MIN_RANK][g] = static_cast<double*>(logfc[g][2].begin());

            ddptrs[scran::differential_analysis::MIN][g] = static_cast<double*>(delta_detected[g][0].begin());
            ddptrs[scran::differential_analysis::MEAN][g] = static_cast<double*>(delta_detected[g][1].begin());
            ddptrs[scran::differential_analysis::MIN_RANK][g] = static_cast<double*>(delta_detected[g][2].begin());
        }
    }

    // Setting up all the damn output vectors (III) - Cohen & AUC
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
            cohens[g].emplace_back(NR);
            aucs[g].emplace_back(NR);
        }

        cptrs[scran::differential_analysis::MIN][g] = static_cast<double*>(cohens[g][0].begin());
        cptrs[scran::differential_analysis::MEAN][g] = static_cast<double*>(cohens[g][1].begin());
        cptrs[scran::differential_analysis::MIN_RANK][g] = static_cast<double*>(cohens[g][2].begin());
        
        aptrs[scran::differential_analysis::MIN][g] = static_cast<double*>(aucs[g][0].begin());
        aptrs[scran::differential_analysis::MEAN][g] = static_cast<double*>(aucs[g][1].begin());
        aptrs[scran::differential_analysis::MIN_RANK][g] = static_cast<double*>(aucs[g][2].begin());
    }

    // Running the marker scoring.
    scran::ScoreMarkers runner;
    runner.set_num_threads(nthreads);
    runner.set_threshold(lfc);
    runner.run_blocked(mat, static_cast<const int*>(groups.begin()), bptr, std::move(mptrs), std::move(dptrs), std::move(cptrs), std::move(aptrs), std::move(lptrs), std::move(ddptrs));

    // Organizing the output.
    Rcpp::List output(ngroups), raw_means(ngroups), raw_detected(ngroups);
    for (int g = 0; g < ngroups; ++g) {
        Rcpp::NumericVector outmean, outdet;
        if (bptr == NULL) {
            outmean = means[g][0];
            outdet = detected[g][0];
        } else {
            outmean = Rcpp::NumericVector(NR); 
            scran::average_vectors(NR, vector_to_pointers<double>(means[g]), static_cast<double*>(outmean.begin()));
            raw_means[g] = Rcpp::List(means[g].begin(), means[g].end());

            outdet = Rcpp::NumericVector(NR);
            scran::average_vectors(NR, vector_to_pointers<double>(detected[g]), static_cast<double*>(outdet.begin()));
            raw_detected[g] = Rcpp::List(detected[g].begin(), detected[g].end());
        }

        if (simple_means_only) {
            output[g] = Rcpp::DataFrame::create(
                Rcpp::Named("mean") = outmean,
                Rcpp::Named("detected") = outdet,
                Rcpp::Named("logFC") = logfc[g][0],
                Rcpp::Named("delta.detected") = delta_detected[g][0],
                Rcpp::Named("cohen.min") = cohens[g][0],
                Rcpp::Named("cohen.mean") = cohens[g][1],
                Rcpp::Named("cohen.rank") = cohens[g][2],
                Rcpp::Named("auc.min") = aucs[g][0],
                Rcpp::Named("auc.mean") = aucs[g][1],
                Rcpp::Named("auc.rank") = aucs[g][2]
            );
        } else {
            output[g] = Rcpp::DataFrame::create(
                Rcpp::Named("mean") = outmean,
                Rcpp::Named("detected") = outdet,
                Rcpp::Named("logFC.min") = logfc[g][0],
                Rcpp::Named("logFC.mean") = logfc[g][1],
                Rcpp::Named("logFC.rank") = logfc[g][2],
                Rcpp::Named("delta.detected.min") = delta_detected[g][0],
                Rcpp::Named("delta.detected.mean") = delta_detected[g][1],
                Rcpp::Named("delta.detected.rank") = delta_detected[g][2],
                Rcpp::Named("cohen.min") = cohens[g][0],
                Rcpp::Named("cohen.mean") = cohens[g][1],
                Rcpp::Named("cohen.rank") = cohens[g][2],
                Rcpp::Named("auc.min") = aucs[g][0],
                Rcpp::Named("auc.mean") = aucs[g][1],
                Rcpp::Named("auc.rank") = aucs[g][2]
            );
        }
    }

    if (bptr == NULL) {
        return Rcpp::List::create(Rcpp::Named("statistics") = output);

    } else {
        return Rcpp::List::create(
            Rcpp::Named("statistics") = output,
            Rcpp::Named("per.batch") = Rcpp::List::create(
                Rcpp::Named("mean") = raw_means,
                Rcpp::Named("detected") = raw_detected
            )
        );
    }
}
