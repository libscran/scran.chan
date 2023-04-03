#include "config.h"

#include "scran/differential_analysis/ScoreMarkers.hpp"
#include "scran/differential_analysis/PairwiseEffects.hpp"
#include "scran/differential_analysis/SummarizeEffects.hpp"
#include "scran/utils/average_vectors.hpp"

#include "Rcpp.h"
#include "tatamize.h"

#include "vector_to_pointers.h"
#include "ResolvedBatch.h"

struct MarkerOutputs {
    MarkerOutputs(int ngr, int nb, size_t nge, bool sm) :
        ngroups(ngr),
        nbatches(nb),
        ngenes(nge),
        simple_means_only(sm),

        means(ngroups), 
        detected(ngroups),
        mptrs(ngroups), 
        dptrs(ngroups),

        logfc(ngroups), 
        delta_detected(ngroups),
        lptrs(scran::differential_analysis::n_summaries), // ddptrs is filled in body.

        cohens(ngroups), 
        aucs(ngroups),
        cptrs(scran::differential_analysis::n_summaries), // aptrs is filled in body.

        output_stats(ngroups),
        output_raw_means(ngroups),
        output_raw_detected(ngroups)
    {
        // Setting up all the damn output vectors (I) - for mean & detected
        for (int g = 0; g < ngroups; ++g) {
            means[g].reserve(nbatches);
            detected[g].reserve(nbatches);
            for (size_t b = 0; b < nbatches; ++b) {
                means[g].emplace_back(ngenes);
                detected[g].emplace_back(ngenes);
            }
        }

        for (int g = 0; g < ngroups; ++g) {
            mptrs[g] = vector_to_pointers<double>(means[g]);
            dptrs[g] = vector_to_pointers<double>(detected[g]);
        }

        // Setting up all the damn output vectors (II) - for log-FC & delta-detected
        lptrs[scran::differential_analysis::MEAN].resize(ngroups);
        if (!simple_means_only) {
            lptrs[scran::differential_analysis::MIN].resize(ngroups);
            lptrs[scran::differential_analysis::MIN_RANK].resize(ngroups);
        }
        ddptrs = lptrs;

        for (int g = 0; g < ngroups; ++g) {
            if (simple_means_only) {
                logfc[g].emplace_back(ngenes);
                lptrs[scran::differential_analysis::MEAN][g] = static_cast<double*>(logfc[g].back().begin());
                delta_detected[g].emplace_back(ngenes);
                ddptrs[scran::differential_analysis::MEAN][g] = static_cast<double*>(delta_detected[g].back().begin());
            } else {
                logfc[g].reserve(3);
                delta_detected[g].reserve(3);
                for (int i = 0; i < 3; ++i) {
                    logfc[g].emplace_back(ngenes);
                    delta_detected[g].emplace_back(ngenes);
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
        cptrs[scran::differential_analysis::MIN].resize(ngroups);
        cptrs[scran::differential_analysis::MEAN].resize(ngroups);
        cptrs[scran::differential_analysis::MIN_RANK].resize(ngroups);
        aptrs = cptrs;

        for (int g = 0; g < ngroups; ++g) {
            cohens[g].reserve(3);
            aucs[g].reserve(3);
            for (int i = 0; i < 3; ++i) {
                cohens[g].emplace_back(ngenes);
                aucs[g].emplace_back(ngenes);
            }

            cptrs[scran::differential_analysis::MIN][g] = static_cast<double*>(cohens[g][0].begin());
            cptrs[scran::differential_analysis::MEAN][g] = static_cast<double*>(cohens[g][1].begin());
            cptrs[scran::differential_analysis::MIN_RANK][g] = static_cast<double*>(cohens[g][2].begin());

            aptrs[scran::differential_analysis::MIN][g] = static_cast<double*>(aucs[g][0].begin());
            aptrs[scran::differential_analysis::MEAN][g] = static_cast<double*>(aucs[g][1].begin());
            aptrs[scran::differential_analysis::MIN_RANK][g] = static_cast<double*>(aucs[g][2].begin());
        }
    }

    void format_output(bool report_batches) {
        for (int g = 0; g < ngroups; ++g) {
            Rcpp::NumericVector outmean, outdet;
            if (!report_batches) {
                outmean = means[g][0];
                outdet = detected[g][0];
            } else {
                outmean = Rcpp::NumericVector(ngenes); 
                scran::average_vectors(ngenes, vector_to_pointers<double>(means[g]), static_cast<double*>(outmean.begin()));
                output_raw_means[g] = Rcpp::List(means[g].begin(), means[g].end());

                outdet = Rcpp::NumericVector(ngenes);
                scran::average_vectors(ngenes, vector_to_pointers<double>(detected[g]), static_cast<double*>(outdet.begin()));
                output_raw_detected[g] = Rcpp::List(detected[g].begin(), detected[g].end());
            }

            if (simple_means_only) {
                output_stats[g] = Rcpp::DataFrame::create(
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
                output_stats[g] = Rcpp::DataFrame::create(
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
    }

public:
    int ngroups;
    int nbatches;
    size_t ngenes;
    bool simple_means_only;

    std::vector<std::vector<Rcpp::NumericVector> > means, detected;
    std::vector<std::vector<double*> > mptrs, dptrs;

    std::vector<std::vector<Rcpp::NumericVector> > logfc, delta_detected, cohens, aucs;
    std::vector<std::vector<double*> > lptrs, ddptrs, cptrs, aptrs; 

    Rcpp::List output_stats;
    Rcpp::List output_raw_means;
    Rcpp::List output_raw_detected;
};

//[[Rcpp::export(rng=false)]]
Rcpp::List score_markers(SEXP x, Rcpp::IntegerVector groups, Rcpp::Nullable<Rcpp::IntegerVector> batch, bool simple_means_only, double lfc, int nthreads) {
    auto mat = extract_NumericMatrix(x);
    size_t NR = mat->nrow();
    const int ngroups = (groups.size() ? *std::max_element(groups.begin(), groups.end()) + 1 : 1);

    auto batch_info = ResolvedBatch(batch);
    auto bptr = batch_info.ptr;
    size_t nbatches = batch_info.number();

    MarkerOutputs store(ngroups, nbatches, NR, simple_means_only);

    // Running the marker scoring.
    scran::ScoreMarkers runner;
    runner.set_num_threads(nthreads);
    runner.set_threshold(lfc);
    runner.run_blocked(mat, 
        static_cast<const int*>(groups.begin()), 
        bptr, 
        std::move(store.mptrs), 
        std::move(store.dptrs), 
        std::move(store.cptrs), 
        std::move(store.aptrs), 
        std::move(store.lptrs), 
        std::move(store.ddptrs)
    );

    // Organizing the output.
    store.format_output(bptr != NULL);
    if (bptr == NULL) {
        return Rcpp::List::create(Rcpp::Named("statistics") = store.output_stats);

    } else {
        return Rcpp::List::create(
            Rcpp::Named("statistics") = store.output_stats,
            Rcpp::Named("per.batch") = Rcpp::List::create(
                Rcpp::Named("mean") = store.output_raw_means,
                Rcpp::Named("detected") = store.output_raw_detected
            )
        );
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::List score_markers_full(SEXP x, Rcpp::IntegerVector groups, Rcpp::Nullable<Rcpp::IntegerVector> batch, bool simple_means_only, double lfc, int nthreads) {
    auto mat = extract_NumericMatrix(x);
    size_t NR = mat->nrow();
    const int ngroups = (groups.size() ? *std::max_element(groups.begin(), groups.end()) + 1 : 1);

    auto batch_info = ResolvedBatch(batch);
    auto bptr = batch_info.ptr;
    size_t nbatches = batch_info.number();

    MarkerOutputs store(ngroups, nbatches, NR, simple_means_only);

    size_t total_size = ngroups * ngroups * NR;
    std::vector<double> pairwise_aucs(total_size), pairwise_cohens(total_size), pairwise_lfcs(total_size), pairwise_deltas(total_size);

    // Computing pairwise effects.
    scran::PairwiseEffects runner;
    runner.set_num_threads(nthreads);
    runner.set_threshold(lfc);
    runner.run_blocked(mat, 
        static_cast<const int*>(groups.begin()), 
        bptr, 
        std::move(store.mptrs), 
        std::move(store.dptrs), 
        pairwise_cohens.data(),
        pairwise_aucs.data(),
        pairwise_lfcs.data(),
        pairwise_deltas.data()
    );

    // Shifting effects into their Rcpp::NumericVectors once summarization is done.
    // This avoids creating two copies of all effects at the same time.
    scran::SummarizeEffects summer;
    summer.set_num_threads(nthreads);

    auto shift = [](const scran::SummarizeEffects& summer, int ngroups, size_t ngenes, std::vector<double>& pairwise) -> auto {
        std::vector<std::vector<Rcpp::NumericVector> > output(ngroups);
        for (int g = 0; g < ngroups; ++g) {
            output[g].resize(ngroups);
            for (int g2 = 0; g2 < ngroups; ++g2) {
                if (g == g2) {
                    continue;
                }

                Rcpp::NumericVector vals(ngenes);
                size_t offset = g * ngroups + g2;
                size_t shift = ngroups * ngroups;
                for (size_t r = 0; r < ngenes; ++r) {
                    vals[r] = pairwise[offset];
                    offset += shift;
                }

                output[g][g2] = vals;
            }
        }

        pairwise.clear();
        pairwise.shrink_to_fit(); // release memory.
        return output;
    };

    summer.run(NR, ngroups, pairwise_cohens.data(), std::move(store.cptrs));
    auto rcpp_cohen = shift(summer, ngroups, NR, pairwise_cohens);

    summer.run(NR, ngroups, pairwise_aucs.data(), std::move(store.aptrs));
    auto rcpp_auc = shift(summer, ngroups, NR, pairwise_aucs);

    summer.run(NR, ngroups, pairwise_lfcs.data(), std::move(store.lptrs));
    auto rcpp_lfc = shift(summer, ngroups, NR, pairwise_lfcs);

    summer.run(NR, ngroups, pairwise_deltas.data(), std::move(store.ddptrs));
    auto rcpp_delta = shift(summer, ngroups, NR, pairwise_deltas);

    // Organizing the output.
    Rcpp::List paired(ngroups);
    for (int g = 0; g < ngroups; ++g) {
        Rcpp::List subpaired(ngroups);
        for (int g2 = 0; g2 < ngroups; ++g2) {
            if (g == g2) {
                continue;
            }

            subpaired[g2] = Rcpp::DataFrame::create(
                Rcpp::Named("logFC") = rcpp_lfc[g][g2],
                Rcpp::Named("delta.detected") = rcpp_delta[g][g2],
                Rcpp::Named("cohen") = rcpp_cohen[g][g2],
                Rcpp::Named("auc") = rcpp_auc[g][g2]
            );
        }

        paired[g] = subpaired;
    }

    store.format_output(bptr != NULL);
    if (bptr == NULL) {
        return Rcpp::List::create(
            Rcpp::Named("statistics") = store.output_stats,
            Rcpp::Named("pairwise") = paired
        );

    } else {
        return Rcpp::List::create(
            Rcpp::Named("statistics") = store.output_stats,
            Rcpp::Named("per.batch") = Rcpp::List::create(
                Rcpp::Named("mean") = store.output_raw_means,
                Rcpp::Named("detected") = store.output_raw_detected
            ),
            Rcpp::Named("pairwise") = paired
        );
    }
}

