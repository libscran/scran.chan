#include "run_tsne.h"
#include "run_umap.h"
#include "cluster_kmeans.h"
#include "cluster_snn_graph.h"
#ifdef _OPENMP
#include "omp.h"
#endif

//[[Rcpp::export(rng=false)]]
Rcpp::List run_all_downstream(SEXP clust_init, SEXP umap_init, SEXP tsne_init, SEXP kmeans_init, int nthreads) {
    std::vector<int> jobs;

    // Defining all the steps to run in parallel.
    SNNGraph* cptr = NULL;
    if (!Rf_isNull(clust_init)) {
        cptr = Rcpp::XPtr<SNNGraph>(clust_init).get();
        jobs.push_back(0);
    }

    InitializedUmap* uptr = NULL;
    if (!Rf_isNull(umap_init)) {
        uptr = Rcpp::XPtr<InitializedUmap>(umap_init).get();
        jobs.push_back(1);
    }

    InitializedTsne* tptr = NULL;
    if (!Rf_isNull(tsne_init)) {
        tptr = Rcpp::XPtr<InitializedTsne>(tsne_init).get();
        jobs.push_back(2);
    }

    std::unique_ptr<Kmeans> kptr;
    if (!Rf_isNull(kmeans_init)) {
        Rcpp::List klist(kmeans_init);
        kptr.reset(new Kmeans(klist[0], Rcpp::IntegerVector(klist[1])[0]));
        jobs.push_back(3);
    }

    // Running the iterations.
    int outer = std::min(nthreads, static_cast<int>(jobs.size()));
    Rcpp::List output(4);

    #pragma omp parallel for num_threads(outer)
    for (size_t j = 0; j < jobs.size(); ++j) {
        auto i = jobs[j];
        if (i == 0) {
            cptr->run();
            #pragma omp critical 
            output[i] = cptr->yield();
        } else if (i == 1) {
            uptr->run();
            #pragma omp critical 
            output[i] = uptr->yield();
        } else if (i == 2) {
#ifdef _OPENMP
            if (nthreads < 3) {
                // I believe this should be okay based on discussions of nested use of OpenMP routines,
                // see https://docs.oracle.com/cd/E19205-01/819-5270/aewbi/index.html.
                omp_set_num_threads(nthreads);
            } else {
                omp_set_num_threads(nthreads - 2);
            }
#endif            
            tptr->run();
            #pragma omp critical 
            output[i] = tptr->yield();
        } else if (i == 3) {
            kptr->run();
            #pragma omp critical
            output[i] = kptr->yield();
        }
    }

    return output;
}
