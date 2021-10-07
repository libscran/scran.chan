#include "run_tsne.h"
#include "run_umap.h"
#include "cluster_snn_graph.h"
#ifdef _OPENMP
#include "omp.h"
#endif

//[[Rcpp::export(rng=false)]]
Rcpp::List run_all_neighbors(SEXP clust_init, SEXP umap_init, SEXP tsne_init, int nthreads) {
    Rcpp::XPtr<SNNGraph> cptr(clust_init);
    Rcpp::XPtr<InitializedUmap> uptr(umap_init);
    Rcpp::XPtr<InitializedTsne> tptr(tsne_init);

    int outer = std::min(nthreads, 3);
    Rcpp::List output(3);

    #pragma omp parallel for num_threads(outer)
    for (int i = 0; i < 3; ++i) {
        if (i == 0) {
            auto clust = cptr->run();
            #pragma omp critical 
            output[i] = cptr->yield(clust);
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
        }
    }

    return output;
}
