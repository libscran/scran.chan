#include "scran/static/bindings.hpp"
#include "scran/quality_control/PerCellQCMetrics.hpp"

void scran::sttc::PerCellQCMetrics(const scran::sttc::Matrix& mat, std::vector<const int*> subsets, double* sums, int* detected, std::vector<double*> subset_proportions) {
    scran::PerCellQCMetrics qc;
    qc.run(mat.get(), std::move(subsets), sums, detected, std::move(subset_proportions));
    return;
}
