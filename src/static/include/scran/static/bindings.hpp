#ifndef SCRAN_STATIC_HPP
#define SCRAN_STATIC_HPP

#include "tatami/tatami.h"
#include <vector>

namespace scran {

namespace sttc {

typedef std::shared_ptr<tatami::NumericMatrix> Matrix;

void PerCellQCMetrics(const Matrix& mat, std::vector<const double*> subsets, double* sums, int* detected, std::vector<double*> subset_proportions);

}

}

#endif
