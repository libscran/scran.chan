#ifndef DIET_SCRAN_TATAMIZE_H
#define DIET_SCRAN_TATAMIZE_H

#include "tatami/tatami.hpp"
#include <memory>
#include "Rcpp.h"

std::shared_ptr<tatami::Matrix<double, int> > tatamize(Rcpp::RObject);

#endif
