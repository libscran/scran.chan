#include "tatami/tatami.h"
#include "Rcpp.h"
#include "tatamize.h"

#include <cstdint>
#include <limits>
#include <type_traits>

template<bool byrow, class Incoming>
SEXP create_integer_matrix(const Incoming& x0, std::vector<int> i_, std::vector<int> p_, int nrow, int ncol) {
    auto maxed = (x0.size() ? *std::max_element(x0.begin(), x0.end()) : 0);
    auto mined = (x0.size() ? *std::min_element(x0.begin(), x0.end()) : 0);

    if (mined < 0) {
        throw std::runtime_error("expression values should be positive");

    } else if (maxed <= std::numeric_limits<uint8_t>::max()) {
        std::vector<uint8_t> x_(x0.begin(), x0.end());
        typedef tatami::CompressedSparseMatrix<byrow, double, int, decltype(x_), decltype(i_), decltype(p_)> SparseMat;
        return new_MatrixChan(new SparseMat(nrow, ncol, std::move(x_), std::move(i_), std::move(p_)));

    } else if (maxed <= std::numeric_limits<uint16_t>::max()) {
        std::vector<uint16_t> x_(x0.begin(), x0.end());
        typedef tatami::CompressedSparseMatrix<byrow, double, int, decltype(x_), decltype(i_), decltype(p_)> SparseMat;
        return new_MatrixChan(new SparseMat(nrow, ncol, std::move(x_), std::move(i_), std::move(p_)));

    } else {
        std::vector<int> x_(x0.begin(), x0.end());
        typedef tatami::CompressedSparseMatrix<byrow, double, int, decltype(x_), decltype(i_), decltype(p_)> SparseMat;
        return new_MatrixChan(new SparseMat(nrow, ncol, std::move(x_), std::move(i_), std::move(p_)));
    }
}

template<bool byrow>
SEXP initialize_with_copy(Rcpp::RObject x, Rcpp::IntegerVector i, Rcpp::IntegerVector p, int nrow, int ncol, bool forced) {
    std::vector<int> i_(i.begin(), i.end());
    std::vector<int> p_(p.begin(), p.end());

    if (x.sexp_type() == INTSXP) {
        Rcpp::IntegerVector x0(x);
        return create_integer_matrix<byrow>(x0, std::move(i_), std::move(p_), nrow, ncol);

    } else if (x.sexp_type() == REALSXP) {
        Rcpp::NumericVector x0(x);

        if (forced) {
            return create_integer_matrix<byrow>(x0, std::move(i_), std::move(p_), nrow, ncol);

        } else {
            std::vector<double> x_(x0.begin(), x0.end());
            typedef tatami::CompressedSparseMatrix<byrow, double, int, decltype(x_), decltype(i_), decltype(p_)> SparseMat;
            return new_MatrixChan(new SparseMat(nrow, ncol, std::move(x_), std::move(i_), std::move(p_)));
        }
    } else {
        throw std::runtime_error("unrecognized type for 'x'");
    }
}

//[[Rcpp::export(rng=false)]]
SEXP initialize_from_CSC(Rcpp::RObject x, Rcpp::IntegerVector i, Rcpp::IntegerVector p, int nrow, int ncol, bool forced) {
    return initialize_with_copy<false>(x, i, p, nrow, ncol, forced);
}

//[[Rcpp::export(rng=false)]]
SEXP initialize_from_CSR(Rcpp::RObject x, Rcpp::IntegerVector i, Rcpp::IntegerVector p, int nrow, int ncol, bool forced) {
    return initialize_with_copy<true>(x, i, p, nrow, ncol, forced);
}



