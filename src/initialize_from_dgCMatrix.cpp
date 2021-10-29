#include "tatami/tatami.h"
#include "Rcpp.h"
#include "tatamize.h"

#include <cstdint>
#include <limits>
#include <type_traits>

template<class V>
struct RcppVectorPlus {
    RcppVectorPlus(V in) : vec(std::move(in)) {}
    typedef typename std::remove_reference<decltype(V()[0])>::type T;

    const T* data() const {
        return static_cast<const T*>(vec.begin());
    }

    auto begin() const {
        return vec.begin();
    }

    auto end() const {
        return vec.end();
    }

    size_t size() const {
        return vec.size();
    }

    T operator[](size_t i) const {
        return vec[i];
    }
private:
    V vec;
};

template<bool byrow>
SEXP initialize_from_Matrix(Rcpp::NumericVector x, Rcpp::IntegerVector i, Rcpp::IntegerVector p, int nrow, int ncol, bool no_copy, bool force_integer) {
    if (no_copy) {
        RcppVectorPlus x_(x);
        RcppVectorPlus i_(i);
        RcppVectorPlus p_(p);
        typedef tatami::CompressedSparseMatrix<byrow, double, int, decltype(x_), decltype(i_), decltype(p_)> SparseMat;
        return new_MatrixChan(new SparseMat(nrow, ncol, std::move(x_), std::move(i_), std::move(p_), false));
    } else {
        std::vector<int> i_(i.begin(), i.end());
        std::vector<int> p_(p.begin(), p.end());

        if (force_integer) {
            auto maxed = (x.size() ? *std::max_element(x.begin(), x.end()) : 0);
            auto mined = (x.size() ? *std::min_element(x.begin(), x.end()) : 0);

            if (mined < 0) {
                throw std::runtime_error("expression values should be positive");

            } else if (maxed <= std::numeric_limits<uint8_t>::max()) {
                std::vector<uint8_t> x_(x.begin(), x.end());
                typedef tatami::CompressedSparseMatrix<byrow, double, int, decltype(x_), decltype(i_), decltype(p_)> SparseMat;
                return new_MatrixChan(new SparseMat(nrow, ncol, std::move(x_), std::move(i_), std::move(p_)));

            } else if (maxed <= std::numeric_limits<uint16_t>::max()) {
                std::vector<uint16_t> x_(x.begin(), x.end());
                typedef tatami::CompressedSparseMatrix<byrow, double, int, decltype(x_), decltype(i_), decltype(p_)> SparseMat;
                return new_MatrixChan(new SparseMat(nrow, ncol, std::move(x_), std::move(i_), std::move(p_)));

            } else {
                std::vector<int> x_(x.begin(), x.end());
                typedef tatami::CompressedSparseMatrix<byrow, double, int, decltype(x_), decltype(i_), decltype(p_)> SparseMat;
                return new_MatrixChan(new SparseMat(nrow, ncol, std::move(x_), std::move(i_), std::move(p_)));
            }

        } else {
            std::vector<double> x_(x.begin(), x.end());
            typedef tatami::CompressedSparseMatrix<byrow, double, int, decltype(x_), decltype(i_), decltype(p_)> SparseMat;
            return new_MatrixChan(new SparseMat(nrow, ncol, std::move(x_), std::move(i_), std::move(p_), false));
        }
    }
}

//[[Rcpp::export(rng=false)]]
SEXP initialize_from_dgCMatrix(Rcpp::NumericVector x, Rcpp::IntegerVector i, Rcpp::IntegerVector p, int nrow, int ncol, bool no_copy, bool force_integer) {
    return initialize_from_Matrix<false>(x, i, p, nrow, ncol, no_copy, force_integer);
}

//[[Rcpp::export(rng=false)]]
SEXP initialize_from_dgRMatrix(Rcpp::NumericVector x, Rcpp::IntegerVector i, Rcpp::IntegerVector p, int nrow, int ncol, bool no_copy, bool force_integer) {
    return initialize_from_Matrix<true>(x, i, p, nrow, ncol, no_copy, force_integer);
}
