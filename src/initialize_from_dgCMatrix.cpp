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

//[[Rcpp::export(rng=false)]]
SEXP initialize_from_dgCMatrix(Rcpp::NumericVector x, Rcpp::IntegerVector i, Rcpp::IntegerVector p, int nrow, int ncol) {
    RcppVectorPlus x_(x);
    RcppVectorPlus i_(i);
    RcppVectorPlus p_(p);
    typedef tatami::CompressedSparseColumnMatrix<double, int, decltype(x_), decltype(i_), decltype(p_)> SparseMat;
    return new_MatrixChan(new SparseMat(nrow, ncol, std::move(x_), std::move(i_), std::move(p_), false));
}

//[[Rcpp::export(rng=false)]]
SEXP initialize_from_dgRMatrix(Rcpp::NumericVector x, Rcpp::IntegerVector i, Rcpp::IntegerVector p, int nrow, int ncol) {
    RcppVectorPlus x_(x);
    RcppVectorPlus i_(i);
    RcppVectorPlus p_(p);
    typedef tatami::CompressedSparseRowMatrix<double, int, decltype(x_), decltype(i_), decltype(p_)> SparseMat;
    return new_MatrixChan(new SparseMat(nrow, ncol, std::move(x_), std::move(i_), std::move(p_), false));
}
