#include "config.h"

#include "tatami/tatami.h"
#include "Rcpp.h"
#include "tatamize.h"

#include <cstdint>
#include <limits>
#include <type_traits>

template<class XVector, class IVector, class PVector>
SEXP create_matrix_copy_byrow(XVector x, IVector i, PVector p, int nrow, int ncol, bool byrow) {
    if (byrow) {
        typedef tatami::CompressedSparseMatrix<true, double, int, XVector, IVector, PVector> SparseMat;
        return new_MatrixChan(new SparseMat(nrow, ncol, std::move(x), std::move(i), std::move(p), false));
    } else {
        typedef tatami::CompressedSparseMatrix<false, double, int, XVector, IVector, PVector> SparseMat;
        return new_MatrixChan(new SparseMat(nrow, ncol, std::move(x), std::move(i), std::move(p), false));
    }
}

template<class Incoming, class RowType>
SEXP create_matrix_copy_x_max(const Incoming& x, std::vector<RowType> i, std::vector<size_t> p, int nrow, int ncol, bool byrow) {
    auto mined = (x.size() ? *std::min_element(x.begin(), x.end()) : 0);
    if (mined < 0) {
        throw std::runtime_error("expression values should be positive");
    }

    auto maxed = (x.size() ? *std::max_element(x.begin(), x.end()) : 0);
    if (maxed <= std::numeric_limits<uint16_t>::max()) {
        return create_matrix_copy_byrow(std::vector<uint16_t>(x.begin(), x.end()), std::move(i), std::move(p), nrow, ncol, byrow);
    } else {
        return create_matrix_copy_byrow(std::vector<int>(x.begin(), x.end()), std::move(i), std::move(p), nrow, ncol, byrow);
    }
}

template<class RowType>
SEXP create_matrix_copy_x_type(const Rcpp::RObject& x, std::vector<RowType> i, std::vector<size_t> p, int nrow, int ncol, bool byrow, bool forced) {
    if (x.sexp_type() == INTSXP) {
        Rcpp::IntegerVector x_(x);
        return create_matrix_copy_x_max(x_, std::move(i), std::move(p), nrow, ncol, byrow);
    } 

    Rcpp::NumericVector x_(x);
    if (forced) {
        return create_matrix_copy_x_max(x_, std::move(i), std::move(p), nrow, ncol, byrow);
    }

    return create_matrix_copy_byrow(std::vector<double>(x_.begin(), x_.end()), std::move(i), std::move(p), nrow, ncol, byrow);
}

std::vector<size_t> transfer_p(const Rcpp::RObject& p) {
    if (p.sexp_type() == INTSXP) {
        Rcpp::IntegerVector p_(p);
        return std::vector<size_t>(p_.begin(), p_.end());
    } else {
        Rcpp::NumericVector p_(p);
        return std::vector<size_t>(p_.begin(), p_.end());
    }
}

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

template<class XVector, class IVector, class PVector>
SEXP create_matrix_nocopy(XVector x, IVector i, PVector p, int nrow, int ncol, bool byrow) {
    RcppVectorPlus<XVector> x_(x);
    RcppVectorPlus<IVector> i_(i);

    if (byrow) {
        typedef tatami::CompressedSparseMatrix<true, double, int, decltype(x_), decltype(i_), decltype(p)> SparseMat;
        return new_MatrixChan(new SparseMat(nrow, ncol, std::move(x_), std::move(i_), std::move(p), false));
    } else {
        typedef tatami::CompressedSparseMatrix<false, double, int, decltype(x_), decltype(i_), decltype(p)> SparseMat;
        return new_MatrixChan(new SparseMat(nrow, ncol, std::move(x_), std::move(i_), std::move(p), false));
    }
}

//[[Rcpp::export(rng=false)]]
SEXP initialize_from_memory(Rcpp::RObject x, Rcpp::RObject i, Rcpp::RObject p, int nrow, int ncol, bool no_copy, bool byrow, bool forced) {
    auto p_ = transfer_p(p);
    if (i.sexp_type() != INTSXP) {
        throw std::runtime_error("'i' vector should be integer");
    }
    Rcpp::IntegerVector i_(i);

    // If no copy is requested, we hold onto the Rcpp vectors.
    if (no_copy) {
        if (x.sexp_type() == INTSXP) {
            Rcpp::IntegerVector x_(x);
            return create_matrix_nocopy(std::move(x_), std::move(i_), std::move(p_), nrow, ncol, byrow);
        } else if (x.sexp_type() != REALSXP) {
            throw std::runtime_error("'x' vector should be integer or real");
        }

        Rcpp::NumericVector x_(x);
        return create_matrix_nocopy(std::move(x_), std::move(i_), std::move(p_), nrow, ncol, byrow);
    }

    // Otherwise, beginning the copying process with the indices.
    if (nrow <= std::numeric_limits<uint16_t>::max()) {
        return create_matrix_copy_x_type(x, std::vector<uint16_t>(i_.begin(), i_.end()), std::move(p_), nrow, ncol, byrow, forced);
    } else {
        return create_matrix_copy_x_type(x, std::vector<int>(i_.begin(), i_.end()), std::move(p_), nrow, ncol, byrow, forced);
    }
}
