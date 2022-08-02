#include "config.h"

#include "tatami/tatami.h"
#include "Rcpp.h"
#include "tatamize.h"

#include <cstdint>
#include <limits>
#include <type_traits>

template<bool byrow, class Incoming, class RowType>
SEXP create_matrix_copy_x_max(const Incoming& x, std::vector<RowType> i, std::vector<size_t> p, int nrow, int ncol) {
    auto maxed = (x.size() ? *std::max_element(x.begin(), x.end()) : 0);
    auto mined = (x.size() ? *std::min_element(x.begin(), x.end()) : 0);

    if (mined < 0) {
        throw std::runtime_error("expression values should be positive");

    } else if (maxed <= std::numeric_limits<uint8_t>::max()) {
        std::vector<uint8_t> x_(x.begin(), x.end());
        typedef tatami::CompressedSparseMatrix<byrow, double, int, decltype(x_), decltype(i), decltype(p)> SparseMat;
        return new_MatrixChan(new SparseMat(nrow, ncol, std::move(x_), std::move(i), std::move(p)));

    } else if (maxed <= std::numeric_limits<uint16_t>::max()) {
        std::vector<uint16_t> x_(x.begin(), x.end());
        typedef tatami::CompressedSparseMatrix<byrow, double, int, decltype(x_), decltype(i), decltype(p)> SparseMat;
        return new_MatrixChan(new SparseMat(nrow, ncol, std::move(x_), std::move(i), std::move(p)));

    } else {
        std::vector<int> x_(x.begin(), x.end());
        typedef tatami::CompressedSparseMatrix<byrow, double, int, decltype(x_), decltype(i), decltype(p)> SparseMat;
        return new_MatrixChan(new SparseMat(nrow, ncol, std::move(x_), std::move(i), std::move(p)));
    }
}

template<bool byrow, class RowType>
SEXP create_matrix_copy_x_type(const Rcpp::RObject& x, std::vector<RowType> i, std::vector<size_t> p, int nrow, int ncol, bool forced) {
    if (x.sexp_type() == INTSXP) {
        Rcpp::IntegerVector x_(x);
        return create_matrix_copy_x_max<byrow>(x_, std::move(i), std::move(p), nrow, ncol);
    } else {
        Rcpp::NumericVector x_(x);
        if (forced) {
            return create_matrix_copy_x_max<byrow>(x_, std::move(i), std::move(p), nrow, ncol);
        } else {
            std::vector<double> xcopy(x_.begin(), x_.end());
            typedef tatami::CompressedSparseMatrix<byrow, double, int, decltype(xcopy), decltype(i), decltype(p)> SparseMat;
            return new_MatrixChan(new SparseMat(nrow, ncol, std::move(xcopy), std::move(i), std::move(p)));
        }
    }
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

template<bool byrow>
SEXP create_matrix_copy(Rcpp::RObject x, Rcpp::RObject i, Rcpp::RObject p, int nrow, int ncol, bool forced) {
    auto p_ = transfer_p(p);

    if (i.sexp_type() != INTSXP) {
        throw std::runtime_error("'i' vector should be integer");
    }
    Rcpp::IntegerVector i_(i);

    if (nrow <= std::numeric_limits<uint16_t>::max()) {
        return create_matrix_copy_x_type<byrow>(x, std::vector<uint16_t>(i_.begin(), i_.end()), std::move(p_), nrow, ncol, forced);
    } else {
        return create_matrix_copy_x_type<byrow>(x, std::vector<int>(i_.begin(), i_.end()), std::move(p_), nrow, ncol, forced);
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

template<bool byrow, class XVector, class IVector, class PVector>
SEXP create_matrix_nocopy_plus(XVector x, IVector i, PVector p, int nrow, int ncol) {
    RcppVectorPlus<XVector> x_(x);
    RcppVectorPlus<IVector> i_(i);
    typedef tatami::CompressedSparseMatrix<byrow, double, int, decltype(x_), decltype(i_), decltype(p)> SparseMat;
    return new_MatrixChan(new SparseMat(nrow, ncol, std::move(x_), std::move(i_), std::move(p), false));
}

template<bool byrow>
SEXP create_matrix_nocopy(Rcpp::RObject x, Rcpp::RObject i, Rcpp::RObject p, int nrow, int ncol) {
    auto p_ = transfer_p(p);

    if (i.sexp_type() != INTSXP) {
        throw std::runtime_error("'i' vector should be integer");
    }
    Rcpp::IntegerVector i_(i);
    
    if (x.sexp_type() == INTSXP) {
        Rcpp::IntegerVector x_(x);
        return create_matrix_nocopy_plus<byrow>(std::move(x_), std::move(i_), std::move(p_), nrow, ncol);
    } else if (x.sexp_type() != REALSXP) {
        throw std::runtime_error("'x' vector should be integer or real");
    }

    Rcpp::NumericVector x_(x);
    return create_matrix_nocopy_plus<byrow>(std::move(x_), std::move(i_), std::move(p_), nrow, ncol);
}

//[[Rcpp::export(rng=false)]]
SEXP initialize_from_CSC(Rcpp::RObject x, Rcpp::RObject i, Rcpp::RObject p, int nrow, int ncol, bool no_copy, bool forced) {
    if (no_copy) {
        return create_matrix_nocopy<false>(x, i, p, nrow, ncol);
    } else {
        return create_matrix_copy<false>(x, i, p, nrow, ncol, forced);
    }
}

//[[Rcpp::export(rng=false)]]
SEXP initialize_from_CSR(Rcpp::RObject x, Rcpp::RObject i, Rcpp::RObject p, int nrow, int ncol, bool no_copy, bool forced) {
    if (no_copy) {
        return create_matrix_nocopy<true>(x, i, p, nrow, ncol);
    } else {
        return create_matrix_copy<true>(x, i, p, nrow, ncol, forced);
    }
}



