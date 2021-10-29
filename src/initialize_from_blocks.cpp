#include "tatami/tatami.h"
#include "Rcpp.h"
#include "tatamize.h"

#include <cstdint>
#include <limits>
#include <type_traits>
#include <deque>
#include <vector>

template<typename I, typename V>
struct accumulated {
    accumulated(int nr, int nc, bool byrow) : nrow(nr), ncol(nc), offsets((byrow ? nr :nc) + 1) {}

    size_t nrow, ncol;

    std::deque<I> indices;
    std::deque<V> values;
    std::vector<size_t> offsets;

    typedef I index_type;
    typedef V value_type;

    int so_far = 0;
};

typedef accumulated<uint16_t, int> short_int; 
typedef accumulated<uint32_t, int> regular_int; 
typedef accumulated<uint16_t, double> short_dbl; 
typedef accumulated<uint32_t, double> regular_dbl; 

Rcpp::List initialize_from_blocks(int nr, int nc, bool is_integer, bool byrow) {
    if (nc <= std::numeric_limits<uint16_t>::max()) {
        if (is_integer) {
            return Rcpp::List::create(
                Rcpp::XPtr<short_int>(new short_int(nr, nc, byrow)),
                Rcpp::IntegerVector::create(0)
            );
        } else {
            return Rcpp::List::create(
                Rcpp::XPtr<short_dbl>(new short_dbl(nr, nc, byrow)),
                Rcpp::IntegerVector::create(1)
            );
        }
    } else {
        if (is_integer) {
            return Rcpp::List::create(
                Rcpp::XPtr<regular_int>(new regular_int(nr, nc, byrow)),
                Rcpp::IntegerVector::create(2)
            );
        } else {
            return Rcpp::List::create(
                Rcpp::XPtr<regular_dbl>(new regular_dbl(nr, nc, byrow)),
                Rcpp::IntegerVector::create(3)
            );
        }
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::List initialize_from_blocks_CSC(int nr, int nc, bool is_integer) {
    return initialize_from_blocks(nr, nc, is_integer, false);    
}

//[[Rcpp::export(rng=false)]]
Rcpp::List initialize_from_blocks_CSR(int nr, int nc, bool is_integer) {
    return initialize_from_blocks(nr, nc, is_integer, true);
}

int fetch_id(Rcpp::IntegerVector x) {
    return x[0];
}

template<class Accumulated, class Vector>
SEXP add_new_block_internal(Rcpp::List ptr0, Rcpp::IntegerVector primary, Rcpp::IntegerVector secondary, Vector values, int nprimary) {
    Rcpp::XPtr<Accumulated> ptr(ptr0[0]);

    ptr->values.insert(ptr->values.end(), values.begin(), values.end());
    ptr->indices.insert(ptr->indices.end(), secondary.begin(), secondary.end());

    int& adjustment = ptr->so_far;
    for (auto x : primary) {
        ++(ptr->offsets[x + adjustment + 1]); // +1 because we want one-after-the-index.
    }

    adjustment += nprimary;
    return R_NilValue;
}

SEXP add_new_block(Rcpp::List ptr0, Rcpp::IntegerVector primary, Rcpp::IntegerVector secondary, SEXP values, int nprimary) {
    int y = fetch_id(ptr0[1]);
    if (y == 0) {
        return add_new_block_internal<short_int, Rcpp::IntegerVector>(ptr0, primary, secondary, values, nprimary);
    } else if (y == 1) {
        return add_new_block_internal<short_dbl, Rcpp::NumericVector>(ptr0, primary, secondary, values, nprimary);
    } else if (y == 2) {
        return add_new_block_internal<regular_int, Rcpp::IntegerVector>(ptr0, primary, secondary, values, nprimary);
    } else if (y == 3) {
        return add_new_block_internal<regular_dbl, Rcpp::NumericVector>(ptr0, primary, secondary, values, nprimary);
    }
}

//[[Rcpp::export(rng=false)]]
SEXP add_new_block_CSC(SEXP ptr0, Rcpp::IntegerVector rows, Rcpp::IntegerVector columns, SEXP values, int ncolumns) {
    return add_new_block(ptr0, columns, rows, values, ncolumns);
}

//[[Rcpp::export(rng=false)]]
SEXP add_new_block_CSR(SEXP ptr0, Rcpp::IntegerVector rows, Rcpp::IntegerVector columns, SEXP values, int nrows) {
    return add_new_block(ptr0, rows, columns, values, nrows);
}

template<bool byrow, class Accumulated>
SEXP finalize_all_blocks_internal(Rcpp::List ptr0) {
    Rcpp::XPtr<Accumulated> ptr(ptr0[0]);

    std::vector<typename Accumulated::index_type> indices(ptr->indices.begin(), ptr->indices.end());
    std::vector<typename Accumulated::value_type> values(ptr->values.begin(), ptr->values.end());

    for (size_t i = 1; i < ptr->offsets.size(); ++i) {
        ptr->offsets[i] += ptr->offsets[i-1];
    }

    typedef tatami::CompressedSparseMatrix<byrow, double, int, decltype(values), decltype(indices), decltype(ptr->offsets)> SparseMat;
    return new_MatrixChan(new SparseMat(ptr->nrow, ptr->ncol, std::move(values), std::move(indices), std::move(ptr->offsets)));
}

template<bool byrow>
SEXP finalize_all_blocks(Rcpp::List ptr0) {
    int y = fetch_id(ptr0[1]);
    if (y == 0) {
        return finalize_all_blocks_internal<byrow, short_int>(ptr0);
    } else if (y == 1) {
        return finalize_all_blocks_internal<byrow, short_dbl>(ptr0);
    } else if (y == 2) {
        return finalize_all_blocks_internal<byrow, regular_int>(ptr0);
    } else if (y == 3) {
        return finalize_all_blocks_internal<byrow, regular_dbl>(ptr0);
    }
}

//[[Rcpp::export(rng=false)]]
SEXP finalize_all_blocks_CSC(SEXP ptr0) {
    return finalize_all_blocks<false>(ptr0);
}

//[[Rcpp::export(rng=false)]]
SEXP finalize_all_blocks_CSR(SEXP ptr0) {
    return finalize_all_blocks<true>(ptr0);
}
