#include "tatami/tatami.h"
#include "Rcpp.h"
#include "tatamize.h"

#include <cstdint>
#include <limits>
#include <type_traits>
#include <deque>
#include <vector>

template<typename I, typename V>
struct accumulated_rows {
    accumulated_rows(int nr, int nc) : nrow(nr), ncol(nc), offsets(nr + 1) {}

    size_t nrow, ncol;

    std::deque<I> columns;
    std::deque<V> values;
    std::vector<size_t> offsets;

    typedef I index_type;
    typedef V value_type;

    int so_far = 0;
};

typedef accumulated_rows<uint16_t, int> short_int; 
typedef accumulated_rows<uint32_t, int> regular_int; 
typedef accumulated_rows<uint16_t, double> short_dbl; 
typedef accumulated_rows<uint32_t, double> regular_dbl; 

//[[Rcpp::export(rng=false)]]
SEXP initialize_from_blocks(int nr, int nc, bool is_integer) {
    if (nc <= std::numeric_limits<uint16_t>::max()) {
        if (is_integer) {
            return Rcpp::XPtr<short_int>(new short_int(nr, nc));
        } else {
            return Rcpp::XPtr<short_dbl>(new short_dbl(nr, nc));
        }
    } else {
        if (is_integer) {
            return Rcpp::XPtr<regular_int>(new regular_int(nr, nc));
        } else {
            return Rcpp::XPtr<regular_dbl>(new regular_dbl(nr, nc));
        }
    }
}

template<class Acc, class Vec>
SEXP add_new_block_internal(SEXP ptr0, Rcpp::IntegerVector rows, Rcpp::IntegerVector columns, Vec values, int nr_block) {
    Rcpp::XPtr<Acc> ptr(ptr0);

    ptr->values.insert(ptr->values.end(), values.begin(), values.end());
    ptr->columns.insert(ptr->columns.end(), columns.begin(), columns.end());

    int& adjustment = ptr->so_far;
    for (auto r : rows) {
        ++(ptr->offsets[r + adjustment + 1]); // +1 because we want one-after-the-index.
    }

    adjustment += nr_block;
    return R_NilValue;
}

//[[Rcpp::export(rng=false)]]
SEXP add_new_block(SEXP ptr0, Rcpp::IntegerVector rows, Rcpp::IntegerVector columns, SEXP values, int nr_block, int nc, bool is_integer) {
    if (nc <= std::numeric_limits<uint16_t>::max()) {
        if (is_integer) {
            add_new_block_internal<short_int, Rcpp::IntegerVector>(ptr0, rows, columns, values, nr_block);
        } else {
            add_new_block_internal<short_dbl, Rcpp::NumericVector>(ptr0, rows, columns, values, nr_block);
        }
    } else {
        if (is_integer) {
            add_new_block_internal<regular_int, Rcpp::IntegerVector>(ptr0, rows, columns, values, nr_block);
        } else {
            add_new_block_internal<regular_dbl, Rcpp::NumericVector>(ptr0, rows, columns, values, nr_block);
        }
    }
    return R_NilValue;
}

template<class Acc>
SEXP finalize_all_blocks_internal(SEXP ptr0) {
    Rcpp::XPtr<Acc> ptr(ptr0);

    std::vector<typename Acc::index_type> indices(ptr->columns.begin(), ptr->columns.end());
    std::vector<typename Acc::value_type> values(ptr->values.begin(), ptr->values.end());

    for (size_t i = 1; i <= ptr->offsets.size(); ++i) {
        ptr->offsets[i] += ptr->offsets[i-1];
    }

    typedef tatami::CompressedSparseRowMatrix<double, int, decltype(values), decltype(indices), decltype(ptr->offsets)> SparseMat;
    return new_MatrixChan(new SparseMat(ptr->nrow, ptr->ncol, std::move(values), std::move(indices), std::move(ptr->offsets)));
}

//[[Rcpp::export(rng=false)]]
SEXP finalize_all_blocks(SEXP ptr0, int nc, bool is_integer) {
    if (nc <= std::numeric_limits<uint16_t>::max()) {
        if (is_integer) {
            return finalize_all_blocks_internal<short_int>(ptr0);
        } else {
            return finalize_all_blocks_internal<short_dbl>(ptr0);
        }
    } else {
        if (is_integer) {
            return finalize_all_blocks_internal<regular_int>(ptr0);
        } else {
            return finalize_all_blocks_internal<regular_dbl>(ptr0);
        }
    }
}
