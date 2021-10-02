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
SEXP initialize_from_sparse(Rcpp::NumericVector x, Rcpp::IntegerVector i, Rcpp::IntegerVector p, int nrow, int ncol) {
    RcppVectorPlus x_(x);
    RcppVectorPlus i_(i);
    RcppVectorPlus p_(p);
    typedef tatami::CompressedSparseColumnMatrix<double, int, decltype(x_), decltype(i_), decltype(p_)> SparseMat;
    return new_MatrixChan(new SparseMat(nrow, ncol, std::move(x_), std::move(i_), std::move(p_), false));
}

template<typename X, typename I, class V>
SEXP initialize_from_blocks_internal(Rcpp::List indmat, const std::vector<V>& vlist, int nrow, Rcpp::IntegerVector ncols, int nthreads) {
    size_t nblocks = vlist.size();

    std::vector<Rcpp::IntegerMatrix> ilist;
    ilist.reserve(nblocks);
    for (size_t i = 0; i < nblocks; ++i) {
        ilist.emplace_back(indmat[i]);
    }

    // Initializing the total space.
    size_t nnzeros = 0;
    for (size_t i = 0; i < nblocks; ++i) {
        nnzeros += vlist[i].size();
    }

    std::vector<X> values(nnzeros);
    std::vector<I> indices(nnzeros);
    int ncol = std::accumulate(ncols.begin(), ncols.end(), 0);
    std::vector<size_t> indptrs(ncol + 1);

    std::vector<int> block_starts(ncols.size());
    for (size_t i = 1; i < nblocks; ++i) {
        block_starts[i] = block_starts[i-1] + vlist[i-1].size();
    }

    // Sorting within each block.
    #pragma omp parallel for num_threads(nthreads)
    for (size_t i = 0; i < nblocks; ++i) {
        const auto& curv = vlist[i];
        size_t curnnz = curv.size(); 
        const auto& curi = ilist[i];
        const int* row = curi.begin(), *col = row + curnnz;

        // Maybe it's already sorted?
        bool not_sorted = false;
        for (size_t z = 1; z < curnnz; ++z) {
            if (col[z] < col[z-1]) {
                not_sorted = true;
                break;
            } else if (col[z] == col[z-1]) {
                if (row[z] < row[z-1]) {
                    not_sorted = true;
                    break;
                }
            }
        }

        if (not_sorted) {
            std::vector<size_t> order(curnnz);
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(), [&](size_t l, size_t r) -> bool {
                if (col[l] == col[r]) {
                    return row[l] < row[r];
                } else {
                    return col[l] < col[r];
                }
            });

            size_t offset = block_starts[i];
            for (size_t z = 0; z < curnnz; ++z, ++offset) {
                auto o = order[z];
                values[offset] = curv[o];
                indices[offset] = row[o];
                ++indptrs[col[o]]; // col is already 1-based, so no need to +1.
            }
        } else {
            std::copy(curv.begin(), curv.end(), values.begin() + block_starts[i]);
            std::copy(row, row + curnnz, indices.begin() + block_starts[i]);
            for (size_t z = 0; z < curnnz; ++z) {
                ++indptrs[col[z]]; // col is already 1-based, so no need to +1.
            }
        }
    }

    // Creating the sparse matrix.
    for (auto& x : indices) {
        --x;
    }

    for (size_t i = 1; i <= ncol; ++i) {
        indptrs[i] += indptrs[i - 1];
    }

    typedef tatami::CompressedSparseColumnMatrix<double, int, decltype(values), decltype(indices), decltype(indptrs)> SparseMat;
    return new_MatrixChan(new SparseMat(nrow, ncol, std::move(values), std::move(indices), std::move(indptrs)));
}

//[[Rcpp::export(rng=false)]]
SEXP initialize_from_blocks(Rcpp::List indices, Rcpp::List values, int nrow, Rcpp::IntegerVector ncols, int nthreads) {
    bool is_integral = true;
    size_t nblocks = values.size();
    if (nblocks) {
        Rcpp::RObject first = values[0];
        if (first.sexp_type() == REALSXP) {
            is_integral = false;
        }
    }

    if (is_integral) {
        std::vector<Rcpp::IntegerVector> vlist;
        vlist.reserve(nblocks);
        for (size_t i = 0; i < nblocks; ++i) {
            vlist.emplace_back(values[i]);
        }
    
        if (nrow <= std::numeric_limits<uint16_t>::max()) {
            return initialize_from_blocks_internal<int, uint16_t>(indices, std::move(vlist), nrow, ncols, nthreads);
        } else {
            return initialize_from_blocks_internal<int, int>(indices, std::move(vlist), nrow, ncols, nthreads);
        }
    } else {
        std::vector<Rcpp::NumericVector> vlist;
        vlist.reserve(nblocks);
        for (size_t i = 0; i < nblocks; ++i) {
            vlist.emplace_back(values[i]);
        }
    
        if (nrow <= std::numeric_limits<uint16_t>::max()) {
            return initialize_from_blocks_internal<double, uint16_t>(indices, std::move(vlist), nrow, ncols, nthreads);
        } else {
            return initialize_from_blocks_internal<double, int>(indices, std::move(vlist), nrow, ncols, nthreads);
        }
    }
}
