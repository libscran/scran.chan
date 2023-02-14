#include "config.h"

#include <array>
#include "tatami/tatami.hpp"
#include "tatami/ext/hdf5/load_hdf5_matrix.hpp"
#include "Rcpp.h"
#include "tatamize.h"
#include "H5Cpp.h"

template<bool ROW, typename Tx, typename Ti>
SEXP initialize_from_hdf5_base(const std::string& file, const std::string& name, size_t nrow, size_t ncol) {
    auto mat = tatami::load_hdf5_compressed_sparse_matrix<ROW, double, int, std::vector<Tx>, std::vector<Ti> >(nrow, ncol, file, name + "/data", name + "/indices", name + "/indptr");
    typedef decltype(mat) Matrix;
    return new_MatrixChan(new Matrix(std::move(mat)));
}

template<typename Tx, typename Ti>
SEXP initialize_from_hdf5_byrow(const std::string& file, const std::string& name, size_t nrow, size_t ncol, bool byrow) {
    if (byrow) {
        return initialize_from_hdf5_base<true, Tx, Ti>(file, name, nrow, ncol);
    } else {
        return initialize_from_hdf5_base<false, Tx, Ti>(file, name, nrow, ncol);
    }
}

template<typename Tx>
SEXP initialize_from_hdf5_i_max(const std::string& file, const std::string& name, size_t nrow, size_t ncol, bool byrow) {
    if ((byrow ? ncol : nrow) <= std::numeric_limits<uint16_t>::max()) {
        return initialize_from_hdf5_byrow<Tx, uint16_t>(file, name, nrow, ncol, byrow);
    } else {
        return initialize_from_hdf5_byrow<Tx, int>(file, name, nrow, ncol, byrow);
    }
}

//[[Rcpp::export(rng=false)]]
SEXP initialize_from_hdf5(std::string file, std::string name, size_t nrow, size_t ncol, bool byrow, bool forced) {
    bool is_float = false;
    bool is_ushort = false;
    {
        H5::H5File handle(file, H5F_ACC_RDONLY);
        auto xhandle = handle.openDataSet(name + "/data");
        auto xtype = xhandle.getDataType();
        is_float = (xtype.getClass() == H5T_FLOAT);
        if (!is_float) {
            H5::IntType xitype(xhandle);
            is_ushort = (xitype.getSize() <= 2 && xitype.getSign() == H5T_SGN_NONE);
        }
    }

    if (is_float && !forced) {
        return initialize_from_hdf5_i_max<double>(file, name, nrow, ncol, byrow);
    } else if (is_ushort) {
        return initialize_from_hdf5_i_max<uint16_t>(file, name, nrow, ncol, byrow);
    } else {
        return initialize_from_hdf5_i_max<int>(file, name, nrow, ncol, byrow);
    }
}
