#include "tatamize.h"
#include <exception>

std::shared_ptr<tatami::Matrix<double, int> > tatamize(Rcpp::RObject x) {
    if (!x.isS4()) {
        // Had better be an ordinary matrix, then.
        if (x.sexp_type() == INTSXP) {
            Rcpp::IntegerMatrix mat(x);
            return std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseColumnMatrix<double, int, Rcpp::IntegerVector>(mat.nrow(), mat.ncol(), Rcpp::IntegerVector(x)));
        } else if (x.sexp_type() == REALSXP) {
            Rcpp::NumericMatrix mat(x);
            return std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseColumnMatrix<double, int, Rcpp::NumericVector>(mat.nrow(), mat.ncol(), Rcpp::NumericVector(x)));
        } else {
            throw std::runtime_error("unknown type for a dense matrix");
        }
    } 
    
    Rcpp::S4 cls(x);
    if (cls.is("dgCMatrix")) {
        Rcpp::IntegerVector dims(cls.slot("Dim"));
        return std::shared_ptr<tatami::NumericMatrix>(new tatami::CompressedSparseColumnMatrix<double, int, Rcpp::NumericVector, Rcpp::IntegerVector, Rcpp::IntegerVector>(
            dims[0], 
            dims[1],
            Rcpp::NumericVector(cls.slot("x")),
            Rcpp::IntegerVector(cls.slot("i")),
            Rcpp::IntegerVector(cls.slot("p"))
        ));
    }

    if (cls.is("dgeMatrix")) {
        Rcpp::IntegerVector dims(cls.slot("Dim"));
        return std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseColumnMatrix<double, int, Rcpp::NumericVector>(
            dims[0], 
            dims[1],
            Rcpp::NumericVector(cls.slot("x"))
        ));
    }

    throw std::runtime_error("unknown matrix type");
    return nullptr;
}
