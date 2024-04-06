/**
 * C++ code for the trace computation (to speed up processing
 * of massive matrices in qtlMSelect()).
 *
 * Code author: Russell Edson, Biometry Hub
 * Date last modified: 06/04/2024
 */

#include <string>
#include <vector>

#include <Rcpp.h>

//' Computes the trace of a given S4 Matrix class
//'
//' @param mat The matrix to compute the trace for
//' @return The trace (i.e. sum of diagonal entries)
//' @noRd
// [[Rcpp::export]]
double computeTrace(const Rcpp::S4 &mat) {
  std::string className = mat.attr("class");
  double trace = 0.0;

  if (className == "dgeMatrix") {
    // S4 structure of a dgeMatrix:
    //   Slot "Dim": IntegerVector of length 2 (for m, n dimensions)
    //   Slot "x": NumericVector of length m*n (values in column-major order)
    Rcpp::IntegerVector dim = mat.slot("Dim");
    Rcpp::NumericVector x = mat.slot("x");

    for (int k = 0; k < dim[1]; k++) {
      trace += x[k*dim[0] + k];
    }
  } else if (className == "dgCMatrix") {
    // S4 structure of a dgCMatrix:
    //   Slot "Dim": IntegerVector of length 2 (for m, n dimensions)
    //   Slot "x": NumericVector of length nz (non-zero entries)
    //   Slot "i": IntegerVector of length nz (row indices, 0-based)
    //   Slot "p": IntegerVector of length n+1 (column non-zero counts)
    //Rcpp::IntegerVector dim = mat.slot("Dim");
    Rcpp::NumericVector x = mat.slot("x");
    Rcpp::IntegerVector i = mat.slot("i");
    Rcpp::IntegerVector p = mat.slot("p");

    int valPtr = 0;
    for (int k = 0; k < p.size() - 1; k++) {
      int nonzero = p[k + 1] - p[k];
      for (int j = 0; j < nonzero; j++) {
        if (i[valPtr + j] == k) {
          trace += x[valPtr + j];
        }
      }
      valPtr += nonzero;
    }
  } else {  //TODO: Other Matrix:: classes?
    Rcpp::stop("computeTrace() not defined for objects of class " + className);
  }
  return trace;
}
