/**
 * C+ code (using the Armadillo libraries via RcppArmadillo) to compute
 * the vqtilde component for the multivariate outlier statistics and
 * speed up processing in the qtlSelect() function.
 *
 * Code author: Russell Edson, Biometry Hub
 * Date last modified: 19/11/2024
 */
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Computes the matrix product vqtilde
//'
//' Note: This function assumes that the matrices passed in are all
//' base R matrices (instead of Matrix::dgeMatrix types, for example).
//'
//' @param trans The matrix product M^T(MM^T)^{-1}, typically many rows
//' @param Ginv The inverse matrix G_a^-
//' @param vatilde The variance matrix var(tilde(a))
//' @param ntrait The number of traits
//' @return The matrix vqtilde
// [[Rcpp::export]]
Rcpp::NumericVector compute_vqtilde(Rcpp::NumericMatrix trans,
    Rcpp::NumericMatrix Ginv, Rcpp::NumericMatrix vatilde, int ntrait) {
  // Map input matrices to Armadillo matrix objects (without copying).
  arma::uword n = trans.nrow(), p = trans.ncol();
  arma::mat trans_mat = arma::mat(trans.begin(), n, p, false);
  arma::mat Ginv_mat = arma::mat(Ginv.begin(), Ginv.nrow(), Ginv.ncol(), false);
  arma::mat vatilde_mat = arma::mat(vatilde.begin(), vatilde.nrow(),
      vatilde.ncol(), false);

  // We augment the return value vqtilde with the expected row names.
  Rcpp::NumericVector vqtilde(n);
  Rcpp::List dimnames = trans.attr("dimnames");
  Rcpp::CharacterVector rownames = dimnames[0];
  vqtilde.attr("names") = rownames;

  // For the computing, we work on the kth slice of trans for each
  // genetic marker separately, and also divvy up vatilde into blocks
  // for the matrix multiplication.
  arma::mat varq(ntrait, ntrait);
  for (arma::uword k = 0; k < n; k++) {
    for (arma::uword i = 0; i < varq.n_rows; i++) {
      for (arma::uword j = 0; j < varq.n_cols; j++) {
        varq.submat(i, j, i, j) = trans_mat.row(k) * 
            vatilde_mat.submat(i*p, j*p, (i+1)*p - 1, (j+1)*p - 1) *
            trans_mat.row(k).t();
      }
    }
    vqtilde[k] = arma::trace(Ginv_mat * varq);
  }
  return vqtilde;
}
