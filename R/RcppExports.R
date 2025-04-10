# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Computes the matrix product vqtilde
#'
#' Note: This function assumes that the matrices passed in are all
#' base R matrices (instead of Matrix::dgeMatrix types, for example).
#'
#' @param trans The matrix product M^T(MM^T)^{-1}, typically many rows
#' @param Ginv The inverse matrix G_a^-
#' @param vatilde The variance matrix var(tilde(a))
#' @param ntrait The number of traits
#' @return The matrix vqtilde
compute_vqtilde <- function(trans, Ginv, vatilde, ntrait) {
    .Call(`_mvwgaim_compute_vqtilde`, trans, Ginv, vatilde, ntrait)
}

