#' Align component labels for K = 2
#'
#' Align labels for K = 2 so that component 1 in the bootstrap sample
#' corresponds to component 1 in the reference estimate (up to permutation).
#'
#' @param coefs_boot A \code{2 x 3} matrix of bootstrap estimates
#'   with columns \code{alpha, beta, gamma}.
#' @param coefs_ref A \code{2 x 3} matrix of reference estimates
#'   with the same column structure.
#'
#' @return A \code{2 x 3} matrix of aligned bootstrap coefficients.
#'
#' @keywords internal
#' @noRd
align_coefs <- function(coefs_boot, coefs_ref) {
  if (!all(dim(coefs_boot) == c(2L, 3L))) {
    stop("align_coefs() is only implemented for K = 2 and 3 parameters.")
  }
  
  err1 <- sum((coefs_boot            - coefs_ref)^2)
  err2 <- sum((coefs_boot[c(2, 1), ] - coefs_ref)^2)
  if (err2 < err1) {
    coefs_boot <- coefs_boot[c(2, 1), ]
  }
  coefs_boot
}
