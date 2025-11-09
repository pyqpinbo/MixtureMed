#' EM estimation for heterogeneous mediation with tree-based \eqn{\pi(X)}
#'
#' This function estimates a K-component heterogeneous mediation model
#' using an EM algorithm where the mixing proportions \eqn{\pi_k(X)}
#' are modeled nonparametrically by a decision tree. Optionally, it
#' computes nonparametric bootstrap confidence intervals for
#' subgroup-specific NIE and NDE.
#'
#' @param K Integer. Number of mixture components (\code{K >= 2}).
#' @param data_model A \code{data.frame} with columns \code{T, X1, X2, M, Y}.
#' @param max_iter Integer. Maximum number of EM iterations.
#' @param tol Numeric. Convergence tolerance for relative log-likelihood
#'   change and coefficient change.
#' @param compute_ci Logical. If \code{TRUE}, perform nonparametric bootstrap.
#' @param sims_ci Integer. Number of bootstrap replications.
#' @param verbose Logical. If \code{TRUE}, print EM and bootstrap progress.
#' @param tree_args A list of arguments passed to the internal
#'   tree-fitting function (e.g. \code{minsplit}, \code{cp},
#'   \code{maxdepth} for \code{rpart}).
#'
#' @return A list extending the output of the internal core EM with:
#' \describe{
#'   \item{coefs}{Estimated subgroup-specific \code{alpha, beta, gamma}.}
#'   \item{piz}{Estimated mixing proportions \eqn{\hat{\pi}_k(X)}.}
#'   \item{w}{Posterior responsibilities.}
#'   \item{sigma_M}{Component-specific residual standard deviations for \code{M}.}
#'   \item{sigma_Y}{Component-specific residual standard deviations for \code{Y}.}
#'   \item{loglik}{Final log-likelihood.}
#'   \item{pi_model}{Fitted tree model(s) for \eqn{\pi_k(X)}.}
#'   \item{n_iter}{Number of EM iterations used.}
#'   \item{ci}{If \code{compute_ci = TRUE}, a list of length \code{K} with
#'     entries \code{NIE_hat}, \code{NIE_ci}, \code{NDE_hat}, \code{NDE_ci}.}
#' }
#'
#' @examples
#' \donttest{
#' set.seed(2025)
#' sim  <- data_generate(200, setting = 1L)
#' fit  <- Med_EM_tree(
#'   K          = 2L,
#'   data_model = sim$data,
#'   max_iter   = 100,
#'   tol        = 1e-4,
#'   compute_ci = FALSE
#' )
#' fit$coefs
#' }
#'
#' @export

Med_EM_tree <- function(
    K          = 2L,
    data_model,
    max_iter   = 200L,
    tol        = 1e-4,
    compute_ci = FALSE,
    sims_ci    = 500L,   # number of bootstrap replications (B)
    verbose    = TRUE,
    tree_args  = list(minsplit = 20L, cp = 0.01, maxdepth = 5L)
) {
  # Wrapper for EM + (optional) nonparametric bootstrap CIs
  # for subgroup-specific NIE and NDE.
  #
  # Args:
  #   K, data_model, max_iter, tol, verbose, tree_args : passed to Med_EM_tree_core()
  #   compute_ci : if TRUE, perform nonparametric bootstrap
  #   sims_ci    : number of bootstrap replications (B)
  #
  # Returns: same as Med_EM_tree_core plus
  #   ci : if compute_ci = TRUE,
  #        list of length K, each with NIE_hat, NIE_ci, NDE_hat, NDE_ci
  
  ## 1) Fit core EM on original data
  fit <- Med_EM_tree_core(
    K          = K,
    data_model = data_model,
    max_iter   = max_iter,
    tol        = tol,
    verbose    = verbose,
    tree_args  = tree_args
  )
  
  ci_list <- NULL
  
  if (compute_ci) {
    n <- nrow(data_model)
    K <- K  # just to emphasize that we use the same K
    
    # Baseline coefficient estimates used to align bootstrap labels
    coefs_ref <- fit$coefs
    
    # Storage for bootstrap NIE / NDE
    NIE_boot <- matrix(NA_real_, nrow = sims_ci, ncol = K)
    NDE_boot <- matrix(NA_real_, nrow = sims_ci, ncol = K)
    
    for (b in 1:sims_ci) {
      idx_b <- sample.int(n, size = n, replace = TRUE)
      data_b <- data_model[idx_b, , drop = FALSE]
      
      # Refit EM on the bootstrap sample (silent)
      fit_b <- Med_EM_tree_core(
        K          = K,
        data_model = data_b,
        max_iter   = max_iter,
        tol        = tol,
        verbose    = FALSE,
        tree_args  = tree_args
      )
      
      # Align labels to reduce label switching issues
      coefs_b <- fit_b$coefs
      if (K == 2L) {
        # For K = 2, use a simple permutation-based alignment
        coefs_b_aligned <- align_coefs(coefs_b, coefs_ref)
      } else {
        # For K > 2 we currently skip explicit label alignment
        coefs_b_aligned <- coefs_b
      }
      
      for (k in 1:K) {
        alpha_b <- coefs_b_aligned[k, "alpha"]
        beta_b  <- coefs_b_aligned[k, "beta"]
        gamma_b <- coefs_b_aligned[k, "gamma"]
        
        NIE_boot[b, k] <- alpha_b * beta_b
        NDE_boot[b, k] <- gamma_b
      }
      
      if (verbose && (b %% 50L == 0L)) {
        cat(sprintf("  Bootstrap %d / %d done\n", b, sims_ci))
      }
    }
    
    # NIE / NDE point estimates and percentile CIs based on the original sample
    ci_list <- vector("list", K)
    for (k in 1:K) {
      alpha_hat <- coefs_ref[k, "alpha"]
      beta_hat  <- coefs_ref[k, "beta"]
      gamma_hat <- coefs_ref[k, "gamma"]
      
      NIE_hat_k <- alpha_hat * beta_hat
      NDE_hat_k <- gamma_hat
      
      NIE_ci_k <- stats::quantile(NIE_boot[, k], probs = c(0.025, 0.975),
                                  na.rm = TRUE)
      NDE_ci_k <- stats::quantile(NDE_boot[, k], probs = c(0.025, 0.975),
                                  na.rm = TRUE)
      
      ci_list[[k]] <- list(
        NIE_hat = NIE_hat_k,
        NIE_ci  = as.numeric(NIE_ci_k),
        NDE_hat = NDE_hat_k,
        NDE_ci  = as.numeric(NDE_ci_k)
      )
    }
  }
  
  fit$ci <- ci_list
  fit
}
