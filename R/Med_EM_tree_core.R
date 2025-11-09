#' Core EM algorithm with tree-based mixing proportions
#'
#' This function implements the core EM algorithm for a K-component
#' heterogeneous mediation model, where the mixing proportions
#' \eqn{\pi_k(X)} are modeled via decision trees.
#'
#' @param K Integer. Number of components (\code{K >= 2}).
#' @param data_model A \code{data.frame} with columns \code{T, X1, X2, M, Y}.
#' @param max_iter Maximum number of EM iterations.
#' @param tol Convergence tolerance for the relative log-likelihood
#'   change and coefficient change.
#' @param verbose Logical. If \code{TRUE}, print EM progress.
#' @param tree_args A list of control parameters passed to
#'   \code{update_pi_tree()}.
#'
#' @return A list with elements:
#' \describe{
#'   \item{piz}{\eqn{n \times K} matrix of mixing proportions.}
#'   \item{w}{\eqn{n \times K} matrix of responsibilities.}
#'   \item{sigma_M}{Length-\code{K} vector of \eqn{\sigma_{Mk}}.}
#'   \item{sigma_Y}{Length-\code{K} vector of \eqn{\sigma_{Yk}}.}
#'   \item{params}{List of fitted regression models per component.}
#'   \item{coefs}{\code{K x 3} matrix of \code{alpha, beta, gamma}.}
#'   \item{loglik}{Final log-likelihood value.}
#'   \item{pi_model}{Tree model(s) used to model \eqn{\pi_k(X)}.}
#'   \item{n_iter}{Number of EM iterations used.}
#' }
#'
#' @keywords internal
#' @noRd

Med_EM_tree_core <- function(
    K          = 2L,
    data_model,
    max_iter   = 200L,
    tol        = 1e-4,
    verbose    = TRUE,
    tree_args  = list(minsplit = 20L, cp = 0.01, maxdepth = 5L)
) {
  # Core EM algorithm for a K-component heterogeneous mediation model
  # with tree-based mixing proportions pi_k(X).
  # This function does NOT compute bootstrap CIs.
  #
  # Args:
  #   K          : number of latent components (K >= 2)
  #   data_model : data.frame with columns T, X1, X2, M, Y
  #   max_iter   : maximum number of EM iterations
  #   tol        : tolerance for convergence (relative loglik + coef change)
  #   verbose    : print EM progress if TRUE
  #   tree_args  : list of control parameters passed to update_pi_tree()
  #
  # Returns:
  #   list(
  #     piz      : n x K matrix of mixing proportions at convergence
  #     w        : n x K matrix of responsibilities at convergence
  #     sigma_M  : length-K vector of sigma for M
  #     sigma_Y  : length-K vector of sigma for Y
  #     params   : list of fitted models per component
  #     coefs    : K x 3 matrix (alpha_k, beta_k, gamma_k)
  #     loglik   : final log-likelihood
  #     pi_model : fitted tree model(s) used for pi(X) (e.g., list of rpart objects)
  #     n_iter   : number of EM iterations actually used
  #   )
  
  if (K < 2L) {
    stop("Med_EM_tree_core() requires K >= 2.")
  }
  
  data_model <- as.data.frame(data_model)
  required_cols <- c("T", "X1", "X2", "M", "Y")
  if (!all(required_cols %in% names(data_model))) {
    stop("data_model must contain columns: T, X1, X2, M, Y.")
  }
  
  n <- nrow(data_model)
  
  # Covariates driving pi(X): here X1 and X2
  Z <- data_model[, c("X1", "X2"), drop = FALSE]
  
  T_vec <- data_model$T
  M_vec <- data_model$M
  Y_vec <- data_model$Y
  
  ## ---------- Initialization ----------
  
  # 1) Initial clustering on (M, Y)
  km <- kmeans(scale(data_model[, c("M", "Y")]), centers = K, nstart = 10L)
  init_group <- km$cluster
  
  # 2) Initial responsibilities: hard labels from k-means
  w <- matrix(0, nrow = n, ncol = K)
  for (k in 1:K) {
    w[init_group == k, k] <- 1
  }
  
  # 3) Initial pi(X) via tree(s), using initial responsibilities
  pi_res  <- do.call(update_pi_tree, c(list(w = w, Z = Z), tree_args))
  piz      <- pi_res$piz
  pi_model <- pi_res$model
  
  # 4) Initial component-specific regression parameters and sigmas
  params   <- vector("list", K)
  sigma_M  <- numeric(K)
  sigma_Y  <- numeric(K)
  coefs    <- matrix(NA_real_, nrow = K, ncol = 3L,
                     dimnames = list(NULL, c("alpha", "beta", "gamma")))
  
  for (k in 1:K) {
    wk <- w[, k]
    eff_n <- sum(wk)
    if (eff_n < 1e-6) {
      # If a component is essentially empty, allocate equal weights
      wk <- rep(1 / n, n)
      eff_n <- 1
    }
    wk_norm <- wk / eff_n
    
    fit_M <- lm(M ~ T + X1 + X2, data = data_model, weights = wk_norm)
    fit_Y <- lm(Y ~ T + M + X1 + X2, data = data_model, weights = wk_norm)
    
    params[[k]] <- list(fit_M = fit_M, fit_Y = fit_Y)
    
    mu_M <- predict(fit_M, newdata = data_model)
    mu_Y <- predict(fit_Y, newdata = data_model)
    
    sigma_M[k] <- sqrt(sum(wk_norm * (M_vec - mu_M)^2))
    sigma_Y[k] <- sqrt(sum(wk_norm * (Y_vec - mu_Y)^2))
    
    sigma_M[k] <- max(sigma_M[k], 1e-3)
    sigma_Y[k] <- max(sigma_Y[k], 1e-3)
    
    coefs[k, "alpha"] <- coef(fit_M)["T"]
    coefs[k, "beta"]  <- coef(fit_Y)["M"]
    coefs[k, "gamma"] <- coef(fit_Y)["T"]
  }
  
  prev_loglik <- -Inf
  prev_coefs  <- coefs
  n_iter_used <- NA_integer_
  
  ## ---------- EM iterations ----------
  
  for (iter in 1:max_iter) {
    ## ---------- E-step ----------
    
    # Log densities for each component
    log_dens <- matrix(NA_real_, nrow = n, ncol = K)
    for (k in 1:K) {
      mu_M <- predict(params[[k]]$fit_M, newdata = data_model)
      mu_Y <- predict(params[[k]]$fit_Y, newdata = data_model)
      
      log_f_M <- dnorm(M_vec, mean = mu_M, sd = sigma_M[k], log = TRUE)
      log_f_Y <- dnorm(Y_vec, mean = mu_Y, sd = sigma_Y[k], log = TRUE)
      
      log_dens[, k] <- log_f_M + log_f_Y
    }
    
    # Combine with log pi(X)
    log_pi <- log(piz)
    log_w  <- log_pi + log_dens
    
    # log-sum-exp stabilization
    row_max      <- apply(log_w, 1L, max)
    log_w_stable <- sweep(log_w, 1L, row_max, FUN = "-")
    w_unnorm     <- exp(log_w_stable)
    
    row_sums <- rowSums(w_unnorm)
    bad <- !is.finite(row_sums) | row_sums <= 0
    if (any(bad)) {
      # Fallback: assign uniform responsibilities
      w_unnorm[bad, ] <- 1 / K
      row_sums[bad]   <- 1
    }
    w <- w_unnorm / row_sums
    
    ## ---------- M-step ----------
    
    params_new  <- vector("list", K)
    sigma_M_new <- numeric(K)
    sigma_Y_new <- numeric(K)
    coefs_new   <- matrix(NA_real_, nrow = K, ncol = 3L,
                          dimnames = list(NULL, c("alpha", "beta", "gamma")))
    
    for (k in 1:K) {
      wk <- w[, k]
      eff_n <- sum(wk)
      if (eff_n < 1e-4) {
        # Component nearly empty: keep previous parameters
        params_new[[k]] <- params[[k]]
        sigma_M_new[k]  <- sigma_M[k]
        sigma_Y_new[k]  <- sigma_Y[k]
        coefs_new[k, ]  <- coefs[k, ]
        next
      }
      
      wk_norm <- wk / eff_n
      
      fit_M <- lm(M ~ T + X1 + X2, data = data_model, weights = wk_norm)
      fit_Y <- lm(Y ~ T + M + X1 + X2, data = data_model, weights = wk_norm)
      
      params_new[[k]] <- list(fit_M = fit_M, fit_Y = fit_Y)
      
      mu_M <- predict(fit_M, newdata = data_model)
      mu_Y <- predict(fit_Y, newdata = data_model)
      
      sigma_M_k <- sqrt(sum(wk_norm * (M_vec - mu_M)^2))
      sigma_Y_k <- sqrt(sum(wk_norm * (Y_vec - mu_Y)^2))
      
      sigma_M_new[k] <- max(sigma_M_k, 1e-3)
      sigma_Y_new[k] <- max(sigma_Y_k, 1e-3)
      
      coefs_new[k, "alpha"] <- coef(fit_M)["T"]
      coefs_new[k, "beta"]  <- coef(fit_Y)["M"]
      coefs_new[k, "gamma"] <- coef(fit_Y)["T"]
    }
    
    params  <- params_new
    sigma_M <- sigma_M_new
    sigma_Y <- sigma_Y_new
    coefs   <- coefs_new
    
    ## ---------- Update pi(X) via tree ----------
    pi_res  <- do.call(update_pi_tree, c(list(w = w, Z = Z), tree_args))
    piz      <- pi_res$piz
    pi_model <- pi_res$model
    
    ## ---------- Compute log-likelihood ----------
    log_dens <- matrix(NA_real_, nrow = n, ncol = K)
    for (k in 1:K) {
      mu_M <- predict(params[[k]]$fit_M, newdata = data_model)
      mu_Y <- predict(params[[k]]$fit_Y, newdata = data_model)
      
      log_f_M <- dnorm(M_vec, mean = mu_M, sd = sigma_M[k], log = TRUE)
      log_f_Y <- dnorm(Y_vec, mean = mu_Y, sd = sigma_Y[k], log = TRUE)
      
      log_dens[, k] <- log_f_M + log_f_Y
    }
    loglik_vec <- log(rowSums(piz * exp(log_dens)))
    loglik     <- sum(loglik_vec)
    
    # Convergence diagnostics
    ll_diff   <- abs((loglik - prev_loglik) / (abs(prev_loglik) + 1e-6))
    coef_diff <- max(abs(coefs - prev_coefs))
    
    if (verbose) {
      cat(sprintf("Iter %3d: loglik = %.4f, rel_ll_diff = %.3e, coef_diff = %.3e\n",
                  iter, loglik, ll_diff, coef_diff))
    }
    
    if (ll_diff < tol && coef_diff < sqrt(tol)) {
      if (verbose) cat("Converged.\n")
      n_iter_used <- iter
      break
    }
    
    prev_loglik <- loglik
    prev_coefs  <- coefs
  }
  
  if (is.na(n_iter_used)) {
    n_iter_used <- max_iter
  }
  
  list(
    piz      = piz,
    w        = w,
    sigma_M  = sigma_M,
    sigma_Y  = sigma_Y,
    params   = params,
    coefs    = coefs,
    loglik   = loglik,
    pi_model = pi_model,
    n_iter   = n_iter_used
  )
}
