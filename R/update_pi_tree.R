#' Update tree-based mixing proportions
#'
#' Internal helper to fit decision tree(s) for \eqn{\pi_k(X)} given
#' current responsibilities.
#'
#' @param w An \eqn{n \times K} matrix of responsibilities.
#' @param Z A \code{data.frame} of covariates used to model \eqn{\pi_k(X)}.
#' @param ... Additional control arguments passed to \code{rpart()}.
#'
#' @return A list with elements \code{piz} (estimated \eqn{\pi_k(X)})
#'   and \code{model} (the fitted tree object(s)).
#'
#' @keywords internal
#' @noRd

update_pi_tree <- function(
    w,
    Z,
    minsplit = 20L,
    cp       = 0.01,
    maxdepth = 3L,
    clip     = 1e-4
) {
  # Update mixing proportions pi_k(X) using regression trees.
  # Supports general K >= 2.
  #
  # Args:
  #   w        : n x K matrix of responsibilities (posterior probabilities)
  #   Z        : data.frame or matrix of covariates driving pi(X)
  #   minsplit : minimum number of observations in a node before split
  #   cp       : complexity parameter for pruning
  #   maxdepth : maximal depth of each tree
  #   clip     : lower/upper bound to avoid probabilities 0 or 1
  #
  # Returns:
  #   list(
  #     piz   : n x K matrix of updated pi_k(X),
  #     model : list of length K containing fitted rpart objects
  #   )
  
  if (!is.matrix(w)) {
    w <- as.matrix(w)
  }
  n <- nrow(w)
  K <- ncol(w)
  
  if (K < 2L) {
    stop("update_pi_tree() requires at least K = 2.")
  }
  
  Z <- as.data.frame(Z)
  
  if (any(!is.finite(w))) {
    stop("Non-finite values in responsibilities w.")
  }
  
  ctrl <- rpart::rpart.control(
    minsplit = minsplit,
    cp       = cp,
    maxdepth = maxdepth
  )
  
  # Fit one regression tree per component to approximate w[, k]
  pi_hat_raw <- matrix(NA_real_, nrow = n, ncol = K)
  model_list <- vector("list", K)
  
  for (k in 1:K) {
    yk <- w[, k]
    df_tree <- data.frame(Z, y = yk)
    
    fit_k <- rpart::rpart(
      formula = y ~ .,
      data    = df_tree,
      method  = "anova",
      control = ctrl
    )
    
    pred_k <- as.numeric(predict(fit_k, newdata = Z))
    pi_hat_raw[, k] <- pred_k
    model_list[[k]] <- fit_k
  }
  
  # Handle non-finite predictions
  bad_pred <- !is.finite(pi_hat_raw)
  if (any(bad_pred)) {
    pi_hat_raw[bad_pred] <- 1 / K
  }
  
  # Clip to avoid exact 0 or 1
  pi_hat_raw <- pmin(pmax(pi_hat_raw, clip), 1 - clip)
  
  # Normalize row-wise to sum to 1
  row_sums <- rowSums(pi_hat_raw)
  bad_row  <- !is.finite(row_sums) | row_sums <= 0
  if (any(bad_row)) {
    pi_hat_raw[bad_row, ] <- 1 / K
    row_sums[bad_row]     <- 1
  }
  piz_hat <- pi_hat_raw / row_sums
  
  list(piz = piz_hat, model = model_list)
}
