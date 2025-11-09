#' Generate simulated data for heterogeneous mediation model
#'
#' This function generates data from a two-component heterogeneous
#' mediation model with different specifications of the mixing
#' proportion \eqn{\pi_1(X_1, X_2)}.
#'
#' @param n Integer. Sample size.
#' @param setting Integer in \code{1:3}.
#'   \itemize{
#'     \item \code{1}: linear logit \eqn{\pi_1(X)}.
#'     \item \code{2}: interaction-based nonlinearity.
#'     \item \code{3}: sinusoidal nonlinearity.
#'   }
#'
#' @return A list with elements:
#' \describe{
#'   \item{data}{A \code{data.frame} with columns \code{T, X1, X2, M, Y}.}
#'   \item{true}{A list containing true parameters and latent structure:
#'     \code{coefs}, \code{NIE}, \code{NDE}, \code{pi}, \code{group},
#'     \code{sigma_M}, \code{sigma_Y}.}
#' }
#'
#' @examples
#' set.seed(123)
#' sim <- data_generate(200, setting = 1L)
#' str(sim$data)
#'
#' @export
data_generate <- function(n, setting = 1L) {
  K <- 2L
  
  ## 1. Generate covariates
  T_var  <- rbinom(n, size = 1L, prob = 0.5)   # binary treatment
  X1_var <- rnorm(n, mean = 0, sd = 1)
  X2_var <- rnorm(n, mean = 0, sd = 1)
  
  ## 2. Mixing proportion π₁(X₁, X₂)
  eta <- switch(
    as.character(setting),
    "1" = 0.5 - 1.0 * X1_var + 1.0 * X2_var,
    "2" = 0.5 - 1.0 * X1_var * X2_var,
    "3" = 0.5 - 1.0 * (sin(X1_var) + 0.5 * cos(X2_var)),
    stop("setting must be 1, 2, or 3")
  )
  
  pi1 <- plogis(eta)
  pi1 <- pmin(pmax(pi1, 0.05), 0.95)  # avoid extreme probabilities
  p_mat <- cbind(pi1, 1 - pi1)
  
  ## 3. Latent group assignment
  true_group <- apply(p_mat, 1L, function(p) sample(1:K, 1L, prob = p))
  
  ## 4. True parameters for each component
  gamma0   <- c(0.4, 0.8)
  gamma_T  <- c(1.0, 1.5)   # αₖ
  gamma_X1 <- c(1.0, 2.0)
  gamma_X2 <- c(-1.0, -2.0)
  
  beta0   <- c(0.3, 0.6)
  beta_T  <- c(1.0, 3.0)    # γₖ (direct effect)
  beta_M  <- c(1.5, -1.5)   # βₖ (mediator effect)
  beta_X1 <- c(2.0, 4.0)
  beta_X2 <- c(-2.0, -4.0)
  
  sigma_M_true <- c(0.6, 0.8)
  sigma_Y_true <- c(0.5, 0.7)
  
  ## 5. Vectorized data generation
  g  <- true_group
  e1 <- rnorm(n, 0, sigma_M_true[g])
  e2 <- rnorm(n, 0, sigma_Y_true[g])
  
  M <- gamma0[g] + gamma_T[g]*T_var + gamma_X1[g]*X1_var + gamma_X2[g]*X2_var + e1
  Y <- beta0[g]  + beta_T[g]*T_var  + beta_M[g]*M      + beta_X1[g]*X1_var + beta_X2[g]*X2_var + e2
  
  ## 6. True mediation effects
  alpha_true <- gamma_T
  beta_true  <- beta_M
  gamma_true <- beta_T
  
  true_coefs <- cbind(alpha_true, beta_true, gamma_true)
  colnames(true_coefs) <- c("alpha", "beta", "gamma")
  
  data <- data.frame(T = T_var, X1 = X1_var, X2 = X2_var, M = M, Y = Y)
  
  list(
    data = data,
    true = list(
      coefs   = true_coefs,
      NIE     = alpha_true * beta_true,
      NDE     = gamma_true,
      pi      = p_mat,
      group   = true_group,
      sigma_M = sigma_M_true,
      sigma_Y = sigma_Y_true
    )
  )
}
