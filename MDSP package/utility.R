# Helper function for scientific notation
scientific <- function(x, digits = 3) {
  sprintf(paste0("%.", digits, "e"), x)
}

# Function to evaluate group assignment accuracy
evaluate_group_accuracy <- function(true_groups, estimated_groups) {
  if (!all(dim(true_groups) == dim(estimated_groups))) {
    stop("true_groups and estimated_groups must have the same dimensions.")
  }

  n <- nrow(true_groups)
  m <- ncol(true_groups)

  time_accuracies <- rep(NA_real_, m)

  has_gtools <- requireNamespace("gtools", quietly = TRUE)

  for (t in 1:m) {
    true_t <- true_groups[, t]
    est_t  <- estimated_groups[, t]

    # keep only non-missing pairs
    valid_idx <- which(!is.na(true_t) & !is.na(est_t))

    if (length(valid_idx) == 0) {
      # no observed labels at this time point
      time_accuracies[t] <- NA_real_
      next
    }

    true_t_obs <- true_t[valid_idx]
    est_t_obs  <- est_t[valid_idx]
    n_obs      <- length(valid_idx)

    K_true <- max(true_t_obs)
    K_est  <- max(est_t_obs)

    # If we have gtools and both K's are small, try all permutations for label matching
    if (has_gtools && K_true == K_est && K_true <= 6) {
      perms <- gtools::permutations(n = K_true, r = K_true)

      best_accuracy <- 0
      for (p in 1:nrow(perms)) {
        perm <- perms[p, ]
        # relabel estimated groups via permutation
        est_perm <- perm[est_t_obs]
        acc <- sum(true_t_obs == est_perm) / n_obs
        if (acc > best_accuracy) best_accuracy <- acc
      }
      time_accuracies[t] <- best_accuracy

    } else {
      # Fallback: no permutation matching (or unequal K), just raw label agreement
      time_accuracies[t] <- sum(true_t_obs == est_t_obs) / n_obs
    }
  }

  list(
    overall_accuracy = mean(time_accuracies, na.rm = TRUE),
    time_accuracies  = time_accuracies
  )
}

# Function to simulate data for the model
simulate_model_data <- function(n_individuals = 100, n_time = 10, n_subgroups = 3,
                                q_controls = 2, p_treatments = 2, d_covariates = 3, ysd = 1) {

  # Set dimensions
  n <- n_individuals
  m <- n_time
  K <- n_subgroups

  # Time points normalized to [0, 1]
  time_points <- seq(0, 1, length.out = m)

  # Simulate control variables Z (homogeneous effects)
  Z <- array(rnorm(n * m * q_controls), dim = c(n, m, q_controls))
  #array(rnorm(n * m * q_controls), dim = c(n, m, q_controls))
  alpha <- matrix(1, nrow = m, ncol = q_controls)
  #matrix(rnorm(m * q_controls, 0, 0.5), nrow = m, ncol = q_controls)

  # Simulate treatment variables X
  X <- array(rnorm(n * m * p_treatments), dim = c(n, m, p_treatments))

  # Simulate time-varying covariates h_t
  h_t <- matrix(0, nrow = m, ncol = d_covariates)
  h_t[, 1] <- time_points  # Linear time trend
  if (d_covariates > 1) {
    h_t[, 2] <- time_points #sin(2 * pi * time_points)  # Sinusoidal component
  }
  if (d_covariates > 2) {
    h_t[, 3:d_covariates] <- time_points
    #+ matrix(rnorm((d_covariates - 2) * m, 0, 1),     nrow = m, ncol = d_covariates - 2)
  }

  # True baseline coefficients g_0k for each subgroup
  g_0k <- matrix(0, nrow = K, ncol = p_treatments)
  g_0k[1,] <- 1
  g_0k[2,] <- -1

  # True time-trend parameters G_k
  G_k <- array(0, dim = c(K, p_treatments, d_covariates))
  G_k[1,,] <- 1
  G_k[2,,] <- -1
  #array(rnorm(K * p_treatments * d_covariates, 0, 0.4),      dim = c(K, p_treatments, d_covariates))

  # Generate parametric subgroup centers: γ_kt = g_0k + G_k * h_t
  gamma_kt <- array(0, dim = c(K, m, p_treatments))
  for (k in 1:K) {
    for (t in 1:m) {
      for (j in 1:p_treatments) {
        gamma_kt[k, t, j] <- g_0k[k, j] + sum(G_k[k, j, ] * h_t[t, ])
      }
    }
  }

  # True subgroup assignments (can change over time)
  true_groups <- array(0, dim = c(n, m))

  # Generate dynamic group assignments and response
  y <- matrix(0, nrow = n, ncol = m)
  beta_it <- array(0, dim = c(n, m, p_treatments))

  for (i in 1:n) {
    # Individual can switch groups over time with some probability
    if(i <= n/2){
      current_group <- 1#sample(1:K, 1)
    }else{
      current_group <- 2
    }
    for (t in 1:m) {
      # Allow group switching with small probability
      if (t > 1 && runif(1) < 0.12) {
        #current_group <- sample(1:K, 1)
      }

      true_groups[i, t] <- current_group
      # Add some noise around the center
      beta_it[i, t, ] <- gamma_kt[current_group, t, ]

      # Response: y_it = Z_it * α_t + X_it * β_it + ε_it
      z_contrib <- sum(Z[i, t, ] * alpha[t, ])
      x_contrib <- sum(X[i, t, ] * beta_it[i, t, ])
      epsilon <- rnorm(1, 0, ysd)

      y[i, t] <- z_contrib + x_contrib + epsilon
    }
  }

  return(list(
    y = y, Z = Z, X = X, h_t = h_t,
    true_groups = true_groups, beta_it = beta_it,
    alpha = alpha, gamma_kt = gamma_kt,
    g_0k = g_0k, G_k = G_k,
    time_points = time_points,
    n = n, m = m, K = K, q = q_controls, p = p_treatments, d = d_covariates
  ))
}



roc_curve <- function(y_true, y_score) {
  stopifnot(length(y_true) == length(y_score))
  y_true <- as.integer(y_true)

  o <- order(-y_score, decreasing = FALSE)  # we’ll reverse below
  o <- order(-y_score)                      # sort by score desc
  y <- y_true[o]
  s <- y_score[o]

  P <- sum(y == 1); N <- sum(y == 0)
  if (P == 0 || N == 0) stop("Need both positives and negatives.")

  # Sweep thresholds at each unique score
  uniq_idx <- c(which(diff(s) != 0), length(s))
  tp <- fp <- numeric(length(uniq_idx))
  cum_pos <- cumsum(y == 1)
  cum_neg <- cumsum(y == 0)

  for (k in seq_along(uniq_idx)) {
    i <- uniq_idx[k]
    tp[k] <- cum_pos[i]
    fp[k] <- cum_neg[i]
  }

  TPR <- tp / P
  FPR <- fp / N

  # Add (0,0) and (1,1) endpoints
  FPR <- c(0, FPR, 1)
  TPR <- c(0, TPR, 1)

  # AUC via trapezoid
  auc <- sum( (FPR[-1] - FPR[-length(FPR)]) * (TPR[-1] + TPR[-length(TPR)]) / 2 )

  list(FPR = FPR, TPR = TPR, thresholds = c(Inf, unique(s), -Inf), AUC = auc)
}


## RMSE between estimated and true gamma, with label-switch fix for K = 2
gamma_rmse <- function(g_hat, g_true) {
  stopifnot(all(dim(g_hat) == dim(g_true)))

  K <- dim(g_true)[1]

  if (K == 1) {
    return(sqrt(mean((g_hat - g_true)^2, na.rm = TRUE)))
  } else if (K == 2) {
    rmse1 <- sqrt(mean((g_hat - g_true)^2, na.rm = TRUE))
    rmse2 <- sqrt(mean((g_hat[c(2, 1), , ] - g_true)^2, na.rm = TRUE))
    return(min(rmse1, rmse2))
  } else {
    # For K > 2, you can later add a Hungarian-matching version if needed
    return(sqrt(mean((g_hat - g_true)^2, na.rm = TRUE)))
  }
}

## From individual-level beta + group labels to subgroup-level gamma
## beta_hat: n x m x p_treat; group_mat: n x m
estimate_gamma_from_beta <- function(beta_hat, group_mat, K) {
  n <- dim(beta_hat)[1]
  m <- dim(beta_hat)[2]
  p_treat <- dim(beta_hat)[3]

  gamma_hat <- array(NA_real_, dim = c(K, m, p_treat))

  for (tt in 1:m) {
    for (kk in 1:K) {
      ind <- which(group_mat[, tt] == kk)
      if (length(ind) == 1) {
        gamma_hat[kk, tt, ] <- beta_hat[ind, tt, ]
      } else if (length(ind) > 1) {
        gamma_hat[kk, tt, ] <- colMeans(beta_hat[ind, tt, , drop = FALSE])
      }
    }
  }
  gamma_hat
}
