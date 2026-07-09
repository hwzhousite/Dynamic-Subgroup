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
