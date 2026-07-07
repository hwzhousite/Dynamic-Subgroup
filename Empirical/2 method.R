##############################################################################
#  02_methods_empirical.R
#  Method wrappers for the rolling-window empirical evaluation.
#
#  Each function signature:
#    fit_*(Y_win, X_win, params) -> list(beta_agg, group_vec, Gamma_hat, K_used, method)
#
#  Arguments:
#    Y_win  : n × T_win  return matrix        (NA = firm inactive that month)
#    X_win  : n × T_win × p  characteristic array
#    params : named list of tuning parameters
#
#  Returns:
#    beta_agg  : n × p   time-aggregated coefficient per firm (used for ranking)
#    group_vec : length-n  integer group membership at end of window
#    Gamma_hat : K × T_win × p  group-level trajectories  (NULL if not applicable)
#    K_used    : integer  number of groups actually used   (NA if not applicable)
#    method    : character label
##############################################################################

suppressPackageStartupMessages({
  library(Matrix)
  library(splines)
})

## ── null-coalescing operator ──────────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

## ── utility: B-spline basis (T points → T × df matrix) ───────────────────────
make_bspline_basis <- function(T, df = 6, degree = 3) {
  bs(seq(0, 1, length.out = T), df = df, degree = degree, intercept = TRUE)
}

## ══════════════════════════════════════════════════════════════════════════════
##  METHOD 0  Homogeneous Factor Model — pooled OLS, time-invariant
## ══════════════════════════════════════════════════════════════════════════════

fit_homogeneous <- function(Y_win, X_win, params = list()) {
  n <- nrow(Y_win); TT <- ncol(Y_win); p <- dim(X_win)[3]

  y_vec <- c(Y_win)
  X_mat <- do.call(cbind, lapply(seq_len(p), function(j) c(X_win[, , j])))
  keep  <- !is.na(y_vec) & rowSums(is.na(X_mat)) == 0

  if (sum(keep) < p + 1)
    stop("fit_homogeneous: insufficient observations after NA removal")

  beta_common <- coef(lm.fit(X_mat[keep, , drop = FALSE], y_vec[keep]))
  beta_agg    <- matrix(rep(beta_common, n), nrow = n, byrow = TRUE)
  colnames(beta_agg) <- dimnames(X_win)[[3]]

  list(beta_agg  = beta_agg,
       group_vec = rep(1L, n),
       Gamma_hat = NULL,
       K_used    = 1L,
       method    = "Homogeneous_FM")
}

## ══════════════════════════════════════════════════════════════════════════════
##  METHOD 1  Static Heterogeneous Regression — per-firm OLS + K-means
##            (Tang et al. 2025 baseline)
## ══════════════════════════════════════════════════════════════════════════════

fit_static_hetero <- function(Y_win, X_win,
                              params = list(K_max = 5)) {
  n <- nrow(Y_win); TT <- ncol(Y_win); p <- dim(X_win)[3]
  K_max <- params$K_max %||% 5L

  # Step 1: per-firm time-pooled OLS
  beta_indiv <- matrix(NA_real_, n, p,
                       dimnames = list(NULL, dimnames(X_win)[[3]]))
  for (i in seq_len(n)) {
    act <- which(!is.na(Y_win[i, ]))
    if (length(act) < p + 1) next
    y_i <- Y_win[i, act]
    X_i <- matrix(X_win[i, act, ], nrow = length(act))
    if (any(!is.finite(X_i))) next
    fit <- tryCatch(lm.fit(X_i, y_i), error = function(e) NULL)
    if (!is.null(fit)) beta_indiv[i, ] <- coef(fit)
  }

  # Step 2: K-means on valid rows
  valid <- which(rowSums(is.na(beta_indiv)) == 0)

  if (length(valid) < 2) {
    beta_agg <- beta_indiv
    beta_agg[is.na(beta_agg)] <- 0
    return(list(beta_agg  = beta_agg,
                group_vec = rep(1L, n),
                Gamma_hat = NULL,
                K_used    = 1L,
                method    = "Static_Hetero"))
  }

  K_use     <- min(K_max, max(2L, floor(length(valid) / 5)))
  km        <- kmeans(beta_indiv[valid, ], centers = K_use,
                      nstart = 20, iter.max = 100)
  group_vec <- rep(1L, n)
  group_vec[valid] <- km$cluster

  # Nearest-centre assignment for firms with partial NAs
  invalid <- setdiff(seq_len(n), valid)
  for (i in invalid) {
    b <- replace(beta_indiv[i, ], is.na(beta_indiv[i, ]), 0)
    dists <- rowSums(sweep(km$centers, 2, b, "-")^2)
    group_vec[i] <- which.min(dists)
  }

  # Step 3: group-mean coefficients
  beta_agg <- matrix(0, n, p, dimnames = list(NULL, dimnames(X_win)[[3]]))
  for (k in seq_len(K_use)) {
    idx <- which(group_vec == k & rowSums(is.na(beta_indiv)) == 0)
    if (length(idx) > 0)
      beta_agg[group_vec == k, ] <-
        matrix(colMeans(beta_indiv[idx, , drop = FALSE]),
               nrow = sum(group_vec == k), byrow = TRUE)
  }

  list(beta_agg  = beta_agg,
       group_vec = group_vec,
       Gamma_hat = NULL,
       K_used    = K_use,
       method    = "Static_Hetero")
}

## ══════════════════════════════════════════════════════════════════════════════
##  METHOD 2  AFEM — dynamic heterogeneous regression (main proposal)
##            Spline group trajectories + adaptive fusion
## ══════════════════════════════════════════════════════════════════════════════

fit_afem_empirical <- function(Y_win, X_win,
                               params = list(K_max    = 5,
                                             df_spline = 6,
                                             lam1     = 0.1,
                                             lam2     = 0.05,
                                             n_agg    = 12,
                                             max_iter = 30,
                                             tol      = 1e-4)) {
  n  <- nrow(Y_win); TT <- ncol(Y_win); p <- dim(X_win)[3]
  K_max    <- params$K_max     %||% 5L
  df_sp    <- params$df_spline %||% 6L
  lam1     <- params$lam1      %||% 0.1
  lam2     <- params$lam2      %||% 0.05
  n_agg    <- min(params$n_agg %||% 12L, TT)
  max_iter <- params$max_iter  %||% 30L
  tol      <- params$tol       %||% 1e-4

  B <- make_bspline_basis(TT, df = df_sp)   # TT × df_sp

  ## Step 1 — ridge shrinkage estimate per (firm, time) → beta_indiv: n×TT×p
  beta_indiv <- array(NA_real_, c(n, TT, p))
  for (i in seq_len(n)) {
    for (t in seq_len(TT)) {
      if (is.na(Y_win[i, t])) next
      x_it <- X_win[i, t, ]
      if (any(is.na(x_it))) next
      beta_indiv[i, t, ] <- as.numeric(
        solve(tcrossprod(matrix(x_it, p, 1)) + diag(lam1, p),
              x_it * Y_win[i, t])
      )
    }
  }

  ## Step 2 — initialise K groups via K-means on time-averaged betas
  beta_mean <- apply(beta_indiv, c(1, 3), mean, na.rm = TRUE)
  valid     <- which(rowSums(is.na(beta_mean)) == 0)
  K_use     <- min(K_max, max(2L, floor(length(valid) / 5)))
  group_vec <- rep(1L, n)
  if (length(valid) >= K_use) {
    km <- kmeans(beta_mean[valid, , drop = FALSE], centers = K_use,
                 nstart = 20, iter.max = 100)
    group_vec[valid] <- km$cluster
  }

  ## Step 3 — alternating: spline fit per group ↔ subject re-assignment
  Gamma <- array(0, c(K_use, TT, p))      # K × TT × p  (smooth trajectories)
  Theta <- array(0, c(K_use, df_sp, p))   # K × df × p  (spline coefficients)

  for (iter in seq_len(max_iter)) {
    Gamma_old <- Gamma

    ## 3a. Spline fit per (group, covariate)
    for (k in seq_len(K_use)) {
      idx_k <- which(group_vec == k)
      if (length(idx_k) == 0) next
      for (j in seq_len(p)) {
        y_kj <- c(beta_indiv[idx_k, , j])            # (n_k · TT) vector
        B_kj <- kronecker(rep(1, length(idx_k)), B)  # (n_k · TT) × df
        keep  <- !is.na(y_kj)
        if (sum(keep) < df_sp) next
        coef_kj <- solve(crossprod(B_kj[keep, ]) + diag(lam2, df_sp),
                         crossprod(B_kj[keep, ], y_kj[keep]))
        Theta[k, , j] <- coef_kj
        Gamma[k, , j] <- B %*% coef_kj
      }
    }

    ## 3b. Adaptive fusion — merge groups whose final-period centres are close
    if (iter > 3 && K_use > 1) {
      centres   <- matrix(Gamma[, TT, ], nrow = K_use)   # K × p
      pair_dist <- as.matrix(dist(centres))
      fuse_thr  <- lam1 * sqrt(p)
      membership <- seq_len(K_use)
      for (k1 in seq_len(K_use - 1))
        for (k2 in (k1 + 1):K_use)
          if (pair_dist[k1, k2] < fuse_thr)
            membership[membership == membership[k2]] <- membership[k1]
      new_labels <- as.integer(factor(membership))
      K_new      <- max(new_labels)
      if (K_new < K_use) {
        Gamma_new <- array(0, c(K_new, TT, p))
        Theta_new <- array(0, c(K_new, df_sp, p))
        for (k in seq_len(K_new)) {
          old_ks <- which(new_labels == k)
          Gamma_new[k, , ] <- apply(Gamma[old_ks, , , drop = FALSE], c(2, 3), mean)
          Theta_new[k, , ] <- apply(Theta[old_ks, , , drop = FALSE], c(2, 3), mean)
        }
        K_use     <- K_new
        group_vec <- new_labels[group_vec]
        Gamma     <- Gamma_new
        Theta     <- Theta_new
      }
    }

    ## 3c. Re-assign each firm to nearest group trajectory
    for (i in seq_len(n)) {
      act <- which(!is.na(beta_indiv[i, , 1]))
      if (length(act) == 0) next
      bi <- beta_indiv[i, act, , drop = FALSE]   # 1 × |act| × p -> act × p
      dists <- vapply(seq_len(K_use), function(k) {
        sum((Gamma[k, act, , drop = FALSE] - bi)^2, na.rm = TRUE)
      }, numeric(1))
      group_vec[i] <- which.min(dists)
    }

    ## Convergence check
    if (max(abs(Gamma - Gamma_old), na.rm = TRUE) < tol && iter > 5) break
  }

  ## Step 4 — aggregate group trajectory over last n_agg months
  t_agg    <- max(1L, TT - n_agg + 1L):TT
  beta_agg <- matrix(NA_real_, n, p, dimnames = list(NULL, dimnames(X_win)[[3]]))
  for (i in seq_len(n))
    beta_agg[i, ] <- colMeans(Gamma[group_vec[i], t_agg, , drop = FALSE])
  beta_agg[is.na(beta_agg)] <- 0

  list(beta_agg  = beta_agg,
       group_vec = group_vec,
       Gamma_hat = Gamma,
       K_used    = K_use,
       method    = "AFEM")
}

## ══════════════════════════════════════════════════════════════════════════════
##  METHOD 3  SplineKM — individual spline fit + K-means (benchmark)
## ══════════════════════════════════════════════════════════════════════════════

fit_splinekm_empirical <- function(Y_win, X_win,
                                   params = list(K_max    = 5,
                                                 df_spline = 6,
                                                 n_agg    = 12)) {
  n  <- nrow(Y_win); TT <- ncol(Y_win); p <- dim(X_win)[3]
  K_max <- params$K_max     %||% 5L
  df_sp <- params$df_spline %||% 6L
  n_agg <- min(params$n_agg %||% 12L, TT)

  B <- make_bspline_basis(TT, df = df_sp)

  # Per-firm, per-covariate spline on y*x projection
  gamma_curves <- array(NA_real_, c(n, TT, p))
  for (i in seq_len(n)) {
    for (j in seq_len(p)) {
      y_ij <- Y_win[i, ] * X_win[i, , j]
      act  <- which(!is.na(y_ij))
      if (length(act) < df_sp) next
      coef_ij <- tryCatch(
        solve(crossprod(B[act, ]) + diag(0.01, df_sp),
              crossprod(B[act, ], y_ij[act])),
        error = function(e) rep(0, df_sp)
      )
      gamma_curves[i, , j] <- B %*% coef_ij
    }
  }

  # K-means on end-of-window values
  feat_mat  <- gamma_curves[, TT, ]                  # n × p
  valid     <- which(rowSums(is.na(feat_mat)) == 0)
  K_use     <- min(K_max, max(2L, floor(length(valid) / 5)))
  group_vec <- rep(1L, n)
  if (length(valid) >= K_use) {
    km <- kmeans(feat_mat[valid, , drop = FALSE], centers = K_use, nstart = 20)
    group_vec[valid] <- km$cluster
  }

  # Time-aggregate over last n_agg months
  t_agg    <- max(1L, TT - n_agg + 1L):TT
  beta_agg <- matrix(0, n, p, dimnames = list(NULL, dimnames(X_win)[[3]]))
  for (i in seq_len(n)) {
    vals <- gamma_curves[i, t_agg, ]
    if (!anyNA(vals)) beta_agg[i, ] <- colMeans(vals)
  }

  list(beta_agg  = beta_agg,
       group_vec = group_vec,
       Gamma_hat = NULL,
       K_used    = K_use,
       method    = "SplineKM")
}
