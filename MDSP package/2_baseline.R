# ============================================================
# baselines.R
#
# Baseline methods for the simulation study:
#   "Heterogeneity-Aware Statistical Learning"
#
# Extended model with HOMOGENEOUS (Z) and HETEROGENEOUS (X) covariates:
#
#   Y_it = z_it' * alpha_t  +  x_it' * beta_t^(k)  +  epsilon_it
#
#   Z      : n x m x r  covariate array, SHARED time-varying coefficient
#   X      : n x m x p  covariate array, GROUP-SPECIFIC spline coefficient
#   alpha  : m x r matrix  — homogeneous coefficient at each time point t
#            (same across all groups; each row alpha[t,] applies at time t)
#
#   When Z = NULL the model reduces to the original heterogeneous-only version.
#
# Strategy applied uniformly across all methods:
#   1. Estimate alpha (m x r) via m independent per-time-point ridge regressions,
#      each pooling all subjects observed at that time (ignoring group structure).
#   2. Residualize: R_it = Y_it - z_it' * alpha_hat[t, ]
#   3. Run the original heterogeneous method on R in place of Y.
#   4. Return alpha_hat (m x r) alongside the standard output fields.
#
# METHOD 1 — FlexMix (mixture-of-regressions)
#   fit_flexmix()
#
# METHOD 2 — Individualized Spline + coefficient k-means (bandwidth window)
#   fit_indiv_spline()
#
# METHOD 3 — Regional Spline + per-region coefficient k-means
#   fit_regional_spline()
#
# METHOD 4 — Individualized Spline + fitted-value k-means (bandwidth window)
#   fit_indiv_spline_fv()
#
# METHOD 5 — Regional Spline + per-region fitted-value k-means
#   fit_regional_spline_fv()
#
# METHOD 7 — Per-time-point heterogeneous regression (Tang et al. style)
# METHOD 8 — Kernel Ridge Pilot + k-means (Steps 1-2 of HetRegL1, no ADMM)
#             L1 fusion penalty solved via ADMM at each t independently;
#             subject labels aggregated across time via majority vote.
#   fit_het_reg_l1()
#
# METHOD 9 — Discrete-time FlexMix (time-wise mixture regression)
#             At each t independently, fit a K-component mixture-of-regressions
#             via FlexMix; align labels across time with sequential Hungarian
#             matching on co-occurrence tables.
#   fit_flex_timewise()
#
# All functions share a common interface:
#   Input  : Y (n x m), X (n x m x p), Z (n x m x r | NULL), K, df, ...
#   Output : list(subject_group, gamma_hat [K x m], group_matrix [n x m],
#                 alpha [m x r] | NULL, time_sec, ...)
# ============================================================

library(splines)

.load_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
.load_pkg("flexmix")
.load_pkg("mclust")

# ============================================================
# SHARED INTERNAL HELPERS — BASIS AND FEATURE CONSTRUCTION
# ============================================================

# Kronecker feature: phi_it = x_it (x) b_t  (length p*q)
.phi <- function(x_it, b_t) as.numeric(outer(x_it, b_t))

# Build long-format data.frame from (residualized) Y and X.
# Rows with any NA in Y or X are dropped.
# Columns: y, phi_1 ... phi_{pq}, subj, tpt
.build_long_df <- function(Y, XZ, B) {
  # XZ is the combined covariate array (n x m x p_total)
  n <- nrow(Y); m <- ncol(Y); p <- dim(XZ)[3]; q <- ncol(B)
  feat_names <- paste0("phi_", seq_len(p * q))
  rows <- vector("list", n * m)
  idx  <- 1L
  for (i in seq_len(n))
    for (t in seq_len(m)) {
      if (is.na(Y[i, t]) || any(is.na(XZ[i, t, ]))) next
      rows[[idx]] <- c(y = Y[i, t],
                       setNames(.phi(XZ[i, t, ], B[t, ]), feat_names),
                       subj = i, tpt = t)
      idx <- idx + 1L
    }
  as.data.frame(do.call(rbind, rows[seq_len(idx - 1L)]),
                stringsAsFactors = FALSE)
}

# B-spline knot attributes fixed from a reference grid.
# Storing knots separately ensures the same spline space is used
# when evaluating at new single points (avoids basis mismatch).
.bs_knots <- function(u_ref, df) {
  b <- bs(u_ref, df = df, intercept = TRUE)
  list(knots          = attr(b, "knots"),
       boundary_knots = attr(b, "Boundary.knots"),
       degree         = attr(b, "degree"))
}

# Evaluate a B-spline at u_new using precomputed knot attributes.
.bs_eval <- function(u_new, bk) {
  bs(u_new,
     knots          = bk$knots,
     Boundary.knots = bk$boundary_knots,
     degree         = bk$degree,
     intercept      = TRUE)
}

# ============================================================
# SHARED INTERNAL HELPERS — Z + X COMBINED COVARIATE STRATEGY
# ============================================================
#
# Baseline estimation strategy: treat Z as heterogeneous
# -----------------------------------------------------------------------
# Rather than estimating a separate homogeneous alpha, we treat Z on an
# equal footing with X and concatenate them into a single (p + r)-dimensional
# heterogeneous covariate array:
#
#   XZ[i, t, ] = c( XZ[i, t, ],  Z[i, t, ] )   # length p + r
#
# Each baseline method then fits its group-specific coefficient model on
# Y and XZ with no separate Z component (Z_arg = NULL).  The group-specific
# coefficients for the first p dimensions correspond to X; the last r
# dimensions correspond to Z (treated as group-specific here).
#
# After fitting, the "alpha" field returned by each method contains the
# group-mean of the Z-column coefficients, averaged across groups weighted
# by group size.  This gives a summary comparable to the true shared alpha,
# though it is no longer constrained to be identical across groups.
#
# When Z is NULL, XZ = X and the methods reduce to the original
# heterogeneous-only version.

#' Concatenate X and Z into a single (p + r)-dimensional covariate array.
#' When Z is NULL, returns X unchanged.
#'
#' @param X  n x m x p array
#' @param Z  n x m x r array  (or NULL)
#' @return   n x m x (p + r) array  (or X when Z is NULL)
.combine_ZX <- function(X, Z) {
  if (is.null(Z)) return(X)
  n <- dim(X)[1]; m <- dim(X)[2]; p <- dim(X)[3]; r <- dim(Z)[3]
  XZ <- array(NA_real_, dim = c(n, m, p + r))
  XZ[, , seq_len(p)]         <- X
  XZ[, , p + seq_len(r)]     <- Z
  XZ
}

#' Extract a group-weighted mean of the Z-column coefficients from gamma_hat.
#' gamma_hat is K x m; the last r rows of the covariate were Z.
#' Returns an m x r matrix comparable to the true alpha (or NULL if Z=NULL).
#'
#' @param gamma_hat  K x m matrix of group trajectory scalars
#'   (NOTE: this is a scalar summary; for the full coefficient extraction
#'    use the Theta matrices directly when available)
#' @param p_x      number of X covariates
#' @param r        number of Z covariates
#' @param pi_vec   length-K group proportions (for weighted average)
#' @return NULL (scalar gamma_hat does not carry per-covariate info)
#'   Full alpha extraction requires access to Theta, handled per-method.
.extract_alpha_from_combined <- function(gamma_hat, p_x, r, pi_vec = NULL) {
  # gamma_hat is a scalar K x m summary (L2 norm or first component),
  # not sufficient to recover the r-dimensional alpha vector.
  # Return NULL and let calling code handle alpha separately if needed.
  NULL
}

#' Standard preamble for all baseline methods.
#' Combines X and Z, returns the combined array XZ and alpha = NULL.
#' Each method receives XZ as its covariate and Z_arg = NULL.
#'
#' @param Y        n x m outcome
#' @param X        n x m x p heterogeneous covariate
#' @param Z        n x m x r homogeneous covariate (or NULL)
#' @return list(XZ = combined n x m x (p+r) array,
#'              p_x = original X dimension,
#'              r   = Z dimension (0 if Z=NULL))
.prep_combined <- function(Y, X, Z) {
  p_x <- dim(X)[3]
  r   <- if (is.null(Z)) 0L else dim(Z)[3]
  XZ  <- .combine_ZX(X, Z)
  list(XZ = XZ, p_x = p_x, r = r)
}

# ============================================================
# SHARED INTERNAL HELPERS — LABEL ALIGNMENT
# ============================================================

#' Align a set of cluster labels cur to a reference ref_col.
#' K=2: bit-flip if it improves agreement.
#' K>2: Hungarian algorithm via clue::solve_LSAP.
.align_to_ref <- function(cur, ref_col, K) {
  if (K == 2L) {
    if (mean(cur == ref_col, na.rm = TRUE) <
        mean((3L - cur) == ref_col, na.rm = TRUE))
      cur <- 3L - cur
    return(cur)
  }
  .load_pkg("clue")
  tab  <- table(factor(ref_col, 1:K), factor(cur, 1:K))
  perm <- clue::solve_LSAP(max(tab) - tab)
  perm[cur]
}

#' Majority-vote aggregation across region label columns (n x n_region matrix).
.majority_vote <- function(aligned_labels, K) {
  apply(aligned_labels, 1L, function(row) {
    tbl <- tabulate(row[!is.na(row)], nbins = K)
    if (all(tbl == 0L)) return(1L)
    as.integer(which.max(tbl))
  })
}

.canonicalise_components <- function(gamma_hat, subject_group,
                                     group_matrix, K) {
  #' Sort K components into a canonical ascending order by mean trajectory.
  #' Permutes gamma_hat rows and remaps subject_group / group_matrix labels
  #' so that component 1 always has the lowest mean trajectory value and
  #' component K the highest.  This makes the output deterministic regardless
  #' of the random-restart ordering inside FlexMix or k-means.
  #'
  #' @param gamma_hat     K x m matrix
  #' @param subject_group integer vector length n
  #' @param group_matrix  n x m integer matrix
  #' @param K             number of components
  #' @return list(gamma_hat, subject_group, group_matrix) — permuted
  if (is.null(gamma_hat) || all(is.na(gamma_hat)))
    return(list(gamma_hat     = gamma_hat,
                subject_group = subject_group,
                group_matrix  = group_matrix,
                perm          = seq_len(K)))

  mean_traj <- rowMeans(gamma_hat, na.rm = TRUE)
  ord       <- order(mean_traj)          # ascending by mean trajectory
  if (identical(ord, seq_len(K))) {      # already canonical — nothing to do
    return(list(gamma_hat     = gamma_hat,
                subject_group = subject_group,
                group_matrix  = group_matrix,
                perm          = seq_len(K)))
  }

  # inv_ord[k] = new label for old component k
  inv_ord <- integer(K)
  inv_ord[ord] <- seq_len(K)

  gamma_hat_new     <- gamma_hat[ord, , drop = FALSE]
  subject_group_new <- inv_ord[subject_group]
  group_matrix_new  <- matrix(inv_ord[as.integer(group_matrix)],
                              nrow = nrow(group_matrix),
                              ncol = ncol(group_matrix))

  list(gamma_hat     = gamma_hat_new,
       subject_group = as.integer(subject_group_new),
       group_matrix  = group_matrix_new,
       perm          = ord)   # ord[k] = which old component maps to new position k
}

# ============================================================
# METHOD 1: FLEXMIX
# ============================================================

#' Fit FlexMix mixture-of-regressions with optional homogeneous Z.
#'
#' When Z is supplied:
#'   1. alpha estimated via pooled ridge-WLS on Z.
#'   2. FlexMix fitted on residual R = Y - Z * alpha_hat.
#'   3. alpha_hat returned in output list.
#'
#' @param Y          n x m outcome matrix  (NAs allowed)
#' @param X          n x m x p heterogeneous covariate array  (NAs allowed)
#' @param Z          n x m x r homogeneous covariate array (NULL = omit)
#' @param K          number of mixture components
#' @param df         B-spline degrees of freedom for time
#' @param lambda_z   ridge penalty for alpha  (default 1e-4)
#' @param iter_max   max EM iterations
#' @param n_start    random restarts (uses stepFlexmix when > 1)
#' @param min_prior  component pruning threshold
#'
#' @return list(subject_group, gamma_hat [K x m], group_matrix [n x m],
#'              alpha [r] | NULL, time_sec, fit)
fit_flexmix <- function(Y, X, Z = NULL, K = 2L, df = 6L,
                        lambda_z  = 1e-4,
                        iter_max  = 200L,
                        n_start   = 3L,
                        min_prior = 0.05) {

  n <- nrow(Y); m <- ncol(Y)

  # Combined XZ covariate
  prep  <- .prep_combined(Y, X, Z)
  XZ    <- prep$XZ        # combined (p + r)-dim covariate
  p_x   <- prep$p_x       # original X columns (1..p_x)
  p     <- dim(XZ)[3]     # total covariate dim = p_x + r
  alpha <- NULL            # Z treated heterogeneously; no shared alpha

  # Global B-spline basis on [0,1]
  t_seq <- seq(0, 1, length.out = m)
  B     <- bs(t_seq, df = df, intercept = TRUE)
  q     <- ncol(B);  pq <- p * q

  df_long    <- .build_long_df(Y, XZ, B)
  feat_names <- paste0("phi_", seq_len(pq))
  fmla       <- as.formula(
    paste("y ~", paste(feat_names, collapse = " + "), "- 1")
  )

  na_result <- list(subject_group = rep(1L, n),
                    gamma_hat     = matrix(NA_real_, K, m),
                    group_matrix  = matrix(1L, n, m),
                    alpha         = alpha,
                    time_sec      = NA_real_,
                    fit           = NULL)

  if (nrow(df_long) < K * (pq + 1L)) {
    warning("fit_flexmix: insufficient observations after NA removal.")
    return(na_result)
  }

  fm_ctrl <- list(iter.max  = iter_max,
                  minprior  = min_prior,
                  tolerance = 1e-5)

  t0  <- proc.time()["elapsed"]
  fit <- tryCatch({
    if (n_start > 1L)
      flexmix::stepFlexmix(fmla, data = df_long, k = K,
                           nrep = n_start, control = fm_ctrl)
    else
      flexmix::flexmix(fmla, data = df_long, k = K, control = fm_ctrl)
  }, error = function(e) {
    warning(sprintf("fit_flexmix: flexmix error — %s", conditionMessage(e)))
    NULL
  })
  t1 <- proc.time()["elapsed"]

  if (is.null(fit)) return(modifyList(na_result, list(time_sec = t1 - t0)))
  if (inherits(fit, "stepFlexmix")) fit <- flexmix::getModel(fit, "BIC")

  # Subject-level group assignment (majority vote over time points)
  cls_long          <- flexmix::clusters(fit)
  df_long$cls       <- cls_long
  subj_tbl          <- tapply(df_long$cls, df_long$subj,
                              function(x) as.integer(names(which.max(table(x)))))
  subject_group     <- rep(1L, n)
  valid_subj        <- as.integer(names(subj_tbl))
  subject_group[valid_subj] <- as.integer(subj_tbl)

  # group_matrix: FlexMix gives subject-level (time-invariant) labels
  group_matrix <- matrix(subject_group, nrow = n, ncol = m)

  # gamma_hat from fitted regression coefficients
  params    <- flexmix::parameters(fit)
  gamma_hat    <- matrix(NA_real_, K, m)
  Gamma_hat <- array(NA_real_, dim = c(K, m, p_x))  # K x m x p_x estimated group trajectories

  for (k in seq_len(K)) {
    raw    <- if (is.list(params)) as.numeric(params[[k]]) else
      as.numeric(params[, k])
    coef_k <- raw[seq_len(min(pq, length(raw)))]
    if (length(coef_k) < pq) {
      warning(sprintf("fit_flexmix: component %d got %d coefs, expected %d",
                      k, length(coef_k), pq))
      next
    }
    Theta_k  <- matrix(coef_k, nrow = p, ncol = q, byrow = FALSE)
    Theta_kx <- Theta_k[seq_len(p_x), , drop = FALSE]  # X-block only
    for (t in seq_len(m)) {
      gvec            <- as.numeric(Theta_kx %*% B[t, ])   # p_x-vector
      gamma_hat[k, t]      <- if (p_x == 1L) gvec[1L] else sqrt(sum(gvec^2))
      Gamma_hat[k, t, ] <- gvec
    }
  }

  # Canonicalise component ordering: sort by ascending mean trajectory
  can <- .canonicalise_components(gamma_hat, subject_group, group_matrix, K)
  Gamma_hat <- Gamma_hat[can$perm, , , drop = FALSE]
  # Broadcast K x m x p -> n x m x p
  beta_hat_arr <- Gamma_hat[can$subject_group, , , drop = FALSE]

  list(subject_group = can$subject_group,
       gamma_hat     = can$gamma_hat,
       Gamma_hat     = Gamma_hat,          # K x m x p_x estimated group trajectories
       beta_hat_arr  = beta_hat_arr,
       group_matrix  = can$group_matrix,
       alpha         = alpha,
       time_sec      = t1 - t0,
       fit           = fit)
}

# ============================================================
# METHOD 2: INDIVIDUALIZED SPLINE + COEFFICIENT k-MEANS
#           (bandwidth windows, global clustering)
# ============================================================

#' Fit individualized B-splines with bandwidth window control, cluster on
#' concatenated window coefficient vectors, with optional homogeneous Z.
#'
#' For each subject i and window center t_c, a ridge-penalized WLS spline
#' is fitted using only observations inside [t_c - delta, t_c + delta].
#' Coefficient vectors are concatenated across windows into a subject feature
#' vector (length n_win * pq) and clustered via k-means.
#'
#' @param Y          n x m outcome matrix  (NAs allowed)
#' @param X          n x m x p heterogeneous covariate array  (NAs allowed)
#' @param Z          n x m x r homogeneous covariate array (NULL = omit)
#' @param K          number of groups
#' @param df         B-spline df per window
#' @param delta      half-bandwidth in [0,1] units  (default 0.20)
#' @param lambda     ridge penalty for het spline  (default 1e-3)
#' @param lambda_z   ridge penalty for alpha  (default 1e-4)
#' @param km_nstart  k-means restarts  (default 20)
#'
#' @return list(subject_group, gamma_hat [K x m], group_matrix [n x m],
#'              alpha [r] | NULL, time_sec, t_centers, delta, n_win)
fit_indiv_spline <- function(Y, X, Z = NULL, K = 2L, df = 6L,
                             lambda    = 1e-3,
                             lambda_z  = 1e-4,
                             km_nstart = 20L) {

  n        <- nrow(Y); m <- ncol(Y)
  t_global <- seq(0, 1, length.out = m)

  prep  <- .prep_combined(Y, X, Z)
  XZ    <- prep$XZ
  p_x   <- prep$p_x
  p     <- dim(XZ)[3]
  alpha <- NULL
  q     <- df; pq <- p * q

  t0_clock <- proc.time()["elapsed"]

  # ── Global reference knots from all observed time points ─────────────
  # Use the full t_global grid so the knot structure is shared across
  # subjects; the ridge penalty regularises sparsely-observed subjects.
  global_bk <- .bs_knots(t_global, df)

  # ── Per-subject global ridge spline ───────────────────────────────────
  # For each subject fit one ridge regression on their observed (t, XZ, Y)
  # triples using the global B-spline basis. Missing (i,t) are skipped.
  coef_mat <- matrix(NA_real_, nrow = n, ncol = pq)

  for (i in seq_len(n)) {
    obs_i <- which(!is.na(Y[i, ]) &
                     rowSums(is.na(matrix(XZ[i, , ], nrow = m, ncol = p))) == 0L)
    if (length(obs_i) < df + 1L) next       # too few observations

    B_i      <- .bs_eval(t_global[obs_i], global_bk)   # |obs_i| x q
    phi_rows <- lapply(obs_i,
                       function(t) .phi(XZ[i, t, ], B_i[match(t, obs_i), ]))
    Des_i    <- do.call(rbind, phi_rows)
    ys_i     <- Y[i, obs_i]
    pen      <- lambda * diag(pq)
    coef_mat[i, ] <- tryCatch(
      as.numeric(solve(crossprod(Des_i) + pen, t(Des_i) %*% ys_i)),
      error = function(e) rep(NA_real_, pq)
    )
  }

  valid <- which(apply(coef_mat, 1, function(r) all(is.finite(r))))
  if (length(valid) < K) {
    warning("fit_indiv_spline: too few valid subjects for k-means.")
    return(list(subject_group = rep(1L, n),
                gamma_hat     = matrix(NA_real_, K, m),
                Gamma_hat     = array(NA_real_, dim = c(K, m, p_x)),
                beta_hat_arr  = array(NA_real_, dim = c(n, m, p_x)),
                group_matrix  = matrix(1L, n, m),
                alpha         = alpha,
                time_sec      = proc.time()["elapsed"] - t0_clock))
  }

  # ── k-means on coefficient vectors ────────────────────────────────────
  km <- tryCatch(
    kmeans(coef_mat[valid, , drop = FALSE],
           centers = K, nstart = km_nstart, iter.max = 100L),
    error = function(e) list(cluster = rep(1L, length(valid)))
  )
  subject_group           <- rep(NA_integer_, n)
  subject_group[valid]    <- as.integer(km$cluster)
  subject_group[is.na(subject_group)] <- 1L
  group_matrix <- matrix(subject_group, nrow = n, ncol = m)

  # ── Group trajectories: mean Theta_k evaluated on global grid ─────────
  B_global <- .bs_eval(t_global, global_bk)   # m x q

  gamma_hat <- matrix(NA_real_, K, m)
  Gamma_hat <- array(NA_real_, dim = c(K, m, p_x))

  for (k in seq_len(K)) {
    idx_k <- which(subject_group == k)
    if (length(idx_k) == 0L) next
    mean_coef <- colMeans(coef_mat[idx_k, , drop = FALSE], na.rm = TRUE)
    if (any(!is.finite(mean_coef))) next
    Theta_k  <- matrix(mean_coef, nrow = p, ncol = q, byrow = FALSE)
    Theta_kx <- Theta_k[seq_len(p_x), , drop = FALSE]
    for (t in seq_len(m)) {
      gvec <- as.numeric(Theta_kx %*% B_global[t, ])
      gamma_hat[k, t]   <- if (p_x == 1L) gvec[1L] else sqrt(sum(gvec^2))
      Gamma_hat[k, t, ] <- gvec
    }
  }

  # Broadcast K x m x p_x -> n x m x p_x
  beta_hat_arr <- Gamma_hat[subject_group, , , drop = FALSE]

  list(subject_group = subject_group,
       gamma_hat     = gamma_hat,
       Gamma_hat     = Gamma_hat,
       beta_hat_arr  = beta_hat_arr,
       group_matrix  = group_matrix,
       alpha         = alpha,
       time_sec      = proc.time()["elapsed"] - t0_clock)
}

# ============================================================
# METHOD 3: REGIONAL SPLINE + PER-REGION COEFFICIENT k-MEANS
# ============================================================

#' Fit individualized B-splines on hard non-overlapping sub-regions,
#' cluster subjects within each region on coefficient vectors, then
#' aggregate via majority vote, with optional homogeneous Z.
#'
#' Algorithm:
#'   1. Partition [0,1] into n_region equal-width intervals.
#'   2. Per region r, per subject i: fit ridge-WLS B-spline on obs in region r.
#'   3. Per region r: k-means on n x pq coefficient matrix -> region label z_r.
#'   4. Align region labels (Hungarian / flip) then majority vote -> subject_group.
#'   5. Group trajectories: per region r, per group k: mean Theta_k^(r) evaluated
#'      on the global time grid.
#'
#' @param Y          n x m outcome matrix  (NAs allowed)
#' @param X          n x m x p heterogeneous covariate array  (NAs allowed)
#' @param Z          n x m x r homogeneous covariate array (NULL = omit)
#' @param K          number of groups
#' @param df         B-spline df per region
#' @param n_region   number of hard sub-regions  (NULL = auto ~5)
#' @param lambda     ridge penalty for het spline  (default 1e-3)
#' @param lambda_z   ridge penalty for alpha  (default 1e-4)
#' @param km_nstart  k-means restarts per region  (default 20)
#' @param min_obs    min valid obs in a region to fit spline  (NULL = df+1)
#'
#' @return list(subject_group, gamma_hat [K x m], group_matrix [n x m],
#'              alpha [r] | NULL, region_labels, region_breaks, time_sec, n_region)
fit_regional_spline <- function(Y, X, Z = NULL,
                                K         = 2L,
                                df        = 6L,
                                n_region  = NULL,
                                lambda    = 1e-3,
                                lambda_z  = 1e-4,
                                km_nstart = 20L,
                                min_obs   = NULL) {

  n <- nrow(Y); m <- ncol(Y)
  t_global <- seq(0, 1, length.out = m)
  if (is.null(min_obs))  min_obs  <- df + 1L
  if (is.null(n_region)) n_region <- max(3L, min(5L, floor(m / (df + 2L))))

  prep  <- .prep_combined(Y, X, Z)
  XZ    <- prep$XZ        # combined (p + r)-dim covariate
  p_x   <- prep$p_x       # original X columns (1..p_x)
  p     <- dim(XZ)[3]     # total covariate dim = p_x + r
  alpha <- NULL            # Z treated heterogeneously; no shared alpha
  q  <- df;  pq <- p * q

  region_breaks <- seq(0, 1, length.out = n_region + 1L)
  region_of_t   <- pmax(1L, pmin(
    findInterval(t_global, region_breaks, rightmost.closed = TRUE), n_region))

  t0_clock <- proc.time()["elapsed"]

  # Precompute per-region B-spline knots
  reg_bk <- vector("list", n_region)
  for (r in seq_len(n_region)) {
    t_in_r <- t_global[region_of_t == r]
    lo_r   <- region_breaks[r]; hi_r <- region_breaks[r + 1L]
    u_ref  <- pmax(pmin((t_in_r - lo_r) / (hi_r - lo_r), 1), 0)
    if (length(unique(u_ref)) < df + 1L)
      u_ref <- seq(0, 1, length.out = df + 2L)
    reg_bk[[r]] <- .bs_knots(u_ref, df)
  }

  # Per-subject, per-region coefficient estimation on residualized R
  coef_array <- array(0, dim = c(n, n_region, pq))

  for (r in seq_len(n_region)) {
    lo_r    <- region_breaks[r]; hi_r <- region_breaks[r + 1L]
    t_idx_r <- which(region_of_t == r)
    for (i in seq_len(n)) {
      obs_r <- t_idx_r[!is.na(Y[i, t_idx_r]) & rowSums(is.na(matrix(XZ[i, t_idx_r, ], nrow=length(t_idx_r), ncol=p))) == 0L]
      if (length(obs_r) < min_obs) next
      u_obs    <- pmax(pmin((t_global[obs_r] - lo_r) / (hi_r - lo_r), 1), 0)
      B_r      <- .bs_eval(u_obs, reg_bk[[r]])
      phi_rows <- lapply(seq_along(obs_r),
                         function(j) .phi(XZ[i, obs_r[j], ], B_r[j, ]))
      ys_r     <- Y[i, obs_r]
      Des_r    <- do.call(rbind, phi_rows)
      pen      <- lambda * diag(pq)
      coef_array[i, r, ] <- tryCatch(
        as.numeric(solve(crossprod(Des_r) + pen, t(Des_r) %*% ys_r)),
        error = function(e) rep(0, pq)
      )
    }
  }

  # Per-region k-means on coefficient vectors
  region_labels <- matrix(NA_integer_, nrow = n, ncol = n_region)
  for (r in seq_len(n_region)) {
    cm_r  <- coef_array[, r, , drop = TRUE]
    valid <- which(rowSums(!is.finite(cm_r)) == 0L)
    if (length(valid) < K) { region_labels[, r] <- 1L; next }
    km_r <- tryCatch(
      kmeans(cm_r[valid, , drop = FALSE], centers = K,
             nstart = km_nstart, iter.max = 200L),
      error = function(e) list(cluster = rep(1L, length(valid)))
    )
    region_labels[valid, r] <- as.integer(km_r$cluster)
    region_labels[is.na(region_labels[, r]), r] <- 1L
  }

  # Align region labels then majority vote
  aligned_labels         <- region_labels
  aligned_labels[, 1L]   <- region_labels[, 1L]
  if (n_region > 1L)
    for (r in 2L:n_region)
      aligned_labels[, r] <- .align_to_ref(region_labels[, r],
                                           region_labels[, 1L], K)

  subject_group <- .majority_vote(aligned_labels, K)

  # Group trajectories: per region, per group -> mean Theta -> evaluate
  gamma_hat    <- matrix(NA_real_, K, m)
  Gamma_hat <- array(NA_real_, dim = c(K, m, p_x))  # K x m x p_x estimated group trajectories
  for (r in seq_len(n_region)) {
    lo_r    <- region_breaks[r]; hi_r <- region_breaks[r + 1L]
    t_idx_r <- which(region_of_t == r)
    for (k in seq_len(K)) {
      idx_k     <- which(subject_group == k)
      if (length(idx_k) == 0L) next
      mean_coef <- colMeans(coef_array[idx_k, r, , drop = FALSE], na.rm = TRUE)
      if (any(!is.finite(mean_coef))) next
      Theta_k  <- matrix(mean_coef, nrow = p, ncol = q, byrow = FALSE)
      Theta_kx <- Theta_k[seq_len(p_x), , drop = FALSE]  # X-block only
      for (t in t_idx_r) {
        u_t      <- pmax(pmin((t_global[t] - lo_r) / (hi_r - lo_r), 1), 0)
        b_t      <- as.numeric(.bs_eval(u_t, reg_bk[[r]]))
        gvec     <- as.numeric(Theta_kx %*% b_t)            # p_x-vector
        gamma_hat[k, t]      <- if (p_x == 1L) gvec[1L] else sqrt(sum(gvec^2))
        Gamma_hat[k, t, ] <- gvec
      }
    }
  }

  # group_matrix: region-specific labels broadcast to time points
  group_matrix <- matrix(NA_integer_, nrow = n, ncol = m)
  for (r in seq_len(n_region))
    group_matrix[, which(region_of_t == r)] <- aligned_labels[, r]

  # Broadcast K x m x p -> n x m x p
  beta_hat_arr <- Gamma_hat[subject_group, , , drop = FALSE]
  list(subject_group  = subject_group,
       gamma_hat      = gamma_hat,
       Gamma_hat      = Gamma_hat,         # K x m x p_x estimated group trajectories
       beta_hat_arr   = beta_hat_arr,
       group_matrix   = group_matrix,
       alpha          = alpha,
       region_labels  = aligned_labels,
       region_breaks  = region_breaks,
       time_sec       = proc.time()["elapsed"] - t0_clock,
       n_region       = n_region)
}

# ============================================================
# METHOD 4: INDIVIDUALIZED SPLINE + FITTED-VALUE k-MEANS
#           (bandwidth windows, global clustering)
# ============================================================

#' Fit per-subject B-splines in bandwidth windows, evaluate fitted values on
#' the global m-point grid, then run a single k-means on the n x m
#' fitted-value matrix, with optional homogeneous Z.
#'
#' Key difference from Method 2: the clustering feature space is the
#' fitted-value curve (length m) rather than concatenated coefficient
#' vectors (length pq). Global B-spline across all observed time points.
#' geometry-consistent — two subjects with similar trajectories but
#' different coefficient representations will be clustered together.
#'
#' @param Y          n x m outcome matrix  (NAs allowed)
#' @param X          n x m x p heterogeneous covariate array  (NAs allowed)
#' @param Z          n x m x r homogeneous covariate array (NULL = omit)
#' @param K          number of groups
#' @param df         B-spline df per window
#' @param delta      half-bandwidth in [0,1] units  (default 0.20)
#' @param lambda     ridge penalty for het spline  (default 1e-3)
#' @param lambda_z   ridge penalty for alpha  (default 1e-4)
#' @param km_nstart  k-means restarts  (default 20)
#'
#' @return list(subject_group, gamma_hat [K x m], group_matrix [n x m],
#'              alpha [r] | NULL, fitted_values [n x m],
#'              time_sec)
fit_indiv_spline_fv <- function(Y, X, Z = NULL, K = 2L, df = 6L,
                                lambda    = 1e-3,
                                lambda_z  = 1e-4,
                                km_nstart = 20L) {

  n        <- nrow(Y); m <- ncol(Y)
  t_global <- seq(0, 1, length.out = m)

  prep  <- .prep_combined(Y, X, Z)
  XZ    <- prep$XZ
  p_x   <- prep$p_x
  p     <- dim(XZ)[3]
  alpha <- NULL
  q     <- df; pq <- p * q

  t0_clock <- proc.time()["elapsed"]

  # ── Global reference knots ────────────────────────────────────────────
  global_bk <- .bs_knots(t_global, df)
  B_global  <- .bs_eval(t_global, global_bk)   # m x q

  # ── Per-subject global ridge spline -> fitted values ──────────────────
  fitted_values <- matrix(NA_real_, nrow = n, ncol = m)
  beta_vec_arr  <- array(NA_real_,  dim  = c(n, m, p_x))

  for (i in seq_len(n)) {
    obs_i <- which(!is.na(Y[i, ]) &
                     rowSums(is.na(matrix(XZ[i, , ], nrow = m, ncol = p))) == 0L)
    if (length(obs_i) < df + 1L) next

    B_i      <- B_global[obs_i, , drop = FALSE]
    phi_rows <- lapply(seq_along(obs_i),
                       function(j) .phi(XZ[i, obs_i[j], ], B_i[j, ]))
    Des_i    <- do.call(rbind, phi_rows)
    ys_i     <- Y[i, obs_i]
    pen      <- lambda * diag(pq)
    coef_i   <- tryCatch(
      as.numeric(solve(crossprod(Des_i) + pen, t(Des_i) %*% ys_i)),
      error = function(e) rep(NA_real_, pq)
    )
    if (!all(is.finite(coef_i))) next

    Theta_i  <- matrix(coef_i, nrow = p, ncol = q, byrow = FALSE)
    # Predict on the full m-point grid (including unobserved times)
    for (t in seq_len(m)) {
      gvec   <- as.numeric(Theta_i %*% B_global[t, ])
      gvec_x <- gvec[seq_len(p_x)]
      fitted_values[i, t]   <- if (p_x == 1L) gvec_x[1L] else sqrt(sum(gvec_x^2))
      beta_vec_arr[i, t, ]  <- gvec_x
    }
  }

  # ── k-means on fitted-value curves ────────────────────────────────────
  valid <- which(rowSums(!is.finite(fitted_values)) == 0L)
  na_return <- list(subject_group = rep(1L, n),
                    gamma_hat     = matrix(NA_real_, K, m),
                    Gamma_hat     = array(NA_real_, dim = c(K, m, p_x)),
                    beta_hat_arr  = array(NA_real_, dim = c(n, m, p_x)),
                    group_matrix  = matrix(1L, n, m),
                    alpha         = alpha,
                    fitted_values = fitted_values,
                    time_sec      = proc.time()["elapsed"] - t0_clock)

  if (length(valid) < K) {
    warning("fit_indiv_spline_fv: too few valid subjects for k-means.")
    return(na_return)
  }

  km <- tryCatch(
    kmeans(fitted_values[valid, , drop = FALSE],
           centers = K, nstart = km_nstart, iter.max = 200L),
    error = function(e) {
      warning(sprintf("fit_indiv_spline_fv: kmeans — %s", conditionMessage(e)))
      list(cluster = rep(1L, length(valid)))
    }
  )
  subject_group        <- rep(1L, n)
  subject_group[valid] <- as.integer(km$cluster)
  group_matrix         <- matrix(subject_group, nrow = n, ncol = m)

  # ── Group trajectories: within-group mean of per-subject beta_vec ─────
  gamma_hat <- matrix(NA_real_, K, m)
  Gamma_hat <- array(NA_real_, dim = c(K, m, p_x))

  for (k in seq_len(K)) {
    idx_k <- which(subject_group == k)
    if (length(idx_k) == 0L) next
    gamma_hat[k, ] <- colMeans(fitted_values[idx_k, , drop = FALSE], na.rm = TRUE)
    for (j in seq_len(p_x))
      Gamma_hat[k, , j] <- colMeans(
        beta_vec_arr[idx_k, , j, drop = FALSE], na.rm = TRUE)
  }

  # Broadcast K x m x p_x -> n x m x p_x
  beta_hat_arr <- Gamma_hat[subject_group, , , drop = FALSE]

  list(subject_group = subject_group,
       gamma_hat     = gamma_hat,
       Gamma_hat     = Gamma_hat,
       beta_hat_arr  = beta_hat_arr,
       group_matrix  = group_matrix,
       alpha         = alpha,
       fitted_values = fitted_values,
       time_sec      = proc.time()["elapsed"] - t0_clock)
}

# ============================================================
# METHOD 5: REGIONAL SPLINE + PER-REGION FITTED-VALUE k-MEANS
# ============================================================

#' Fit per-subject B-splines on hard non-overlapping sub-regions, cluster on
#' per-region fitted-value submatrices, then aggregate via majority vote,
#' with optional homogeneous Z.
#'
#' Key difference from Method 3: per-region k-means uses the n x |t_r|
#' fitted-value submatrix (|t_r| = time points in region r) rather than the
#' n x pq coefficient matrix. Group trajectories are within-group means of
#' per-region fitted values stitched across regions.
#'
#' @param Y          n x m outcome matrix  (NAs allowed)
#' @param X          n x m x p heterogeneous covariate array  (NAs allowed)
#' @param Z          n x m x r homogeneous covariate array (NULL = omit)
#' @param K          number of groups
#' @param df         B-spline df per region
#' @param n_region   number of hard sub-regions  (NULL = auto ~5)
#' @param lambda     ridge penalty for het spline  (default 1e-3)
#' @param lambda_z   ridge penalty for alpha  (default 1e-4)
#' @param km_nstart  k-means restarts per region  (default 20)
#' @param min_obs    min valid obs in a region to fit spline  (NULL = df+1)
#'
#' @return list(subject_group, gamma_hat [K x m], group_matrix [n x m],
#'              alpha [r] | NULL, fitted_values [n x m],
#'              region_labels, region_breaks, time_sec, n_region)
fit_regional_spline_fv <- function(Y, X, Z = NULL,
                                   K         = 2L,
                                   df        = 6L,
                                   n_region  = NULL,
                                   lambda    = 1e-3,
                                   lambda_z  = 1e-4,
                                   km_nstart = 20L,
                                   min_obs   = NULL) {

  n <- nrow(Y); m <- ncol(Y)
  t_global <- seq(0, 1, length.out = m)
  if (is.null(min_obs))  min_obs  <- df + 1L
  if (is.null(n_region)) n_region <- max(3L, min(5L, floor(m / (df + 2L))))

  prep  <- .prep_combined(Y, X, Z)
  XZ    <- prep$XZ        # combined (p + r)-dim covariate
  p_x   <- prep$p_x       # original X columns (1..p_x)
  p     <- dim(XZ)[3]     # total covariate dim = p_x + r
  alpha <- NULL            # Z treated heterogeneously; no shared alpha
  q  <- df;  pq <- p * q  # must come AFTER p is set

  region_breaks <- seq(0, 1, length.out = n_region + 1L)
  region_of_t   <- pmax(1L, pmin(
    findInterval(t_global, region_breaks, rightmost.closed = TRUE), n_region))

  t0_clock <- proc.time()["elapsed"]

  # Precompute per-region B-spline knots
  reg_bk <- vector("list", n_region)
  for (r in seq_len(n_region)) {
    t_in_r <- t_global[region_of_t == r]
    lo_r   <- region_breaks[r]; hi_r <- region_breaks[r + 1L]
    u_ref  <- pmax(pmin((t_in_r - lo_r) / (hi_r - lo_r), 1), 0)
    if (length(unique(u_ref)) < df + 1L)
      u_ref <- seq(0, 1, length.out = df + 2L)
    reg_bk[[r]] <- .bs_knots(u_ref, df)
  }

  # Per-subject, per-region coefficient estimation on residualized R
  coef_array <- array(0, dim = c(n, n_region, pq))

  for (r in seq_len(n_region)) {
    lo_r    <- region_breaks[r]; hi_r <- region_breaks[r + 1L]
    t_idx_r <- which(region_of_t == r)
    for (i in seq_len(n)) {
      obs_r <- t_idx_r[!is.na(Y[i, t_idx_r]) & rowSums(is.na(matrix(XZ[i, t_idx_r, ], nrow=length(t_idx_r), ncol=p))) == 0L]
      if (length(obs_r) < min_obs) next
      u_obs    <- pmax(pmin((t_global[obs_r] - lo_r) / (hi_r - lo_r), 1), 0)
      B_r      <- .bs_eval(u_obs, reg_bk[[r]])
      phi_rows <- lapply(seq_along(obs_r),
                         function(j) .phi(XZ[i, obs_r[j], ], B_r[j, ]))
      ys_r     <- Y[i, obs_r]
      Des_r    <- do.call(rbind, phi_rows)
      pen      <- lambda * diag(pq)
      coef_array[i, r, ] <- tryCatch(
        as.numeric(solve(crossprod(Des_r) + pen, t(Des_r) %*% ys_r)),
        error = function(e) rep(0, pq)
      )
    }
  }

  # Evaluate per-subject fitted values on the full m-point grid
  fitted_values <- matrix(NA_real_, nrow = n, ncol = m)
  beta_vec_arr  <- array(NA_real_, dim = c(n, m, p_x))  # per-subject per-cov
  for (r in seq_len(n_region)) {
    lo_r    <- region_breaks[r]; hi_r <- region_breaks[r + 1L]
    t_idx_r <- which(region_of_t == r)
    for (i in seq_len(n)) {
      coef_r <- coef_array[i, r, ]
      if (!all(is.finite(coef_r))) next
      Theta_i <- matrix(coef_r, nrow = p, ncol = q, byrow = FALSE)
      for (t in t_idx_r) {
        u_t  <- pmax(pmin((t_global[t] - lo_r) / (hi_r - lo_r), 1), 0)
        b_t  <- as.numeric(.bs_eval(u_t, reg_bk[[r]]))
        gvec <- as.numeric(Theta_i %*% b_t)
        gvec_x <- gvec[seq_len(p_x)]    # X-block only -> correct scale
        fitted_values[i, t] <- if (p_x == 1L) gvec_x[1L] else sqrt(sum(gvec_x^2))
        beta_vec_arr[i, t, ] <- gvec_x   # store full p_x-vector
      }
    }
  }

  # Per-region k-means on the n x |t_r| fitted-value submatrix
  region_labels <- matrix(NA_integer_, nrow = n, ncol = n_region)
  for (r in seq_len(n_region)) {
    t_idx_r <- which(region_of_t == r)
    fv_r    <- fitted_values[, t_idx_r, drop = FALSE]   # n x |t_r|
    valid   <- which(rowSums(!is.finite(fv_r)) == 0L)
    if (length(valid) < K) { region_labels[, r] <- 1L; next }
    km_r <- tryCatch(
      kmeans(fv_r[valid, , drop = FALSE], centers = K,
             nstart = km_nstart, iter.max = 200L),
      error = function(e) {
        warning(sprintf("fit_regional_spline_fv: kmeans region %d — %s",
                        r, conditionMessage(e)))
        list(cluster = rep(1L, length(valid)))
      }
    )
    region_labels[valid, r] <- as.integer(km_r$cluster)
    region_labels[is.na(region_labels[, r]), r] <- 1L
  }

  # Align region labels then majority vote
  aligned_labels       <- region_labels
  aligned_labels[, 1L] <- region_labels[, 1L]
  if (n_region > 1L)
    for (r in 2L:n_region)
      aligned_labels[, r] <- .align_to_ref(region_labels[, r],
                                           region_labels[, 1L], K)

  subject_group <- .majority_vote(aligned_labels, K)

  # Group trajectories: within-group means of per-region fitted values
  gamma_hat <- matrix(NA_real_, K, m)
  Gamma_hat <- array(NA_real_, dim = c(K, m, p_x))  # K x m x p_x
  for (r in seq_len(n_region)) {
    t_idx_r <- which(region_of_t == r)
    for (k in seq_len(K)) {
      idx_k <- which(subject_group == k)
      if (length(idx_k) == 0L) next
      gm <- colMeans(fitted_values[idx_k, t_idx_r, drop = FALSE], na.rm = TRUE)
      if (any(!is.finite(gm))) next
      gamma_hat[k, t_idx_r] <- gm
      for (j in seq_len(p_x))
        Gamma_hat[k, t_idx_r, j] <- colMeans(
          matrix(beta_vec_arr[idx_k, t_idx_r, j],
                 nrow = length(idx_k)), na.rm = TRUE)
    }
  }

  # group_matrix: region-specific labels broadcast to time points
  group_matrix <- matrix(NA_integer_, nrow = n, ncol = m)
  for (r in seq_len(n_region))
    group_matrix[, which(region_of_t == r)] <- aligned_labels[, r]

  # Broadcast K x m x p -> n x m x p
  beta_hat_arr <- Gamma_hat[subject_group, , , drop = FALSE]
  list(subject_group  = subject_group,
       gamma_hat      = gamma_hat,
       Gamma_hat      = Gamma_hat,         # K x m x p_x estimated group trajectories
       beta_hat_arr   = beta_hat_arr,
       group_matrix   = group_matrix,
       alpha          = alpha,
       fitted_values  = fitted_values,
       region_labels  = aligned_labels,
       region_breaks  = region_breaks,
       time_sec       = proc.time()["elapsed"] - t0_clock,
       n_region       = n_region)
}



# ============================================================
# METHOD 7: CONCAVE PAIRWISE FUSION (Ma & Huang 2015)
#           Longitudinal adaptation for time-varying heterogeneous regression
# ============================================================

#' Concave pairwise fusion ADMM for time-varying subgroup analysis.
#'
#' Direct longitudinal adaptation of Ma & Huang (2015) "A concave pairwise
#' fusion approach to subgroup analysis."
#'
#' ORIGINAL PAPER model (eq. 1):
#'   y_i = mu_i + x_i' beta + eps_i
#'   mu_i in {alpha_1,...,alpha_K}  (subject-specific intercepts, K groups)
#'   ADMM on pairwise differences eta_{ij} = mu_i - mu_j with MCP/SCAD penalty.
#'
#' OUR ADAPTATION for time-varying coefficients:
#'   Y_it = x_it' b_i(t) + z_it' alpha_t + eps_it
#'
#' Two-step algorithm:
#'
#'   Step 1. Residualise Y on Z (shared alpha_t, pooled per-time ridge).
#'
#'   Step 2. For each subject i and time t, get a pilot estimate v_i(t) via
#'     local kernel-weighted ridge regression (Shujie Ma nonparametric estimator):
#'       A_it = sum_s K_h(s,t) x_is x_is' + lambda_r I     [p x p]
#'       c_it = sum_s K_h(s,t) x_is R_is                   [p]
#'       v_i(t) = A_it^{-1} c_it
#'
#'   Step 3. At each time t, apply Ma & Huang (2015) ADMM fusion on the
#'     n pilot vectors {v_i(t)} in R^p:
#'
#'       min_{b_1,...,b_n} (1/2) sum_i ||b_i - v_i||^2
#'                         + lambda_f * sum_{i<j} p_gamma(||b_i - b_j||_2)
#'
#'     where p_gamma is MCP or SCAD applied to the L2 norm of pairwise
#'     differences (group-lasso extension of Ma & Huang eq. 3).
#'
#'     ADMM variables (following Section 3.1 of the paper):
#'       eta_{ij} = b_i - b_j  in R^p  (auxiliary; one per pair i<j)
#'       upsilon_{ij} in R^p           (dual variable)
#'
#'     Updates at each ADMM iteration:
#'
#'     (a) b-update, EXACT closed form.  The b-step FOC is
#'           (I + rho * Delta'Delta) b = v + rho * Delta'(eta - upsilon/rho)
#'         and Delta'Delta = n_v * I - 11', so by Sherman-Morrison
#'           (I + rho*Delta'Delta)^{-1} = (I + rho * 11') / (1 + rho * n_v).
#'         With w_i = v_i + rho * A_row_i,
#'           A_row_i =  sum_{j>i}(eta_{ij} - upsilon_{ij}/rho)
#'                    - sum_{j<i}(eta_{ji} - upsilon_{ji}/rho)
#'         and noting sum_i A_row_i = 0 (each pair contributes +/-), we get
#'           b_i = ( v_i + rho * A_row_i + rho * sum_j v_j ) / (1 + rho * n_v).
#'
#'     (b) eta-update (eq. 7/8/9, extended to vector norm):
#'         delta_{ij} = b_i - b_j + upsilon_{ij}/rho
#'         d = ||delta_{ij}||_2   (scalar norm of the p-vector)
#'
#'         L1:   eta_{ij} = (1 - lambda_f/(rho*d))_+ * delta_{ij}
#'               (group soft threshold)
#'
#'         MCP:  if d <= gamma*lambda_f:
#'                 eta_{ij} = (1 - lambda_f/(rho*d))_+ * delta_{ij} / (1 - 1/(gamma*rho))
#'               else:
#'                 eta_{ij} = delta_{ij}
#'
#'         SCAD: if d <= lambda_f + lambda_f/rho:
#'                 eta_{ij} = (1 - lambda_f/(rho*d))_+ * delta_{ij}
#'               elif d <= gamma*lambda_f:
#'                 eta_{ij} = (1 - gamma*lambda_f/((gamma-1)*rho*d))_+ * delta_{ij}
#'                             / (1 - 1/((gamma-1)*rho))
#'               else:
#'                 eta_{ij} = delta_{ij}
#'
#'     (c) upsilon-update:
#'         upsilon_{ij} <- upsilon_{ij} + rho*(b_i - b_j - eta_{ij})
#'
#'     Convergence: primal residual ||Delta*b - eta||_F < tol_admm.
#'
#'   Step 4. Group assignment: k-means (K supplied) on fused b_hat at each t,
#'     rolling Hungarian alignment across t; majority vote -> subject_group.
#'     The fused pairwise structure (||eta_{ij}|| < tol_fuse) is additionally
#'     summarised as a per-time diagnostic K_fuse_vec via transitive closure
#'     (union-find) of below-threshold pairs.
#'
#' @param Y         n x m  (NAs allowed)
#' @param X         n x m x p_x  heterogeneous covariates
#' @param Z         n x m x r    homogeneous covariates (NULL = omit)
#' @param K         number of groups (used for k-means post-processing)
#' @param h         kernel bandwidth on [0,1] for pilot estimate (NULL = 0.25)
#' @param lambda_r  ridge penalty for pilot kernel regression  (default 1e-3)
#' @param lambda_f  fusion penalty lambda                      (default 0.1)
#' @param penalty   one of "MCP", "SCAD", "L1"                (default "MCP")
#' @param gamma     concavity parameter for MCP/SCAD           (default 3.0)
#' @param rho       ADMM augmented Lagrangian step size        (default 1.0)
#' @param max_admm  max ADMM iterations                        (default 200)
#' @param tol_admm  ADMM primal residual convergence tolerance (default 1e-3)
#' @param tol_fuse  threshold below which ||eta_ij|| is treated as zero (default 1e-2)
#' @param lambda_z  ridge penalty for Z residualisation        (default 1e-4)
#' @param km_nstart k-means random restarts                    (default 20)
#' @param min_obs_t min observed subjects at t for clustering  (default K+1)
fit_het_reg_l1 <- function(Y, X, Z = NULL, K = 2L,
                           h         = NULL,
                           lambda_r  = 1e-3,
                           lambda_f  = 0.1,
                           penalty   = c("MCP", "SCAD", "L1"),
                           gamma     = 3.0,
                           rho       = 1.0,
                           max_admm  = 200L,
                           tol_admm  = 1e-3,
                           tol_fuse  = 1e-2,
                           lambda_z  = 1e-4,
                           km_nstart = 20L,
                           min_obs_t = NULL) {

  penalty <- match.arg(penalty)
  n <- nrow(Y); m <- ncol(Y)
  t_global <- seq(0, 1, length.out = m)
  p_x <- dim(X)[3]
  if (is.null(min_obs_t)) min_obs_t <- K + 1L
  if (is.null(h))         h         <- 0.25

  # Concave-penalty proximal operators require these for positive denominators
  if (penalty == "MCP"  && gamma * rho <= 1)
    stop("MCP requires gamma * rho > 1")
  if (penalty == "SCAD" && (gamma - 1) * rho <= 1)
    stop("SCAD requires (gamma - 1) * rho > 1")

  t0_clock <- proc.time()["elapsed"]

  # ── Vectorised group-threshold operator ───────────────────────────────
  # Delta: n_pairs x p matrix. Returns the shrunk matrix (Ma & Huang
  # eq. 7/8/9 extended to the L2 norm of each row).
  .thresh_mat <- function(Delta, lf, rho_, gamma_, penalty_) {
    d <- sqrt(rowSums(Delta^2))
    s <- numeric(length(d))
    tiny <- d < 1e-12
    dp <- pmax(d, 1e-12)                       # safe divisor
    if (penalty_ == "L1") {
      s <- pmax(1 - lf / (rho_ * dp), 0)
    } else if (penalty_ == "MCP") {
      inside <- d <= gamma_ * lf
      s <- ifelse(inside,
                  pmax(1 - lf / (rho_ * dp), 0) / (1 - 1 / (gamma_ * rho_)),
                  1)
    } else {  # SCAD
      r1 <- d <= lf + lf / rho_
      r2 <- !r1 & d <= gamma_ * lf
      s <- rep(1, length(d))
      s[r1] <- pmax(1 - lf / (rho_ * dp[r1]), 0)
      s[r2] <- pmax(1 - gamma_ * lf / ((gamma_ - 1) * rho_ * dp[r2]), 0) /
        (1 - 1 / ((gamma_ - 1) * rho_))
    }
    s[tiny] <- 0
    Delta * s
  }

  # ── Step 1: residualise Y on Z ────────────────────────────────────────
  R     <- Y
  alpha <- NULL
  if (!is.null(Z)) {
    r_z   <- dim(Z)[3]
    alpha <- matrix(NA_real_, m, r_z)
    for (t in seq_len(m)) {
      obs_t <- which(!is.na(Y[, t]) &
                       rowSums(is.na(matrix(Z[, t, ], nrow=n, ncol=r_z))) == 0L)
      if (length(obs_t) < r_z + 1L) next
      Z_t  <- matrix(Z[obs_t, t, ], nrow=length(obs_t), ncol=r_z)
      y_t  <- Y[obs_t, t]
      a_t  <- tryCatch(
        as.numeric(solve(crossprod(Z_t) + lambda_z*diag(r_z), t(Z_t) %*% y_t)),
        error = function(e) rep(0, r_z))
      alpha[t, ]  <- a_t
      R[obs_t, t] <- y_t - as.numeric(Z_t %*% a_t)
    }
  }

  # ── Step 2: per-subject local kernel ridge pilot estimates ────────────
  # v_hat[i, t, ] = A_it^{-1} c_it  (Shujie Ma nonparametric estimator)
  pen_r <- lambda_r * diag(p_x)
  v_hat <- array(NA_real_, dim = c(n, m, p_x))

  for (i in seq_len(n)) {
    obs_i <- which(!is.na(R[i, ]) &
                     rowSums(is.na(matrix(X[i, , ], nrow=m, ncol=p_x))) == 0L)
    if (length(obs_i) < 2L) next
    x_obs <- matrix(X[i, obs_i, ], nrow=length(obs_i), ncol=p_x)
    r_obs <- R[i, obs_i]
    t_obs <- t_global[obs_i]
    for (t in seq_len(m)) {
      w    <- exp(-(t_obs - t_global[t])^2 / (2*h^2))
      if (sum(w) < 1e-10) next
      Xw   <- x_obs * w
      A_it <- crossprod(Xw, x_obs) + pen_r
      c_it <- as.numeric(t(Xw) %*% r_obs)
      v_hat[i, t, ] <- tryCatch(as.numeric(solve(A_it, c_it)),
                                error = function(e) rep(NA_real_, p_x))
    }
  }

  # ── Step 3: per-time ADMM fusion (Ma & Huang 2015, Section 3.1) ──────
  # At each t, fuse {v_i(t)} with concave pairwise penalty.
  # Denoising form: min (1/2)||b-v||^2 + lf * sum_{i<j} p_gamma(||b_i-b_j||)
  b_fused    <- array(NA_real_, dim = c(n, m, p_x))
  K_fuse_vec <- rep(NA_integer_, m)     # per-t fused-structure diagnostic

  for (t in seq_len(m)) {
    v_t   <- matrix(v_hat[, t, ], nrow=n, ncol=p_x)
    valid <- which(rowSums(is.na(v_t)) == 0L)
    n_v   <- length(valid)
    if (n_v < min_obs_t) next
    v_mat <- v_t[valid, , drop=FALSE]   # n_v x p_x

    if (n_v == 1L || lambda_f == 0) {
      b_fused[valid, t, ] <- v_mat
      K_fuse_vec[t] <- n_v
      next
    }

    # All (i < j) pairs, vectorised (column-major upper.tri gives i < j)
    iu <- which(upper.tri(diag(n_v)), arr.ind = TRUE)
    ei <- iu[, 1L]; ej <- iu[, 2L]
    n_pairs <- length(ei)
    pair_grp <- c(ei, ej)               # for rowsum() in the b-update

    # Initialise: warm-start eta at the observed pairwise differences
    # (eta = 0 would pull the first iterations toward total fusion)
    b       <- v_mat                                          # n_v x p_x
    eta     <- v_mat[ei, , drop=FALSE] - v_mat[ej, , drop=FALSE]
    upsilon <- matrix(0, nrow=n_pairs, ncol=p_x)              # dual

    sum_v <- matrix(colSums(v_mat), n_v, p_x, byrow = TRUE)
    denom <- 1 + rho * n_v

    for (iter in seq_len(max_admm)) {

      # (a) b-update: EXACT closed form (see header).  A_row via signed
      # aggregation of ZU = eta - upsilon/rho over pair endpoints; rowsum
      # returns rows in sorted group order = 1..n_v (all present).
      ZU    <- eta - upsilon / rho
      A_mat <- rowsum(rbind(ZU, -ZU), group = pair_grp)
      b     <- (v_mat + rho * A_mat + rho * sum_v) / denom

      # (b) eta-update: group threshold at delta = Delta b + upsilon/rho
      Bdiff <- b[ei, , drop=FALSE] - b[ej, , drop=FALSE]
      eta   <- .thresh_mat(Bdiff + upsilon / rho,
                           lambda_f, rho, gamma, penalty)

      # (c) dual update + primal residual
      Rmat    <- Bdiff - eta
      upsilon <- upsilon + rho * Rmat
      if (sqrt(sum(Rmat^2)) < tol_admm) break
    }

    b_fused[valid, t, ] <- b

    # Fused-structure diagnostic: number of connected components under
    # ||eta_{ij}|| < tol_fuse (union-find / transitive closure)
    d_eta  <- sqrt(rowSums(eta^2))
    parent <- seq_len(n_v)
    find_root <- function(x) { while (parent[x] != x) x <- parent[x]; x }
    for (e in which(d_eta < tol_fuse)) {
      ra <- find_root(ei[e]); rb <- find_root(ej[e])
      if (ra != rb) parent[rb] <- ra
    }
    K_fuse_vec[t] <- length(unique(vapply(seq_len(n_v), find_root, integer(1))))
  }

  # ── Step 4: extract groups ────────────────────────────────────────────
  # K-means on fused estimates (K supplied); rolling Hungarian alignment
  # across t; majority vote -> subject_group.
  time_labels <- matrix(NA_integer_, nrow=n, ncol=m)

  for (t in seq_len(m)) {
    b_t   <- matrix(b_fused[, t, ], nrow=n, ncol=p_x)
    valid <- which(rowSums(is.na(b_t)) == 0L)
    if (length(valid) < min_obs_t) next

    km_t <- tryCatch(
      kmeans(b_t[valid, , drop=FALSE],
             centers=K, nstart=km_nstart, iter.max=200L),
      error = function(e) list(cluster=rep(1L, length(valid)))
    )
    lbl_t        <- rep(NA_integer_, n)
    lbl_t[valid] <- as.integer(km_t$cluster)

    if (t == 1L || all(is.na(time_labels[, t-1L]))) {
      time_labels[, t] <- lbl_t
    } else {
      time_labels[, t] <- .align_to_ref(lbl_t, time_labels[, t-1L], K)
    }
  }

  # Fill gaps; majority vote -> subject_group
  for (i in seq_len(n)) {
    na_t    <- which(is.na(time_labels[i, ]))
    valid_t <- which(!is.na(time_labels[i, ]))
    if (length(valid_t) == 0L) { time_labels[i, ] <- 1L; next }
    for (t in na_t)
      time_labels[i, t] <- time_labels[i, valid_t[which.min(abs(valid_t - t))]]
  }

  subject_group <- .majority_vote(time_labels, K)
  group_matrix  <- time_labels
  for (i in seq_len(n))
    group_matrix[i, is.na(group_matrix[i, ])] <- subject_group[i]

  # ── Group trajectories ────────────────────────────────────────────────
  gamma_hat <- matrix(NA_real_, K, m)
  Gamma_hat <- array(NA_real_, dim=c(K, m, p_x))

  for (k in seq_len(K)) {
    idx_k <- which(subject_group == k)
    if (length(idx_k) == 0L) next
    for (t in seq_len(m)) {
      b_kt <- matrix(b_fused[idx_k, t, ], nrow=length(idx_k), ncol=p_x)
      ok   <- which(rowSums(is.na(b_kt)) == 0L)
      if (length(ok) == 0L) next
      mean_b            <- colMeans(b_kt[ok, , drop=FALSE])
      gamma_hat[k, t]   <- if (p_x == 1L) mean_b[1L] else sqrt(sum(mean_b^2))
      Gamma_hat[k, t, ] <- mean_b
    }
  }

  # beta_hat_arr: the method's genuine individual-level estimate is the
  # fused b_i(t), NOT the group-mean broadcast.  NAs (times where fusion
  # was skipped) are backfilled from the assigned group's trajectory.
  beta_hat_arr <- b_fused
  for (i in seq_len(n)) {
    na_t <- which(is.na(beta_hat_arr[i, , 1L]))
    if (length(na_t) > 0L && !is.na(subject_group[i]))
      beta_hat_arr[i, na_t, ] <- Gamma_hat[subject_group[i], na_t, ]
  }

  list(subject_group = subject_group,
       gamma_hat     = gamma_hat,
       Gamma_hat     = Gamma_hat,
       Gamma         = Gamma_hat,   # alias: runner / .adapt_fit expect $Gamma
       beta_hat_arr  = beta_hat_arr,
       group_matrix  = group_matrix,
       K_fuse_vec    = K_fuse_vec,  # per-t fused-K diagnostic (tol_fuse)
       alpha         = alpha,
       b_hat         = v_hat,       # pilot (pre-fusion) estimates
       b_fused       = b_fused,     # post-fusion estimates
       time_sec      = proc.time()["elapsed"] - t0_clock)
}

# ============================================================
# METHOD 8: KERNEL RIDGE PILOT + K-MEANS
#           (Steps 1-2 of HetRegL1 without ADMM fusion)
# ============================================================
#
#' Per-subject local kernel ridge pilot estimates, clustered via k-means.
#'
#' This is a standalone baseline consisting of exactly the first two steps
#' of fit_het_reg_l1 (Ma & Huang 2015 adaptation), without the ADMM fusion
#' step. It isolates what a pointwise nonparametric estimator can achieve:
#'
#'   Step 1. Residualise Y on Z (pooled per-time ridge, same as HetRegL1).
#'
#'   Step 2. Per-subject local kernel ridge pilot:
#'             A_it = sum_s K_h(s,t) x_is x_is' + lambda_r * I
#'             c_it = sum_s K_h(s,t) x_is R_is
#'             v_hat[i,t,] = solve(A_it, c_it)
#'
#'   Step 3. At each time t, run k-means on v_hat[valid,t,] with K clusters.
#'           Rolling Hungarian alignment + majority vote -> subject_group.
#'
#' Use as a reference to measure the contribution of fusion (ADMM) above
#' the raw kernel pilot. HetRegL1 should always beat this on merging settings;
#' if it does not, the ADMM fusion is adding noise rather than signal.
#'
#' @param Y         n x m  (NAs allowed)
#' @param X         n x m x p_x
#' @param Z         n x m x r  or NULL
#' @param K         number of groups
#' @param h         kernel bandwidth (default 0.25)
#' @param lambda_r  ridge penalty for pilot regression (default 1e-3)
#' @param lambda_z  ridge penalty for Z residualisation (default 1e-4)
#' @param km_nstart k-means random restarts (default 20)
#' @param min_obs_t minimum valid subjects at t to run k-means (default K+1)
fit_kernel_pilot <- function(Y, X, Z = NULL, K = 2L,
                             h         = 0.25,
                             lambda_r  = 1e-3,
                             lambda_z  = 1e-4,
                             km_nstart = 20L,
                             min_obs_t = NULL) {

  n <- nrow(Y); m <- ncol(Y)
  t_global <- seq(0, 1, length.out = m)
  p_x <- dim(X)[3]
  if (is.null(min_obs_t)) min_obs_t <- K + 1L

  t0_clock <- proc.time()["elapsed"]

  # ── Step 1: residualise Y on Z ────────────────────────────────────────
  R     <- Y
  alpha <- NULL
  if (!is.null(Z)) {
    r_z   <- dim(Z)[3]
    alpha <- matrix(NA_real_, m, r_z)
    for (t in seq_len(m)) {
      obs_t <- which(!is.na(Y[, t]) &
                       rowSums(is.na(matrix(Z[, t, ], nrow=n, ncol=r_z))) == 0L)
      if (length(obs_t) < r_z + 1L) next
      Z_t   <- matrix(Z[obs_t, t, ], nrow=length(obs_t), ncol=r_z)
      y_t   <- Y[obs_t, t]
      a_t   <- tryCatch(
        as.numeric(solve(crossprod(Z_t) + lambda_z * diag(r_z), t(Z_t) %*% y_t)),
        error = function(e) rep(0, r_z))
      alpha[t, ]  <- a_t
      R[obs_t, t] <- y_t - as.numeric(Z_t %*% a_t)
    }
  }

  # ── Step 2: per-subject local kernel ridge pilot estimates ────────────
  # v_hat[i, t, ] = A_it^{-1} c_it  (Shujie Ma nonparametric estimator)
  pen_r <- lambda_r * diag(p_x)
  v_hat <- array(NA_real_, dim = c(n, m, p_x))

  for (i in seq_len(n)) {
    obs_i <- which(!is.na(R[i, ]) &
                     rowSums(is.na(matrix(X[i, , ], nrow=m, ncol=p_x))) == 0L)
    if (length(obs_i) < 2L) next
    x_obs <- matrix(X[i, obs_i, ], nrow=length(obs_i), ncol=p_x)
    r_obs <- R[i, obs_i]
    t_obs <- t_global[obs_i]
    for (t in seq_len(m)) {
      w    <- exp(-(t_obs - t_global[t])^2 / (2 * h^2))
      if (sum(w) < 1e-10) next
      Xw   <- x_obs * w
      A_it <- crossprod(Xw, x_obs) + pen_r
      c_it <- as.numeric(t(Xw) %*% r_obs)
      v_hat[i, t, ] <- tryCatch(
        as.numeric(solve(A_it, c_it)),
        error = function(e) rep(NA_real_, p_x))
    }
  }

  # ── Step 3: per-time k-means on pilot estimates ───────────────────────
  time_labels <- matrix(NA_integer_, nrow = n, ncol = m)

  for (t in seq_len(m)) {
    v_t   <- matrix(v_hat[, t, ], nrow = n, ncol = p_x)
    valid <- which(rowSums(is.na(v_t)) == 0L)
    if (length(valid) < min_obs_t) next

    km_t <- tryCatch(
      kmeans(v_t[valid, , drop=FALSE],
             centers = K, nstart = km_nstart, iter.max = 200L),
      error = function(e) list(cluster = rep(1L, length(valid))))

    lbl_t        <- rep(NA_integer_, n)
    lbl_t[valid] <- as.integer(km_t$cluster)

    if (t == 1L || all(is.na(time_labels[, t - 1L]))) {
      time_labels[, t] <- lbl_t
    } else {
      time_labels[, t] <- .align_to_ref(lbl_t, time_labels[, t - 1L], K)
    }
  }

  # Fill gaps; majority vote -> subject_group
  for (i in seq_len(n)) {
    na_t    <- which(is.na(time_labels[i, ]))
    valid_t <- which(!is.na(time_labels[i, ]))
    if (length(valid_t) == 0L) { time_labels[i, ] <- 1L; next }
    for (t in na_t)
      time_labels[i, t] <- time_labels[i, valid_t[which.min(abs(valid_t - t))]]
  }

  subject_group <- .majority_vote(time_labels, K)
  group_matrix  <- time_labels
  for (i in seq_len(n))
    group_matrix[i, is.na(group_matrix[i, ])] <- subject_group[i]

  # ── Group trajectories ────────────────────────────────────────────────
  gamma_hat <- matrix(NA_real_, K, m)
  Gamma_hat <- array(NA_real_, dim = c(K, m, p_x))

  for (k in seq_len(K)) {
    idx_k <- which(subject_group == k)
    if (length(idx_k) == 0L) next
    for (t in seq_len(m)) {
      v_kt <- matrix(v_hat[idx_k, t, ], nrow = length(idx_k), ncol = p_x)
      ok   <- which(rowSums(is.na(v_kt)) == 0L)
      if (length(ok) == 0L) next
      mean_v            <- colMeans(v_kt[ok, , drop = FALSE])
      gamma_hat[k, t]   <- if (p_x == 1L) mean_v[1L] else sqrt(sum(mean_v^2))
      Gamma_hat[k, t, ] <- mean_v
    }
  }

  beta_hat_arr <- Gamma_hat[subject_group, , , drop = FALSE]

  list(subject_group = subject_group,
       gamma_hat     = gamma_hat,
       Gamma_hat     = Gamma_hat,
       beta_hat_arr  = beta_hat_arr,
       group_matrix  = group_matrix,
       alpha         = alpha,
       v_hat         = v_hat,       # pilot estimates (n x m x p_x)
       time_sec      = proc.time()["elapsed"] - t0_clock)
}


# ============================================================
# METHOD 9 — DISCRETE-TIME FLEXMIX
# ============================================================

#' Fit a K-component mixture-of-regressions independently at each time point t,
#' then align component labels sequentially across time via Hungarian matching.
#'
#' Algorithm:
#'   For each t = 1, ..., m:
#'     1. Build the cross-sectional design matrix W_t = [Z_t | X_t] (n_obs x (r+p))
#'        for all subjects observed at t.  When Z = NULL, W_t = X_t.
#'     2. Fit a K-component FlexMix mixture-of-regressions:
#'          Y_it ~ N( w_it' beta_k, sigma_k^2 )  with mixing weights pi_k
#'        using the no-intercept formula  y ~ W - 1.
#'     3. Record hard cluster assignments and per-component coefficient vectors.
#'     4. Align component labels to the previous time point via Hungarian
#'        matching on the K x K co-occurrence table (rolling reference).
#'        Simultaneously permute the coefficient matrix to stay consistent.
#'
#'   After the time loop:
#'     5. Forward/backward-fill NA labels (time points with too few obs or
#'        degenerate fits) using each subject's nearest observed time.
#'     6. Majority-vote across time -> subject_group.
#'     7. gamma_hat[k, t] = L2 norm of the X-block of the group-k coefficient
#'        at time t (consistent with all other baselines that use X only for
#'        the trajectory scalar; Z coefficients are returned separately as alpha).
#'
#' Degeneracy handling: if FlexMix collapses to one component at time t,
#' labels are carried forward from t-1 (or randomly seeded at t=1).
#' The degenerate flag vector records which time points were affected.
#'
#' @param Y        n x m outcome matrix  (NAs allowed)
#' @param X        n x m x p heterogeneous covariate array
#' @param Z        n x m x r homogeneous covariate array  (NULL = omit)
#' @param K        number of mixture components
#' @param lambda_z unused; kept for interface compatibility
#' @param iter_max max FlexMix EM iterations per time point  (default 200)
#' @param min_prior minimum component prior weight; components below this are
#'                  pruned by FlexMix  (default 0.05)
#' @param min_obs_t minimum observed subjects at t to attempt fitting (default K+1)
#'
#' @return list(subject_group [n], gamma_hat [K x m], group_matrix [n x m],
#'              alpha [m x r] | NULL, degenerate [m], time_sec)
fit_flex_timewise <- function(Y, X, Z = NULL, K = 2L,
                              lambda_z  = 1e-4,   # unused; interface compat.
                              iter_max  = 200L,
                              min_prior = 0.05,
                              min_obs_t = NULL,
                              ...) {

  n <- nrow(Y); m <- ncol(Y)
  p <- dim(X)[3L]
  r <- if (!is.null(Z)) dim(Z)[3L] else 0L
  if (is.null(min_obs_t)) min_obs_t <- K + 1L

  t0          <- proc.time()[["elapsed"]]
  group_mat   <- matrix(NA_integer_, n, m)   # per-subject per-time labels
  coef_arr    <- array(NA_real_, c(K, m, p + r))  # K x m x (p+r) coefficients
  degenerate  <- logical(m)

  # BUG FIX 1: S4 constructor — no namespace prefix on new()
  fm_ctrl <- new("FLXcontrol",
                 iter.max  = iter_max,
                 minprior  = min_prior,
                 verbose   = 0L)

  # ── Pass 1: fit mixture regression at each t independently ───────────────
  for (t in seq_len(m)) {

    y_t <- Y[, t]

    # Build W_t = [Z_t | X_t]  (n x (r+p)); select complete cases
    W_t   <- if (r > 0L) cbind(matrix(Z[, t, ], n, r),
                               matrix(X[, t, ], n, p))
    else         matrix(X[, t, ], n, p)
    obs   <- which(!is.na(y_t) & rowSums(is.na(W_t)) == 0L)
    n_obs <- length(obs)

    if (n_obs < min_obs_t) {
      degenerate[t] <- TRUE
      next
    }

    y_obs <- y_t[obs]
    W_obs <- W_t[obs, , drop = FALSE]

    fit_t <- tryCatch(
      flexmix::flexmix(y_obs ~ W_obs - 1, k = K, control = fm_ctrl),
      error = function(e) NULL
    )

    # BUG FIX 2: degeneracy check — use @k slot, not flexmix::length()
    is_degen <- is.null(fit_t) || fit_t@k < K

    if (is_degen) {
      degenerate[t] <- TRUE
      if (t > 1L && any(!is.na(group_mat[obs, t - 1L]))) {
        group_mat[obs, t] <- group_mat[obs, t - 1L]
        coef_arr[, t, ]   <- coef_arr[, t - 1L, ]
      } else {
        group_mat[obs, t] <- as.integer(
          sample(rep(seq_len(K), length.out = n_obs))
        )
        coef_arr[, t, ] <- 0
      }
      next
    }

    # ── Non-degenerate: extract assignments and coefficients ─────────────
    degenerate[t]     <- FALSE
    group_mat[obs, t] <- as.integer(flexmix::clusters(fit_t))

    # BUG FIX 3: flexmix::parameters() appends sigma at the end of each
    # component vector. Strip to the first (p + r) entries (covariate coefs).
    coef_t <- sapply(seq_len(K), function(k) {
      cv <- as.numeric(flexmix::parameters(fit_t, component = k))
      cv[seq_len(p + r)]   # drop trailing sigma
    })
    # coef_t is (p+r) x K; transpose to K x (p+r) for coef_arr
    coef_arr[, t, ] <- t(coef_t)

    # ── Rolling Hungarian alignment to t-1 ───────────────────────────────
    if (t > 1L && any(!is.na(group_mat[obs, t - 1L]))) {
      prev <- group_mat[obs, t - 1L]
      cur  <- group_mat[obs, t]
      ok   <- !is.na(prev) & !is.na(cur)

      if (any(ok)) {
        tab <- matrix(0L, K, K)
        for (i in which(ok))
          tab[prev[i], cur[i]] <- tab[prev[i], cur[i]] + 1L

        perm <- as.integer(clue::solve_LSAP(tab, maximum = TRUE))

        # BUG FIX 4: remap labels without in-place aliasing
        remapped <- rep(NA_integer_, length(cur))
        for (k in seq_len(K))
          remapped[which(cur == perm[k])] <- k
        group_mat[obs, t] <- remapped

        # BUG FIX 5: permute coefficient rows via temp copy (avoid aliasing)
        tmp <- coef_arr[perm, t, , drop = FALSE]
        coef_arr[, t, ] <- tmp
      }
    }
  }

  # ── Pass 2: fill NA time points (degenerate or insufficient obs) ─────────
  for (i in seq_len(n)) {
    na_t    <- which(is.na(group_mat[i, ]))
    valid_t <- which(!is.na(group_mat[i, ]))
    if (length(valid_t) == 0L) { group_mat[i, ] <- 1L; next }
    for (t in na_t)
      group_mat[i, t] <- group_mat[i, valid_t[which.min(abs(valid_t - t))]]
  }

  # ── Pass 3: majority-vote -> subject_group ───────────────────────────────
  subject_group <- .majority_vote(group_mat, K)

  # ── gamma_hat: L2 norm of X-block of group-mean coefficient at each t ────
  # Coefficient layout in coef_arr[k, t, ]: first r entries = Z block,
  # next p entries = X block.
  x_idx     <- r + seq_len(p)           # column indices for X covariates
  gamma_hat    <- matrix(NA_real_, K, m)
  # Build beta_hat_arr as n x m x p directly (no intermediate K x m x p array)
  # subject i gets the coefficient vector of their assigned group at each t
  beta_hat_arr  <- array(NA_real_, dim = c(n, m, p))  # n x m x p
  Gamma_hat     <- array(NA_real_, dim = c(K, m, p))  # K x m x p estimated group trajectories
  for (k in seq_len(K))
    for (t in seq_len(m)) {
      cv <- coef_arr[k, t, x_idx]
      if (any(is.na(cv))) next
      gamma_hat[k, t]     <- if (p == 1L) cv[1L] else sqrt(sum(cv^2))
      Gamma_hat[k, t, ] <- cv
    }

  # ── alpha: m x r matrix of mean Z-block coefficients across components ───
  alpha <- NULL
  if (r > 0L) {
    z_idx <- seq_len(r)
    alpha <- matrix(NA_real_, m, r)
    for (t in seq_len(m)) {
      # Weight by group sizes at this time point
      sz  <- tabulate(group_mat[, t], nbins = K)
      wts <- sz / max(sum(sz), 1L)
      a_t <- numeric(r)
      for (k in seq_len(K))
        a_t <- a_t + wts[k] * coef_arr[k, t, z_idx]
      if (all(is.finite(a_t))) alpha[t, ] <- a_t
    }
  }

  # Canonicalise component ordering: sort by ascending mean trajectory
  can <- .canonicalise_components(gamma_hat, subject_group, group_mat, K)

  # Broadcast: permute Gamma_hat to canonical order, then expand to n x m x p
  gamma_k_can <- Gamma_hat[can$perm, , , drop = FALSE]  # K x m x p, canonical order
  for (i in seq_len(n))
    beta_hat_arr[i, , ] <- gamma_k_can[can$subject_group[i], , ]

  list(
    subject_group = can$subject_group,
    gamma_hat     = can$gamma_hat,
    Gamma_hat     = gamma_k_can,     # K x m x p estimated group trajectories
    beta_hat_arr  = beta_hat_arr,
    group_matrix  = can$group_matrix,
    alpha         = alpha,
    degenerate    = degenerate,
    time_sec      = proc.time()[["elapsed"]] - t0
  )
}


# Timewise label alignment for per-t mixture fits (MoR TC).
# At each t, exhaustively search the K! permutations matching fitted rows to
# TRUE gamma rows in the full p-dimensional space; among permutations within
# rel_tol of the per-t minimum cost, prefer the one used at t-1 (continuity
# through merged/crossing phases, where the assignment is otherwise arbitrary).
#
# Direction convention: perm[a] = b means fitted row a matches truth row b.
#   - Reordering the fitted ARRAY into truth order uses order(perm):
#       Gamma_new[, t, ] = Gamma[order(perm), t, ]
#   - Relabelling fitted LABELS uses perm directly:  lbl_new = perm[lbl]
align_timewise_gamma <- function(fit, dat, rel_tol = 1e-8) {
  Gt <- if (!is.null(dat$Gamma_true)) dat$Gamma_true else dat$G
  Gm <- fit$Gamma
  if (is.null(Gt) || is.null(Gm)) return(fit)
  if (length(dim(Gt)) == 2L) Gt <- array(Gt, c(dim(Gt), 1L))   # K x M -> K x M x 1

  p_use <- min(dim(Gm)[3L], dim(Gt)[3L])
  perms <- .all_perms(K)
  prev  <- seq_len(K)
  perm_by_t <- matrix(NA_integer_, M, K)

  for (tt in seq_len(M)) {
    # cost of assigning fitted row a to truth row b at this t (full p-space)
    cost <- matrix(NA_real_, K, K)
    for (a in seq_len(K))
      for (b in seq_len(K))
        cost[a, b] <- sum((Gm[a, tt, seq_len(p_use)] -
                             Gt[b, tt, seq_len(p_use)])^2, na.rm = TRUE)

    vals <- vapply(perms, function(pp)
      sum(cost[cbind(seq_len(K), pp)]), numeric(1))
    v_min <- min(vals)
    tied  <- which(vals <= v_min * (1 + rel_tol) + rel_tol)

    # among (near-)ties, pick the permutation most consistent with t-1
    pick <- tied[which.max(vapply(tied, function(j)
      sum(perms[[j]] == prev), numeric(1)))]
    perm <- perms[[pick]]

    perm_by_t[tt, ] <- perm
    inv <- order(perm)                       # array reindexing: order(perm)
    Gm[, tt, ] <- Gm[inv, tt, , drop = FALSE]

    # relabel per-t soft/hard assignments if the fit carries them
    if (!is.null(fit$tau) && length(dim(fit$tau)) == 3L)
      fit$tau[, tt, ] <- fit$tau[, tt, inv, drop = FALSE]
    if (!is.null(fit$labels_t) && is.matrix(fit$labels_t))
      fit$labels_t[, tt] <- perm[fit$labels_t[, tt]]

    prev <- perm
  }

  fit$Gamma     <- Gm
  fit$perm_by_t <- perm_by_t                 # keep for diagnostics
  fit
}




# ============================================================
# METHOD 9: ROLLING WINDOW ESTIMATION
#           Local-constant windowed regression for time-varying
#           heterogeneous coefficients
# ============================================================

#' Rolling-window estimation of individual and group coefficient trajectories.
#'
#' Model:  Y_it = x_it' b_i(t) + z_it' alpha_t + eps_it
#'
#' Local-constant approximation: within a window W(t) of time points around t,
#' coefficients are treated as constant.  This is the boxcar-kernel analogue
#' of the local kernel pilot (KernelPilot), and the simplest honest
#' time-varying baseline.
#'
#' Algorithm, for each t = 1..m:
#'
#'   Step 0. Residualise Y on Z per time point (pooled per-time ridge),
#'     identical to the other baselines: R_it = Y_it - z_it' alpha_hat_t.
#'
#'   Step 1. INDIVIDUAL estimate (beta): per-subject ridge on the subject's
#'     observations pooled over W(t):
#'       b_i(t) = ( sum_{s in W(t)} x_is x_is' + lambda_r I )^{-1}
#'                  sum_{s in W(t)} x_is R_is
#'
#'   Step 2. GROUPING: k-means (K supplied) on { b_i(t) } at each t,
#'     rolling Hungarian alignment across t; majority vote -> subject_group.
#'
#'   Step 3. GROUP estimate (gamma): pooled within-group ridge over W(t),
#'     using all observations of subjects assigned to group k:
#'       gamma_k(t) = ( sum_{i in G_k} sum_{s in W(t)} x_is x_is' + lambda_r I )^{-1}
#'                      sum_{i in G_k} sum_{s in W(t)} x_is R_is
#'     Pooling gives gamma the full sqrt(n_k * |W|) variance reduction;
#'     beta stays genuinely individual-level (Step 1), consistent with the
#'     convention used for the concave-fusion baseline.
#'
#' Window definition (indices, inclusive):
#'   side = "centered":  W(t) = [t - w, t + w]   (clipped to [1, m])
#'   side = "trailing":  W(t) = [t - w + 1, t]   (clipped; strict
#'                        no-lookahead -- use this for the CRSP protocol)
#'
#' @param Y         n x m  (NAs allowed)
#' @param X         n x m x p_x  heterogeneous covariates
#' @param Z         n x m x r    homogeneous covariates (NULL = omit)
#' @param K         number of groups (k-means)
#' @param w         window half-width (centered) or full width (trailing),
#'                  in time-point units (default: max(3, round(m/8)))
#' @param side      "centered" or "trailing"                  (default centered)
#' @param lambda_r  ridge penalty for windowed regressions    (default 1e-3)
#' @param lambda_z  ridge penalty for Z residualisation       (default 1e-4)
#' @param km_nstart k-means random restarts                   (default 20)
#' @param min_obs_i min obs per subject in window to fit b_i  (default p_x + 1)
#' @param min_obs_t min subjects with valid b_i at t to cluster (default K + 1)
fit_rolling_window <- function(Y, X, Z = NULL, K = 2L,
                               w         = NULL,
                               side      = c("centered", "trailing"),
                               lambda_r  = 1e-3,
                               lambda_z  = 1e-4,
                               km_nstart = 20L,
                               min_obs_i = NULL,
                               min_obs_t = NULL) {

  side <- match.arg(side)
  n <- nrow(Y); m <- ncol(Y)
  p_x <- dim(X)[3]
  if (is.null(w))         w         <- max(3L, round(m / 8))
  if (is.null(min_obs_i)) min_obs_i <- p_x + 1L
  if (is.null(min_obs_t)) min_obs_t <- K + 1L

  t0_clock <- proc.time()["elapsed"]
  pen_r <- lambda_r * diag(p_x)

  # Window index sets, precomputed once
  win_idx <- lapply(seq_len(m), function(t) {
    if (side == "centered") seq(max(1L, t - w), min(m, t + w))
    else                    seq(max(1L, t - w + 1L), t)
  })

  # ── Step 0: residualise Y on Z (per-time pooled ridge) ────────────────
  R     <- Y
  alpha <- NULL
  if (!is.null(Z)) {
    r_z   <- dim(Z)[3]
    alpha <- matrix(NA_real_, m, r_z)
    for (t in seq_len(m)) {
      obs_t <- which(!is.na(Y[, t]) &
                       rowSums(is.na(matrix(Z[, t, ], nrow=n, ncol=r_z))) == 0L)
      if (length(obs_t) < r_z + 1L) next
      Z_t <- matrix(Z[obs_t, t, ], nrow=length(obs_t), ncol=r_z)
      y_t <- Y[obs_t, t]
      a_t <- tryCatch(
        as.numeric(solve(crossprod(Z_t) + lambda_z*diag(r_z), t(Z_t) %*% y_t)),
        error = function(e) rep(0, r_z))
      alpha[t, ]  <- a_t
      R[obs_t, t] <- y_t - as.numeric(Z_t %*% a_t)
    }
  }

  # ── Step 1: per-subject rolling-window ridge -> beta_hat_arr ─────────
  beta_hat_arr <- array(NA_real_, dim = c(n, m, p_x))

  for (i in seq_len(n)) {
    X_i   <- matrix(X[i, , ], nrow = m, ncol = p_x)
    ok_i  <- which(!is.na(R[i, ]) & rowSums(is.na(X_i)) == 0L)
    if (length(ok_i) < min_obs_i) next
    for (t in seq_len(m)) {
      s <- intersect(win_idx[[t]], ok_i)
      if (length(s) < min_obs_i) next
      X_w <- X_i[s, , drop = FALSE]
      r_w <- R[i, s]
      beta_hat_arr[i, t, ] <- tryCatch(
        as.numeric(solve(crossprod(X_w) + pen_r, crossprod(X_w, r_w))),
        error = function(e) rep(NA_real_, p_x))
    }
  }

  # ── Step 2: per-time clustering + rolling alignment ───────────────────
  time_labels <- matrix(NA_integer_, nrow = n, ncol = m)

  for (t in seq_len(m)) {
    b_t   <- matrix(beta_hat_arr[, t, ], nrow = n, ncol = p_x)
    valid <- which(rowSums(is.na(b_t)) == 0L)
    if (length(valid) < min_obs_t) next

    km_t <- tryCatch(
      kmeans(b_t[valid, , drop = FALSE],
             centers = K, nstart = km_nstart, iter.max = 200L),
      error = function(e) list(cluster = rep(1L, length(valid)))
    )
    lbl_t        <- rep(NA_integer_, n)
    lbl_t[valid] <- as.integer(km_t$cluster)

    if (t == 1L || all(is.na(time_labels[, t-1L]))) {
      time_labels[, t] <- lbl_t
    } else {
      time_labels[, t] <- .align_to_ref(lbl_t, time_labels[, t-1L], K)
    }
  }

  # Fill gaps.  Centered: nearest valid label (either direction).
  # Trailing (causal): most recent PAST label only -- filling from the
  # future would leak information into the backtest.
  for (i in seq_len(n)) {
    na_t    <- which(is.na(time_labels[i, ]))
    valid_t <- which(!is.na(time_labels[i, ]))
    if (length(valid_t) == 0L) { time_labels[i, ] <- 1L; next }
    for (t in na_t) {
      if (side == "trailing") {
        past <- valid_t[valid_t < t]
        time_labels[i, t] <- if (length(past) > 0L) time_labels[i, max(past)]
        else time_labels[i, min(valid_t)]  # burn-in only
      } else {
        time_labels[i, t] <- time_labels[i, valid_t[which.min(abs(valid_t - t))]]
      }
    }
  }

  subject_group <- .majority_vote(time_labels, K)   # full-sample summary label
  group_matrix  <- time_labels

  # Per-t group composition used for the pooled gamma step and backfill:
  #   centered: static full-sample majority label (standard convention)
  #   trailing: RUNNING majority vote over labels[1:t] (causal)
  grp_at_t <- matrix(NA_integer_, n, m)
  if (side == "trailing") {
    for (i in seq_len(n)) {
      cnt <- integer(K)
      for (t in seq_len(m)) {
        l <- time_labels[i, t]
        if (!is.na(l)) cnt[l] <- cnt[l] + 1L
        grp_at_t[i, t] <- if (sum(cnt) == 0L) NA_integer_ else which.max(cnt)
      }
    }
  } else {
    grp_at_t[] <- subject_group
  }

  # ── Step 3: pooled within-group rolling-window ridge -> Gamma ────────
  gamma_hat <- matrix(NA_real_, K, m)
  Gamma_hat <- array(NA_real_, dim = c(K, m, p_x))

  for (k in seq_len(K)) {
    for (t in seq_len(m)) {
      idx_k <- which(grp_at_t[, t] == k)
      if (length(idx_k) == 0L) next
      s_win <- win_idx[[t]]
      # stack all (i in G_k, s in W(t)) observations
      XtX <- matrix(0, p_x, p_x); Xty <- numeric(p_x); n_obs <- 0L
      for (i in idx_k) {
        X_i <- matrix(X[i, s_win, ], nrow = length(s_win), ncol = p_x)
        r_i <- R[i, s_win]
        ok  <- which(!is.na(r_i) & rowSums(is.na(X_i)) == 0L)
        if (length(ok) == 0L) next
        Xo  <- X_i[ok, , drop = FALSE]
        XtX <- XtX + crossprod(Xo)
        Xty <- Xty + as.numeric(crossprod(Xo, r_i[ok]))
        n_obs <- n_obs + length(ok)
      }
      if (n_obs < p_x + 1L) next
      g_kt <- tryCatch(as.numeric(solve(XtX + pen_r, Xty)),
                       error = function(e) rep(NA_real_, p_x))
      Gamma_hat[k, t, ] <- g_kt
      gamma_hat[k, t]   <- if (p_x == 1L) g_kt[1L] else sqrt(sum(g_kt^2))
    }
  }

  # Backfill beta NAs (subjects/times where the individual window was too
  # thin) from the per-t assigned group's pooled trajectory.  In trailing
  # mode grp_at_t is the running (causal) label, so backfilled values also
  # respect no-lookahead.
  for (i in seq_len(n)) {
    na_t <- which(is.na(beta_hat_arr[i, , 1L]))
    for (t in na_t) {
      g <- grp_at_t[i, t]
      if (!is.na(g)) beta_hat_arr[i, t, ] <- Gamma_hat[g, t, ]
    }
  }

  list(subject_group = subject_group,
       gamma_hat     = gamma_hat,
       Gamma_hat     = Gamma_hat,
       Gamma         = Gamma_hat,   # alias: runner / .adapt_fit expect $Gamma
       beta_hat_arr  = beta_hat_arr,
       group_matrix  = group_matrix,
       alpha         = alpha,
       w             = w,
       side          = side,
       time_sec      = proc.time()["elapsed"] - t0_clock)
}

# All permutations of 1..K_ as a list of integer vectors (K_! total).
.all_perms <- function(K_) {
  if (K_ == 1L) return(list(1L))
  out <- list()
  for (i in seq_len(K_)) {
    sub  <- .all_perms(K_ - 1L)
    rest <- setdiff(seq_len(K_), i)
    for (s in sub) out[[length(out) + 1L]] <- c(i, rest[s])
  }
  out
}