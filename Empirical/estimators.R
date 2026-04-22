##############################################################################
#  estimators.R  —  Beta estimator functions for the backtest framework
#
#  Each function follows the contract:
#    f(Y_hist, X_hist, t_end, ...) -> list($beta)
#
#  Arguments:
#    Y_hist    : n × W matrix — returns in window [t-WIDTH, t-1]
#    X_hist    : n × W × p array — characteristics in same window
#    t_end     : integer — prediction target month (NOT included in Y/X_hist)
#    MV_hist   : n × W matrix — log market-cap (optional, for segmentation)
#    ...       : estimator-specific parameters
#
#  Return value:
#    list($beta) where $beta is either:
#      (p+1) vector  — shared beta, broadcast to all n firms
#      n × (p+1) matrix — firm-specific beta
#
#  Available estimators:
#    est_homogeneous   pooled OLS + causal rolling mean
#    est_seg_spline    segmented OLS (K size deciles) + rolling-mean + spline
#
#  To add a new estimator, define a function here with the same signature
#  and source() this file alongside backtest.R.
##############################################################################

## ── est_homogeneous ──────────────────────────────────────────────────────────
#  Pooled cross-sectional OLS at each window month w, then causal rolling mean.
#  Returns a single (p+1) vector broadcast to all firms.
est_homogeneous <- function(Y_hist, X_hist, t_end, ...) {
  W <- ncol(Y_hist); p <- dim(X_hist)[3]
  beta_lm <- matrix(NA_real_, W, p + 1L)

  for (w in seq_len(W)) {
    y_w <- Y_hist[, w];  X_w <- X_hist[, w, ]
    ok  <- !is.na(y_w) & rowSums(is.na(X_w)) == 0
    if (sum(ok) < p + 2L) next
    fit <- tryCatch(lm.fit(cbind(1, X_w[ok, ]), y_w[ok]), error = function(e) NULL)
    if (!is.null(fit) && all(is.finite(coef(fit))))
      beta_lm[w, ] <- coef(fit)
  }

  valid <- beta_lm[!is.na(beta_lm[, 1]), , drop = FALSE]
  if (nrow(valid) == 0) return(list(beta = rep(NA_real_, p + 1L)))
  list(beta = colMeans(valid))
}

## ── est_seg_spline ───────────────────────────────────────────────────────────
#  Segmented OLS (K NYSE size deciles) at each window month w,
#  then causal rolling-mean history, then smoothing spline evaluated at
#  the last K_AVG window points (averaged).
#  Returns an n × (p+1) matrix — each firm gets its size-segment's beta.
#
#  Parameters (passed via ...):
#    MV_hist      n × W log-MV matrix           [required for segmentation]
#    SIZE_BREAKS  numeric vector of quantile cuts [default: deciles 0.1..0.9]
#    SPLINE_DF    integer or NULL (GCV)           [default: NULL]
#    SPLINE_MIN_OBS integer                       [default: 24L]
#    K_AVG        integer                         [default: 12L]
#    exchcd_hist  n × W exchange-code matrix      [optional, for NYSE filter]
est_seg_spline <- function(Y_hist, X_hist, t_end,
                           MV_hist        = NULL,
                           SIZE_BREAKS    = seq(0.10, 0.90, by = 0.10),
                           SPLINE_DF      = NULL,
                           SPLINE_MIN_OBS = 24L,
                           K_AVG          = 12L,
                           exchcd_hist    = NULL,
                           ...) {
  n_firms <- nrow(Y_hist)
  W       <- ncol(Y_hist)
  p       <- dim(X_hist)[3]
  N_SEG   <- length(SIZE_BREAKS) + 1L

  ## Step 1 — per-segment OLS at each window month ───────────────────────────
  beta_lm_seg <- array(NA_real_, dim = c(N_SEG, W, p + 1L))
  seg_hist    <- matrix(NA_integer_, n_firms, W)

  for (w in seq_len(W)) {
    y_w  <- Y_hist[, w];  X_w <- X_hist[, w, ]
    mv_w <- if (!is.null(MV_hist)) MV_hist[, w] else NULL
    ok   <- !is.na(y_w) & rowSums(is.na(X_w)) == 0
    if (!is.null(mv_w)) ok <- ok & !is.na(mv_w)
    if (sum(ok) < N_SEG * (p + 2L)) next

    firm_idx <- which(ok)
    mv_ok    <- if (!is.null(mv_w)) mv_w[firm_idx] else NULL

    ## NYSE or all-stock breakpoints
    if (!is.null(mv_ok)) {
      if (!is.null(exchcd_hist)) {
        ex_w   <- exchcd_hist[firm_idx, w]
        ref_mv <- mv_ok[!is.na(ex_w) & ex_w == 1L]
        if (length(ref_mv) < 10L) ref_mv <- mv_ok
      } else {
        ref_mv <- mv_ok
      }
      brk_w <- quantile(ref_mv, probs = SIZE_BREAKS, na.rm = TRUE)
      seg_w <- as.integer(cut(mv_ok, breaks = c(-Inf, brk_w, Inf),
                               labels = FALSE, right = TRUE))
    } else {
      seg_w <- rep(1L, length(firm_idx))
    }
    seg_hist[firm_idx, w] <- seg_w

    for (s in seq_len(N_SEG)) {
      idx_s <- firm_idx[seg_w == s]
      if (length(idx_s) < p + 2L) next
      fit_s <- tryCatch(
        lm.fit(cbind(1, X_w[idx_s, ]), y_w[idx_s]),
        error = function(e) NULL
      )
      if (!is.null(fit_s) && all(is.finite(coef(fit_s))))
        beta_lm_seg[s, w, ] <- coef(fit_s)
    }
  }

  ## Step 2 — rolling-mean history + spline at last K_AVG points ────────────
  beta_sp_seg <- matrix(NA_real_, N_SEG, p + 1L)

  for (s in seq_len(N_SEG)) {
    for (j in seq_len(p + 1L)) {
      raw_s <- beta_lm_seg[s, , j]   # W-vector of monthly OLS betas

      ## Causal rolling-mean history: rw_hist[w] = mean(raw_s[1..w])
      rw_hist <- rep(NA_real_, W)
      for (w in seq_len(W)) {
        v <- raw_s[seq_len(w)]
        v <- v[!is.na(v) & is.finite(v)]
        if (length(v) > 0L) rw_hist[w] <- mean(v)
      }

      ok_h <- which(!is.na(rw_hist) & is.finite(rw_hist))
      if (length(ok_h) < SPLINE_MIN_OBS) next

      sp_fit <- tryCatch({
        if (is.null(SPLINE_DF))
          smooth.spline(x = ok_h, y = rw_hist[ok_h], cv = FALSE)
        else
          smooth.spline(x = ok_h, y = rw_hist[ok_h],
                        df = min(SPLINE_DF, length(ok_h) - 1L))
      }, error = function(e) NULL)
      if (is.null(sp_fit)) next

      avg_pts <- seq(max(1L, W - K_AVG + 1L), W)
      pv <- tryCatch(predict(sp_fit, x = avg_pts)$y, error = function(e) NULL)
      if (!is.null(pv) && all(is.finite(pv)))
        beta_sp_seg[s, j] <- mean(pv)
    }
  }

  ## Step 3 — assign each firm its segment's beta (from last window month) ───
  seg_last <- seg_hist[, W]
  beta_out <- matrix(NA_real_, n_firms, p + 1L)
  for (i in seq_len(n_firms)) {
    s_i <- seg_last[i]
    if (!is.na(s_i) && s_i >= 1L && s_i <= N_SEG)
      beta_out[i, ] <- beta_sp_seg[s_i, ]
  }
  list(beta = beta_out)
}
