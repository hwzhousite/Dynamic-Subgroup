##############################################################################
#  03_rolling_window.R
#  Rolling-window out-of-sample evaluation engine
#
#  For each prediction month t in the evaluation period:
#    1. Extract training window [t - WIN_SIZE, t-1]
#    2. Fit each method independently on the training window
#    3. Predict E[ret] at t via X[i,t,] · beta_agg[i,]
#    4. Form H-L portfolio: long top TOP_PCTILE, short bottom BOT_PCTILE
#    5. Record portfolio return and diagnostics
#
#  Each method saves its own files immediately after finishing:
#    results/{Method}_results.rds   results/{Method}_returns.csv
#  Combined files written at the end:
#    results/rolling_results.rds    results/portfolio_returns.csv
##############################################################################

suppressPackageStartupMessages({
  library(dplyr)
})

source("2 method.R")

## ── config ────────────────────────────────────────────────────────────────────
WIN_SIZE    <- 120L   # training window length (months)
N_AGG       <- 36L    # AFEM/SplineKM aggregation window (months)
PRED_MONTHS <- 12L    # number of out-of-sample prediction months
TOP_PCTILE  <- 0.10   # long  : top 10%
BOT_PCTILE  <- 0.10   # short : bottom 10%
MIN_OBS_WIN <- 12L    # min active months a firm needs in window to be included
MIN_ACTIVE  <- 20L    # min active firms per prediction month

dir.create("results", showWarnings = FALSE)

## ── load wide panel ───────────────────────────────────────────────────────────
message("[rolling] loading data ...")
wide <- readRDS("/Users/haowen/Dropbox/Dynamic\ Subgroup/DSHR/weekly/Dec28\ Empirical/crsp_wide.rds")
Y    <- wide$Y          # n × T
X    <- wide$X          # n × T × p
n    <- nrow(Y); TT <- ncol(Y); p <- dim(X)[3]
message(sprintf("[rolling] panel: n=%d firms | T=%d months | p=%d chars", n, TT, p))

## ── method registry ───────────────────────────────────────────────────────────
# Method names use underscores — safe as filename components
METHODS <- list(
  list(name   = "Homogeneous_FM",
       fn     = fit_homogeneous,
       params = list()),

  list(name   = "Static_Hetero",
       fn     = fit_static_hetero,
       params = list(K_max = 5)),

  list(name   = "AFEM",
       fn     = fit_afem_empirical,
       params = list(K_max = 5, df_spline = 6, lam1 = 0.08, lam2 = 0.04,
                     n_agg = N_AGG, max_iter = 40, tol = 1e-4)),

  list(name   = "SplineKM",
       fn     = fit_splinekm_empirical,
       params = list(K_max = 5, df_spline = 6, n_agg = N_AGG))
)

# Prediction months: last PRED_MONTHS months of the panel
pred_idx <- seq(TT - PRED_MONTHS + 1L, TT)
message(sprintf("[rolling] prediction period: months %d-%d  (%s to %s)",
                pred_idx[1], pred_idx[length(pred_idx)],
                format(wide$dates[pred_idx[1]], "%Y-%m"),
                format(wide$dates[pred_idx[length(pred_idx)]], "%Y-%m")))

## ── helper: evaluate one (method × month) ────────────────────────────────────
evaluate_one_month <- function(mth, t) {
  t_win  <- max(1L, t - WIN_SIZE):(t - 1L)
  if (length(t_win) < 24) return(NULL)

  # Active firms: observed at t AND have >= MIN_OBS_WIN obs in the window
  active <- which(!is.na(Y[, t]) &
                    rowSums(!is.na(Y[, t_win])) >= MIN_OBS_WIN)
  if (length(active) < MIN_ACTIVE) return(NULL)

  Y_win <- Y[active, t_win, drop = FALSE]
  X_win <- X[active, t_win, , drop = FALSE]
  Y_t   <- Y[active, t]

  fit      <- mth$fn(Y_win, X_win, mth$params)
  beta_agg <- fit$beta_agg                       # n_active × p
  X_t      <- X[active, t, , drop = FALSE]       # n_active × p

  # Predicted return: dot product row-wise
  pred_r <- rowSums(X_t * beta_agg, na.rm = FALSE)
  pred_r[rowSums(is.na(X_t)) > 0] <- NA          # NA out if any char missing

  valid <- !is.na(pred_r) & !is.na(Y_t)
  if (sum(valid) < MIN_ACTIVE) return(NULL)

  ranks     <- rank(pred_r[valid], ties.method = "average")
  long_idx  <- which(ranks >= quantile(ranks, 1 - TOP_PCTILE))
  short_idx <- which(ranks <= quantile(ranks, BOT_PCTILE))
  if (length(long_idx) == 0 || length(short_idx) == 0) return(NULL)

  data.frame(
    method    = mth$name,
    t_idx     = t,
    date      = format(wide$dates[t], "%Y-%m"),
    hl_ret    = mean(Y_t[valid][long_idx]) - mean(Y_t[valid][short_idx]),
    long_ret  = mean(Y_t[valid][long_idx]),
    short_ret = mean(Y_t[valid][short_idx]),
    n_long    = length(long_idx),
    n_short   = length(short_idx),
    n_active  = length(active),
    K_used    = fit$K_used %||% NA_integer_,
    stringsAsFactors = FALSE
  )
}

## ── run one method over all prediction months ─────────────────────────────────
# pred_idx is passed explicitly — no hidden global dependency
run_method <- function(mth, pred_idx) {
  label  <- mth$name
  n_pred <- length(pred_idx)
  message(sprintf("\n[rolling] ══ START: %-15s  (%d prediction months) ══════",
                  label, n_pred))
  t0      <- proc.time()["elapsed"]
  records <- list()

  for (i in seq_along(pred_idx)) {
    t <- pred_idx[i]
    rec <- tryCatch(
      evaluate_one_month(mth, t),
      error = function(e) {
        message(sprintf("  [!] %s | t=%d (%s): %s",
                        label, t,
                        format(wide$dates[t], "%Y-%m"),
                        conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(rec)) records[[length(records) + 1]] <- rec

    # Progress every 6 iterations and on final step
    if (i %% 6 == 0 || i == n_pred)
      message(sprintf("  [%s] iter %d/%d | t=%d (%s) | recorded=%d",
                      label, i, n_pred, t,
                      format(wide$dates[t], "%Y-%m"),
                      length(records)))
  }

  elapsed <- proc.time()["elapsed"] - t0

  if (length(records) == 0) {
    message(sprintf("[rolling] %s: no valid records produced (%.1fs)", label, elapsed))
    return(invisible(NULL))
  }

  df <- do.call(rbind, records)

  # Per-method save
  rds_path <- sprintf("results/%s_results.rds", label)
  csv_path <- sprintf("results/%s_returns.csv", label)
  saveRDS(df, rds_path)
  write.csv(df, csv_path, row.names = FALSE)

  # Annualised Sharpe (guard: n>=2 and sd>0)
  sr <- if (nrow(df) >= 2 && sd(df$hl_ret) > 0)
    round((mean(df$hl_ret) / sd(df$hl_ret)) * sqrt(12), 3)
  else NA_real_

  message(sprintf(
    "[rolling] %s: %d months | mean HL=%.4f | Sharpe=%s | -> %s  (%.1fs)",
    label, nrow(df), mean(df$hl_ret),
    if (is.na(sr)) "NA" else sprintf("%.3f", sr),
    rds_path, elapsed))
  df
}

## ── execute all methods independently ────────────────────────────────────────
results_list <- list()
for (mth in METHODS) {
  results_list[[mth$name]] <- run_method(mth, pred_idx)
}

## ── combined save ─────────────────────────────────────────────────────────────
valid_results <- Filter(Negate(is.null), results_list)

if (length(valid_results) > 0) {
  ret_df <- do.call(rbind, valid_results)
  saveRDS(ret_df, "results/rolling_results.rds")
  write.csv(ret_df, "results/portfolio_returns.csv", row.names = FALSE)
  message(sprintf(
    "\n[rolling] combined save: %d records | %d methods | %d pred months",
    nrow(ret_df), length(valid_results), PRED_MONTHS))
} else {
  warning("[rolling] no results produced by any method")
}

## ── final summary ─────────────────────────────────────────────────────────────
message("\n[rolling] === Per-method summary ===")
for (nm in names(results_list)) {
  df <- results_list[[nm]]
  if (is.null(df)) {
    message(sprintf("  %-18s  NO OUTPUT", nm))
  } else {
    sr <- if (nrow(df) >= 2 && sd(df$hl_ret) > 0)
      round((mean(df$hl_ret) / sd(df$hl_ret)) * sqrt(12), 3) else NA_real_
    message(sprintf("  %-18s  %2d months | Sharpe = %s  [results/%s_results.rds]",
                    nm, nrow(df),
                    if (is.na(sr)) "  NA  " else sprintf("%6.3f", sr),
                    nm))
  }
}
