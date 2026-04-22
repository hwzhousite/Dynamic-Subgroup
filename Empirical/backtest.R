##############################################################################
#  backtest.R  —  General rolling-window backtesting framework
#
#  Core design:
#    backtest(Y, X, MV, WIDTH, beta_estimator, ...)
#      ├─ For each t = WIDTH+1 … T:
#      │    beta_hat_t  <- beta_estimator(Y[, 1:(t-1)], X[,1:(t-1),], t, ...)
#      │    y_hat_t     <- X[, t, ] %*% beta_hat_t  (firm-level prediction)
#      │    return_t    <- portfolio(y_hat_t, Y[, t], MV[, t])
#      └─ Aggregate into H-L decile returns, Sharpe, cumulative log-return
#
#  beta_estimator is a user-supplied function with signature:
#    f(Y_win, X_win, t_end, ...) -> named list:
#      $beta  : n × (p+1) matrix  — per-firm predicted beta at t_end
#               (or p+1 vector if all firms share the same beta)
#      (any other elements are ignored)
#
#  Estimators are defined in estimators.R (sourced automatically).
#  Add new estimators there without modifying this file.
##############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

## Load beta estimators — place estimators.R in the same directory
source("estimators.R")

## ══════════════════════════════════════════════════════════════════════════════
##  SECTION 1 — Portfolio utilities (exact originals)
## ══════════════════════════════════════════════════════════════════════════════

portfolio <- function(et, rt, wt) {
  if (length(et) != length(rt)) stop("portfolio: length(et) != length(rt)")
  if (missing(wt)) {
    dt  <- data.frame(et=et, rt=rt)
    if (any(is.na(dt))) dt <- dt[complete.cases(dt), ]
    rs  <- dt$rt[order(dt$et, decreasing=TRUE)]
    N   <- length(rs); dec <- round(N / 10)
    c(colMeans(matrix(rs[1:(dec*9)], dec, 9)), mean(rs[-(1:(dec*9))]))
  } else {
    dt  <- data.frame(et=et, rt=rt, wt=wt)
    if (any(is.na(dt))) dt <- dt[complete.cases(dt), ]
    dt  <- dt[order(dt$et, decreasing=TRUE), ]
    rw  <- dt$rt * dt$wt
    N   <- length(rw); dec <- round(N / 10)
    sr  <- c(colSums(matrix(rw[1:(dec*9)], dec, 9)),  sum(rw[-(1:(dec*9))]))
    sw  <- c(colSums(matrix(dt$wt[1:(dec*9)], dec, 9)), sum(dt$wt[-(1:(dec*9))]))
    sr / sw
  }
}

sharpe_ann <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2 || sd(x) == 0) return(NA_real_)
  mean(x) / sd(x) * sqrt(12)
}

## ══════════════════════════════════════════════════════════════════════════════
##  SECTION 2 — Core backtest engine
## ══════════════════════════════════════════════════════════════════════════════

#' @param Y        n × T return matrix
#' @param X        n × T × p characteristic array (NO intercept column)
#' @param MV       n × T log market-cap matrix (for VW weights), or NULL
#' @param WIDTH    integer — rolling window width (months)
#' @param beta_estimator  function(Y_hist, X_hist, t_end, ...) -> list($beta)
#'                        beta can be (p+1)-vector (shared) or n×(p+1) matrix
#' @param N_DECILES integer — number of portfolio deciles (default 10)
#' @param ROLL_K    integer — rolling Sharpe window in months (default 12)
#' @param verbose   logical — print progress (default TRUE)
#' @param ...       extra arguments passed to beta_estimator
#'
#' @return list with:
#'   $pf_eq, $pf_vw   : T × N_DECILES return matrices
#'   $HL_eq, $HL_vw   : H-L return vectors (length T - WIDTH)
#'   $window_records  : data.frame with dates, H-L returns, rolling Sharpe
#'   $sharpe_eq, $sharpe_vw : final expanding Sharpe on H-L
#'   $beta_arr        : n × T × (p+1) array of predicted betas (may be large)

backtest <- function(Y, X, MV        = NULL,
                     WIDTH            = 120L,
                     beta_estimator   = est_homogeneous,
                     N_DECILES        = 10L,
                     ROLL_K           = 12L,
                     dates            = NULL,
                     verbose          = TRUE,
                     ...) {

  n  <- nrow(Y); TT <- ncol(Y); p <- dim(X)[3]
  T_start <- WIDTH + 1L

  if (T_start > TT)
    stop(sprintf("backtest: WIDTH=%d >= TT=%d", WIDTH, TT))

  pf_eq   <- matrix(NA_real_, TT, N_DECILES)
  pf_vw   <- matrix(NA_real_, TT, N_DECILES)
  beta_arr <- array(NA_real_, dim=c(n, TT, p+1L))  # store predicted betas

  if (verbose)
    message(sprintf("[backtest] T_start=%d | T_end=%d | WIDTH=%d | n=%d | p=%d",
                    T_start, TT, WIDTH, n, p))

  pb <- if (verbose) txtProgressBar(min=T_start, max=TT, style=3) else NULL

  for (t in seq(T_start, TT)) {

    ## ── Extract causal window [t-WIDTH, t-1] ──────────────────────────────
    t_win <- seq(t - WIDTH, t - 1L)
    Y_win <- Y[, t_win, drop=FALSE]                    # n × WIDTH
    X_win <- X[, t_win, , drop=FALSE]                  # n × WIDTH × p
    MV_win <- if (!is.null(MV)) MV[, t_win, drop=FALSE] else NULL

    ## ── Estimate beta for month t ──────────────────────────────────────────
    est <- tryCatch(
      beta_estimator(Y_win, X_win, t, MV_hist=MV_win, ...),
      error = function(e) { message(sprintf("  [t=%d] estimator error: %s", t, e$message)); NULL }
    )
    if (is.null(est) || is.null(est$beta)) {
      if (!is.null(pb)) setTxtProgressBar(pb, t)
      next
    }

    beta_t <- est$beta   # either (p+1) vector or n × (p+1) matrix

    ## ── Broadcast scalar beta to all firms ────────────────────────────────
    if (is.vector(beta_t)) {
      if (length(beta_t) != p+1L) { if (!is.null(pb)) setTxtProgressBar(pb, t); next }
      beta_mat <- matrix(beta_t, nrow=n, ncol=p+1L, byrow=TRUE)
    } else {
      if (!identical(dim(beta_t), c(n, p+1L))) {
        if (!is.null(pb)) setTxtProgressBar(pb, t); next
      }
      beta_mat <- beta_t
    }
    beta_arr[, t, ] <- beta_mat

    ## ── Predict y_hat at t ────────────────────────────────────────────────
    y_t   <- Y[, t]
    X_t   <- X[, t, ]                                  # n × p

    ok <- !is.na(y_t) & rowSums(is.na(X_t)) == 0 &
          rowSums(is.na(beta_mat)) == 0 &
          rowSums(!is.finite(beta_mat)) == 0
    if (sum(ok) < N_DECILES * 2L) {
      if (!is.null(pb)) setTxtProgressBar(pb, t); next
    }

    fi   <- which(ok)
    Xd   <- cbind(1, X_t[fi, ])                        # n_ok × (p+1)
    yhat <- rowSums(Xd * beta_mat[fi, ])                # predicted return

    pf_eq[t, ] <- portfolio(et=yhat, rt=y_t[fi])

    if (!is.null(MV)) {
      mv_t <- exp(MV[fi, t])
      mv_t[!is.finite(mv_t) | mv_t <= 0] <- NA
      pf_vw[t, ] <- portfolio(et=yhat, rt=y_t[fi], wt=mv_t)
    }

    if (!is.null(pb)) setTxtProgressBar(pb, t)
  }
  if (!is.null(pb)) close(pb)

  ## ── Aggregate results ─────────────────────────────────────────────────────
  eval_idx <- seq(T_start, TT)
  HL_eq    <- pf_eq[eval_idx, 1] - pf_eq[eval_idx, N_DECILES]
  HL_vw    <- pf_vw[eval_idx, 1] - pf_vw[eval_idx, N_DECILES]

  dt_labels <- if (!is.null(dates)) format(dates[eval_idx], "%Y-%m") else
               as.character(eval_idx)

  sharpe_exp_eq <- sharpe_roll_eq <- sharpe_exp_vw <- sharpe_roll_vw <-
    rep(NA_real_, length(eval_idx))

  for (i in seq_along(eval_idx)) {
    sharpe_exp_eq[i]  <- sharpe_ann(HL_eq[1:i])
    sharpe_exp_vw[i]  <- sharpe_ann(HL_vw[1:i])
    ri <- max(1L, i-ROLL_K+1L):i
    sharpe_roll_eq[i] <- sharpe_ann(HL_eq[ri])
    sharpe_roll_vw[i] <- sharpe_ann(HL_vw[ri])
  }

  window_records <- data.frame(
    t_idx          = eval_idx,
    date           = dt_labels,
    hl_ret_eq      = HL_eq,
    hl_ret_vw      = HL_vw,
    sharpe_exp_eq  = sharpe_exp_eq,
    sharpe_roll_eq = sharpe_roll_eq,
    sharpe_exp_vw  = sharpe_exp_vw,
    sharpe_roll_vw = sharpe_roll_vw,
    stringsAsFactors = FALSE
  )

  list(
    pf_eq          = pf_eq,
    pf_vw          = pf_vw,
    HL_eq          = HL_eq,
    HL_vw          = HL_vw,
    window_records = window_records,
    sharpe_eq      = sharpe_ann(HL_eq),
    sharpe_vw      = sharpe_ann(HL_vw),
    beta_arr       = beta_arr,
    eval_idx       = eval_idx
  )
}

## ══════════════════════════════════════════════════════════════════════════════
##  SECTION 3 — Summary and plot helpers
## ══════════════════════════════════════════════════════════════════════════════

#' Print a compact Sharpe + decile summary for a backtest result
summarise_backtest <- function(res, label = "Method", N_DECILES = 10L) {
  wr  <- res$window_records
  idx <- res$eval_idx
  cat(sprintf("\n══ %s ══\n", label))
  cat(sprintf("  EW H-L Sharpe : %.4f\n", res$sharpe_eq))
  cat(sprintf("  VW H-L Sharpe : %.4f\n", res$sharpe_vw))
  cat(sprintf("  Eval period   : %s to %s (%d months)\n",
              wr$date[1], wr$date[nrow(wr)], nrow(wr)))

  mk_tbl <- function(pf) {
    avg <- colMeans(pf[idx, ], na.rm=TRUE)*100
    sdv <- apply(pf[idx, ], 2, sd, na.rm=TRUE)*100
    sr  <- round(avg/sdv*sqrt(12), 3)
    cbind(AVG=round(rev(avg),3), SD=round(rev(sdv),3), SR=rev(sr))
  }
  HL_eq <- res$HL_eq; HL_vw <- res$HL_vw
  r_eq <- rbind(mk_tbl(res$pf_eq),
                `H-L`=round(c(mean(HL_eq,na.rm=T)*100, sd(HL_eq,na.rm=T)*100,
                               res$sharpe_eq),3))
  r_vw <- rbind(mk_tbl(res$pf_vw),
                `H-L`=round(c(mean(HL_vw,na.rm=T)*100, sd(HL_vw,na.rm=T)*100,
                               res$sharpe_vw),3))
  rownames(r_eq) <- rownames(r_vw) <- c("H", as.character((N_DECILES-1):2), "L", "H-L")
  cat("\n── Equal-Weighted ──\n"); print(r_eq)
  cat("\n── Value-Weighted ──\n"); print(r_vw)
  invisible(list(eq=r_eq, vw=r_vw))
}

#' Cumulative log-return plot for one or more backtest results
plot_backtest <- function(results_list, dates = NULL, N_DECILES = 10L,
                          title = "Cumulative Log Returns") {
  cumlog <- function(x) cumsum(log(1 + ifelse(is.na(x), 0, x)))

  df_all <- purrr::map_dfr(names(results_list), function(nm) {
    res      <- results_list[[nm]]
    idx      <- res$eval_idx
    dt_vals  <- if (!is.null(dates)) as.Date(paste0(format(dates[idx], "%Y-%m"), "-01")) else
                as.Date("1970-01-01") + idx
    rbind(
      data.frame(date=dt_vals, cum_ret=cumlog(res$pf_eq[idx, 1]),
                 method=nm, position="Long",  weighting="EW"),
      data.frame(date=dt_vals, cum_ret=cumlog(res$pf_eq[idx, N_DECILES]),
                 method=nm, position="Short", weighting="EW"),
      data.frame(date=dt_vals, cum_ret=cumlog(res$pf_vw[idx, 1]),
                 method=nm, position="Long",  weighting="VW"),
      data.frame(date=dt_vals, cum_ret=cumlog(res$pf_vw[idx, N_DECILES]),
                 method=nm, position="Short", weighting="VW")
    )
  })

  ggplot(df_all, aes(date, cum_ret, colour=method, linetype=position)) +
    geom_line(linewidth=0.7) +
    geom_hline(yintercept=0, colour="black", linewidth=0.3) +
    facet_wrap(~weighting, ncol=2) +
    scale_linetype_manual(values=c(Long="solid", Short="dotted")) +
    scale_x_date(date_breaks="10 years", date_labels="%Y") +
    labs(title=title, x=NULL, y="Cumulative log return",
         colour="Method", linetype="Position") +
    theme_bw(base_size=11) +
    theme(plot.title=element_text(face="bold"), panel.grid.minor=element_blank(),
          legend.position="bottom")
}

## ══════════════════════════════════════════════════════════════════════════════
##  SECTION 4 — Example usage
## ══════════════════════════════════════════════════════════════════════════════
if (FALSE) {  # set to TRUE to run

  wide <- readRDS("/Users/haowen/Dropbox/Dynamic Subgroup/DSHR/weekly/Dec28 Empirical/crsp_wide.rds")
  Y      <- wide$Y
  X      <- wide$X
  MV     <- wide$X[, , which(wide$chars == "MV")]
  dates  <- wide$dates

  ## ── Run homogeneous baseline ───────────────────────────────────────────────
  res_homo <- backtest(
    Y = Y, X = X, MV = MV, WIDTH = 120L,
    beta_estimator = est_homogeneous,
    dates = dates
  )
  summarise_backtest(res_homo, label = "Homogeneous OLS")

  ## ── Run segmented spline ───────────────────────────────────────────────────
  res_seg <- backtest(
    Y = Y, X = X, MV = MV, WIDTH = 120L,
    beta_estimator = est_seg_spline,
    dates = dates,
    # extra args forwarded to est_seg_spline:
    SIZE_BREAKS    = seq(0.10, 0.90, by=0.10),
    SPLINE_DF      = NULL,
    SPLINE_MIN_OBS = 24L,
    K_AVG          = 12L,
    exchcd_hist    = wide$exchcd   # NULL if not available
  )
  summarise_backtest(res_seg, label = "Segmented Spline (K=10)")

  ## ── Compare ───────────────────────────────────────────────────────────────
  p_comp <- plot_backtest(
    list("Homogeneous" = res_homo, "Seg-Spline K=10" = res_seg),
    dates = dates,
    title = "Backtest Comparison — CRSP 1964-2022"
  )
  ggsave("plots/backtest_comparison.png", p_comp, width=12, height=5, dpi=150)

  ## ── Save ──────────────────────────────────────────────────────────────────
  saveRDS(list(homo = res_homo, seg = res_seg), "results/backtest_results.rds")
}
