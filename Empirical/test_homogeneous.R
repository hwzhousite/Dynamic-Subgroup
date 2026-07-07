##############################################################################
#  test_homogeneous.R
#  Homogeneous linear model — replicates Reg15_v2_Dec2022.R logic
#  adapted to the crsp_wide.rds panel format.
#
#  Step 1 — For every month t, fit cross-sectional OLS: ret ~ X
#            → store beta.lm  (T × (p+1))
#
#  Step 2 — Rolling average of betas over WIDTH months
#            → beta.rw  (T × (p+1))
#
#  Step 3 — For t = WIDTH+1 … T, predict with beta.rw[t-1,]
#            Form equal-weighted AND value-weighted decile portfolios
#            Record decile returns and H-L each month
#
#  Step 4 — Compute per-rolling-window Sharpe (expanding & rolling-K)
#            and final summary table matching the original output format
#
#  Prerequisites:
#    data1964_2022 in workspace  OR  data/crsp_wide.rds already exists
##############################################################################

suppressPackageStartupMessages({
  library(dplyr)
})

## ── config ────────────────────────────────────────────────────────────────────
WIDTH         <- 60L   # rolling average window (= estimation window)
N_DECILES     <- 10L    # number of portfolio deciles
ROLL_SHARPE_K <- 12L    # window for rolling Sharpe display

dir.create("results", showWarnings = FALSE)

## ── load panel ────────────────────────────────────────────────────────────────

Y      <- wide$Y          # n × T  returns
X      <- wide$X          # n × T × p  characteristics
MV_raw <- wide$X[, , which(wide$chars == "MV")]  # n × T  log(MV) for VW weights
n      <- nrow(Y); TT <- ncol(Y); p <- dim(X)[3]
message(sprintf("[test] panel: n=%d firms | T=%d months | p=%d chars", n, TT, p))
message(sprintf("[test] date range: %s to %s",
                format(wide$dates[1], "%Y-%m"),
                format(wide$dates[TT], "%Y-%m")))

## ── portfolio functions (exact originals from empirical functions script) ──────

portfolio <- function(et, rt, wt) {
  if (length(et) != length(rt)) stop("WARNING: lenghth different!")
  if (missing(wt)) {
    dt <- data.frame(et = et, rt = rt)
    if (sum(is.na(dt)) > 0) {
      dt2 <- dt[-which(is.na(rowSums(dt))), ]
    } else { dt2 <- dt }
    rt.sort <- (dt2$rt)[order(dt2$et, decreasing = TRUE)]
    rt.n    <- length(rt.sort)
    dec     <- round(rt.n / 10)
    rt.pf   <- c(colMeans(matrix(rt.sort[1:(dec * 9)], dec, 9)),
                 mean(rt.sort[-(1:(dec * 9))]))
  } else {
    dt <- data.frame(et = et, rt = rt, wt = wt)
    if (sum(is.na(dt)) > 0) {
      dt2 <- dt[-which(is.na(rowSums(dt))), ]
    } else { dt2 <- dt }
    dt.sort    <- dt2[order(dt2$et, decreasing = TRUE), ]
    rt.w.sort  <- dt.sort$rt * dt.sort$wt
    rt.n       <- length(rt.w.sort)
    dec        <- round(rt.n / 10)
    srt.pf     <- c(colSums(matrix(rt.w.sort[1:(dec * 9)], dec, 9)),
                    sum(rt.w.sort[-(1:(dec * 9))]))
    wt.pf      <- c(colSums(matrix(dt.sort$wt[1:(dec * 9)], dec, 9)),
                    sum(dt.sort$wt[-(1:(dec * 9))]))
    rt.pf      <- srt.pf / wt.pf
  }
  return(rt.pf)
}

portfolio.nyse <- function(et, rt, nyse, wt) {
  if (length(et) != length(rt)) stop("WARNING: lenghth different!")
  et.nyse  <- et[which(nyse == 1)]
  qt.nyse  <- quantile(et.nyse, prob = c(1:9) / 10, na.rm = TRUE)
  rt.pf    <- rep(NA, 10)
  if (missing(wt)) {
    rt.pf[1]  <- mean(rt[which(et < qt.nyse[1])], na.rm = TRUE)
    rt.pf[10] <- mean(rt[which(et > qt.nyse[9])], na.rm = TRUE)
    for (q in 2:9)
      rt.pf[q] <- mean(rt[which(et > qt.nyse[q - 1] & et < qt.nyse[q])],
                       na.rm = TRUE)
  } else {
    wt2       <- wt * as.numeric(!is.na(rt))
    rt.pf[1]  <- sum((rt * wt2)[which(et < qt.nyse[1])],  na.rm = TRUE) /
                 sum(wt2[which(et < qt.nyse[1])],          na.rm = TRUE)
    rt.pf[10] <- sum((rt * wt2)[which(et > qt.nyse[9])],  na.rm = TRUE) /
                 sum(wt2[which(et > qt.nyse[9])],          na.rm = TRUE)
    for (q in 2:9) {
      id.qt    <- which(et > qt.nyse[q - 1] & et < qt.nyse[q])
      rt.pf[q] <- sum((rt * wt2)[id.qt], na.rm = TRUE) /
                  sum(wt2[id.qt],         na.rm = TRUE)
    }
  }
  return(rt.pf)
}

## ══════════════════════════════════════════════════════════════════════════════
##  STEP 1 — cross-sectional OLS for every month t
## ══════════════════════════════════════════════════════════════════════════════
message("\n[Step 1] Cross-sectional OLS  (T=", TT, " months) ...")
beta.lm <- matrix(NA_real_, TT, p + 1L)   # intercept + p slopes

pb <- txtProgressBar(min = 1, max = TT, style = 3)
for (t in seq_len(TT)) {
  y_t <- Y[, t]
  X_t <- X[, t, ]                                  # n × p
  ok  <- !is.na(y_t) & rowSums(is.na(X_t)) == 0
  if (sum(ok) < p + 2L) { setTxtProgressBar(pb, t); next }

  fit <- tryCatch(
    lm.fit(cbind(1, X_t[ok, ]), y_t[ok]),
    error = function(e) NULL
  )
  if (!is.null(fit)) beta.lm[t, ] <- coef(fit)
  setTxtProgressBar(pb, t)
}
close(pb)
message(sprintf("[Step 1] done — %d / %d months with valid beta",
                sum(!is.na(beta.lm[, 1])), TT))

## ══════════════════════════════════════════════════════════════════════════════
##  STEP 2 — rolling average of betas (WIDTH-month window)
##           beta.rw[t,] = mean of beta.lm[(t-WIDTH+1):t, ]  (same as filter)
## ══════════════════════════════════════════════════════════════════════════════
message(sprintf("\n[Step 2] Rolling average (width=%d) ...", WIDTH))
conv.w  <- rep(1 / WIDTH, WIDTH)
beta.rw <- matrix(NA_real_, TT, p + 1L)
for (j in seq_len(p + 1L))
  beta.rw[, j] <- as.numeric(stats::filter(beta.lm[, j], sides = 1, filter = conv.w))
message("[Step 2] done")

## ══════════════════════════════════════════════════════════════════════════════
##  STEP 3 — rolling prediction and portfolio formation
##           predict at t using beta.rw[t-1,]  (no look-ahead)
## ══════════════════════════════════════════════════════════════════════════════
message(sprintf("\n[Step 3] Portfolio construction  (t = %d to %d) ...",
                WIDTH + 1L, TT))

pf.eq <- matrix(NA_real_, TT, N_DECILES)   # equal-weighted decile returns
pf.vw <- matrix(NA_real_, TT, N_DECILES)   # value-weighted decile returns

# For per-window Sharpe tracking
hl_eq_series <- rep(NA_real_, TT)
hl_vw_series <- rep(NA_real_, TT)

pb <- txtProgressBar(min = WIDTH + 1L, max = TT, style = 3)

for (t in seq(WIDTH + 1L, TT)) {
  if (any(is.na(beta.rw[t - 1L, ]))) { setTxtProgressBar(pb, t); next }

  y_t  <- Y[, t]
  X_t  <- X[, t, ]                                 # n × p
  ok   <- !is.na(y_t) & rowSums(is.na(X_t)) == 0
  if (sum(ok) < N_DECILES * 2L) { setTxtProgressBar(pb, t); next }

  y_hat <- as.numeric(cbind(1, X_t[ok, ]) %*% beta.rw[t - 1L, ])

  # Equal-weighted
  pf.eq[t, ] <- portfolio(et = y_hat, rt = y_t[ok])          # EW: wt missing

  # Value-weighted (weight = exp(log(MV)) = MV; log(MV) stored after prep)
  mv_t <- exp(MV_raw[ok, t])             # raw MV (exponentiate log MV stored in X)
  mv_t[!is.finite(mv_t) | mv_t <= 0] <- NA
  pf.vw[t, ] <- portfolio(et = y_hat, rt = y_t[ok], wt = mv_t)  # VW

  setTxtProgressBar(pb, t)
}
close(pb)
message("[Step 3] done")

## ══════════════════════════════════════════════════════════════════════════════
##  STEP 4 — H-L returns, per-window Sharpe, and summary table
## ══════════════════════════════════════════════════════════════════════════════
message("\n[Step 4] Computing Sharpe ratios ...")

# Evaluation period: drop first WIDTH months (no prediction there)
eval_idx <- seq(WIDTH + 1L, TT)

HL_eq <- pf.eq[eval_idx, 1] - pf.eq[eval_idx, N_DECILES]  # High - Low
HL_vw <- pf.vw[eval_idx, 1] - pf.vw[eval_idx, N_DECILES]

sharpe_ann <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2 || sd(x) == 0) return(NA_real_)
  (mean(x) / sd(x)) * sqrt(12)
}

# Build per-window record table (one row per prediction month)
window_records <- data.frame(
  window           = seq_along(eval_idx),
  t_idx            = eval_idx,
  date             = format(wide$dates[eval_idx], "%Y-%m"),
  hl_ret_eq        = HL_eq,
  hl_ret_vw        = HL_vw,
  stringsAsFactors = FALSE
)

# Expanding and rolling Sharpe computed after each window
sharpe_exp_eq  <- rep(NA_real_, nrow(window_records))
sharpe_roll_eq <- rep(NA_real_, nrow(window_records))
sharpe_exp_vw  <- rep(NA_real_, nrow(window_records))
sharpe_roll_vw <- rep(NA_real_, nrow(window_records))

for (i in seq_len(nrow(window_records))) {
  sharpe_exp_eq[i]  <- sharpe_ann(window_records$hl_ret_eq[1:i])
  sharpe_exp_vw[i]  <- sharpe_ann(window_records$hl_ret_vw[1:i])
  roll_idx          <- max(1L, i - ROLL_SHARPE_K + 1L):i
  sharpe_roll_eq[i] <- sharpe_ann(window_records$hl_ret_eq[roll_idx])
  sharpe_roll_vw[i] <- sharpe_ann(window_records$hl_ret_vw[roll_idx])
}

window_records$sharpe_exp_eq  <- sharpe_exp_eq
window_records$sharpe_roll_eq <- sharpe_roll_eq
window_records$sharpe_exp_vw  <- sharpe_exp_vw
window_records$sharpe_roll_vw <- sharpe_roll_vw

## ── summary table (mirrors original Reg15 output) ────────────────────────────
# Equal-weighted
avg_eq <- colMeans(pf.eq[eval_idx, ], na.rm = TRUE) * 100
sd_eq  <- apply(pf.eq[eval_idx, ], 2, sd, na.rm = TRUE) * 100
sr_eq  <- round(avg_eq / sd_eq * sqrt(12), 3)
res_eq <- cbind(AVG = round(rev(avg_eq), 3),
                SD  = round(rev(sd_eq),  3),
                SR  = rev(sr_eq))
hl_row_eq <- round(c(mean(HL_eq, na.rm = TRUE) * 100,
                     sd(HL_eq,   na.rm = TRUE) * 100,
                     sharpe_ann(HL_eq)), 3)
result_eq  <- rbind(res_eq, `H-L` = hl_row_eq)
rownames(result_eq) <- c("H", as.character(9:2), "L", "H-L")

# Value-weighted
avg_vw <- colMeans(pf.vw[eval_idx, ], na.rm = TRUE) * 100
sd_vw  <- apply(pf.vw[eval_idx, ], 2, sd, na.rm = TRUE) * 100
sr_vw  <- round(avg_vw / sd_vw * sqrt(12), 3)
res_vw <- cbind(AVG = round(rev(avg_vw), 3),
                SD  = round(rev(sd_vw),  3),
                SR  = rev(sr_vw))
hl_row_vw <- round(c(mean(HL_vw, na.rm = TRUE) * 100,
                     sd(HL_vw,   na.rm = TRUE) * 100,
                     sharpe_ann(HL_vw)), 3)
result_vw  <- rbind(res_vw, `H-L` = hl_row_vw)
rownames(result_vw) <- c("H", as.character(9:2), "L", "H-L")

## ── print results ─────────────────────────────────────────────────────────────
cat("\n")
cat("══════════════════════════════════════════════════════════════════\n")
cat("  Homogeneous Linear Model — Decile Portfolio Results\n")
cat(sprintf("  Evaluation: %s to %s  (%d months)\n",
            window_records$date[1],
            window_records$date[nrow(window_records)],
            nrow(window_records)))
cat("══════════════════════════════════════════════════════════════════\n\n")

cat("── Equal-Weighted ──────────────────────────────────────────────\n")
print(result_eq)

cat("\n── Value-Weighted ──────────────────────────────────────────────\n")
print(result_vw)

cat("\n── Combined (EW | VW) ──────────────────────────────────────────\n")
colnames(result_eq) <- paste0(colnames(result_eq), "_EW")
colnames(result_vw) <- paste0(colnames(result_vw), "_VW")
print(cbind(result_eq, result_vw))

cat("\n── Per-Window H-L Return & Sharpe ─────────────────────────────\n")
cat(sprintf("  %-10s  %+10s  %+10s  %+10s  %+10s  %+10s\n",
            "Date", "HL_EW", "SR_exp_EW", "SR_roll_EW", "HL_VW", "SR_exp_VW"))
cat(sprintf("  %s\n", strrep("-", 65)))
for (i in seq_len(nrow(window_records))) {
  r <- window_records[i, ]
  cat(sprintf("  %-10s  %+10.4f  %+10.4f  %+10.4f  %+10.4f  %+10.4f\n",
              r$date,
              ifelse(is.na(r$hl_ret_eq), 0, r$hl_ret_eq),
              ifelse(is.na(r$sharpe_exp_eq), 0, r$sharpe_exp_eq),
              ifelse(is.na(r$sharpe_roll_eq), 0, r$sharpe_roll_eq),
              ifelse(is.na(r$hl_ret_vw), 0, r$hl_ret_vw),
              ifelse(is.na(r$sharpe_exp_vw), 0, r$sharpe_exp_vw)))
}

## ── save ──────────────────────────────────────────────────────────────────────
saveRDS(list(window_records = window_records,
             result_eq      = result_eq,
             result_vw      = result_vw,
             pf_eq          = pf.eq,
             pf_vw          = pf.vw,
             beta_lm        = beta.lm,
             beta_rw        = beta.rw),
        "results/Homogeneous_FM_results.rds")
write.csv(window_records, "results/Homogeneous_FM_returns.csv", row.names = FALSE)
message("\n[test] saved results/Homogeneous_FM_results.rds")
message("[test] saved results/Homogeneous_FM_returns.csv")

## ══════════════════════════════════════════════════════════════════════════════
##  STEP 5 — Final Sharpe ratio summary + cumulative return plot
## ══════════════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages({ library(ggplot2); library(tidyr); library(dplyr) })
dir.create("plots", showWarnings = FALSE)

## ── 5a. Final Sharpe ratios ───────────────────────────────────────────────────
sr_HL_eq <- sharpe_ann(window_records$hl_ret_eq)   # annualised Sharpe H-L EW
sr_HL_vw <- sharpe_ann(window_records$hl_ret_vw)   # annualised Sharpe H-L VW

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║   FINAL SHARPE RATIOS (Annualised, H-L Portfolio)        ║\n")
cat("╠══════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Equal-Weighted  H-L Sharpe : %+.4f                    ║\n", sr_HL_eq))
cat(sprintf("║  Value-Weighted  H-L Sharpe : %+.4f                    ║\n", sr_HL_vw))
cat("╠══════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  EW  Mean H-L (ann. %%)      : %+.3f                     ║\n",
            mean(window_records$hl_ret_eq, na.rm = TRUE) * 12 * 100))
cat(sprintf("║  VW  Mean H-L (ann. %%)      : %+.3f                     ║\n",
            mean(window_records$hl_ret_vw, na.rm = TRUE) * 12 * 100))
cat(sprintf("║  EW  SD   H-L (ann. %%)      : %+.3f                     ║\n",
            sd(window_records$hl_ret_eq, na.rm = TRUE) * sqrt(12) * 100))
cat(sprintf("║  VW  SD   H-L (ann. %%)      : %+.3f                     ║\n",
            sd(window_records$hl_ret_vw, na.rm = TRUE) * sqrt(12) * 100))
cat(sprintf("║  Eval. period                : %s to %s           ║\n",
            window_records$date[1],
            window_records$date[nrow(window_records)]))
cat(sprintf("║  N months                    : %-4d                      ║\n",
            sum(!is.na(window_records$hl_ret_eq))))
cat("╚══════════════════════════════════════════════════════════╝\n")

## ── 5b. Build cumulative log return series ────────────────────────────────────
# Matches Figure 2: cumulative log returns of the Long (decile 1, highest pred)
# and Short (decile 10, lowest pred) legs plotted separately.
# portfolio() returns deciles sorted high→low, so:
#   col 1  = highest predicted return  = Long
#   col 10 = lowest  predicted return  = Short

wr <- window_records %>%
  mutate(date = as.Date(paste0(date, "-01"))) %>%
  arrange(date)

# Long leg = decile 1 (col 1),  Short leg = decile 10 (col 10)
# Use eval_idx rows of pf.eq / pf.vw (drop first WIDTH burn-in rows)
long_eq  <- pf.eq[eval_idx, 1]
short_eq <- pf.eq[eval_idx, 10]
long_vw  <- pf.vw[eval_idx, 1]
short_vw <- pf.vw[eval_idx, 10]

# Cumulative log return: sum of log(1 + r_t)
cumlog <- function(x) cumsum(log(1 + ifelse(is.na(x), 0, x)))

cum_long_eq  <- cumlog(long_eq)
cum_short_eq <- cumlog(short_eq)
cum_long_vw  <- cumlog(long_vw)
cum_short_vw <- cumlog(short_vw)

## ── 5c. Two-panel plot matching Figure 2 format ───────────────────────────────
make_panel <- function(cum_long, cum_short, dates, title_panel, sr_val) {
  df <- data.frame(
    date  = dates,
    Long  = cum_long,
    Short = cum_short
  ) %>%
    pivot_longer(-date, names_to = "Position", values_to = "cum_log_ret")

  ggplot(df, aes(date, cum_log_ret, colour = Position, linetype = Position)) +
    geom_line(linewidth = 0.7) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
    scale_colour_manual(
      values   = c("Long" = "red", "Short" = "blue"),
      labels   = c("Long" = "Homogeneous Long", "Short" = "Homogeneous Short")
    ) +
    scale_linetype_manual(
      values   = c("Long" = "solid", "Short" = "dotted"),
      labels   = c("Long" = "Long", "Short" = "Short")
    ) +
    scale_x_date(date_breaks = "10 years", date_labels = "%Y") +
    annotate("text",
             x = max(dates) - 365 * 2,
             y = max(cum_long, na.rm = TRUE) * 0.95,
             label = sprintf("Sharpe = %.3f", sr_val),
             hjust = 1, size = 3.2, colour = "black") +
    labs(
      title    = title_panel,
      x        = "date",
      y        = "Cumulative log return",
      colour   = NULL,
      linetype = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position   = c(0.18, 0.82),
      legend.background = element_rect(fill = "white", colour = "grey70", linewidth = 0.3),
      legend.key.width  = unit(1.5, "cm"),
      legend.text       = element_text(size = 8),
      plot.title        = element_text(face = "bold", size = 10),
      panel.grid.minor  = element_blank(),
      axis.title.y      = element_text(size = 8)
    )
}

pA <- make_panel(cum_long_eq, cum_short_eq, wr$date,
                 "Panel A: Equal-Weighted",   sr_HL_eq)
pB <- make_panel(cum_long_vw, cum_short_vw, wr$date,
                 "Panel B: Value-Weighted",   sr_HL_vw)

# Combine side-by-side with patchwork or gridExtra
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p_combined <- pA + pB +
    plot_annotation(
      title   = "Cumulative Log Portfolio Returns — Homogeneous Linear Model",
      caption = sprintf(
        "Rolling %d-month OLS beta. Decile portfolios sorted by predicted return.\nEval: %s to %s  |  EW Sharpe = %.3f  |  VW Sharpe = %.3f",
        WIDTH,
        format(min(wr$date), "%Y-%m"), format(max(wr$date), "%Y-%m"),
        sr_HL_eq, sr_HL_vw),
      theme = theme(
        plot.title   = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.caption = element_text(size = 8, hjust = 0)
      )
    )
  ggsave("plots/Homogeneous_FM_cumulative.png", p_combined,
         width = 12, height = 5, dpi = 150)
} else {
  # Fallback: save panels individually
  ggsave("plots/Homogeneous_FM_cumulative_EW.png", pA, width = 6, height = 5, dpi = 150)
  ggsave("plots/Homogeneous_FM_cumulative_VW.png", pB, width = 6, height = 5, dpi = 150)
  message("[plot] patchwork not installed — saved panels separately")
  message("[plot]   install.packages('patchwork') to get combined figure")
}
message("[plot] saved plots/Homogeneous_FM_cumulative*.png")
