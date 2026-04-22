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
