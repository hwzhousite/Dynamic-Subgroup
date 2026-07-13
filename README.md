# Dynamic-Subgroup

**Heterogeneous regression with time-varying latent subgroups.**

Panel data are often assumed to be either homogeneous (one coefficient vector for
everyone) or *statically* heterogeneous (each unit belongs to a fixed latent group). Both
assumptions are strong. In many applications — equity return prediction being the
motivating example here — the number of distinct behavioural regimes and the membership
of each unit both change over time, and regimes that were once distinct can converge and
merge.

This repository implements **MDSP** (merge-detecting subgroup penalty) and its empirical
counterpart **AFEM** (adaptive-fusion EM): a framework that jointly

1. assigns each unit `i` to a latent subgroup at each time `t`,
2. estimates a smooth coefficient trajectory `γ_k(t)` for each subgroup, and
3. **detects merges** — via an adaptive-lasso fusion penalty on pairwise trajectory
   differences, so groups whose coefficients converge are collapsed into one rather than
   being held apart by a fixed `K`.

The repo contains a simulation study with baselines, and an out-of-sample long–short
portfolio evaluation on CRSP (1964–2022).

---

## Model

For unit `i = 1..n` at time `t = 1..m`, with covariates `X[i,t,·] ∈ R^p` and optional
homogeneous controls `Z[i,t,·] ∈ R^q`:

```
y_it = Z[i,t,·]ᵀ α(t) + X[i,t,·]ᵀ β_it,     β_it = γ_{k(i,t)}(t) + ε_it
```

- `k(i,t) ∈ {1,…,K}` is the unobserved, time-varying subgroup label.
- `γ_k(·)` is a smooth trajectory, modelled with B-splines / penalised splines.
- A fusion penalty on `‖γ_k(t) − γ_l(t)‖` with **adaptive weights** drives exact merges:
  where two trajectories coincide, their difference is shrunk to exactly zero, so the
  effective number of groups is data-determined and may vary over time.

Estimation alternates between a classification step and a trajectory-fitting step
(ADMM / proximal updates in the nonparametric version).

---

## Repository layout

```
MDSP package/            Method + simulation study
  1_data_simulation.R      Data-generating processes (incl. market-realistic noise)
  2_baseline.R             Competing methods, common interface
  3_MDSP(nonpara).R        Main estimator — splines + adaptive-lasso fusion (ADMM)
  3_MDSP(para).R           Parametric variant (basis-coefficient EM)
  utility.R                Group accuracy, ROC, gamma-RMSE helpers
  case1.Rmd / case1.html   Worked example: recovering an engineered merge point
  tra_lagret.rds           Real trajectory used to seed the case study

Empirical/               CRSP application
  1 data load.R            Load, winsorise, cross-sectionally rank-standardise
  2 method.R               Homogeneous / Static-Hetero / AFEM / SplineKM
  3 Rolling Window.R       Rolling out-of-sample refit + long–short portfolios
  4 metrics.R              Sharpe, t-stat, drawdown, avg K → CSV + LaTeX table
  5 visualize.R            Cumulative returns, rolling Sharpe, detected K, γ̂ paths
  test_homogeneous.R       Sanity checks under the homogeneous null
```

---

## Simulation study

### Data generation (`1_data_simulation.R`)

Trajectories are a smooth **backbone** plus optional market-realistic fluctuations:

- a **common-factor** AR(1) shock shared by all groups (mimics market beta), and
- **per-group idiosyncratic** AR(1) shocks.

A tapered variant (`.add_market_fluct_tapered`) zeroes the idiosyncratic noise wherever
backbones coincide, so groups that are truly merged are *exactly* identical in the
realised data — which is what makes merge recovery a well-posed target. Setting
`fluct = NULL` recovers the smooth, noise-free curves exactly (backward compatible).

### Main estimator

```r
source("3_MDSP(nonpara).R")

fit <- MDSP_nonpara(
  Y, X, Z = NULL,
  time     = time_points,
  n = n, m = m, K = 2,
  l1       = 0.5,          # fusion penalty
  kp_vec   = rep(0.1, m),
  beta.int = int           # warm start, n × m × p
)

fit$gamma   # K × m × p estimated subgroup trajectories
```

`MDSP_nonpara` is sensitive to initialisation. `case1.Rmd` demonstrates the recommended
warm start: draw many random binary partitions, keep the one that **maximises** the
separation `Σ_t ‖γ₁(t) − γ₂(t)‖²` between the two implied group fits, then classify each
unit by minimum absolute residual.

### Baselines (`2_baseline.R`)

All share the interface `fit_*(Y, X, Z, K, df, ...)` and are made comparable by
residualising out the homogeneous `Z` component before fitting.

| Function | Idea |
|---|---|
| `fit_flexmix` | Mixture of regressions on spline-expanded covariates |
| `fit_indiv_spline` | Per-unit spline fit → K-means on coefficients |
| `fit_regional_spline` | Region-wise spline fits |
| `fit_kernel_pilot` | Kernel-smoothed pilot estimate |
| `fit_flex_timewise` | Independent mixture fit at each time point |
| `fit_rolling_window` | Rolling-window refit |

`utility.R` supplies the evaluation metrics: group-recovery accuracy (with label
alignment), ROC curves, and RMSE of `γ̂` against the truth.

---

## Empirical study (CRSP, 1964–2022)

### Data

Monthly CRSP returns with **p = 15** firm characteristics:

```
ACC, ROA, LogAG, MV, LagRet, LogIssues, DY, LogRet,
LogIssues_short, ATO, BM, DP, SP, SD, betamkt
```

Preprocessing, per cross-section: winsorise at [1%, 99%], then rank-standardise to
`(−0.5, 0.5)`. Reshaped into arrays `Y (n × T)` and `X (n × T × p)`.

> **Data access.** CRSP is licensed, so the raw panel is not distributed here.
> `1 data load.R` and `3 Rolling Window.R` currently reference absolute paths on the
> author's machine (`/Users/haowen/Dropbox/...`) and expect `HG1964_2022.RData` and
> `crsp_wide.rds`. Point these at your own copy before running.

### Methods compared (`2 method.R`)

| Method | Description |
|---|---|
| **Homogeneous** | Pooled OLS; one time-invariant coefficient vector |
| **Static-Hetero** | Per-firm OLS + K-means on coefficients (static-subgroup baseline) |
| **AFEM** | *Proposed.* Ridge pilot → spline group trajectories → adaptive fusion; `K` chosen up to `K_max` and allowed to vary over time |
| **SplineKM** | Per-firm spline fit + K-means (dynamic-but-unfused benchmark) |

AFEM defaults: `K_max = 5`, `df_spline = 6`, `lam1 = 0.1` (fusion), `lam2 = 0.05`
(smoothness), `n_agg = 12`, `max_iter = 30`, `tol = 1e-4`.

### Evaluation protocol (`3 Rolling Window.R`)

For each out-of-sample month `t`:

1. Train on the window `[t − 120, t − 1]` (`WIN_SIZE = 120` months).
2. Refit every method independently on that window.
3. Predict `E[r_it] = X[i,t,·]ᵀ β̂_i`.
4. Form a long–short portfolio: long the top decile, short the bottom decile.
5. Record the realised return and diagnostics.

Filters: a firm needs ≥ 12 active months in the window; a month needs ≥ 20 active firms.

### Metrics (`4 metrics.R`)

Annualised return, annualised Sharpe, t-statistic, win rate, maximum drawdown, and — for
AFEM — the mean number of detected groups. Written to `results/perf_table.csv` and
`results/perf_table.tex` (booktabs).

---

## Getting started

### Dependencies

```r
install.packages(c(
  "splines", "Matrix", "abind",                    # core
  "mixtools", "flexmix", "genlasso", "mgcv",       # estimation
  "locfit", "cluster", "igraph", "clue",
  "dplyr", "tidyr", "ggplot2"                      # data + plots
))
```

### Reproduce the simulation

```r
setwd("MDSP package")
rmarkdown::render("case1.Rmd")
```

This builds two trajectories with a known merge point, simulates a panel from them, fits
`MDSP_nonpara`, and plots true vs. fitted trajectories side by side.

### Reproduce the empirical study

```r
setwd("Empirical")
source("1 data load.R")        # → crsp_wide.rds        (requires CRSP access)
source("3 Rolling Window.R")   # → results/rolling_results.rds
source("4 metrics.R")          # → results/perf_table.{csv,tex}
source("5 visualize.R")        # → figures
```

---

## Status

Research code, actively developed. The simulation study and the empirical pipeline are
both functional; the CRSP file paths need to be parameterised before the empirical
scripts will run on another machine.

## Citation

If you use this code, please cite the accompanying manuscript (citation to be added on
publication).
