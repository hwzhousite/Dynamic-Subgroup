suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})


## ── helper: winsorise at [plo, phi] per cross-section ─────────────────────────
winsorise <- function(x, plo = 0.01, phi = 0.99) {
  lo <- quantile(x, plo, na.rm = TRUE)
  hi <- quantile(x, phi, na.rm = TRUE)
  pmax(pmin(x, hi), lo)
}

## ── helper: cross-sectional rank-standardise → (-0.5, 0.5) ───────────────────
rank_std <- function(x) {
  r <- rank(x, na.last = "keep", ties.method = "average")
  n <- sum(!is.na(x))
  if (n == 0) return(x)
  (r - (n + 1) / 2) / n
}

## ── MAIN: load or simulate ────────────────────────────────────────────────────


load_crsp <- function(start_year = 1964L,
                      end_year   = 2022L) {
  df <- data1964_2022

  # Validate required columns
  required <- c("permno", "year", "month", "ret", CHARS)
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0)
    stop("[data] missing columns: ", paste(missing_cols, collapse = ", "))

  # Coerce types
  df$permno <- as.integer(df$permno)
  df$year   <- as.integer(df$year)
  df$month  <- as.integer(df$month)
  df$ret    <- as.numeric(df$ret)
  for (ch in CHARS) df[[ch]] <- as.numeric(df[[ch]])

  # Build proper Date column from year + month
  df$date <- as.Date(sprintf("%04d-%02d-01", df$year, df$month))

  df <- df %>%
    filter(!is.na(ret), !is.na(permno)) %>%
    filter(year >= start_year, year <= end_year)

  if (nrow(df) == 0)
    stop(sprintf("[data] no rows after filtering to %d-%d", start_year, end_year))

  message(sprintf("[data] loaded %d rows | %d firms | %d-%d | date range %s to %s",
                  nrow(df), n_distinct(df$permno),
                  start_year, end_year,
                  min(df$date), max(df$date)))
  return(df)
}

## ── preprocessing ─────────────────────────────────────────────────────────────

preprocess_crsp <- function(df) {
  message("[data] preprocessing: log(MV), winsorise, rank-std per cross-section ...")

  # MV is raw market cap -> take log before standardising
  if ("MV" %in% names(df)) {
    df$MV <- ifelse(df$MV > 0, log(df$MV), NA_real_)
    message("[data] MV -> log(MV) applied")
  }

  df <- df %>%
    group_by(date) %>%
    mutate(
      # 1. Winsorise return and each characteristic
      ret  = winsorise(ret),
      across(all_of(CHARS), winsorise),
      # 2. Cross-sectional rank-standardise
      across(all_of(CHARS), rank_std),
      # 3. Count missing covariates per row
      n_missing = rowSums(is.na(across(all_of(CHARS))))
    ) %>%
    # Keep rows with at most 3 missing covariates; zero-fill the rest
    filter(n_missing <= 3) %>%
    mutate(across(all_of(CHARS), ~ ifelse(is.na(.), 0, .))) %>%
    select(-n_missing) %>%
    ungroup()

  message(sprintf("[data] after cleaning: %d stock-month obs | %d firms | %d months",
                  nrow(df), n_distinct(df$permno), n_distinct(df$date)))

  # Per-characteristic coverage report
  cov_tbl <- df %>%
    summarise(across(all_of(CHARS),
                     ~ sprintf("%.1f%%", 100 * mean(!is.na(.))))) %>%
    pivot_longer(everything(), names_to = "char", values_to = "coverage")
  message("[data] covariate coverage (post-clean):")
  for (i in seq_len(nrow(cov_tbl)))
    message(sprintf("         %-20s %s", cov_tbl$char[i], cov_tbl$coverage[i]))

  return(df)
}

## ── reshape to wide arrays ────────────────────────────────────────────────────
# Returns list:
#   $Y       : n x T  matrix          (returns; NA = inactive)
#   $X       : n x T x p array        (characteristics; NA = inactive)
#   $permnos : length-n integer
#   $dates   : length-T Date (sorted)
#   $chars   : length-p character (= CHARS)
#   $years   : length-T integer
#   $months  : length-T integer

panel_to_arrays <- function(df) {
  message("[data] reshaping to wide arrays ...")
  permnos <- sort(unique(df$permno))
  dates   <- sort(unique(df$date))
  n <- length(permnos); TT <- length(dates); p <- length(CHARS)

  Y <- matrix(NA_real_, n, TT,
              dimnames = list(as.character(permnos),
                              format(dates, "%Y-%m")))
  X <- array(NA_real_, dim = c(n, TT, p),
             dimnames = list(as.character(permnos),
                             format(dates, "%Y-%m"),
                             CHARS))

  pm <- match(df$permno, permnos)
  dm <- match(df$date,   dates)
  Y[cbind(pm, dm)] <- df$ret
  for (j in seq_along(CHARS))
    X[cbind(pm, dm, j)] <- df[[CHARS[j]]]

  list(
    Y       = Y,
    X       = X,
    permnos = permnos,
    dates   = dates,
    chars   = CHARS,
    years   = as.integer(format(dates, "%Y")),
    months  = as.integer(format(dates, "%m"))
  )
}

## ── data quality summary ──────────────────────────────────────────────────────
summarise_panel <- function(wide) {
  Y <- wide$Y
  n_obs    <- sum(!is.na(Y))
  n_firms  <- nrow(Y)
  n_months <- ncol(Y)
  avg_act  <- mean(colSums(!is.na(Y)))   # avg active firms per month
  miss_pct <- 100 * mean(is.na(Y))

  message(sprintf("[data] === Panel Summary ==="))
  message(sprintf("[data]   Firms        : %d", n_firms))
  message(sprintf("[data]   Months       : %d  (%s to %s)",
                  n_months,
                  colnames(Y)[1], colnames(Y)[n_months]))
  message(sprintf("[data]   Obs (non-NA) : %d", n_obs))
  message(sprintf("[data]   Avg firms/mo : %.0f", avg_act))
  message(sprintf("[data]   Missingness  : %.1f%%", miss_pct))
  message(sprintf("[data]   Covariates   : %d  (%s)",
                  length(wide$chars),
                  paste(wide$chars, collapse = ", ")))
}
##################################################################################
# Main
##################################################################################
#OUTPUT_DIR = "/Users/haowen/Dropbox/Dynamic\ Subgroup/DSHR/weekly/Dec28\ Empirical/"
## ── covariate names (p = 15) ──────────────────────────────────────────────────
CHARS <- c("ACC", "ROA", "LogAG", "MV", "LagRet",
           "LogIssues", "DY", "LogRet",
           "LogIssues_short", "ATO", "BM", "DP", "SP", "SD", "betamkt")


## ── run ───────────────────────────────────────────────────────────────────────
load("/Users/haowen/Dropbox/Dynamic\ Subgroup/DSHR/weekly/Dec28\ Empirical/HG1964_2022.RData")
df_raw   <- load_crsp()
df_clean <- preprocess_crsp(df_raw)
wide     <- panel_to_arrays(df_clean)
summarise_panel(wide)

saveRDS(df_clean, file.path(OUTPUT_DIR, "crsp_clean.rds"))
saveRDS(wide,     file.path(OUTPUT_DIR, "crsp_wide.rds"))
message(sprintf("[data] saved crsp_clean.rds and crsp_wide.rds  (n=%d, T=%d, p=%d)",
                nrow(wide$Y), ncol(wide$Y), dim(wide$X)[3]))
