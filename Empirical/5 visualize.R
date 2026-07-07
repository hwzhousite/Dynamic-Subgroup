##############################################################################
#  05_visualise.R
#  Empirical analysis plots:
#    (A) Cumulative H-L portfolio return by method
#    (B) Monthly H-L return bar chart with rolling Sharpe
#    (C) AFEM detected K over time
#    (D) Representative coefficient trajectories (Gamma_hat)
#         — requires one saved fit from rolling window
##############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

dir.create("plots", showWarnings = FALSE)

METHOD_COLORS <- c(
  "Homogeneous FM" = "#999999",
  "Static Hetero"  = "#E69F00",
  "SplineKM"       = "#56B4E9",
  "AFEM"           = "#D55E00"
)
METHOD_ORDER <- c("Homogeneous FM", "Static Hetero", "SplineKM", "AFEM")

ret_df <- tryCatch(readRDS("results/rolling_results.rds"),
                   error = function(e) stop("[plot] run 03_rolling_window.R first"))

ret_df <- ret_df %>%
  mutate(method = factor(method, levels = METHOD_ORDER),
         date   = as.Date(paste0(date, "-01")))

## ── (A) Cumulative return ────────────────────────────────────────────────────
cum_df <- ret_df %>%
  arrange(method, date) %>%
  group_by(method) %>%
  mutate(cum_ret = cumprod(1 + hl_ret) - 1) %>%
  ungroup()

pA <- ggplot(cum_df, aes(date, cum_ret * 100, colour = method, linetype = method)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  scale_colour_manual(values = METHOD_COLORS, name = "Method") +
  scale_linetype_manual(values = c("dashed","dotdash","dotted","solid"), name = "Method") +
  labs(title    = "Cumulative H-L Portfolio Return",
       subtitle = "Out-of-sample rolling window (equal-weighted)",
       x = "Date", y = "Cumulative Return (%)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave("plots/cumulative_return.png", pA, width = 8, height = 4.5, dpi = 150)
message("[plot] saved plots/cumulative_return.png")

## ── (B) Monthly bar + rolling Sharpe ────────────────────────────────────────
roll_sharpe <- function(x, k = 12) {
  n <- length(x)
  sapply(seq_len(n), function(i) {
    idx <- max(1, i - k + 1):i
    if (length(idx) < 4) return(NA_real_)
    (mean(x[idx]) / sd(x[idx])) * sqrt(12)
  })
}

sharpe_df <- ret_df %>%
  arrange(method, date) %>%
  group_by(method) %>%
  mutate(rolling_sharpe = roll_sharpe(hl_ret)) %>%
  ungroup()

pB <- ggplot(sharpe_df %>% filter(method == "AFEM"),
             aes(date, rolling_sharpe)) +
  geom_line(colour = METHOD_COLORS["AFEM"], linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1.95, linetype = "dotted", colour = "steelblue",
             linewidth = 0.8) +
  annotate("text", x = min(sharpe_df$date, na.rm = TRUE), y = 2.05,
           label = "Target (1.95)", hjust = 0, size = 3.2, colour = "steelblue") +
  labs(title = "AFEM: 12-month Rolling Annualised Sharpe Ratio",
       x = "Date", y = "Rolling Sharpe") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave("plots/rolling_sharpe_afem.png", pB, width = 8, height = 3.5, dpi = 150)
message("[plot] saved plots/rolling_sharpe_afem.png")

## ── (C) K over time (AFEM) ──────────────────────────────────────────────────
k_df <- ret_df %>%
  filter(method == "AFEM", !is.na(K_used)) %>%
  arrange(date)

if (nrow(k_df) > 0) {
  pC <- ggplot(k_df, aes(date, K_used)) +
    geom_step(colour = METHOD_COLORS["AFEM"], linewidth = 1) +
    geom_point(colour = METHOD_COLORS["AFEM"], size = 1.5) +
    scale_y_continuous(breaks = 1:10) +
    labs(title    = "AFEM: Number of Detected Subgroups over Time",
         subtitle = "Dynamic group structure detected per rolling window",
         x = "Prediction month", y = expression(hat(K))) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  ggsave("plots/K_over_time.png", pC, width = 8, height = 3.5, dpi = 150)
  message("[plot] saved plots/K_over_time.png")
}

## ── (D) All-methods monthly H-L return comparison ───────────────────────────
pD <- ggplot(ret_df, aes(date, hl_ret * 100, fill = method)) +
  geom_col(position = "dodge", width = 20) +
  scale_fill_manual(values = METHOD_COLORS) +
  facet_wrap(~method, ncol = 2) +
  labs(title = "Monthly H-L Portfolio Return by Method (%)",
       x = "Date", y = "H-L Return (%)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

ggsave("plots/monthly_hl_return.png", pD, width = 10, height = 6, dpi = 150)
message("[plot] saved plots/monthly_hl_return.png")

## ── (E) Sharpe ratio comparison bar ─────────────────────────────────────────
perf <- tryCatch(read.csv("results/perf_table.csv", stringsAsFactors = FALSE),
                 error = function(e) NULL)

if (!is.null(perf)) {
  perf$method <- factor(perf$method, levels = METHOD_ORDER)
  pE <- ggplot(perf, aes(method, sharpe_ann, fill = method)) +
    geom_col(width = 0.6) +
    geom_hline(yintercept = 1.95, linetype = "dashed", colour = "grey30") +
    annotate("text", x = 0.6, y = 2.05, label = "Target (1.95)",
             hjust = 0, size = 3.5, colour = "grey30") +
    scale_fill_manual(values = METHOD_COLORS) +
    labs(title = "Annualised Sharpe Ratio — Method Comparison",
         x = NULL, y = "Annualised Sharpe Ratio") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))

  ggsave("plots/sharpe_comparison.png", pE, width = 6, height = 4, dpi = 150)
  message("[plot] saved plots/sharpe_comparison.png")
}

message("[plot] all plots done")
