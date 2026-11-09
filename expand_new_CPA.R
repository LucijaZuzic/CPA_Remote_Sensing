# Clear environment
rm(list = ls())

## Expanded Statistical & Physical Analysis for Žužić et al. manuscript
## Builds directly on the outputs of make_time_series_Darwin_error.R, new_CPA.R, and print_tables_RS.R

library(tidyverse)
library(data.table)
library(lubridate)
library(readr)
library(stringr)
library(fitdistrplus)
library(gridExtra)
library(extrafont)

# Avoid namespace conflicts
filter <- dplyr::filter
select <- dplyr::select

getCurrentFileLocation <- function() {
  this_file <- commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col = value, into = c("key", "value"), sep = "=", fill = "right") %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file) == 0) {
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      this_file <- rstudioapi::getSourceEditorContext()$path
    } else {
      return(getwd())
    }
  }
  return(dirname(this_file))
}

setwd(getCurrentFileLocation())
loadfonts(device = "pdf", quiet = TRUE)

# Create output directory for expanded analysis
exp_dir <- "CPD_Darwin_2014_Expanded_Analysis"
dir.create(exp_dir, showWarnings = FALSE)

sink(file.path(exp_dir, "Expanded_CPA_Analysis_Log.txt"), split = TRUE)

# Shared publication theme matching make_time_series_Darwin_error.R
bigger_theme <- theme_classic(base_size = 18) +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 18, color = "black", face = "bold"),
    axis.title   = element_text(size = 18, color = "black"),
    axis.text    = element_text(size = 16, color = "black"),
    axis.ticks   = element_line(color = "black", linewidth = 0.8),
    axis.line    = element_line(color = "black", linewidth = 0.8),
    legend.text  = element_text(size = 15, color = "black"),
    legend.title = element_text(size = 16, color = "black", face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    strip.text   = element_text(size = 16, color = "black", face = "bold")
  )

# ==============================================================================
# Helper Functions (from print_tables_RS.R and new_CPA.R)
# ==============================================================================
format_latex_scientific <- function(val, threshold = 2.2e-16) {
  val <- as.numeric(val)
  if (is.na(val)) return("NA")
  if (val < threshold) {
    return("$< 2.2 \\times 10^{-16}$")
  } else if (val < 0.001) {
    sci <- sprintf("%.2e", val)
    parts <- unlist(strsplit(sci, "e"))
    base <- parts[1]
    exp <- as.integer(parts[2])
    return(paste0("$", base, " \\times 10^{", exp, "}$"))
  } else {
    formatted_num <- formatC(val, format = "f", digits = 4, big.mark = ",")
    return(paste0("$", formatted_num, "$"))
  }
}

df_to_custom_tabularx <- function(df, caption = "Table Caption", label = "tab:my_label",
                                  digits = 4, col_formats = NULL, file_out = NULL) {
  nc <- ncol(df)
  col_spec <- paste0(rep("C", nc), collapse = "")
  headers <- paste0("\\textbf{", str_replace_all(colnames(df), "_", " "), "}", collapse = " & ")
  
  rows <- apply(df, 1, function(row) {
    formatted_cells <- sapply(seq_along(row), function(i) {
      col_name <- colnames(df)[i]
      val <- row[i]
      
      if (!is.null(col_formats) && col_name %in% names(col_formats)) {
        rule <- col_formats[[col_name]]
        if (is.function(rule)) {
          return(rule(val))
        } else if (is.character(rule)) {
          if (rule == "string") {
            return(str_replace_all(val, "_", " "))
          } else if (rule == "integer") {
            formatted_num <- format(as.integer(val), big.mark = ",", scientific = FALSE)
            return(paste0("$", formatted_num, "$"))
          } else {
            digits_val <- suppressWarnings(as.integer(sub("%.([0-9]+)f.*", "\\1", rule)))
            if (is.na(digits_val)) digits_val <- digits
            formatted_num <- formatC(as.numeric(val), format = "f", digits = digits_val, big.mark = ",")
            return(paste0("$", formatted_num, "$"))
          }
        }
      }
      
      if (!is.na(suppressWarnings(as.numeric(val)))) {
        num_val <- as.numeric(val)
        if (num_val == round(num_val) && abs(num_val) < 1e6) {
          formatted_num <- format(as.integer(num_val), big.mark = ",")
          paste0("$", formatted_num, "$")
        } else {
          formatted_num <- formatC(num_val, format = "f", digits = digits, big.mark = ",")
          paste0("$", formatted_num, "$")
        }
      } else {
        as.character(val)
      }
    })
    paste0(paste0(formatted_cells, collapse = " & "), " \\\\")
  })
  row_string <- paste0(rows, collapse = "\n")
  
  latex_output <- sprintf(
    "\\vspace{-3pt}
\\begin{table}[H]
\\caption{%s\\label{%s}}
\\begin{tabularx}{\\textwidth}{%s}
\\toprule
%s \\\\
\\midrule
%s
\\bottomrule
\\end{tabularx}
\\end{table}
\\vspace{-3pt}",
    caption, label, col_spec, headers, row_string
  )
  
  if (!is.null(file_out)) {
    writeLines(latex_output, file_out)
    cat(sprintf("LaTeX table written to: %s\n", file_out))
  }
  return(latex_output)
}

link_cause_effect_cps <- function(cause_cps_sec, effect_cps_sec, max_lag_sec) {
  cause_cps_sec  <- sort(na.omit(unique(as.numeric(cause_cps_sec))))
  effect_cps_sec <- sort(na.omit(unique(as.numeric(effect_cps_sec))))
  
  if (length(cause_cps_sec) == 0 || length(effect_cps_sec) == 0) return(tibble())
  
  candidates <- expand_grid(
    cause_idx  = seq_along(cause_cps_sec),
    effect_idx = seq_along(effect_cps_sec)
  ) %>%
    mutate(
      t_cause_B_sec    = cause_cps_sec[cause_idx],
      t_effect_GPS_sec = effect_cps_sec[effect_idx],
      delay_sec        = t_effect_GPS_sec - t_cause_B_sec
    ) %>%
    filter(delay_sec >= 0 & delay_sec <= max_lag_sec) %>%
    arrange(delay_sec, t_cause_B_sec)
  
  if (nrow(candidates) == 0) return(tibble())
  
  used_cause   <- rep(FALSE, length(cause_cps_sec))
  used_effect  <- rep(FALSE, length(effect_cps_sec))
  matched_rows <- list()
  
  for (i in seq_len(nrow(candidates))) {
    c_i <- candidates$cause_idx[i]
    e_i <- candidates$effect_idx[i]
    if (!used_cause[c_i] && !used_effect[e_i]) {
      used_cause[c_i]  <- TRUE
      used_effect[e_i] <- TRUE
      matched_rows[[length(matched_rows) + 1]] <- candidates[i, ]
    }
  }
  
  bind_rows(matched_rows) %>%
    arrange(t_cause_B_sec) %>%
    mutate(
      pair_id     = row_number(),
      DOY_approx  = floor(t_cause_B_sec / 86400) + 1,
      delay_min   = delay_sec / 60,
      delay_hours = delay_sec / 3600
    )
}

# ==============================================================================
# 1. Expanded Descriptive Statistics & Correlation Analysis (including TEC & Dst)
# ==============================================================================
cat("\n==============================================================================\n")
cat("1. Expanded Descriptive Statistics & Correlation Matrix (GPS, Bx, By, Bz, TEC, Dst)\n")
cat("==============================================================================\n")

epoch_bin_sec <- 300
cause_cols <- c("Bx", "By", "Bz", "TEC", "Dst")

kp_csv_path <- file.path(paste0("CPD_Darwin_2014_Bx_", epoch_bin_sec), "kp_file.csv")
merged_csv_path <- file.path(paste0("CPD_Darwin_2014_Bx_", epoch_bin_sec), "Darwin_2014_Merged_Unfiltered_GPS_Bxyz.csv")

if (file.exists(kp_csv_path)) {
  full_merged_df <- read_csv(kp_csv_path, show_col_types = FALSE)
} else if (file.exists(merged_csv_path)) {
  full_merged_df <- read_csv(merged_csv_path, show_col_types = FALSE)
} else {
  stop("Merged dataset not found. Please run new_CPA.R first.")
}

vars_all <- intersect(c("horizontal_error", "Bx", "By", "Bz", "TEC", "Dst", "Kp"), colnames(full_merged_df))

desc_stats_all <- map_dfr(vars_all, function(v) {
  x <- na.omit(as.numeric(full_merged_df[[v]]))
  n <- length(x)
  m <- mean(x)
  s <- sd(x)
  skew <- sum((x - m)^3 / s^3) / n
  kurt <- sum((x - m)^4 / s^4) / n - 3
  
  tibble(
    Variable = v,
    N        = n,
    Min      = min(x),
    Q1       = quantile(x, 0.25),
    Median   = median(x),
    Mean     = m,
    Q3       = quantile(x, 0.75),
    Max      = max(x),
    Variance = var(x),
    SD       = s,
    IQR      = IQR(x),
    MAD      = mad(x),
    Skewness = skew,
    Kurtosis = kurt
  )
})

print(desc_stats_all)
write_csv(desc_stats_all, file.path(exp_dir, "Descriptive_Statistics_All_Variables.csv"))

df_to_custom_tabularx(
  df = desc_stats_all %>% select(-N),
  caption = "Descriptive statistics for the unfiltered 2014 Darwin dataset: GPS horizontal positioning error [m], geomagnetic field components ($B_x, B_y, B_z$ [nT]), ionospheric TEC [TECU], ring current index Dst [nT], and planetary $Kp$ index.",
  label = "tab:Descriptive_Statistics_All_Variables",
  col_formats = list(
    "Variable" = "string", "Min" = "%.2f", "Q1" = "%.2f", "Median" = "%.2f",
    "Mean" = "%.2f", "Q3" = "%.2f", "Max" = "%.2f", "Variance" = "%.2f",
    "SD" = "%.2f", "IQR" = "%.2f", "MAD" = "%.2f", "Skewness" = "%.2f", "Kurtosis" = "%.2f"
  ),
  file_out = file.path(exp_dir, "Descriptive_Statistics_All_Variables.tex")
)

# Pearson and Spearman Correlation Matrices
cor_pearson  <- cor(full_merged_df[, vars_all], use = "complete.obs", method = "pearson")
cor_spearman <- cor(full_merged_df[, vars_all], use = "complete.obs", method = "spearman")

cat("\n--- Pearson Correlation Matrix (All Variables) ---\n")
print(round(cor_pearson, 4))
cat("\n--- Spearman Rank Correlation Matrix (All Variables) ---\n")
print(round(cor_spearman, 4))

write_csv(as_tibble(cor_pearson, rownames = "Variable"), file.path(exp_dir, "Correlation_Matrix_Pearson.csv"))
write_csv(as_tibble(cor_spearman, rownames = "Variable"), file.path(exp_dir, "Correlation_Matrix_Spearman.csv"))

# ==============================================================================
# 2. Load All Change Points & Linked Pairs Across All 5 Cause Components
# ==============================================================================
cat("\n==============================================================================\n")
cat("2. Loading Detected Change Points & Linked Pairs for Bx, By, Bz, TEC, Dst\n")
cat("==============================================================================\n")

cp_gps_df <- read_csv(
  file.path(paste0("CPD_Darwin_2014_Bx_", epoch_bin_sec),
            paste0("ChangePoints_GPS_Horizontal_Error_", epoch_bin_sec, ".csv")),
  show_col_types = FALSE
)

cp_cause_list     <- list()
linked_pairs_list <- list()

for (comp in cause_cols) {
  comp_dir <- paste0("CPD_Darwin_2014_", comp, "_", epoch_bin_sec)
  cp_file  <- file.path(comp_dir, paste0("ChangePoints_Geomagnetic_", comp, "_", epoch_bin_sec, ".csv"))
  lnk_file <- file.path(comp_dir, paste0("Linked_CP_Pairs_", comp, "_to_GPS_", epoch_bin_sec, ".csv"))
  
  if (file.exists(cp_file)) {
    cp_cause_list[[comp]] <- read_csv(cp_file, show_col_types = FALSE) %>%
      mutate(Cause_Component = comp)
  }
  if (file.exists(lnk_file)) {
    linked_pairs_list[[comp]] <- read_csv(lnk_file, show_col_types = FALSE) %>%
      mutate(Cause_Component = comp)
  }
}

all_cause_cps_df    <- bind_rows(cp_cause_list)
all_linked_pairs_df <- bind_rows(linked_pairs_list) %>%
  mutate(Cause_Component = factor(Cause_Component, levels = cause_cols))

# ==============================================================================
# 3. Statistical Hypothesis Testing & Monte Carlo Null Model for CP Linking
# ==============================================================================
cat("\n==============================================================================\n")
cat("3. Hypothesis Testing & Monte Carlo Permutation Test for Cause-Effect CP Linking\n")
cat("==============================================================================\n")

# Global omnibus tests across the 5 cause components
kw_test    <- kruskal.test(delay_sec ~ Cause_Component, data = all_linked_pairs_df)
anova_test <- summary(aov(delay_sec ~ Cause_Component, data = all_linked_pairs_df))

cat(sprintf("Kruskal-Wallis test across 5 components: chi-squared = %.4f, df = %d, p-value = %.4g\n",
            kw_test$statistic, kw_test$parameter, kw_test$p.value))
print(anova_test)

# Monte Carlo Null Model:
# Within each active day (364 days, excluding DOY 300), 2 CPs are placed uniformly at random
# in [0, 86400) for the cause series, and matched to the observed GPS CPs using the same greedy algorithm.
set.seed(42)
n_mc <- 500
max_lag_default <- 14400
valid_doys <- unique(floor(cp_gps_df$cp_time_sec / 86400))

mc_null_results <- replicate(n_mc, {
  sim_cause_cps <- as.vector(sapply(valid_doys, function(d) d * 86400 + runif(2, 0, 86400)))
  sim_linked    <- link_cause_effect_cps(sim_cause_cps, cp_gps_df$cp_time_sec, max_lag_default)
  c(n_pairs = nrow(sim_linked), mean_delay = mean(sim_linked$delay_sec, na.rm = TRUE))
})

mc_null_df <- as_tibble(t(mc_null_results))
cat(sprintf("Monte Carlo Null Model (%d runs): Expected N_pairs = %.2f +/- %.2f, Expected Mean Delay = %.2f +/- %.2f s\n",
            n_mc, mean(mc_null_df$n_pairs), sd(mc_null_df$n_pairs),
            mean(mc_null_df$mean_delay), sd(mc_null_df$mean_delay)))

# Component-wise KS test against Uniform(0, 14400) and Monte Carlo p-values
hypothesis_table <- map_dfr(cause_cols, function(comp) {
  sub_df <- all_linked_pairs_df %>% filter(Cause_Component == comp)
  delays <- sub_df$delay_sec
  n_obs  <- length(delays)
  n_tot  <- nrow(all_cause_cps_df %>% filter(Cause_Component == comp))
  
  ks_res   <- ks.test(delays, "punif", 0, max_lag_default)
  wilc_res <- wilcox.test(delays, mu = max_lag_default / 2)
  
  # Empirical p-value for N_pairs exceeding random uniform CPs
  p_mc_n    <- (sum(mc_null_df$n_pairs >= n_obs) + 1) / (n_mc + 1)
  z_score_n <- (n_obs - mean(mc_null_df$n_pairs)) / sd(mc_null_df$n_pairs)
  
  # Parametric distribution fits (Weibull, Gamma, Normal)
  scaled_delays <- pmax(delays, 1) # ensure strictly positive for Weibull/Gamma
  fit_weib <- tryCatch(fitdist(scaled_delays, "weibull"), error = function(e) NULL)
  fit_gam  <- tryCatch(fitdist(scaled_delays, "gamma"), error = function(e) NULL)
  fit_norm <- tryCatch(fitdist(scaled_delays, "norm"), error = function(e) NULL)
  
  aics <- c(
    Weibull = if (!is.null(fit_weib)) fit_weib$aic else Inf,
    Gamma   = if (!is.null(fit_gam))  fit_gam$aic  else Inf,
    Normal  = if (!is.null(fit_norm)) fit_norm$aic else Inf
  )
  best_dist <- names(which.min(aics))
  
  tibble(
    Cause_Component = comp,
    Total_CPs       = n_tot,
    Linked_Pairs    = n_obs,
    Match_Rate_Pct  = 100 * n_obs / n_tot,
    MC_Null_Pairs   = mean(mc_null_df$n_pairs),
    Z_Score_Pairs   = z_score_n,
    p_val_MC_Pairs  = p_mc_n,
    KS_Stat_Unif    = as.numeric(ks_res$statistic),
    KS_p_val        = ks_res$p.value,
    Wilcox_p_val    = wilc_res$p.value,
    Best_Fit_Dist   = best_dist
  )
})

print(hypothesis_table)
write_csv(hypothesis_table, file.path(exp_dir, "Hypothesis_Testing_and_MC_Null_Summary.csv"))

df_to_custom_tabularx(
  df = hypothesis_table,
  caption = "Monte Carlo significance testing and goodness-of-fit analysis of linked Cause--Effect Change Point (CP) pairs ($x = 14,400$~s).",
  label = "tab:Hypothesis_Testing_and_MC_Null",
  col_formats = list(
    "Cause_Component" = "string",
    "Total_CPs"       = "integer",
    "Linked_Pairs"    = "integer",
    "Match_Rate_Pct"  = "%.2f",
    "MC_Null_Pairs"   = "%.2f",
    "Z_Score_Pairs"   = "%.2f",
    "p_val_MC_Pairs"  = format_latex_scientific,
    "KS_Stat_Unif"    = "%.4f",
    "KS_p_val"        = format_latex_scientific,
    "Wilcox_p_val"    = format_latex_scientific,
    "Best_Fit_Dist"   = "string"
  ),
  file_out = file.path(exp_dir, "Hypothesis_Testing_and_MC_Null_Summary.tex")
)

# ==============================================================================
# 4. Sensitivity Analysis to Maximum Allowed Lag Window x (1h, 2h, 3h, 4h, 6h)
# ==============================================================================
cat("\n==============================================================================\n")
cat("4. Sensitivity Analysis Across Maximum Allowed Lag Thresholds x\n")
cat("==============================================================================\n")

lag_grid_sec <- c(3600, 7200, 10800, 14400, 21600) # 1h, 2h, 3h, 4h, 6h

sensitivity_df <- expand_grid(
  Cause_Component = cause_cols,
  Max_Lag_s       = lag_grid_sec
) %>%
  pmap_dfr(function(Cause_Component, Max_Lag_s) {
    c_cps <- all_cause_cps_df %>% filter(Cause_Component == !!Cause_Component) %>% pull(cp_time_sec)
    lnk   <- link_cause_effect_cps(c_cps, cp_gps_df$cp_time_sec, Max_Lag_s)
    
    tibble(
      Cause_Component = Cause_Component,
      Max_Lag_Hours   = Max_Lag_s / 3600,
      Max_Lag_s       = Max_Lag_s,
      N_Pairs         = nrow(lnk),
      Match_Rate_Pct  = 100 * nrow(lnk) / length(c_cps),
      Mean_Delay_s    = mean(lnk$delay_sec, na.rm = TRUE),
      SD_Delay_s      = sd(lnk$delay_sec, na.rm = TRUE),
      Median_Delay_s  = median(lnk$delay_sec, na.rm = TRUE),
      IQR_Delay_s     = IQR(lnk$delay_sec, na.rm = TRUE)
    )
  })

print(sensitivity_df)
write_csv(sensitivity_df, file.path(exp_dir, "Sensitivity_Analysis_Max_Lag.csv"))

df_to_custom_tabularx(
  df = sensitivity_df,
  caption = "Sensitivity analysis of linked CP pairs and response delay statistics as a function of maximum allowed lag $x \\in \\{1, 2, 3, 4, 6\\}$ hours.",
  label = "tab:Sensitivity_Analysis_Max_Lag",
  col_formats = list(
    "Cause_Component" = "string",
    "Max_Lag_Hours"   = "integer",
    "Max_Lag_s"       = "integer",
    "N_Pairs"         = "integer",
    "Match_Rate_Pct"  = "%.2f",
    "Mean_Delay_s"    = "%.2f",
    "SD_Delay_s"      = "%.2f",
    "Median_Delay_s"  = "%.2f",
    "IQR_Delay_s"     = "%.2f"
  ),
  file_out = file.path(exp_dir, "Sensitivity_Analysis_Max_Lag.tex")
)

# ==============================================================================
# 5. Stratification by Geomagnetic Activity Regime (Quiet vs. Active/Storm Days)
# ==============================================================================
cat("\n==============================================================================\n")
cat("5. Stratification by GFZ Kp Geomagnetic Regime (Active Kp_max >= 4 vs Quiet Kp_max < 4)\n")
cat("==============================================================================\n")

if ("Kp_daily_max" %in% colnames(full_merged_df)) {
  daily_kp_lookup <- full_merged_df %>%
    group_by(DOY_approx = DOY + 1) %>%
    summarise(
      Kp_daily_max = max(Kp_daily_max, na.rm = TRUE),
      Dst_daily_min = min(Dst, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Regime = ifelse(Kp_daily_max >= 4, "Active/Storm (Kp >= 4)", "Quiet (Kp < 4)")
    )
  
  stratified_pairs_df <- all_linked_pairs_df %>%
    left_join(daily_kp_lookup, by = "DOY_approx") %>%
    filter(!is.na(Regime))
  
  regime_summary_df <- stratified_pairs_df %>%
    group_by(Cause_Component, Regime) %>%
    summarise(
      N_Pairs        = n(),
      Mean_Delay_s   = mean(delay_sec, na.rm = TRUE),
      SD_Delay_s     = sd(delay_sec, na.rm = TRUE),
      Median_Delay_s = median(delay_sec, na.rm = TRUE),
      IQR_Delay_s    = IQR(delay_sec, na.rm = TRUE),
      .groups        = "drop"
    )
  
  # Wilcoxon test between Active/Storm and Quiet regimes for each component
  regime_pvals <- map_dfr(cause_cols, function(comp) {
    sub_df <- stratified_pairs_df %>% filter(Cause_Component == comp)
    w_res  <- tryCatch(wilcox.test(delay_sec ~ Regime, data = sub_df), error = function(e) list(p.value = NA))
    tibble(Cause_Component = comp, Wilcoxon_p_val = w_res$p.value)
  })
  
  regime_wide_df <- regime_summary_df %>%
    pivot_wider(
      names_from  = Regime,
      values_from = c(N_Pairs, Mean_Delay_s, SD_Delay_s, Median_Delay_s)
    ) %>%
    left_join(regime_pvals, by = "Cause_Component")
  
  print(regime_summary_df)
  write_csv(regime_summary_df, file.path(exp_dir, "Stratified_Delay_By_Kp_Regime.csv"))
  
  df_to_custom_tabularx(
    df = regime_summary_df,
    caption = "Response delay statistics stratified by geomagnetic activity regime: Active/Storm days ($Kp_{\\max} \\ge 4$, 61 days) versus Quiet days ($Kp_{\\max} < 4$, 304 days).",
    label = "tab:Stratified_Delay_By_Kp_Regime",
    col_formats = list(
      "Cause_Component" = "string",
      "Regime"          = "string",
      "N_Pairs"         = "integer",
      "Mean_Delay_s"    = "%.2f",
      "SD_Delay_s"      = "%.2f",
      "Median_Delay_s"  = "%.2f",
      "IQR_Delay_s"     = "%.2f"
    ),
    file_out = file.path(exp_dir, "Stratified_Delay_By_Kp_Regime.tex")
  )
}

# ==============================================================================
# 6. Publication-Ready Comparative Figures
# ==============================================================================
cat("\n==============================================================================\n")
cat("6. Generating Publication-Ready Multi-Component Figures\n")
cat("==============================================================================\n")

# Figure 1: Combined Violin + Boxplot + Mean Comparison of Response Delays Across All 5 Cause Components
p_violin <- ggplot(all_linked_pairs_df, aes(x = Cause_Component, y = delay_sec / 60, fill = Cause_Component)) +
  geom_violin(alpha = 0.35, color = "black", trim = FALSE) +
  geom_boxplot(width = 0.22, outlier.shape = 16, outlier.size = 2, alpha = 0.85, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 4, color = "red") +
  scale_y_continuous(limits = c(0, 260), breaks = seq(0, 240, by = 30)) +
  labs(
    title = "Distribution of GPS Positioning Error Response Delays\nAcross Cause Components (Darwin, 2014)",
    x     = "Cause Component (Change Point Source)",
    y     = "Response Delay [minutes]"
  ) +
  bigger_theme +
  theme(legend.position = "none")

ggsave(file.path(exp_dir, "Combined_Delay_Violin_Boxplot_All_Components.pdf"),
       plot = p_violin, width = 10, height = 6, device = cairo_pdf)
ggsave(file.path(exp_dir, "Combined_Delay_Violin_Boxplot_All_Components.png"),
       plot = p_violin, width = 10, height = 6, dpi = 300)

# Figure 2: Overlaid Experimental Probability Density Functions (PDFs) Across All 5 Components
p_dens <- ggplot(all_linked_pairs_df, aes(x = delay_sec, color = Cause_Component, fill = Cause_Component)) +
  geom_density(alpha = 0.15, linewidth = 1.2, bounds = c(0, max_lag_default)) +
  geom_hline(yintercept = 1 / max_lag_default, linetype = "dashed", color = "black", linewidth = 0.9) +
  scale_x_continuous(breaks = seq(0, 14400, by = 3600),
                     labels = c("0 (0h)", "3,600 (1h)", "7,200 (2h)", "10,800 (3h)", "14,400 (4h)")) +
  labs(
    title = "Experimental PDFs of GPS Positioning Error Response Delays\nvs. Uniform Null Baseline (Dashed Line)",
    x     = "Response Delay [seconds (hours)]",
    y     = "Probability Density",
    color = "Cause",
    fill  = "Cause"
  ) +
  bigger_theme

ggsave(file.path(exp_dir, "Combined_Delay_PDFs_All_Components.pdf"),
       plot = p_dens, width = 10.5, height = 6, device = cairo_pdf)
ggsave(file.path(exp_dir, "Combined_Delay_PDFs_All_Components.png"),
       plot = p_dens, width = 10.5, height = 6, dpi = 300)

# Figure 3: Sensitivity Curves (Matching Rate % and Mean Delay vs Max Allowed Lag x)
p_sens1 <- ggplot(sensitivity_df, aes(x = Max_Lag_Hours, y = Match_Rate_Pct, color = Cause_Component, group = Cause_Component)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 6)) +
  labs(
    title = "Sensitivity of CP Linking Rate to Max Allowed Lag (x)",
    x     = "Maximum Allowed Lag x [hours]",
    y     = "Linked CP Pairs [% of Total CPs]",
    color = "Cause"
  ) +
  bigger_theme

p_sens2 <- ggplot(sensitivity_df, aes(x = Max_Lag_Hours, y = Mean_Delay_s / 60, color = Cause_Component, group = Cause_Component)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 6)) +
  labs(
    title = "Sensitivity of Mean Response Delay to Max Allowed Lag (x)",
    x     = "Maximum Allowed Lag x [hours]",
    y     = "Arithmetic Mean Delay [minutes]",
    color = "Cause"
  ) +
  bigger_theme

p_sens_combined <- arrangeGrob(p_sens1, p_sens2, ncol = 2)
ggsave(file.path(exp_dir, "Sensitivity_Analysis_Max_Lag_Plot.pdf"),
       plot = p_sens_combined, width = 15, height = 6, device = cairo_pdf)
ggsave(file.path(exp_dir, "Sensitivity_Analysis_Max_Lag_Plot.png"),
       plot = p_sens_combined, width = 15, height = 6, dpi = 300)

# Figure 4: Stratified Boxplot (Active/Storm vs Quiet Geomagnetic Conditions)
if (exists("stratified_pairs_df")) {
  p_strat <- ggplot(stratified_pairs_df, aes(x = Cause_Component, y = delay_sec / 60, fill = Regime)) +
    geom_boxplot(position = position_dodge(0.8), width = 0.65, alpha = 0.85, color = "black") +
    labs(
      title = "GPS Positioning Error Response Delays by Geomagnetic Activity Regime",
      x     = "Cause Component",
      y     = "Response Delay [minutes]",
      fill  = "GFZ Kp Regime"
    ) +
    bigger_theme +
    theme(legend.position = "top")
  
  ggsave(file.path(exp_dir, "Stratified_Delay_Kp_Regime_Boxplot.pdf"),
         plot = p_strat, width = 11, height = 6.5, device = cairo_pdf)
  ggsave(file.path(exp_dir, "Stratified_Delay_Kp_Regime_Boxplot.png"),
         plot = p_strat, width = 11, height = 6.5, dpi = 300)
}

# Figure 5: Diurnal (Time-of-Day UT) Distribution of Cause vs Effect Change Points
all_cps_diurnal <- bind_rows(
  cp_gps_df %>% mutate(Component = "GPS Error") %>% select(Component, cp_time_sec),
  all_cause_cps_df %>% mutate(Component = Cause_Component) %>% select(Component, cp_time_sec)
) %>%
  mutate(
    UT_Hour   = (cp_time_sec %% 86400) / 3600,
    Component = factor(Component, levels = c("GPS Error", cause_cols))
  )

p_diurnal <- ggplot(all_cps_diurnal, aes(x = UT_Hour, color = Component, fill = Component)) +
  geom_density(alpha = 0.12, linewidth = 1.1, bounds = c(0, 24)) +
  scale_x_continuous(breaks = seq(0, 24, by = 4), limits = c(0, 24)) +
  labs(
    title = "Diurnal (UT Hour) Density of Detected Change Points (Darwin, 2014)",
    x     = "Universal Time [Hours of Day]",
    y     = "Density",
    color = "Series",
    fill  = "Series"
  ) +
  bigger_theme

ggsave(file.path(exp_dir, "Diurnal_Distribution_ChangePoints_UT.pdf"),
       plot = p_diurnal, width = 10.5, height = 6, device = cairo_pdf)
ggsave(file.path(exp_dir, "Diurnal_Distribution_ChangePoints_UT.png"),
       plot = p_diurnal, width = 10.5, height = 6, dpi = 300)

cat("\nExpanded CPA analysis completed successfully. All tables and figures saved in:", exp_dir, "\n")
sink()