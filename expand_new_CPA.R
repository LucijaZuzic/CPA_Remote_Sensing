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
library(scales)
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

# Optional: if Definitions/ folder exists for LaTeX compilation, update plots there too
def_dir <- "Definitions"
has_def_dir <- dir.exists(def_dir)

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
# Helper Functions (Math Italic Labels, LaTeX Tables, and CP Linking)
# ==============================================================================

# Returns R plotmath expressions with LaTeX-style math italic for variables
math_label_vec <- c(
  "GPS Error"        = "italic(e)[italic(H)]",
  "horizontal_error" = "italic(e)[italic(H)]",
  "Bx"               = "italic(B)[italic(x)]",
  "By"               = "italic(B)[italic(y)]",
  "Bz"               = "italic(B)[italic(z)]",
  "TEC"              = "TEC",
  "Dst"              = "Dst",
  "Kp"               = "italic(Kp)"
)

parse_math_labels <- function(x) {
  mapped <- ifelse(x %in% names(math_label_vec), math_label_vec[as.character(x)], as.character(x))
  parse(text = mapped)
}

get_comp_expr <- function(comp) {
  switch(as.character(comp),
         "Bx"  = expression(italic(B)[italic(x)]),
         "By"  = expression(italic(B)[italic(y)]),
         "Bz"  = expression(italic(B)[italic(z)]),
         "TEC" = expression(TEC),
         "Dst" = expression(Dst),
         expression(italic(B))
  )
}

format_latex_var <- function(val) {
  val_chr <- as.character(val)
  dplyr::case_when(
    val_chr %in% c("horizontal_error", "horizontal error") ~ "$e_H$",
    val_chr == "Bx"  ~ "$B_x$",
    val_chr == "By"  ~ "$B_y$",
    val_chr == "Bz"  ~ "$B_z$",
    val_chr == "Kp"  ~ "$Kp$",
    TRUE             ~ str_replace_all(val_chr, "_", " ")
  )
}

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
            return(format_latex_var(val))
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

# Helper to draw base-R axes with comma-formatted thousands (e.g. 1,000 instead of 1000)
draw_comma_axes <- function(side_vec = c(1, 2)) {
  for (s in side_vec) {
    ticks <- axTicks(s)
    labs  <- format(ticks, big.mark = ",", scientific = FALSE, trim = TRUE)
    axis(s, at = ticks, labels = labs)
  }
}

copy_to_definitions <- function(file_path) {
  if (has_def_dir && file.exists(file_path)) {
    file.copy(file_path, file.path(def_dir, basename(file_path)), overwrite = TRUE)
  }
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
  stop("Merged dataset not found. Please ensure new_CPA.R has been run at least once.")
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
  caption = "Descriptive statistics for the unfiltered 2014 Darwin dataset: GPS horizontal positioning error $e_H$ [m], geomagnetic field components ($B_x, B_y, B_z$ [nT]), ionospheric TEC [TECU], ring current index Dst [nT], and planetary $Kp$ index.",
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
# 1B. Regenerate Figures 1, 2, and 3 (Unfiltered Time Series, Density, Boxplots)
#     with LaTeX Math Italics and Comma-Formatted Thousands Labels
# ==============================================================================
cat("\n==============================================================================\n")
cat("1B. Regenerating Figures 1-3 (Unfiltered Time Series, Density, Boxplots)\n")
cat("==============================================================================\n")

# Figure 1: Time Series
ts_pdf_path <- file.path(exp_dir, "GPS_Bxyz_DOY_Time_Series_Unfiltered.pdf")
cairo_pdf(file = ts_pdf_path, width = 12, height = 8)
par(mfrow = c(2, 2), mar = c(4.5, 5, 3.5, 1.5) + 0.1)

plot(full_merged_df$time, full_merged_df$horizontal_error, type = "l", col = "blue", axes = FALSE,
     main = expression(bold("GPS Horizontal Error " * italic(e)[italic(H)] * " [m] in 2014 (Unfiltered)")),
     xlab = "Time [days in 2014]", ylab = expression("Horizontal error " * italic(e)[italic(H)] * " [m]"))
draw_comma_axes(c(1, 2)); box()

plot(full_merged_df$time, full_merged_df$Bx, type = "l", col = "red", axes = FALSE,
     main = expression(bold(italic(B)[italic(x)] * " [nT] for days in 2014")),
     xlab = "Time [days in 2014]", ylab = expression(italic(B)[italic(x)] * " [nT]"))
draw_comma_axes(c(1, 2)); box()

plot(full_merged_df$time, full_merged_df$By, type = "l", col = "darkgreen", axes = FALSE,
     main = expression(bold(italic(B)[italic(y)] * " [nT] for days in 2014")),
     xlab = "Time [days in 2014]", ylab = expression(italic(B)[italic(y)] * " [nT]"))
draw_comma_axes(c(1, 2)); box()

plot(full_merged_df$time, full_merged_df$Bz, type = "l", col = "purple", axes = FALSE,
     main = expression(bold(italic(B)[italic(z)] * " [nT] for days in 2014")),
     xlab = "Time [days in 2014]", ylab = expression(italic(B)[italic(z)] * " [nT]"))
draw_comma_axes(c(1, 2)); box()

par(mfrow = c(1, 1))
dev.off()
copy_to_definitions(ts_pdf_path)

# Figure 2: Density Plots
dens_pdf_path <- file.path(exp_dir, "GPS_Bxyz_Density_Plots_Unfiltered.pdf")
cairo_pdf(file = dens_pdf_path, width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4.5, 5, 3.5, 1.5) + 0.1)

plot(density(na.omit(full_merged_df$horizontal_error)), axes = FALSE,
     main = expression(bold("Probability density: Horizontal Error " * italic(e)[italic(H)] * " [m]")),
     xlab = expression("Horizontal Error " * italic(e)[italic(H)] * " [m]"), ylab = "Density")
draw_comma_axes(c(1, 2)); box()

plot(density(na.omit(full_merged_df$Bx)), axes = FALSE,
     main = expression(bold("Probability density: " * italic(B)[italic(x)] * " [nT]")),
     xlab = expression(italic(B)[italic(x)] * " [nT]"), ylab = "Density")
draw_comma_axes(c(1, 2)); box()

plot(density(na.omit(full_merged_df$By)), axes = FALSE,
     main = expression(bold("Probability density: " * italic(B)[italic(y)] * " [nT]")),
     xlab = expression(italic(B)[italic(y)] * " [nT]"), ylab = "Density")
draw_comma_axes(c(1, 2)); box()

plot(density(na.omit(full_merged_df$Bz)), axes = FALSE,
     main = expression(bold("Probability density: " * italic(B)[italic(z)] * " [nT]")),
     xlab = expression(italic(B)[italic(z)] * " [nT]"), ylab = "Density")
draw_comma_axes(c(1, 2)); box()

par(mfrow = c(1, 1))
dev.off()
copy_to_definitions(dens_pdf_path)

# Figure 3: Boxplots
box_pdf_path <- file.path(exp_dir, "GPS_Bxyz_Boxplots_Unfiltered.pdf")
cairo_pdf(file = box_pdf_path, width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4.5, 5, 3.5, 1.5) + 0.1)

boxplot(full_merged_df$horizontal_error, axes = FALSE,
        main = expression(bold("Boxplot: Horizontal Error " * italic(e)[italic(H)] * " [m]")),
        ylab = expression("Horizontal Error " * italic(e)[italic(H)] * " [m]"))
draw_comma_axes(2); box()

boxplot(full_merged_df$Bx, axes = FALSE,
        main = expression(bold("Boxplot: " * italic(B)[italic(x)] * " [nT]")),
        ylab = expression(italic(B)[italic(x)] * " [nT]"))
draw_comma_axes(2); box()

boxplot(full_merged_df$By, axes = FALSE,
        main = expression(bold("Boxplot: " * italic(B)[italic(y)] * " [nT]")),
        ylab = expression(italic(B)[italic(y)] * " [nT]"))
draw_comma_axes(2); box()

boxplot(full_merged_df$Bz, axes = FALSE,
        main = expression(bold("Boxplot: " * italic(B)[italic(z)] * " [nT]")),
        ylab = expression(italic(B)[italic(z)] * " [nT]"))
draw_comma_axes(2); box()

par(mfrow = c(1, 1))
dev.off()
copy_to_definitions(box_pdf_path)

# ==============================================================================
# 2. Load All Change Points & Linked Pairs Across All 5 Cause Components
#    and Regenerate Component-Specific Figures 5-9 in exp_dir
# ==============================================================================
cat("\n==============================================================================\n")
cat("2. Loading Detected Change Points & Regenerating Component Figures 5-9\n")
cat("==============================================================================\n")

cp_gps_df <- read_csv(
  file.path(paste0("CPD_Darwin_2014_Bx_", epoch_bin_sec),
            paste0("ChangePoints_GPS_Horizontal_Error_", epoch_bin_sec, ".csv")),
  show_col_types = FALSE
)

cp_cause_list     <- list()
linked_pairs_list <- list()
max_lag_default   <- 14400

render_component_delay_plot <- function(delays_sec, delay_density, delay_mean, delay_sd, comp, max_lag_sec) {
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 2) + 0.1)
  
  max_lag_fmt <- format(max_lag_sec, big.mark = ",", scientific = FALSE, trim = TRUE)
  mean_fmt    <- format(round(delay_mean, 1), big.mark = ",", nsmall = 1, scientific = FALSE, trim = TRUE)
  sd_fmt      <- format(round(delay_sd, 1), big.mark = ",", nsmall = 1, scientific = FALSE, trim = TRUE)
  
  if (comp == "Bx") {
    main_box <- bquote(bold("Boxplot of GPS Response Delays"))
    sub_box  <- bquote("(Lag behind " * italic(B)[italic(x)] * ", max = " * .(max_lag_fmt) * " s)")
    xlab_box <- expression(italic(B)[italic(x)] * " " %->% " GPS Horizontal Error " * italic(e)[italic(H)])
  } else if (comp == "By") {
    main_box <- bquote(bold("Boxplot of GPS Response Delays"))
    sub_box  <- bquote("(Lag behind " * italic(B)[italic(y)] * ", max = " * .(max_lag_fmt) * " s)")
    xlab_box <- expression(italic(B)[italic(y)] * " " %->% " GPS Horizontal Error " * italic(e)[italic(H)])
  } else if (comp == "Bz") {
    main_box <- bquote(bold("Boxplot of GPS Response Delays"))
    sub_box  <- bquote("(Lag behind " * italic(B)[italic(z)] * ", max = " * .(max_lag_fmt) * " s)")
    xlab_box <- expression(italic(B)[italic(z)] * " " %->% " GPS Horizontal Error " * italic(e)[italic(H)])
  } else if (comp == "TEC") {
    main_box <- bquote(bold("Boxplot of GPS Response Delays"))
    sub_box  <- bquote("(Lag behind TEC, max = " * .(max_lag_fmt) * " s)")
    xlab_box <- expression("TEC " %->% " GPS Horizontal Error " * italic(e)[italic(H)])
  } else {
    main_box <- bquote(bold("Boxplot of GPS Response Delays"))
    sub_box  <- bquote("(Lag behind Dst, max = " * .(max_lag_fmt) * " s)")
    xlab_box <- expression("Dst " %->% " GPS Horizontal Error " * italic(e)[italic(H)])
  }
  
  # 1. Boxplot of response delays with comma-formatted Y axis (e.g. 2,000; 4,000; ...; 14,000)
  boxplot(
    delays_sec,
    main   = "",
    ylab   = "Response Delay [seconds]",
    xlab   = xlab_box,
    col    = "lightblue",
    border = "darkblue",
    axes   = FALSE,
    ylim   = c(0, max_lag_sec)
  )
  title(main = main_box, line = 2.2)
  mtext(sub_box, side = 3, line = 0.7, cex = 0.95, font = 2)
  y_ticks <- seq(0, 14000, by = 2000)
  axis(2, at = y_ticks, labels = format(y_ticks, big.mark = ",", scientific = FALSE, trim = TRUE))
  box()
  points(1, delay_mean, col = "red", pch = 19, cex = 1.3)
  legend("topright", legend = c("Arithmetic Mean"), col = c("red"), pch = 19, bty = "n")
  
  # 2. Experimental probability density function with comma-formatted X axis
  plot(
    delay_density,
    main = "",
    xlab = "Response Delay [seconds]",
    ylab = "Density",
    col  = "red",
    lwd  = 2,
    axes = FALSE
  )
  title(main = bquote(bold("Experimental PDF of Response Delays")), line = 2.2)
  mtext(bquote("(Mean = " * .(mean_fmt) * " s, SD = " * .(sd_fmt) * " s)"),
        side = 3, line = 0.7, cex = 0.95, font = 2)
  x_ticks <- seq(0, 14000, by = 2000)
  axis(1, at = x_ticks, labels = format(x_ticks, big.mark = ",", scientific = FALSE, trim = TRUE))
  axis(2, at = axTicks(2), labels = format(axTicks(2), scientific = FALSE, trim = TRUE))
  box()
  polygon(delay_density, col = rgb(1, 0, 0, 0.15), border = "red")
  abline(v = delay_mean, col = "blue", lwd = 2, lty = 2)
  legend("topright", legend = c("Experimental PDF", "Arithmetic Mean"),
         col = c("red", "blue"), lwd = 2, lty = c(1, 2), bty = "n")
  
  par(mfrow = c(1, 1))
}

for (comp in cause_cols) {
  comp_dir <- paste0("CPD_Darwin_2014_", comp, "_", epoch_bin_sec)
  cp_file  <- file.path(comp_dir, paste0("ChangePoints_Geomagnetic_", comp, "_", epoch_bin_sec, ".csv"))
  lnk_file <- file.path(comp_dir, paste0("Linked_CP_Pairs_", comp, "_to_GPS_", epoch_bin_sec, ".csv"))
  
  if (file.exists(cp_file)) {
    cp_cause_list[[comp]] <- read_csv(cp_file, show_col_types = FALSE) %>%
      mutate(Cause_Component = comp)
  }
  if (file.exists(lnk_file)) {
    lnk_df <- read_csv(lnk_file, show_col_types = FALSE) %>%
      mutate(Cause_Component = comp)
    linked_pairs_list[[comp]] <- lnk_df
    
    # Regenerate Delay_Boxplot_and_Density_<comp>_to_GPS_300.pdf & .png inside exp_dir
    delays_sec    <- lnk_df$delay_sec
    delay_mean    <- mean(delays_sec, na.rm = TRUE)
    delay_sd      <- sd(delays_sec, na.rm = TRUE)
    delay_density <- density(delays_sec, from = 0, to = max_lag_default)
    
    pdf_out <- file.path(exp_dir, paste0("Delay_Boxplot_and_Density_", comp, "_to_GPS_", epoch_bin_sec, ".pdf"))
    png_out <- file.path(exp_dir, paste0("Delay_Boxplot_and_Density_", comp, "_to_GPS_", epoch_bin_sec, ".png"))
    
    cairo_pdf(file = pdf_out, width = 11, height = 5.5)
    render_component_delay_plot(delays_sec, delay_density, delay_mean, delay_sd, comp, max_lag_default)
    dev.off()
    copy_to_definitions(pdf_out)
    
    png(filename = png_out, width = 11, height = 5.5, units = "in", res = 300)
    render_component_delay_plot(delays_sec, delay_density, delay_mean, delay_sd, comp, max_lag_default)
    dev.off()
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

kw_test    <- kruskal.test(delay_sec ~ Cause_Component, data = all_linked_pairs_df)
anova_test <- summary(aov(delay_sec ~ Cause_Component, data = all_linked_pairs_df))

cat(sprintf("Kruskal-Wallis test across 5 components: chi-squared = %.4f, df = %d, p-value = %.4g\n",
            kw_test$statistic, kw_test$parameter, kw_test$p.value))
print(anova_test)

set.seed(42)
n_mc <- 500
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

hypothesis_table <- map_dfr(cause_cols, function(comp) {
  sub_df <- all_linked_pairs_df %>% filter(Cause_Component == comp)
  delays <- sub_df$delay_sec
  n_obs  <- length(delays)
  n_tot  <- nrow(all_cause_cps_df %>% filter(Cause_Component == comp))
  
  ks_res   <- ks.test(delays, "punif", 0, max_lag_default)
  wilc_res <- wilcox.test(delays, mu = max_lag_default / 2)
  
  p_mc_n    <- (sum(mc_null_df$n_pairs >= n_obs) + 1) / (n_mc + 1)
  z_score_n <- (n_obs - mean(mc_null_df$n_pairs)) / sd(mc_null_df$n_pairs)
  
  scaled_delays <- pmax(delays, 1)
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

lag_grid_sec <- c(3600, 7200, 10800, 14400, 21600)

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
  }) %>%
  mutate(Cause_Component = factor(Cause_Component, levels = cause_cols))

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
      Kp_daily_max  = max(Kp_daily_max, na.rm = TRUE),
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
  
  print(regime_summary_df)
  write_csv(regime_summary_df, file.path(exp_dir, "Stratified_Delay_By_Kp_Regime.csv"))
  
  df_to_custom_tabularx(
    df = regime_summary_df,
    caption = "Response delay statistics stratified by geomagnetic activity regime: Active/Storm days ($Kp_{\\max} \\ge 4$, 61 days) versus Quiet days ($Kp_{\\max} < 4$, 303 active days).",
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
# 6. Publication-Ready Comparative Figures (with LaTeX Math Italics & Commas)
# ==============================================================================
cat("\n==============================================================================\n")
cat("6. Generating Publication-Ready Multi-Component Figures\n")
cat("==============================================================================\n")

# Figure 10: Combined Violin + Boxplot + Mean Comparison of Response Delays
p_violin <- ggplot(all_linked_pairs_df, aes(x = Cause_Component, y = delay_sec / 60, fill = Cause_Component)) +
  geom_violin(alpha = 0.35, color = "black", trim = FALSE) +
  geom_boxplot(width = 0.22, outlier.shape = 16, outlier.size = 2, alpha = 0.85, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 4, color = "red") +
  scale_x_discrete(labels = parse_math_labels) +
  scale_y_continuous(limits = c(0, 260), breaks = seq(0, 240, by = 30), labels = label_comma()) +
  labs(
    title = expression(bold("Distribution of GPS Positioning Error (" * italic(e)[italic(H)] * ") Response Delays")),
    subtitle = "Across Cause Components (Darwin, 2014)",
    x       = "Cause Component (Change-Point Source)",
    y       = "Response Delay [minutes]"
  ) +
  bigger_theme +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(b = 15)),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

violin_pdf <- file.path(exp_dir, "Combined_Delay_Violin_Boxplot_All_Components.pdf")
ggsave(violin_pdf, plot = p_violin, width = 11, height = 7, device = cairo_pdf)
ggsave(file.path(exp_dir, "Combined_Delay_Violin_Boxplot_All_Components.png"),
       plot = p_violin, width = 11, height = 7, dpi = 300)
copy_to_definitions(violin_pdf)

# Figure 11: Overlaid Experimental Probability Density Functions (PDFs)
p_dens <- ggplot(all_linked_pairs_df, aes(x = delay_sec, color = Cause_Component, fill = Cause_Component)) +
  geom_density(alpha = 0.15, linewidth = 1.2, bounds = c(0, max_lag_default)) +
  geom_hline(yintercept = 1 / max_lag_default, linetype = "dashed", color = "black", linewidth = 0.9) +
  scale_x_continuous(breaks = seq(0, 14400, by = 3600),
                     labels = c("0 (0h)", "3,600 (1h)", "7,200 (2h)", "10,800 (3h)", "14,400 (4h)")) +
  scale_color_discrete(labels = parse_math_labels) +
  scale_fill_discrete(guide = "none") + # Hides the redundant fill legend
  labs(
    title = expression(bold("Experimental PDFs of GPS Positioning Error (" * italic(e)[italic(H)] * ") Response Delays")),
    subtitle = "vs. Uniform Null Baseline (Dashed Line)",
    x       = "Response Delay [seconds (hours)]",
    y       = "Probability Density",
    color   = "Cause Component"
  ) +
  bigger_theme +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(b = 15)),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

pdfs_pdf <- file.path(exp_dir, "Combined_Delay_PDFs_All_Components.pdf")
ggsave(pdfs_pdf, plot = p_dens, width = 11.5, height = 7, device = cairo_pdf)
ggsave(file.path(exp_dir, "Combined_Delay_PDFs_All_Components.png"),
       plot = p_dens, width = 11.5, height = 7, dpi = 300)
copy_to_definitions(pdfs_pdf)

# Figure 12: Sensitivity Curves (Matching Rate % and Mean Delay vs Max Allowed Lag x)
p_sens1 <- ggplot(sensitivity_df, aes(x = Max_Lag_Hours, y = Match_Rate_Pct, color = Cause_Component, group = Cause_Component)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 6), labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  scale_color_discrete(labels = parse_math_labels) +
  labs(
    title = expression(bold("Sensitivity of CP Linking Rate to Max Allowed Lag (" * italic(x) * ")")),
    x     = expression("Maximum Allowed Lag " * italic(x) * " [hours]"),
    y     = "Linked CP Pairs [% of Total CPs]",
    color = "Cause"
  ) +
  bigger_theme

p_sens2 <- ggplot(sensitivity_df, aes(x = Max_Lag_Hours, y = Mean_Delay_s / 60, color = Cause_Component, group = Cause_Component)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 6), labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  scale_color_discrete(labels = parse_math_labels) +
  labs(
    title = expression(bold("Sensitivity of Mean Response Delay to Max Allowed Lag (" * italic(x) * ")")),
    x     = expression("Maximum Allowed Lag " * italic(x) * " [hours]"),
    y     = "Arithmetic Mean Delay [minutes]",
    color = "Cause"
  ) +
  bigger_theme

p_sens_combined <- arrangeGrob(p_sens1, p_sens2, ncol = 2)
sens_pdf <- file.path(exp_dir, "Sensitivity_Analysis_Max_Lag_Plot.pdf")
ggsave(sens_pdf, plot = p_sens_combined, width = 15, height = 6, device = cairo_pdf)
ggsave(file.path(exp_dir, "Sensitivity_Analysis_Max_Lag_Plot.png"),
       plot = p_sens_combined, width = 15, height = 6, dpi = 300)
copy_to_definitions(sens_pdf)

# Figure 13: Stratified Boxplot (Active/Storm vs Quiet Geomagnetic Conditions)
if (exists("stratified_pairs_df")) {
  p_strat <- ggplot(stratified_pairs_df, aes(x = Cause_Component, y = delay_sec / 60, fill = Regime)) +
    geom_boxplot(position = position_dodge(0.8), width = 0.65, alpha = 0.85, color = "black") +
    scale_x_discrete(labels = parse_math_labels) +
    scale_y_continuous(labels = label_comma()) +
    scale_fill_discrete(
      labels = c(
        "Active/Storm (Kp >= 4)" = expression("Active/Storm (" * italic(Kp)[max] >= 4 * ")"),
        "Quiet (Kp < 4)"         = expression("Quiet (" * italic(Kp)[max] < 4 * ")")
      )
    ) +
    labs(
      title = expression(bold("GPS Positioning Error (" * italic(e)[italic(H)] * ") Response Delays by Geomagnetic Activity Regime")),
      x     = "Cause Component",
      y     = "Response Delay [minutes]",
      fill  = expression(bold("GFZ " * italic(Kp) * " Regime"))
    ) +
    bigger_theme +
    theme(legend.position = "top")
  
  strat_pdf <- file.path(exp_dir, "Stratified_Delay_Kp_Regime_Boxplot.pdf")
  ggsave(strat_pdf, plot = p_strat, width = 11, height = 6.5, device = cairo_pdf)
  ggsave(file.path(exp_dir, "Stratified_Delay_Kp_Regime_Boxplot.png"),
         plot = p_strat, width = 11, height = 6.5, dpi = 300)
  copy_to_definitions(strat_pdf)
}

# Figure 4: Diurnal (Time-of-Day UT) Distribution of Cause vs Effect Change Points
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
  scale_x_continuous(breaks = seq(0, 24, by = 4), limits = c(0, 24), labels = label_comma()) +
  scale_color_discrete(labels = parse_math_labels) +
  scale_fill_discrete(guide = "none") + # Hides the redundant fill legend
  labs(
    title = "Diurnal (UT Hour) Density of Detected Change Points",
    subtitle = "Darwin, 2014",
    x       = "Universal Time [Hours of Day]",
    y       = "Density",
    color   = "Series"
  ) +
  bigger_theme +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(b = 15)),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

diurnal_pdf <- file.path(exp_dir, "Diurnal_Distribution_ChangePoints_UT.pdf")
ggsave(diurnal_pdf, plot = p_diurnal, width = 11, height = 7, device = cairo_pdf)
ggsave(file.path(exp_dir, "Diurnal_Distribution_ChangePoints_UT.png"),
       plot = p_diurnal, width = 11, height = 7, dpi = 300)
copy_to_definitions(diurnal_pdf)

cat("\nExpanded CPA analysis completed successfully. All tables and figures saved in:", exp_dir, "\n")
sink()