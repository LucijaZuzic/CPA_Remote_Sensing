# Clear environment
rm(list=ls())

## The R script for Žužić et al. manuscript

library(tidyverse)

# Avoid namespace conflicts
filter <- dplyr::filter
select <- dplyr::select

getCurrentFileLocation<-function()
{
  this_file<-commandArgs() %>% 
    tibble::enframe(name=NULL) %>%
    tidyr::separate(col=value,into=c('key','value'),sep='=',fill='right') %>%
    dplyr::filter(key=='--file') %>%
    dplyr::pull(value)
  if (length(this_file)==0)
  {      this_file<-rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}

# Declaration of the R libraries utilised
library(data.table)
library(lubridate)
library(geosphere)
library(readr)
library(stringr)
library(mcp)
library(posterior)
library(tidybayes)
library(loo)
library(coda)
library(DataExplorer)
library(fitdistrplus)
library(gridExtra)

# Set path to the working directory

setwd(getCurrentFileLocation())

# To make full-year 2014 Bayesian MCMC computationally tractable while preserving
# high temporal resolution (seconds), we aggregate to 5-minute (300 s) epochs within each day
# (288 points/day) without filtering any outliers, and run `mcp` across the days of 2014.
epoch_bin_sec <- 300

for (selected_B_col in c("Bx", "By", "Bz", "TEC", "Dst")) {
  
  # Default to Bx (horizontal northward component, most sensitive to equatorial ring current/storms),
  # or automatically pick the component with largest standardized jump / kurtosis if desired:
  
  # Create output directory for Change Point Detection results and plots
  out_dir <- paste0("CPD_Darwin_2014_", selected_B_col, "_", epoch_bin_sec)
  dir.create(file.path(out_dir), showWarnings = FALSE)
  
  # Start logging console output
  sink()
  closeAllConnections()
  sink(file.path(out_dir, paste0("CPD_analysis_output_Darwin_2014_", selected_B_col, "_", epoch_bin_sec, ".txt")), split = TRUE)
  
  oldw <- getOption("warn")
  options(warn = -1)
  
  # Lock RNG seed for reproducibility
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(42)
  
  # ==============================================================================
  # 1. GPS Horizontal Positioning Error Time Series (Darwin, 2014)
  # ==============================================================================
  cat("\n==============================================================================\n")
  cat("1. Loading Darwin 2014 GPS Horizontal Positioning Error Time Series\n")
  cat("==============================================================================\n")
  
  # Darwin reference coordinates
  ref_coords <- list(
    lon    = 131.132744,
    lat    = -12.843697,
    height = 125.1
  )
  
  station_name   <- "Darwin"
  station_abbrev <- "darw"
  station_long   <- "Darwin, NT, Australia"
  numcol         <- 8
  num_days       <- 365
  
  # ==============================================================================
  # 1. Pre-generated unfiltered CSV from make_time_series_Darwin_error.R already exists
  # ==============================================================================
  
  cat("Reading pre-compiled unfiltered Darwin POS CSV...\n")
  
  file_path <- paste0("afullframe_", station_name, "_new_365.csv")
  
  gps_darwin_df <- read_csv(file_path, show_col_types = FALSE)
  
  cat("Total unfiltered GPS epochs loaded for Darwin:", nrow(gps_darwin_df), "\n")
  
  # ==============================================================================
  # 2. Load Geomagnetic Field Components (Bx, By, Bz) and Merge
  # ==============================================================================
  cat("\n==============================================================================\n")
  cat("2. Loading Geomagnetic Field Data (Bx, By, Bz) & Merging with GPS Data\n")
  cat("==============================================================================\n")
  
  darwin_data_dir <- "darwin_data"
  
  geom1               <- list()
  iono1               <- list()
  all_previous_months <- 0
  month_names         <- c("January", "February", "March", "April", "May", "June",
                           "July", "August", "September", "October", "November", "December")
  
  for (j in 1:12) {
    month_file_name <- file.path(darwin_data_dir, paste0("data_all_", j, ".csv"))
    if (!file.exists(month_file_name)) {
      cat("Warning: Missing file", month_file_name, "\n")
      next
    }
    
    # Read full monthly dataset WITHOUT filtering TEC >= limited or any outliers
    data_month <- as.data.frame(read.csv(month_file_name, header = TRUE, sep = ",", skip = 0))
    
    day_value <- all_previous_months + data_month[["d"]]
    timer     <- round(day_value + data_month[["total.seconds"]] / 86400, digits = 4)
    
    geomi <- data.frame(
      time               = timer,
      DOY                = as.numeric(day_value - 1),
      total_seconds_day  = as.numeric(data_month[["total.seconds"]]),
      total_seconds      = (as.numeric(day_value) - 1) * 86400 + as.numeric(data_month[["total.seconds"]]),
      month              = as.numeric(data_month[["m"]]),
      Bx                 = as.numeric(data_month[["Bx"]]),
      By                 = as.numeric(data_month[["By"]]),
      Bz                 = as.numeric(data_month[["Bz"]]),
      TEC                = as.numeric(data_month[["TEC"]]),
      Dst                = as.numeric(data_month[["Dst"]])
    )
    
    # If horizontal(m) is also inside data_all_j.csv, retain it for verification
    if ("horizontal.m." %in% colnames(data_month)) {
      geomi$horizontal_csv <- as.numeric(data_month[["horizontal.m."]])
    }
    
    geom1[[j]] <- geomi
    all_previous_months <- all_previous_months + length(unique(data_month[["d"]]))
    cat(sprintf("Loaded %s: %d rows (Unfiltered)\n", month_names[j], nrow(geomi)))
  }
  
  geom_df <- bind_rows(geom1)
  
  # Merge Darwin GPS horizontal error (`horizontal(deg)`) with Geomagnetic Bx, By, Bz
  # Primary merge key: DOY and total_seconds_day
  merged_df <- gps_darwin_df %>%
    rename(
      DOY               = doy,
      total_seconds_day = `total seconds day`,
      horizontal_error  = `horizontal(deg)`
    ) %>%
    inner_join(
      geom_df %>% select(DOY, total_seconds_day, time, total_seconds, month, Bx, By, Bz, TEC, Dst),
      by = c("DOY", "total_seconds_day")
    ) %>%
    arrange(total_seconds)
  
  cat("Total merged unfiltered observations (Darwin GPS + Kakadu Bx, By, Bz):", nrow(merged_df), "\n")
  write_csv(merged_df, file.path(out_dir, "Darwin_2014_Merged_Unfiltered_GPS_Bxyz.csv"))
  
  # ==============================================================================
  # 3. Standard Descriptive Statistical Analysis & Component Selection
  # ==============================================================================
  cat("\n==============================================================================\n")
  cat("3. Descriptive Statistical Analysis of Unfiltered Time Series\n")
  cat("==============================================================================\n")
  
  vars_to_analyse <- c("horizontal_error", "Bx", "By", "Bz")
  
  desc_stats_table <- map_dfr(vars_to_analyse, function(v) {
    x <- na.omit(merged_df[[v]])
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
  
  print(desc_stats_table)
  write_csv(desc_stats_table, file.path(out_dir, "Descriptive_Statistics_GPS_Bxyz.csv"))
  
  cat("\n--- Standard R summary() Output ---\n")
  print(summary(merged_df[, vars_to_analyse]))
  
  cat("\n--- Pearson Correlation Matrix ---\n")
  print(cor(merged_df[, vars_to_analyse], use = "complete.obs"))
  
  # Diagnostic plots (Time Series, Density, Boxplots)
  pdf(file = file.path(out_dir, "GPS_Bxyz_DOY_Time_Series_Unfiltered.pdf"), width = 12, height = 8)
  par(mfrow = c(2, 2))
  plot(merged_df$time, merged_df$horizontal_error, type = "l", col = "blue",
       main = "GPS Horizontal Error [m] for days in 2014 (Unfiltered)",
       xlab = "Time [days in 2014]", ylab = "Horizontal positioning error [m]")
  plot(merged_df$time, merged_df$Bx, type = "l", col = "red",
       main = "Bx [nT] for days in 2014", xlab = "Time [days in 2014]", ylab = "Bx [nT]")
  plot(merged_df$time, merged_df$By, type = "l", col = "darkgreen",
       main = "By [nT] for days in 2014", xlab = "Time [days in 2014]", ylab = "By [nT]")
  plot(merged_df$time, merged_df$Bz, type = "l", col = "purple",
       main = "Bz [nT] for days in 2014", xlab = "Time [days in 2014]", ylab = "Bz [nT]")
  par(mfrow = c(1, 1))
  dev.off()
  
  pdf(file = file.path(out_dir, "GPS_Bxyz_Density_Plots_Unfiltered.pdf"), width = 10, height = 8)
  par(mfrow = c(2, 2))
  plot(density(na.omit(merged_df$horizontal_error)),
       main = "Probability density: Horizontal Error [m]", xlab = "Horizontal Error [m]", ylab = "Density")
  plot(density(na.omit(merged_df$Bx)),
       main = "Probability density: Bx [nT]", xlab = "Bx [nT]", ylab = "Density")
  plot(density(na.omit(merged_df$By)),
       main = "Probability density: By [nT]", xlab = "By [nT]", ylab = "Density")
  plot(density(na.omit(merged_df$Bz)),
       main = "Probability density: Bz [nT]", xlab = "Bz [nT]", ylab = "Density")
  par(mfrow = c(1, 1))
  dev.off()
  
  pdf(file = file.path(out_dir, "GPS_Bxyz_Boxplots_Unfiltered.pdf"), width = 10, height = 8)
  par(mfrow = c(2, 2))
  boxplot(merged_df$horizontal_error, main = "Boxplot: Horizontal Error [m]", ylab = "Horizontal Error [m]")
  boxplot(merged_df$Bx, main = "Boxplot: Bx [nT]", ylab = "Bx [nT]")
  boxplot(merged_df$By, main = "Boxplot: By [nT]", ylab = "By [nT]")
  boxplot(merged_df$Bz, main = "Boxplot: Bz [nT]", ylab = "Bz [nT]")
  par(mfrow = c(1, 1))
  dev.off()
  
  # Select the geomagnetic component (Bx, By, or Bz) where special events are best visible.
  # In sub-equatorial regions (Kakadu/Darwin), geomagnetic storms and Dst depressions
  # are most strongly manifested in the northward horizontal component Bx (or highest excess kurtosis / SD).
  b_stats <- desc_stats_table %>% filter(Variable %in% c("Bx", "By", "Bz"))
  cat("\nGeomagnetic Component Event Visibility Metrics (SD and Excess Kurtosis):\n")
  print(b_stats %>% select(Variable, Mean, SD, IQR, Skewness, Kurtosis))
  
  cat("\nSelected geomagnetic component for Change Point Detection:", selected_B_col, "\n")
  
  # ==============================================================================
  # 4. Change Point Detection using `mcp`
  # ==============================================================================
  cat("\n==============================================================================\n")
  cat("4. Running Bayesian Change Point Detection (mcp) on Darwin 2014 Time Series\n")
  cat("==============================================================================\n")
  
  # Define the 2-changepoint model:
  # Segment 1: baseline plateau
  # Segment 2: joined slope at cp_1 (onset of disturbance/event)
  # Segment 3: disjoined slope at cp_2 (recovery / new regime)
  mcp_model_2cp <- list(
    response ~ 1,  # plateau (Intercept_1)
    ~ 0 + time,    # joined slope (time_2) at cp_1
    ~ 1 + time     # disjoined slope (Intercept_3, time_3) at cp_2
  )
  
  # Define the 1-changepoint null model for LOO comparison:
  mcp_model_null <- list(
    response ~ 1 + time, # intercept (Intercept_1) and slope (time_1)
    ~ 1 + time           # disjoined slope (Intercept_2, time_2) at cp_1
  )
  
  # Helper function to fit mcp on a 1-day (or window) time series and extract CPs
  run_mcp_window <- function(df_window, response_col, time_col = "total_seconds",
                             window_label = "", plot_fit = FALSE, compare_null = FALSE) {
    df_sub <- df_window %>%
      select(time = all_of(time_col), response = all_of(response_col)) %>%
      drop_na() %>%
      arrange(time)
    
    if (nrow(df_sub) < 10 || sd(df_sub$response) == 0) {
      cat("Skipped due to low rows or zero variance in:", window_label, "\n")
      return(tibble())
    }
    
    # Fit main mcp model
    fit <- tryCatch(
      mcp(mcp_model_2cp, data = df_sub, par_x = "time",
          sample = "post", adapt = 1000, iter = 2000, chains = 3),
      error = function(e) {
        cat("MCP Error in", window_label, ":", e$message, "\n")
        NULL
      }
    )
    
    if (is.null(fit)) return(tibble())
    
    # Extract summary table of parameters
    sum_fit <- summary(fit)
    
    # Optional LOO comparison against null model
    if (compare_null) {
      fit_null <- tryCatch(
        mcp(mcp_model_null, data = df_sub, par_x = "time",
            sample = "post", adapt = 1000, iter = 2000, chains = 3),
        error = function(e) NULL
      )
      if (!is.null(fit_null)) {
        fit$loo      <- tryCatch(loo(fit), error = function(e) NULL)
        fit_null$loo <- tryCatch(loo(fit_null), error = function(e) NULL)
        if (!is.null(fit$loo) && !is.null(fit_null$loo)) {
          cat("\n--- LOO Comparison for", window_label, response_col, "---\n")
          print(loo::loo_compare(fit$loo, fit_null$loo))
        }
      }
    }
    
    # Save diagnostic plots for representative windows
    if (plot_fit) {
      p_fit  <- plot(fit) + ggtitle(paste("MCP Fit:", response_col, "-", window_label))
      p_pars <- plot_pars(fit, regex_pars = "cp_")
      ggsave(file.path(out_dir, paste0("mcp_fit_", response_col, "_", window_label, "_", epoch_bin_sec, ".pdf")),
             plot = p_fit, width = 8, height = 5)
      ggsave(file.path(out_dir, paste0("mcp_pars_", response_col, "_", window_label, "_", epoch_bin_sec, ".pdf")),
             plot = p_pars, width = 8, height = 5)
      
      # Demonstrate tidybayes & coda posterior extraction
      draws_df <- tidybayes::tidy_draws(fit$mcmc_post)
      mcmc_obj <- coda::mcmc(do.call(rbind, fit$mcmc_post))
    }
    
    # Extract change points (cp_1, cp_2, ...) safely
    sum_df <- as.data.frame(sum_fit)
    
    # If 'name' column doesn't exist, check if parameters are in rownames or another column
    if (!"name" %in% colnames(sum_df)) {
      if ("parameter" %in% colnames(sum_df)) {
        sum_df <- sum_df %>% rename(name = parameter)
      } else {
        sum_df <- sum_df %>% tibble::rownames_to_column("name")
      }
    }
    
    cp_rows <- sum_df %>%
      as_tibble() %>%
      filter(grepl("^cp_", name)) %>%
      mutate(
        series       = response_col,
        window_label = window_label,
        cp_time_sec  = mean,
        cp_lower_sec = lower,
        cp_upper_sec = upper
      )
    
    return(cp_rows)
  }
  
  # Cross-reference with GFZ Potsdam Kp index (https://kp.gfz.de/en/data)
  kp_file <- "kp_2014.txt"
  if (!file.exists(kp_file)) {
    tryCatch({
      download.file("https://kp.gfz.de/app/files/Kp_ap_Ap_SN_F107_since_1932.txt",
                    kp_file, quiet = TRUE)
    }, error = function(e) {
      message("Could not auto-download GFZ Kp file; place Kp_ap_Ap_SN_F107_since_1932.txt as kp_2014.txt")
    })
  }
  
  if (file.exists(kp_file)) {
    kp_lines <- readLines(kp_file)
    kp_lines <- kp_lines[!grepl("^#", kp_lines)]
    kp_raw <- read.table(text = kp_lines, header = FALSE)
    kp_2014 <- kp_raw %>%
      filter(V1 == 2014) %>%
      select(y = V1, m = V2, d = V3,
             Kp0 = V8, Kp3 = V9, Kp6 = V10, Kp9 = V11,
             Kp12 = V12, Kp15 = V13, Kp18 = V14, Kp21 = V15, Ap = V24) %>%
      pivot_longer(cols = starts_with("Kp"), names_to = "block_str", values_to = "Kp") %>%
      mutate(kp_block = as.numeric(sub("Kp", "", block_str))) %>%
      group_by(y, m, d) %>%
      mutate(Kp_daily_max = max(Kp, na.rm = TRUE)) %>%
      ungroup()
    
    # Add kp_block to merged_df based on total_seconds_day (each block spans 3 hours = 10800 seconds)
    # Blocks are at 0, 3, 6, 9, 12, 15, 18, 21 hours
    out_binned_df <- merged_df %>%
      mutate(kp_block = floor(total_seconds_day / 10800) * 3) %>%
      # Join the specific 3-hour Kp index
      left_join(kp_2014, by = c("y", "m", "d", "kp_block"))
    
    n_out <- length(out_binned_df$DOY)
    n_storm <- length(out_binned_df[out_binned_df$Kp >= 5, ]$DOY)
    n_active <- length(out_binned_df[out_binned_df$Kp_daily_max >= 4, ]$DOY)
    n_quiet <- length(out_binned_df[out_binned_df$Kp_daily_max < 4, ]$DOY)
    cat("\n--- GFZ Kp Verification of Extreme Outliers ---\n")
    cat(sprintf("Total epochs: %d\n", n_out))
    cat(sprintf("Storm conditions (Kp >= 5): %d epochs (%.2f%%)\n", n_storm, 100 * n_storm / n_out))
    cat(sprintf("Active/Substorm conditions (Kp_daily_max >= 4): %d epochs (%.2f%%)\n", n_active, 100 * n_active / n_out))
    cat(sprintf("Quiet conditions (Kp_daily_max < 4): %d epochs (%.2f%%)\n", n_quiet, 100 * n_quiet / n_out))
    write_csv(out_binned_df, file.path(out_dir, "kp_file.csv"))
    merged_df <- out_binned_df
  }
  
  if (epoch_bin_sec == 1) {
    merged_binned_df <- merged_df
  } else {
    merged_binned_df <- merged_df %>%
      mutate(time_bin_sec = floor(total_seconds / epoch_bin_sec) * epoch_bin_sec) %>%
      group_by(DOY, time_bin_sec) %>%
      summarise(
        total_seconds    = mean(total_seconds, na.rm = TRUE),
        horizontal_error = mean(horizontal_error, na.rm = TRUE),
        Kp               = mean(Kp, na.rm = TRUE),
        Kp_daily_max     = mean(Kp_daily_max, na.rm = TRUE),
        Bx               = mean(Bx, na.rm = TRUE),
        By               = mean(By, na.rm = TRUE),
        Bz               = mean(Bz, na.rm = TRUE),
        Dst              = mean(Dst, na.rm = TRUE),
        TEC              = mean(TEC, na.rm = TRUE),
        .groups          = "drop"
      ) %>%
      arrange(total_seconds)
  }
  
  merged_binned_df_filter <- merged_binned_df[merged_binned_df$Kp_daily_max >= 4, ]
  
  available_days <- sort(unique(merged_binned_df$DOY))
  cat("Total available DOYs in merged dataset:", length(available_days), "\n")
  cat("Total available rows in merged dataset:", length(merged_binned_df$DOY), "\n")
  
  available_filter_days <- sort(unique(merged_binned_df_filter$DOY))
  cat("Total available DOYs in merged filter dataset:", length(available_filter_days), "\n")
  cat("Total available rows in merged filter dataset:", length(merged_binned_df_filter$DOY), "\n")
  
  # Run across all available days (or select most disturbed days if running interactively)
  days_to_process <- available_days
  
  cp_gps_list <- list()
  cp_B_list   <- list()
  
  for ( idx in seq_along(days_to_process) ) {
    doy_val <- days_to_process[idx]
    df_day  <- merged_binned_df %>% filter(DOY == doy_val)
    
    # Plot & run LOO comparison as demonstration
    do_diag <- TRUE
    
    cat(sprintf("Processing DOY %d (%d / %d)...\n", doy_val + 1, idx, length(days_to_process)))
    
    # 4a. Detect Change Points in GPS Horizontal Positioning Error
    cp_gps_day <- run_mcp_window(
      df_window    = df_day,
      response_col = "horizontal_error",
      time_col     = "total_seconds",
      window_label = paste0("DOY_", doy_val + 1),
      plot_fit     = do_diag,
      compare_null = do_diag
    )
    if (nrow(cp_gps_day) > 0) cp_gps_list[[length(cp_gps_list) + 1]] <- cp_gps_day
    
    # 4b. Detect Change Points in Selected Geomagnetic Component (Bx / By / Bz)
    cp_B_day <- run_mcp_window(
      df_window    = df_day,
      response_col = selected_B_col,
      time_col     = "total_seconds",
      window_label = paste0("DOY_", doy_val + 1),
      plot_fit     = do_diag,
      compare_null = do_diag
    )
    if (nrow(cp_B_day) > 0) cp_B_list[[length(cp_B_list) + 1]] <- cp_B_day
  }
  
  cp_gps_df <- bind_rows(cp_gps_list) %>% arrange(cp_time_sec)
  cp_B_df   <- bind_rows(cp_B_list)   %>% arrange(cp_time_sec)
  
  write_csv(cp_gps_df, file.path(out_dir, paste0("ChangePoints_GPS_Horizontal_Error_", epoch_bin_sec, ".csv")))
  write_csv(cp_B_df,   file.path(out_dir, paste0("ChangePoints_Geomagnetic_", selected_B_col, "_", epoch_bin_sec, ".csv")))
  
  cat("\nTotal Change Points detected in GPS Horizontal Error:", nrow(cp_gps_df), "\n")
  cat("Total Change Points detected in Geomagnetic", selected_B_col, ":", nrow(cp_B_df), "\n")
  
  # ==============================================================================
  # 5. Cause-and-Effect Linking Algorithm (Geomagnetic CP -> GPS Error CP)
  # ==============================================================================
  cat("\n==============================================================================\n")
  cat("5. Linking Cause (", selected_B_col, "CP) and Effect (GPS Horizontal Error CP)\n")
  cat("==============================================================================\n")
  
  # Maximum allowed response delay `x` in seconds:
  # A pair of CPs (one in B, one in GPS horizontal error) is linked ONLY if:
  #   0 <= (t_GPS_CP - t_B_CP) <= max_lag_seconds
  # Otherwise, they are considered NOT connected.
  max_lag_seconds <- 14400 # 4 hours (14,400 seconds); adjust x as needed (e.g., 3600 to 21600 s)
  
  link_cause_effect_cps <- function(cause_cps_sec, effect_cps_sec, max_lag_sec) {
    cause_cps_sec  <- sort(na.omit(unique(as.numeric(cause_cps_sec))))
    effect_cps_sec <- sort(na.omit(unique(as.numeric(effect_cps_sec))))
    
    if (length(cause_cps_sec) == 0 || length(effect_cps_sec) == 0) {
      return(tibble())
    }
    
    # Build all candidate cause-effect pairs where GPS CP lags behind B CP by [0, max_lag_sec]
    candidates <- expand_grid(
      cause_idx = seq_along(cause_cps_sec),
      effect_idx = seq_along(effect_cps_sec)
    ) %>%
      mutate(
        t_cause_B_sec   = cause_cps_sec[cause_idx],
        t_effect_GPS_sec = effect_cps_sec[effect_idx],
        delay_sec       = t_effect_GPS_sec - t_cause_B_sec
      ) %>%
      filter(delay_sec >= 0 & delay_sec <= max_lag_sec) %>%
      arrange(delay_sec, t_cause_B_sec)
    
    if (nrow(candidates) == 0) {
      return(tibble())
    }
    
    # Greedy 1-to-1 matching (each B CP and each GPS CP is matched at most once,
    # prioritizing the closest valid causal pair in time)
    used_cause  <- rep(FALSE, length(cause_cps_sec))
    used_effect <- rep(FALSE, length(effect_cps_sec))
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
    
    linked_df <- bind_rows(matched_rows) %>%
      arrange(t_cause_B_sec) %>%
      mutate(
        pair_id     = row_number(),
        DOY_approx  = floor(t_cause_B_sec / 86400) + 1,
        delay_min   = delay_sec / 60,
        delay_hours = delay_sec / 3600
      )
    
    return(linked_df)
  }
  
  linked_pairs_df <- link_cause_effect_cps(
    cause_cps_sec  = cp_B_df$cp_time_sec,
    effect_cps_sec = cp_gps_df$cp_time_sec,
    max_lag_sec    = max_lag_seconds
  )
  
  write_csv(linked_pairs_df, file.path(out_dir, paste0("Linked_CP_Pairs_", selected_B_col, "_to_GPS_", epoch_bin_sec, ".csv")))
  
  cat(sprintf("Maximum allowed lag (x): %d seconds (%.2f hours)\n", max_lag_seconds, max_lag_seconds / 3600))
  cat("Number of linked Cause-Effect CP pairs identified:", nrow(linked_pairs_df), "\n")
  print(head(linked_pairs_df, 20))
  
  # ==============================================================================
  # 6. Statistical Analysis of Delay Values for Selected Pairs of CPs
  #    (Arithmetic Mean, Variance, Standard Deviation, Boxplot, and density())
  # ==============================================================================
  cat("\n==============================================================================\n")
  cat("6. Statistical Analysis of Response Delays (GPS CP lag behind ", selected_B_col, " CP)\n")
  cat("==============================================================================\n")
  
  if (nrow(linked_pairs_df) > 0) {
    delays_sec <- linked_pairs_df$delay_sec
    
    delay_mean <- mean(delays_sec, na.rm = TRUE)
    delay_var  <- var(delays_sec, na.rm = TRUE)
    delay_sd   <- sd(delays_sec, na.rm = TRUE)
    delay_med  <- median(delays_sec, na.rm = TRUE)
    delay_iqr  <- IQR(delays_sec, na.rm = TRUE)
    
    cat(sprintf("Number of matched CP pairs (N) : %d\n", length(delays_sec)))
    cat(sprintf("Arithmetic Mean Delay [s]      : %.4f s (%.2f min)\n", delay_mean, delay_mean / 60))
    cat(sprintf("Variance of Delay [s^2]        : %.4f s^2\n", delay_var))
    cat(sprintf("Standard Deviation [s]         : %.4f s (%.2f min)\n", delay_sd, delay_sd / 60))
    cat(sprintf("Response (Mean +/- SD) [s]     : %.2f +/- %.2f s\n", delay_mean, delay_sd))
    cat(sprintf("Median (IQR) Delay [s]         : %.2f (%.2f) s\n", delay_med, delay_iqr))
    cat("\n--- Summary of Delays [seconds] ---\n")
    print(summary(delays_sec))
    
    delay_summary_df <- tibble(
      Cause_Component     = selected_B_col,
      Effect_Component    = "GPS_Horizontal_Error",
      Max_Allowed_Lag_s   = max_lag_seconds,
      N_Pairs             = length(delays_sec),
      Arithmetic_Mean_s   = delay_mean,
      Variance_s2         = delay_var,
      Std_Deviation_s     = delay_sd,
      Min_s               = min(delays_sec),
      Q1_s                = quantile(delays_sec, 0.25),
      Median_s            = delay_med,
      Q3_s                = quantile(delays_sec, 0.75),
      Max_s               = max(delays_sec),
      IQR_s               = delay_iqr
    )
    write_csv(delay_summary_df, file.path(out_dir, paste0("Response_Delay_Summary_Statistics_", epoch_bin_sec, ".csv")))
    
    # Compute experimental probability density function using density()
    if (length(delays_sec) >= 2) {
      delay_density <- density(delays_sec, from = 0, to = max_lag_seconds)
      
      # Save Boxplot and Experimental Probability Density Function (density()) side-by-side
      pdf(file = file.path(out_dir, paste0("Delay_Boxplot_and_Density_", selected_B_col, "_to_GPS_", epoch_bin_sec, ".pdf")),
          width = 11, height = 5.5)
      par(mfrow = c(1, 2), mar = c(5, 5, 4, 2) + 0.1)
      
      # 1. Boxplot of response delays
      boxplot(
        delays_sec,
        main = paste0("Boxplot of GPS Response Delays\n(Lag behind ", selected_B_col, ", max = ", max_lag_seconds, " s)"),
        ylab = "Response Delay [seconds]",
        xlab = paste0(selected_B_col, " -> GPS Horizontal Error"),
        col  = "lightblue",
        border = "darkblue"
      )
      points(1, delay_mean, col = "red", pch = 19, cex = 1.3)
      legend("topright", legend = c("Arithmetic Mean"), col = c("red"), pch = 19, bty = "n")
      
      # 2. Experimental probability density function (density())
      plot(
        delay_density,
        main = paste0("Experimental PDF of Response Delays\n(Mean = ", round(delay_mean, 1),
                      " s, SD = ", round(delay_sd, 1), " s)"),
        xlab = "Response Delay [seconds]",
        ylab = "Density",
        col  = "red",
        lwd  = 2
      )
      polygon(delay_density, col = rgb(1, 0, 0, 0.15), border = "red")
      abline(v = delay_mean, col = "blue", lwd = 2, lty = 2)
      legend("topright", legend = c("Experimental PDF", "Arithmetic Mean"),
             col = c("red", "blue"), lwd = 2, lty = c(1, 2), bty = "n")
      
      par(mfrow = c(1, 1))
      dev.off()
      
      # Also save high-resolution PNG
      png(filename = file.path(out_dir, paste0("Delay_Boxplot_and_Density_", selected_B_col, "_to_GPS_", epoch_bin_sec, ".png")),
          width = 11, height = 5.5, units = "in", res = 300)
      par(mfrow = c(1, 2), mar = c(5, 5, 4, 2) + 0.1)
      boxplot(
        delays_sec,
        main = paste0("Boxplot of GPS Response Delays\n(Lag behind ", selected_B_col, ", max = ", max_lag_seconds, " s)"),
        ylab = "Response Delay [seconds]",
        xlab = paste0(selected_B_col, " -> GPS Horizontal Error"),
        col  = "lightblue",
        border = "darkblue"
      )
      points(1, delay_mean, col = "red", pch = 19, cex = 1.3)
      plot(
        delay_density,
        main = paste0("Experimental PDF of Response Delays\n(Mean = ", round(delay_mean, 1),
                      " s, SD = ", round(delay_sd, 1), " s)"),
        xlab = "Response Delay [seconds]",
        ylab = "Density",
        col  = "red",
        lwd  = 2
      )
      polygon(delay_density, col = rgb(1, 0, 0, 0.15), border = "red")
      abline(v = delay_mean, col = "blue", lwd = 2, lty = 2)
      par(mfrow = c(1, 1))
      dev.off()
    }
  } else {
    cat("No CP pairs found within the specified maximum lag of", max_lag_seconds, "seconds.\n")
  }
  
  # Restore options and close sink
  options(warn = oldw)
  sink()
  closeAllConnections()
}