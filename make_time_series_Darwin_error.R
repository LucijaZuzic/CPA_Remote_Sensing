# Clear environment
rm(list=ls())

## The R script for Žužić et al. manuscript

# Set path to the working directory containing prepared data set

library(tidyverse)
getCurrentFileLocation <-  function()
{
  this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col = value, into = c("key", "value"), sep = "=", fill = 'right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file) == 0) {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}

setwd(getCurrentFileLocation())

library(dplyr)
library(tidyr)
library(lubridate)
library(geosphere)
library(ggplot2)
library(readr)
library(stringr)
library(extrafont) # Added extrafont

# (Optional: run loadfonts() if your PDF output isn't rendering the font correctly)
loadfonts(device = "pdf", quiet = TRUE)

# 1. Initialization and reference coordinates

ref_lons <- c(131.132744)
ref_lats <- c(-12.843697)
ref_heights <- c(125.1)

# Shared theme
bigger_theme <- theme_classic(base_size = 20, base_family = "DejaVu Sans") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 20, color = "black"),       # Pure black text
    axis.ticks = element_line(color = "black", linewidth = 1),  # Pure black tick marks
    axis.line = element_line(color = "black", linewidth = 1),   # Pure black axis lines
    legend.text = element_text(size = 20, color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

# 2. File parsing and distance computation (vectorized)
num_days <- 365
stormnum <- 1
station_names <- c("Darwin")
station_abbrevs <- c("darw")
station_longs <- c("Darwin, NT, Australia")
numcols <- c(8)
makeplot <- TRUE
for (numloc in 1:3) {
  full_data <- list()
  case_data <- list()
  date_data <- list()
  ref_coords <- list(
    lon = ref_lons[numloc],
    lat = ref_lats[numloc],
    height = ref_heights[numloc]
  )
  station_name <- station_names[numloc]
  station_abbrev <- station_abbrevs[numloc]
  station_long <- station_longs[numloc]
  numcol <- numcols[numloc]
  storm_dir <- paste0(station_name,"_2014_POS")
  case_data[[as.character(stormnum)]] <- list()
  date_data[[as.character(stormnum)]] <- list()
  
  for (shortfile in 1:num_days) {
    file_path <- paste0(storm_dir,"/",station_abbrev,shortfile,"0.pos")
    
    if (file.exists(file_path)) {
      filename <- basename(file_path)
      
      # Read the file header to extract column names
      raw_lines <- readLines(file_path, n = 20)
      header_line <- raw_lines[numcol] # skip = numcol
      
      # 1. Clean the header and remove the "%" and "GPST" keywords
      header_line <- gsub("%|GPST", "", header_line) 
      # 2. Trim excess spaces
      header_line <- trimws(gsub("\\s+", " ", header_line)) 
      
      # 3. Split into remaining column names and drop any accidental empty strings
      data_cols <- unlist(strsplit(header_line, " "))
      data_cols <- data_cols[data_cols != ""] 
      
      # 4. Prepend date and time to perfectly match the 15 data columns
      final_col_names <- c("date", "time", data_cols)
      
      # Read the data block using the perfectly aligned names
      df_raw <- read_table(file_path, skip = 15, col_names = final_col_names, show_col_types = FALSE)
      
      if (nrow(df_raw) == 0) next
      
      mvals_full <- c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365)
      mvals <- c()
      for (m in mvals_full) {
        if (m <= num_days) {
          mvals <- c(mvals, m)
        }
      }
  
      # Vectorized timestamp parsing
      df <- df_raw %>%
        mutate(
          datetime = ymd_hms(paste(gsub("/", "-", date), time)),
          d = day(datetime),
          m = month(datetime),
          y = year(datetime),
          hr = hour(datetime),
          min = minute(datetime),
          sec = second(datetime),
          doy = mvals[m] + d - 1,
          `total seconds day` = hr * 3600 + min * 60 + sec,
          `total seconds` = doy * 24 * 3600 + `total seconds day`
        )
      
      # Vectorized geodesic distance calculations
      df <- df %>%
        mutate(
          # Point 1 is a simple vector, Point 2 uses cbind() to create an N x 2 matrix
          dy = distGeo(c(ref_coords$lon, ref_coords$lat), 
                        cbind(ref_coords$lon, `latitude(deg)`)),
          
          dx = distGeo(c(ref_coords$lon, ref_coords$lat), 
                        cbind(`longitude(deg)`, ref_coords$lat)),
          
          `absolute height(m)` = abs(ref_coords$height - `height(m)`),
          `absolute latitude(deg)` = abs(dy),
          `absolute longitude(deg)` = abs(dx),
          
          `horizontal(deg)` = distGeo(c(ref_coords$lon, ref_coords$lat), 
                                      cbind(`longitude(deg)`, `latitude(deg)`)),
          
          ang = atan2(dy, dx) / pi * 180,
          `horizontal angle(deg)` = ang,
          `absolute horizontal angle(deg)` = abs(ang),
          `minimum horizontal angle(deg)` = ifelse(abs(ang) > 90, 180 - abs(ang), abs(ang))
        ) %>%
        select(-date, -time, -datetime, -dy, -dx, -ang)
      
      # Store in memory hierarchies
      full_data[[length(full_data) + 1]] <- df
      case_data[[as.character(stormnum)]][[length(case_data[[as.character(stormnum)]]) + 1]] <- df
      date_data[[as.character(stormnum)]][[shortfile]] <- df
      
      print(paste(stormnum, filename, ncol(df), nrow(df)))
    } else {
      print(paste0("Missing",file_path))
    }
  }
  
  # Combine lists into dataframes
  full_df <- bind_rows(full_data)
  write_csv(full_df, paste0("afullframe_",station_name,"_new_",num_days,".csv"))
  
  if ("GPST" %in% names(full_df)) {
    print(length(full_df$GPST))
  }
  
  rdim <- list(hr = 0, len = 0, Q1 = 3, Q2 = 3, avg = 2, Q3 = 2, max = 3)
  nks <- c("hr", "len", "Q1", "Q2", "avg", "Q3", "max")
  
  # 3. Storm day generation (statistics and plots)
  for (stormnum in names(case_data)) {
    storm_df <- bind_rows(case_data[[stormnum]])
    
    # Plotting variables
    maxofs <- num_days

    # Make plot (year)
    if (makeplot) {
      secticks_full <- c(0, 31 * 24 * 3600, 59 * 24 * 3600, 90 * 24 * 3600, 120 * 24 * 3600, 151 * 24 * 3600, 181 * 24 * 3600, 212 * 24 * 3600, 243 * 24 * 3600, 273 * 24 * 3600, 304 * 24 * 3600, 334 * 24 * 3600)
      seclabs_full <- c("1", "32", "60", "91", "121", "152", "182", "213", "244", "274", "305", "335")
      secticks <- c()
      seclabs <- c()
      for (ml in 1:12) {
        if (secticks_full[ml] <= maxofs * 24 * 3600) {
          secticks <- c(secticks, secticks_full[ml])
          seclabs <- c(seclabs, seclabs_full[ml])
        }
      }
      
      filtval <- 10
      tffilt <- storm_df[storm_df$`horizontal(deg)` < filtval,]
      cat(length(storm_df$`horizontal(deg)`), length(tffilt$`horizontal(deg)`), "\n")
      
      plot_title <- paste0("Time series of horizontal positioning errors [m]\nfor ",station_long," (2014)")
      
      p <- ggplot(storm_df, aes(x = `total seconds`, y = `horizontal(deg)`)) +
        geom_line(color = "blue", linewidth = 1) +
        scale_x_continuous(breaks = secticks, labels = seclabs, limits = c(0, min(c((maxofs + 1) * 24 * 3600, max(secticks) + 31 * 24 * 3600))), expand = c(0, 0)) +
        labs(title = plot_title, x = "DOY", y = "Horizontal\npositioning\nerrors [m]") +
        bigger_theme
      
      ggsave(paste0("afullframe_time_",station_name,"_new_",num_days,"_R.pdf"), plot = p, width = 10, height = 4.5, device = cairo_pdf)
      ggsave(paste0("afullframe_time_",station_name,"_new_",num_days,"_R.png"), plot = p, width = 10, height = 4.5)
      
      write_csv(tffilt, paste0("afullframe_filter_",station_name,"_new_",num_days,".csv"))
      
      dff <- read_csv(paste0("afullframe_filter_",station_name,"_new_",num_days,".csv"), show_col_types = FALSE)
      
      plot_title <- paste0("Time series of filtered horizontal positioning errors [m]\n(< ", filtval, ") for ",station_long," (2014)")
      
      p <- ggplot(dff, aes(x = `total seconds`, y = `horizontal(deg)`)) +
        geom_line(color = "blue", linewidth = 1) +
        scale_x_continuous(breaks = secticks, labels = seclabs, limits = c(0, min(c((maxofs + 1) * 24 * 3600, max(secticks) + 31 * 24 * 3600))), expand = c(0, 0)) +
        labs(title = plot_title, x = "DOY", y = "Horizontal\npositioning\nerrors [m]") +
        bigger_theme
      
      ggsave(paste0("afullframe_filter_",station_name,"_new_",num_days,"_R.pdf"), plot = p, width = 10, height = 4.5, device = cairo_pdf)
      ggsave(paste0("afullframe_filter_",station_name,"_new_",num_days,"_R.png"), plot = p, width = 10, height = 4.5)
    }
  }
}