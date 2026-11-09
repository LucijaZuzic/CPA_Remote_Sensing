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
  {     this_file<-rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}

# Set path to the working directory

setwd(getCurrentFileLocation())

format_latex_scientific <- function(val, threshold = 2.2e-16) {
  val <- as.numeric(val)
  if (is.na(val)) return("NA")
  
  if (val < threshold) {
    # Handles extremely small values
    return("$< 2.2 \\times 10^{-16}$")
  } else if (val < 0.001) {
    # Convert standard R scientific notation (e.g., 1.45e-04) to LaTeX math mode
    sci <- sprintf("%.2e", val)
    parts <- unlist(strsplit(sci, "e"))
    base <- parts[1]
    exp <- as.integer(parts[2]) # Strips unnecessary leading zeros in the exponent
    return(paste0("$", base, " \\times 10^{", exp, "}$"))
  } else {
    # Standard decimal formatting for normal values with thousands separator
    formatted_num <- formatC(val, format = "f", digits = 4, big.mark = ",")
    return(paste0("$", formatted_num, "$"))
  }
}

df_to_custom_tabularx <- function(df, caption = "Table Caption", label = "tab:my_label", 
                                  digits = 4, col_formats = NULL, file_out = NULL) {
  
  nc <- ncol(df)
  col_spec <- paste0(rep("C", nc), collapse = "")
  headers <- paste0("\\textbf{", str_replace_all(colnames(df), "_", " "), "}", collapse = " & ")
  
  # Process each row cell-by-cell
  rows <- apply(df, 1, function(row) {
    formatted_cells <- sapply(seq_along(row), function(i) {
      col_name <- colnames(df)[i]
      val <- row[i]
      
      # 1. Check if a custom format rule is defined for this specific column
      if (!is.null(col_formats) && col_name %in% names(col_formats)) {
        rule <- col_formats[[col_name]]
        
        if (is.function(rule)) {
          return(rule(val)) # Custom function returning LaTeX string
        } else if (is.character(rule)) {
          if (rule == "string") {
            return(str_replace_all(val, "_", " "))
          } else {
            if (rule == "integer") {
              formatted_num <- format(as.integer(val), big.mark = ",", scientific = FALSE)
              return(paste0("$", formatted_num, "$"))
            } else {
              # Treat as a sprintf format string like "%.2f" or "%.4f", extracting decimal digits and applying big.mark
              digits_val <- suppressWarnings(as.integer(sub("%.([0-9]+)f.*", "\\1", rule)))
              if (is.na(digits_val)) digits_val <- digits
              formatted_num <- formatC(as.numeric(val), format = "f", digits = digits_val, big.mark = ",")
              return(paste0("$", formatted_num, "$"))
            }
          }
        }
      }
      
      # 2. Default global formatting fallback with thousands separator
      if (is.numeric(val)) {
        if (val == round(val) && abs(val) < 1e6) {
          formatted_num <- format(as.integer(val), big.mark = ",")
          paste0("$", formatted_num, "$")
        } else {
          formatted_num <- formatC(val, format = "f", digits = digits, big.mark = ",")
          paste0("$", formatted_num, "$")
        }
      } else {
        as.character(val)
      }
    })
    paste0(paste0(formatted_cells, collapse = " & "), " \\\\")
  })
  row_string <- paste0(rows, collapse = "\n")
  
  # 2. Assemble the complete LaTeX block
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
    cat(sprintf("\n\nLaTeX table successfully written to: %s\n", file_out))
  }
  
  return(latex_output)
}

epoch_bin_sec <- 300
Response_Delay_Summary_Statistics_df_merged <- c()
for (selected_B_col in c("Bx", "By", "Bz", "TEC", "Dst")) {
  out_dir <- paste0("CPD_Darwin_2014_", selected_B_col, "_", epoch_bin_sec)
  Response_Delay_Summary_Statistics_df <- read_csv(file.path(out_dir, paste0("Response_Delay_Summary_Statistics_", epoch_bin_sec, ".csv")), show_col_types = FALSE)
  Response_Delay_Summary_Statistics_df <- Response_Delay_Summary_Statistics_df[, !(colnames(Response_Delay_Summary_Statistics_df) %in% c("Effect_Component", "Max_Allowed_Lag_s"))]
  Response_Delay_Summary_Statistics_df_merged <- rbind(Response_Delay_Summary_Statistics_df_merged, Response_Delay_Summary_Statistics_df)
  # Generate table using the custom tabularx function
  df_to_custom_tabularx(
    df = Response_Delay_Summary_Statistics_df,
    caption = paste0("Statistical analysis of response delays (GPS CP lag behind ", selected_B_col, " CP)."),
    label = paste0("tab:Response_Delay_Summary_Statistics_", epoch_bin_sec, "_", selected_B_col),
    col_formats = list(
      "Cause_Component" = "string",
      "N_Pairs" = "integer",
      "Arithmetic_Mean_s" = "%.2f",
      "Variance_s2" = "%.2f",
      "Std_Deviation_s" = "%.2f",
      "Min_s" = "%.2f",
      "Q1_s" = "%.2f",
      "Median_s" = "%.2f",
      "Q3_s" = "%.2f",
      "Max_s" = "%.2f",
      "IQR_s" = "%.2f"
    ),
    file_out = file.path(out_dir, paste0("Response_Delay_Summary_Statistics_", epoch_bin_sec, ".tex"))
  )
}

df_to_custom_tabularx(
  df = Response_Delay_Summary_Statistics_df_merged,
  caption = paste0("Statistical analysis of response delays (GPS CP lag behind cause component CP)."),
  label = paste0("tab:Response_Delay_Summary_Statistics_", epoch_bin_sec),
  col_formats = list(
    "Cause_Component" = "string",
    "N_Pairs" = "integer",
    "Arithmetic_Mean_s" = "%.2f",
    "Variance_s2" = "%.2f",
    "Std_Deviation_s" = "%.2f",
    "Min_s" = "%.2f",
    "Q1_s" = "%.2f",
    "Median_s" = "%.2f",
    "Q3_s" = "%.2f",
    "Max_s" = "%.2f",
    "IQR_s" = "%.2f"
  ),
  file_out = "Response_Delay_Summary_Statistics.tex"
)

df1 <- read_csv(file.path(out_dir, paste0("Descriptive_Statistics_GPS_Bxyz.csv")), show_col_types = FALSE)
df1 <- df1[, !(colnames(df1) %in% c("N"))]

df_to_custom_tabularx(
  df = df1,
  caption = paste0("Descriptive statistics for the GPS horizontal positioning error, Bx, By, and Bz."),
  label = paste0("tab:Descriptive_Statistics_GPS_Bxyz"),
  col_formats = list(
    "Variable" = "string",
    "Min" = "%.2f",
    "Q1" = "%.2f",
    "Median" = "%.2f",
    "Mean" = "%.2f",
    "Q3" = "%.2f",
    "Max" = "%.2f",
    "Variance" = "%.2f",
    "SD" = "%.2f",
    "IQR" = "%.2f",
    "MAD" = "%.2f",
    "Skewness" = "%.2f",
    "Kurtosis" = "%.2f"
  ),
  file_out = "Descriptive_Statistics_GPS_Bxyz.tex"
)