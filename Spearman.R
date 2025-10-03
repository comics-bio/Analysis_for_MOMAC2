#!/usr/bin/env Rscript
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript script.R <cancer_type> <file1_path> <file2_path>", call. = FALSE)
}

cancer_type <- args[1]
file1_path <- args[2]
file2_path <- args[3]

file1 <- tryCatch(
  read.table(file1_path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE),
  error = function(e) stop("Error reading ", file1_path, ": ", e$message)
)

file2 <- tryCatch(
  read.table(file2_path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE),
  error = function(e) stop("Error reading ", file2_path, ": ", e$message)
)

common_samples <- intersect(rownames(file1), rownames(file2))
common_cols <- intersect(colnames(file1), colnames(file2))

if (length(common_samples) == 0) {
  stop("No common samples found between the two files")
}
if (length(common_cols) == 0) {
  stop("No common species/columns found between the two files")
}
file1_common <- file1[common_samples, common_cols, drop = FALSE]
file2_common <- file2[common_samples, common_cols, drop = FALSE]

results <- lapply(common_cols, function(col) {
  x <- file1_common[[col]]
  y <- file2_common[[col]]
  
  complete_cases <- complete.cases(x, y)
  if (sum(complete_cases) < 3) {
    return(data.frame(
      Cancer_Type = cancer_type,
      Species = col,
      Spearman_rho = NA,
      P_value = NA,
      N_samples = sum(complete_cases),
      stringsAsFactors = FALSE
    ))
  }
  
  test_result <- suppressWarnings(
    cor.test(x, y, method = "spearman", exact = FALSE, use = "complete.obs")
  )
  
  data.frame(
    Cancer_Type = cancer_type,
    Species = col,
    Spearman_rho = test_result$estimate,
    P_value = test_result$p.value,
    N_samples = sum(complete_cases),
    stringsAsFactors = FALSE
  )
})

result_df <- bind_rows(results) %>%
  mutate(
    P_adj_BH = p.adjust(P_value, method = "BH"),
    Significance = case_when(
      is.na(P_adj_BH) ~ "NA",
      P_adj_BH < 0.001 ~ "***",
      P_adj_BH < 0.01 ~ "**",
      P_adj_BH < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(P_adj_BH)

output_file <- paste0("spearman_results_", cancer_type, "_", format(Sys.Date(), "%Y%m%d"), ".tsv")
write.table(result_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)