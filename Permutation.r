suppressPackageStartupMessages({
  library(tidyverse)
})


generate_derangement <- function(n) {
  if (n == 1) return(1)  
  
  while (TRUE) {
    perm <- sample(n)
    if (!any(perm == 1:n)) {  
      return(perm)
    }
  }
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript shuffle_data_derangement.R <cancer_type> <input_file>", call. = FALSE)
}

cancer_type <- args[1]
input_file <- args[2]

original_data <- read_tsv(input_file, col_types = cols(.default = "c"))
cols_to_shuffle <- 2:min(269, ncol(original_data))

for (i in 1:100) {
  shuffled_data <- original_data
  
  shuffled_data[, cols_to_shuffle] <- lapply(
    original_data[, cols_to_shuffle],
    function(col) {
      col[generate_derangement(length(col))]
    }
  )
  
  output_file <- paste0(cancer_type, "_", i, ".tsv")
  write_tsv(shuffled_data, output_file)
}