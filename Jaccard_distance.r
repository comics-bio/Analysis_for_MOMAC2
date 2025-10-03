#!/usr/bin/env Rscript
library(tidyverse)
library(vegan)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript cross_jaccard.R <project_name> <table1.tsv> <table2.tsv>\nTables format: Rows=samples, Columns=microbes", call. = FALSE)
}

project_name <- args[1]
table1_file <- args[2]
table2_file <- args[3]

read_microbe_table <- function(file) {
  read_tsv(file, col_types = cols()) %>%
    column_to_rownames(var = names(.[1])) %>%  
    t() %>%                                   
    {ifelse(. > 0, 1, 0)}                    
}


table1 <- read_microbe_table(table1_file)
table2 <- read_microbe_table(table2_file)

common_microbes <- intersect(rownames(table1), rownames(table2))
if (length(common_microbes) == 0) stop("No common microbes between tables!")


table1_common <- table1[common_microbes, , drop = FALSE]
table2_common <- table2[common_microbes, , drop = FALSE]

calculate_cross_jaccard <- function(mat1, mat2) {

  combined <- cbind(mat1, mat2)
  

  dist_mat <- vegdist(t(combined), method = "jaccard", binary = TRUE) %>%
    as.matrix()
  

  table1_samples <- colnames(mat1)
  table2_samples <- colnames(mat2)
  

  expand_grid(
    Table1_Sample = table1_samples,
    Table2_Sample = table2_samples
  ) %>%
    mutate(
      Jaccard_Distance = map2_dbl(
        Table1_Sample, 
        Table2_Sample,
        ~dist_mat[.x, .y]
      )
    )
}

result <- calculate_cross_jaccard(table1_common, table2_common)


output_file <- paste0(project_name, "_cross_sample_jaccard.tsv")
write_tsv(result, output_file)
message(paste("Results saved to:", output_file))
message(paste("Total pairs:", nrow(result)))