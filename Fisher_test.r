#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 3) {
  stop("Usage: Rscript column_fisher_test.R <project_name> <file1.tsv> <file2.tsv>")
}

project <- args[1]
file1 <- args[2]
file2 <- args[3]

file1_name <- basename(file1)

read_data <- function(file) {
  as.matrix(read.delim(file, row.names = 1, check.names = FALSE))
}

data1 <- read_data(file1)
data2 <- read_data(file2)

common_microbes <- intersect(colnames(data1), colnames(data2))
common_samples <- intersect(rownames(data1), rownames(data2))

cat(paste("Common samples:", length(common_samples), "\n"))
cat(paste("Common microbes:", length(common_microbes), "\n"))

results <- data.frame(
  Microbe = character(),
  File1_Name = character(),    
  Co_occurrence = numeric(),   
  Only_File1 = numeric(),      
  Only_File2 = numeric(),      
  Co_absence = numeric(),      
  Odds_Ratio = numeric(),
  P_value = numeric(),
  Significant_Cooccurrence = character(),  
  stringsAsFactors = FALSE
)

for(microbe in common_microbes) {
  vec1 <- ifelse(data1[common_samples, microbe] > 0, 1, 0)
  vec2 <- ifelse(data2[common_samples, microbe] > 0, 1, 0)
  
  a <- sum(vec1 == 1 & vec2 == 1)  
  b <- sum(vec1 == 1 & vec2 == 0)  
  c <- sum(vec1 == 0 & vec2 == 1)  
  d <- sum(vec1 == 0 & vec2 == 0)  
  cont_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  
  if(sum(cont_table) < 1) next
  
  test <- try(fisher.test(cont_table), silent = TRUE)
  
  if(!inherits(test, "try-error")) {
    odds_ratio <- ifelse(is.null(test$estimate), NA, test$estimate)
    p_value <- test$p.value
    
    is_significant <- ifelse(p_value < 0.05 & odds_ratio > 1, "Yes", "No")
    
    results[nrow(results)+1, ] <- list(
      Microbe = microbe,
      File1_Name = file1_name,  
      Co_occurrence = a,
      Only_File1 = b,
      Only_File2 = c,
      Co_absence = d,
      Odds_Ratio = odds_ratio,
      P_value = p_value,
      Significant_Cooccurrence = is_significant
    )
  }
}


results$Adj_P <- p.adjust(results$P_value, method = "BH")


results <- results[order(results$Adj_P), ]

output_file <- paste0(project, "_cooccurrence_results.tsv")
write.table(results, output_file, sep = "\t", quote = FALSE, row.names = FALSE)