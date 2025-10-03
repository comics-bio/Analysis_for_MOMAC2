#install.packages("foreach")  
#install.packages("doParallel")
library(data.table)  
library(foreach)  

args <- commandArgs(trailingOnly = TRUE)  
project_code <- args[1]  
df <- read.table(args[2], sep = '\t', header = TRUE, check.names = FALSE)  
cell_types <- read.table(args[3], header = FALSE, sep = "\t", stringsAsFactors = FALSE)[, 1]
final_path_out <- args[4]  # 输出文件路径 



find_nonzero_columns <- function(data) {
  # 检查每列是否包含"Nonzero"值
  is_nonzero_col <- sapply(data, function(col) {
    if(is.character(col) | is.factor(col)) {
      any(col == "Nonzero", na.rm = TRUE)
    } else {
      FALSE
    }
  })
  
  # 返回列名
  names(data)[is_nonzero_col]
}

groups <- find_nonzero_columns(df)

  
results <- data.frame(  
  CellType = character(),  
  Group = character(),  
  AvgZero = numeric(),  
  AvgNonzero = numeric(),  
  MedZero = numeric(),  
  MedNonzero = numeric(),  
  PValue = numeric(),
  FoldChange = numeric(),  # 新增FoldChange列
  stringsAsFactors = FALSE  
)  

# 修改后的处理函数，添加fold change计算
process_and_cell_type <- function(group, cell_type) {  
  selected_A_values <- df[[cell_type]][df[[group]] == "Zero"]  
  selected_B_values <- df[[cell_type]][df[[group]] == "Nonzero"]  
    
  if (length(selected_A_values) >= 3 && length(selected_B_values) >= 3) {  
    test_result <- wilcox.test(selected_A_values, selected_B_values)  
    p_value <- test_result$p.value  
      
    avg_zero <- mean(selected_A_values, na.rm = TRUE)  
    avg_nonzero <- mean(selected_B_values, na.rm = TRUE)  
    med_zero <- median(selected_A_values, na.rm = TRUE)  
    med_nonzero <- median(selected_B_values, na.rm = TRUE)  
    
    # 计算FoldChange (Nonzero/Zero)，处理除零情况
    if (avg_zero == 0) {
      if (avg_nonzero == 0) {
        fold_change <- 1  # 两者都为0时设为1
      } else {
        fold_change <- Inf  # Zero为0，Nonzero不为0
      }
    } else {
      fold_change <- avg_nonzero / avg_zero
    }
      
    return(data.frame(  
      CellType = cell_type,  
      Group = group,  
      AvgZero = avg_zero,  
      AvgNonzero = avg_nonzero,  
      MedZero = med_zero,  
      MedNonzero = med_nonzero,  
      PValue = p_value,
      FoldChange = fold_change,  # 添加fold change
      stringsAsFactors = FALSE  
    ))  
  } else {  
    return(data.frame(  
      CellType = cell_type,  
      Group = group,  
      AvgZero = NA,  
      AvgNonzero = NA,  
      MedZero = NA,  
      MedNonzero = NA,  
      PValue = NA,
      FoldChange = NA,  # 添加fold change
      stringsAsFactors = FALSE  
    ))
  }  
}  

# 使用foreach并行处理每个组和细胞类型的组合  
combined_results <- foreach(group = groups, .combine = rbind, .inorder = TRUE, .packages = c("stats"))%do%{  
  cell_type_results <- foreach(cell_type = cell_types, .combine = rbind, .inorder = TRUE) %do% {  
    process_and_cell_type(group, cell_type)  
  }  
  return(cell_type_results)  
}  
 

# 合并结果并计算调整后的p值  
results <- combined_results  
results$AdjustedPValue <- p.adjust(results$PValue, method = "BH")  

# 根据FoldChange确定富集组
results$enriched <- ifelse(results$FoldChange > 1, "Nonzero", 
                                ifelse(results$FoldChange < 1, "Zero", "Equal"))

# 写入结果文件  
write.table(results, file = final_path_out, sep = "\t", row.names = FALSE, quote = FALSE)