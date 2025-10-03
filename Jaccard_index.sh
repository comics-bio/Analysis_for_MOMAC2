#!/bin/bash

cancer_list="24-cancer.txt"
mkdir -p comparison_results

cat "$cancer_list" | while read cancer; do
    result_file="comparison_results/${cancer}_all_comparisons.txt"
    > "$result_file"
    
    # 获取rob和sal文件列表
    rob_files=($(ls /${cancer}-*-filter.tsv 2>/dev/null))
    sal_files=($(ls /${cancer}-*-filter.tsv 2>/dev/null))
    
    if [ ${#rob_files[@]} -eq 0 ] || [ ${#sal_files[@]} -eq 0 ]; then
        echo "Warning: No valid file pairs found for $cancer" >&2
        continue
    fi
    
    for rob_file in "${rob_files[@]}"; do
        rob_base=$(basename "$rob_file" -filter.tsv)
        awk -F'\t' '{print $1"-"$2"-"$10}' "$rob_file" > "tmp_rob.asso"
        sort "tmp_rob.asso" > "tmp_rob.sorted"
        
        for sal_file in "${sal_files[@]}"; do
            sal_base=$(basename "$sal_file" -filter.tsv)
            awk -F'\t' '{print $1"-"$2"-"$10}' "$sal_file" > "tmp_sal.asso"
            sort "tmp_sal.asso" > "tmp_sal.sorted"
            
            common_lines=$(comm -12 "tmp_rob.sorted" "tmp_sal.sorted")
            
            if [ -n "$common_lines" ]; then
                echo "=== ${rob_base} vs ${sal_base} ===" >> "$result_file"
                echo "$common_lines" | awk -v rob="$rob_base" -v sal="$sal_base" '{
                    print rob"|"sal"|"$0  # 输出格式：rob文件|sal文件|共同内容
                }' >> "$result_file"
                echo "" >> "$result_file"
            fi
            
            rm "tmp_sal.asso" "tmp_sal.sorted"
        done
        
        rm "tmp_rob.asso" "tmp_rob.sorted"
    done
    
    total_common=$(grep -v '^===.*===$' "$result_file" | grep -v '^$' | wc -l)
    echo "${cancer}: Found $total_common common associations" >> "comparison_results/summary.txt"
done

echo "Processing complete. Results saved in comparison_results directory."