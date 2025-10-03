import pandas as pd
import numpy as np
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Kaplan-Meier生存分析脚本')
    parser.add_argument('--input_data', required=True, help='输入数据文件路径')
    parser.add_argument('--cancer_list', required=True, help='癌症类型列表文件路径')
    parser.add_argument('--taxa_list', required=True, help='微生物列表文件路径')
    parser.add_argument('--output_dir', required=True, help='输出目录路径')
    parser.add_argument('--os_output', default='OS_result.txt', help='OS结果文件名')
    parser.add_argument('--pfs_output', default='PFS_result.txt', help='PFS结果文件名')
    
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    OS_data = pd.read_table(args.input_data)

    OS_result = open(os.path.join(args.output_dir, args.os_output), "w")
    PFS_result = open(os.path.join(args.output_dir, args.pfs_output), "w")

    with open(args.cancer_list, "r") as cancer_list:
        for cancer_name in cancer_list:
            cancer_name = cancer_name.rstrip("\n")
            cancer_data = OS_data.loc[OS_data["cancer"] == cancer_name]
            
            with open(args.taxa_list, "r") as taxa_list:
                for taxa_name in taxa_list:
                    taxa_name = taxa_name.rstrip("\n")
                    cancer_data[taxa_name + "_group"] = np.where(
                        cancer_data[taxa_name] == 0, "Zero", "Nonzero"
                    )
                    
                    enriched_group_os = "Undetermined"
                    enriched_group_pfs = "Undetermined"
                    
                    OS_T = cancer_data["OS_time"]
                    OS_E = cancer_data["OS"]
                    PFS_T = cancer_data["PFI_time"]
                    PFS_E = cancer_data["PFI"]
                    
                    dem = (cancer_data[taxa_name + "_group"] == "Nonzero")
                    
                    if len(dem[dem]) > 0 and len(dem[~dem]) > 0:
                        kmf = KaplanMeierFitter()
                        
                        kmf.fit(OS_T[dem], event_observed=OS_E[dem])
                        median_nonzero = kmf.median_survival_time_
                        
                        kmf.fit(OS_T[~dem], event_observed=OS_E[~dem])
                        median_zero = kmf.median_survival_time_
                        
                        if not np.isnan(median_nonzero) and not np.isnan(median_zero):
                            if median_nonzero > median_zero:
                                enriched_group_os = "Nonzero"
                            else:
                                enriched_group_os = "Zero"
                        elif not np.isnan(median_nonzero):
                            enriched_group_os = "Nonzero"
                        elif not np.isnan(median_zero):
                            enriched_group_os = "Zero"
                       
                        kmf.fit(PFS_T[dem], event_observed=PFS_E[dem])
                        median_nonzero_pfs = kmf.median_survival_time_
                        
                        kmf.fit(PFS_T[~dem], event_observed=PFS_E[~dem])
                        median_zero_pfs = kmf.median_survival_time_
                        
                        if not np.isnan(median_nonzero_pfs) and not np.isnan(median_zero_pfs):
                            if median_nonzero_pfs > median_zero_pfs:
                                enriched_group_pfs = "Nonzero"
                            else:
                                enriched_group_pfs = "Zero"
                        elif not np.isnan(median_nonzero_pfs):
                            enriched_group_pfs = "Nonzero"
                        elif not np.isnan(median_zero_pfs):
                            enriched_group_pfs = "Zero"
                    os_survival = logrank_test(
                        OS_T[dem], OS_T[~dem], 
                        event_observed_A=OS_E[dem], 
                        event_observed_B=OS_E[~dem]
                    )
                    
                    pfs_survival = logrank_test(
                        PFS_T[dem], PFS_T[~dem],
                        event_observed_A=PFS_E[dem],
                        event_observed_B=PFS_E[~dem]
                    )
                    
                    print(f"{cancer_name},{taxa_name},OS_p={os_survival.p_value},Enriched={enriched_group_os}", file=OS_result)
                    print(f"{cancer_name},{taxa_name},PFS_p={pfs_survival.p_value},Enriched={enriched_group_pfs}", file=PFS_result)

                    group_output = os.path.join(args.output_dir, f"{cancer_name}_group.tsv")
                    cancer_data.to_csv(group_output, sep='\t', index=False)
    OS_result.close()
    PFS_result.close()

if __name__ == "__main__":
    main()