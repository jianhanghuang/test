#!/usr/bin/env python
import pandas as pd

def calculate_frequency(file_path, output_result_path):
    #file_path：上一步计算出的密码子和氨基酸的数目表
    #output_result_path：结果输出路径

    # 读取数据文件
    df = pd.read_csv(file_path)

    # 根据密码子分组，计算密码子的总数
    total_codon_counts = df['codon_counts'].sum()
    # 计算每个 protein1 的出现总次数
    total_protein_counts = df.groupby('protein1')['codon_counts'].sum().reset_index()
    total_protein_counts = total_protein_counts.rename(columns={'codon_counts': 'protein1_counts'})

    #合并数据
    merged_df = pd.merge(df, total_protein_counts, on='protein1')

    # 计算每个密码子的频率
    merged_df['codon_frequency'] = merged_df['codon_counts'] / total_codon_counts
    #计算同义密码子频率
    merged_df['protein1_frequency'] = merged_df['codon_counts'] / merged_df['protein1_counts']
    #计算每个氨基酸的频率
    merged_df['acid_frequency'] = merged_df['protein1_counts'] / total_codon_counts
    #输出
    merged_df.to_csv(output_result_path, index=False)
    return merged_df

calculate_frequency('../Result/Counts/Arabidopsis_counts.csv', '../Result/Frequency/Arabidopsis_freq.csv')
calculate_frequency('../Result/Counts/C_elegants_counts.csv', '../Result/FrequencyC_elegants_freq.csv')
calculate_frequency('../Result/Counts/Drosophila_counts.csv', '../Result/FrequencyDrosophila_freq.csv')
calculate_frequency('../Result/Counts/Human_counts.csv', '../Result/FrequencyHuman_freq.csv')
calculate_frequency('../Result/Counts/Mouse_counts.csv', '../Result/FrequencyMouse_freq.csv')
calculate_frequency('../Result/Counts/Yeast_counts.csv', '../Result/FrequencyYeast_freq.csv')


