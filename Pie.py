import matplotlib.pyplot as plt
import pandas as pd


def draw_codon_pie(file_path):
    # 读取数据文件
    df1 = pd.read_csv(file_path)
    df1 = df1.sort_values(by='codon_frequency', ascending=False)

    # 绘制密码子的饼图
    # 选取 codon_frequency 大于2%的数据
    selected_df = df1[df1['codon_frequency'] > 0.02]

    # 计算其他的频率总和，用于生成“其他”类别
    other_frequency = df1[df1['codon_frequency'] <= 0.02]['codon_frequency'].sum()

    # 将“其他”类别加入到选取的数据中
    selected_df = selected_df.append({'codon': 'Other', 'codon_frequency': other_frequency}, ignore_index=True)
    return selected_df


def draw_protein_pie(file_path):
    df2 = pd.read_csv(file_path)
    df2 = df2.sort_values(by='acid_frequency', ascending=False)
    df_cleaned2 = df2.drop_duplicates(subset='protein3', keep='first')
    return df_cleaned2

# 绘制饼图
#拟南芥
plt.rcParams['font.size'] = 12
plt.figure(figsize=(5, 5))
df_Arabidopsis = draw_codon_pie('../../Result/Arabidopsis_freq_拟南芥.csv')
plt.pie(df_Arabidopsis['codon_frequency'], labels=df_Arabidopsis['codon'], autopct='%1.1f%%', startangle=90, counterclock=False, textprops={'fontsize': 10})
plt.title('Arabidopsis Codon Frequency Pie Chart')
plt.savefig('Arabidopsis_Codon_frequency_pie_chart.pdf')

plt.rcParams['font.size'] = 12
plt.figure(figsize=(5, 5))
df_cleaned_Arabidopsis = draw_protein_pie('../../Result/Arabidopsis_freq_拟南芥.csv')
plt.pie(df_cleaned_Arabidopsis['acid_frequency'], labels=df_cleaned_Arabidopsis['protein3'], autopct='%1.1f%%', startangle=90, counterclock=False, textprops={'fontsize': 10})
plt.title('Arabidopsis Protein Frequency Pie Chart')
plt.savefig('Arabidopsis_Protein_frequency_pie_chart.pdf')

#线虫
plt.rcParams['font.size'] = 12
plt.figure(figsize=(5, 5))
df_elegants = draw_codon_pie('../../Result/C_elegants_freq_线虫.csv')
plt.pie(df_elegants['codon_frequency'], labels=df_elegants['codon'], autopct='%1.1f%%', startangle=90, counterclock=False, textprops={'fontsize': 10})
plt.title('C_elegants Codon Frequency Pie Chart')
plt.savefig('C_elegants_Codon_frequency_pie_chart.pdf')

plt.rcParams['font.size'] = 12
plt.figure(figsize=(5, 5))
df_cleaned_elegants = draw_protein_pie('../../Result/C_elegants_freq_线虫.csv')
plt.pie(df_cleaned_elegants['acid_frequency'], labels=df_cleaned_elegants['protein3'], autopct='%1.1f%%', startangle=90, counterclock=False, textprops={'fontsize': 10})
plt.title('C_elegants Protein Frequency Pie Chart')
plt.savefig('C_elegants_Protein_frequency_pie_chart.pdf')

#果蝇
plt.rcParams['font.size'] = 12
plt.figure(figsize=(5, 5))
df_Drosophila = draw_codon_pie('../Result/Frequency/Drosophila_freq.csv')
plt.pie(df_Drosophila['codon_frequency'], labels=df_Drosophila['codon'], autopct='%1.1f%%', startangle=90, counterclock=False, textprops={'fontsize': 10})
plt.title('Drosophila Codon Frequency Pie Chart')
plt.savefig('Drosophila_Codon_frequency_pie_chart.pdf')

plt.rcParams['font.size'] = 12
plt.figure(figsize=(5, 5))
df_cleaned_Drosophila = draw_protein_pie('../Result/Frequency/Drosophila_freq.csv')
plt.pie(df_cleaned_Drosophila['acid_frequency'], labels=df_cleaned_Drosophila['protein3'], autopct='%1.1f%%', startangle=90, counterclock=False, textprops={'fontsize': 10})
plt.title('Drosophila Protein Frequency Pie Chart')
plt.savefig('Drosophila_Protein_frequency_pie_chart.pdf')

#人
plt.rcParams['font.size'] = 12
plt.figure(figsize=(5, 5))
df_Human = draw_codon_pie('../Result/Frequency/Human_freq.csv')
plt.pie(df_Human['codon_frequency'], labels=df_Human['codon'], autopct='%1.1f%%', startangle=90, counterclock=False, textprops={'fontsize': 10})
plt.title('Human Codon Frequency Pie Chart')
plt.savefig('Human_Codon_frequency_pie_chart.pdf')

plt.rcParams['font.size'] = 12
plt.figure(figsize=(5, 5))
df_cleaned_Human = draw_protein_pie('../Result/Frequency/Human_freq.csv')
plt.pie(df_cleaned_Human['acid_frequency'], labels=df_cleaned_Human['protein3'], autopct='%1.1f%%', startangle=90, counterclock=False, textprops={'fontsize': 10})
plt.title('Human Protein Frequency Pie Chart')
plt.savefig('Human_Protein_frequency_pie_chart.pdf')

#小鼠
plt.rcParams['font.size'] = 12
plt.figure(figsize=(5, 5))
df_Mouse = draw_codon_pie('../Result/Frequency/Mouse_freq.csv')
plt.pie(df_Mouse['codon_frequency'], labels=df_Mouse['codon'], autopct='%1.1f%%', startangle=90, counterclock=False, textprops={'fontsize': 10})
plt.title('Mouse Codon Frequency Pie Chart')
plt.savefig('Mouse_Codon_frequency_pie_chart.pdf')

plt.rcParams['font.size'] = 12
plt.figure(figsize=(5, 5))
df_cleaned_Mouse = draw_protein_pie('../Result/Frequency/Mouse_freq.csv')
plt.pie(df_cleaned_Mouse['acid_frequency'], labels=df_cleaned_Mouse['protein3'], autopct='%1.1f%%', startangle=90, counterclock=False, textprops={'fontsize': 10})
plt.title('Mouse Protein Frequency Pie Chart')
plt.savefig('Mouse_Protein_frequency_pie_chart.pdf')

#酵母
plt.rcParams['font.size'] = 12
plt.figure(figsize=(5, 5))
df_Yeast = draw_codon_pie('../Result/Frequency/Yeast_freq.csv')
plt.pie(df_Yeast['codon_frequency'], labels=df_Yeast['codon'], autopct='%1.1f%%', startangle=90, counterclock=False, textprops={'fontsize': 10})
plt.title('Yeast Codon Frequency Pie Chart')
plt.savefig('Yeast_Codon_frequency_pie_chart.pdf')

plt.rcParams['font.size'] = 12
plt.figure(figsize=(5, 5))
df_cleaned_Yeast = draw_protein_pie('../Result/Frequency/Yeast_freq.csv')
plt.pie(df_cleaned_Yeast['acid_frequency'], labels=df_cleaned_Yeast['protein3'], autopct='%1.1f%%', startangle=90, counterclock=False, textprops={'fontsize': 10})
plt.title('Yeast Protein Frequency Pie Chart')
plt.savefig('Yeast_Protein_frequency_pie_chart.pdf')
