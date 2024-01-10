import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA


# 假设你的8个物种的数据存储在8个文件中，文件名分别为 species1.csv, species2.csv, ..., species8.csv
# 可以根据实际情况修改文件名和路径
file_paths = ['../Result/Frequency/Arabidopsis_freq.csv','../Result/Frequency/Drosophila_freq.csv',
              '../Result/Frequency/C_elegants_freq.csv','../Result/Frequency/Human_freq.csv','../Result/Frequency/Yeast_freq.csv',
              '../Result/Frequency/Mouse_freq.csv']

# 读取数据
dfs = [pd.read_csv(file_path) for file_path in file_paths]

# 合并数据，假设数据中的列名为 'codon' 和 'codon_frequency'
merged_df = pd.concat(dfs, keys=['Arabidopsis', 'Drosophila', 'Nematode', 'Human', 'Yeast', 'Mouse'], names=['Species'])

X = merged_df[['acid_frequency', 'codon_frequency']]
#Kmeans聚类
kmeans = KMeans(n_clusters=6)
kmeans.fit(X)
merged_df['cluster'] = kmeans.labels_


sns.set(style="whitegrid")
plt.rcParams['font.size'] = 14
plt.figure(figsize=(6, 5))
sns.scatterplot(data=merged_df, x='acid_frequency', y='codon_frequency', palette='Set2', hue='Species', edgecolor="black", linewidth=0.01, s=40)
plt.title('KMeans Clustering of Species based on Codon and Acid Frequency')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.legend(title='Species', bbox_to_anchor=(1, 1), loc='upper left')
plt.savefig('../Visualization/C/Synonymous_codon_cor_scatter_plot.pdf')

plt.figure(figsize=(6, 5))
sns.scatterplot(data=merged_df, x='acid_frequency', y='codon_frequency', palette='Set2', hue='cluster', edgecolor="black", linewidth=0.01, s=40)
plt.title('KMeans Clustering of Species based on Codon and Acid Frequency')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.legend(title='Cluster', bbox_to_anchor=(1, 1), loc='upper left')
plt.savefig('../Visualization/C/Synonymous_codon_cor_scatter_plot_cluster.pdf')
# X = merged_df[['protein1_frequency', 'codon_frequency']]
# #Kmeans聚类
# kmeans = KMeans(n_clusters=6)
# kmeans.fit(X)
# merged_df['cluster'] = kmeans.labels_
#
# plt.figure(figsize=(6, 5))
# sns.scatterplot(data=merged_df, x='protein1_frequency', y='codon_frequency', palette='Set2', hue='Species', edgecolor="black", linewidth=0.01, s=40)
# plt.title('KMeans Clustering of Species based on Codon and Synonym Frequency')
# plt.xlabel('Principal Component 1')
# plt.ylabel('Principal Component 2')
# plt.legend(title='Species', bbox_to_anchor=(1, 1), loc='upper left')
# plt.savefig('../Visualization/C/Synonymous_codon_cor_scatter_plot1.pdf')