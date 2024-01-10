#图D 跟公共数据比对Fraction和Frequency
import pandas as pd
import matplotlib.pyplot as plt
import re
from matplotlib.backends.backend_pdf import PdfPages

# Get a list of all CSV files in the directory
file_list = ['../Result/Fraction/public_data_F_Arabidopsis.csv', '../Result/Fraction/public_data_F_C_elegants.csv',
              '../Result/Fraction/public_data_F_Drosophila.csv', '../Result/Fraction/public_data_F_Human.csv',
             '../Result/Fraction/public_data_F_Mouse.csv', '../Result/Fraction/public_data_F_Yeast.csv']

# Define colors for Fraction and Frequency
colors = {'Fraction': 'blue', 'Frequency': 'orange'}

plt.rcParams['font.size'] = 14
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(23, 10))
i = 0
j = 0

# Loop through each file
# with PdfPages('../Visualization/D/Comparison_Fraction_dot.pdf') as pdf:
#     for file_path in file_list:
#         Species = re.search(r'_(.*?)\.', file_path)
#
#         # Read the CSV file
#         data = pd.read_csv(file_path, index_col='codon')
#         data = data.drop(['protein1', 'protein3', 'codon_counts', 'protein_counts', 'protein1_counts'], axis=1)
#
#         # Calculate the difference between "Fraction" and "Fraction_in_the_DataFrame_Self"
#         data['Difference'] = data['Fraction'] - data['Fraction_self_count']
#         data['Difference_'] = data['Frequency/ Thousand'] - data['Frequency_self_count/ Thousand']
#
#         # Plot a scatter plot for Fraction
#         axes[i, j].scatter(data.index, data['Difference'], label=f'{Species.group(1)} - Fraction', color=colors['Fraction'])
#
#         if i < 1:
#             i += 1
#         elif i == 1 and j < 2:
#             i = 0
#             j += 1
#
#     # Add an empty subplot for the 21st amino acid
#     if j == 2 and i == 1:
#         axes[i, j].scatter([], [], color='white')
#
#     # Customize the plot
#     for ax in axes.flat:
#         ax.set_xlabel('Codons')
#         ax.set_ylabel('Values')
#         ax.legend(loc = 1)
#         ax.set_ylim(-0.05, 0.07)
#         ax.set_title('Difference in Fraction from public data')
#         ax.set_xticks([])
#
#     pdf.savefig()
#     plt.close()


with PdfPages('../Visualization/D/Comparison_Frequency_dot.pdf') as pdf:
    for file_path in file_list:
        Species = re.search(r'_(.*?)\.', file_path)

        # Read the CSV file
        data = pd.read_csv(file_path, index_col='codon')
        data = data.drop(['protein1', 'protein3', 'codon_counts', 'protein_counts', 'protein1_counts'], axis=1)

        # Calculate the difference between "Fraction" and "Fraction_in_the_DataFrame_Self"
        data['Difference'] = data['Fraction'] - data['Fraction_self_count']
        data['Difference_'] = data['Frequency/ Thousand'] - data['Frequency_self_count/ Thousand']

        # Plot a scatter plot for Frequency
        axes[i, j].scatter(data.index, data['Difference_'], label=f'{Species.group(1)} - Frequency', color=colors['Frequency'])

        if i < 1:
            i += 1
        elif i == 1 and j < 2:
            i = 0
            j += 1

    # Add an empty subplot for the 21st amino acid
    if j == 2 and i == 1:
        axes[i, j].scatter([], [], color='white')

    # Customize the plot
    for ax in axes.flat:
        ax.set_xlabel('Codons')
        ax.set_ylabel('Values')
        ax.legend(loc=1)
        ax.set_ylim(-3, 3.5)
        ax.set_title('Difference in Frequency/thousand from public data')
        ax.set_xticks([])

    pdf.savefig()
    plt.close()