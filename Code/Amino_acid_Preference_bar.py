#图B 同一氨基酸偏好性（6个物种6个图）
import pandas as pd
import matplotlib.pyplot as plt
import re
from matplotlib.backends.backend_pdf import PdfPages


# Get a list of all CSV files in the directory
file_list = ['../Result/Frequency/Arabidopsis_freq.csv','../Result/Frequency/Drosophila_freq.csv',
              '../Result/Frequency/C_elegants_freq.csv','../Result/Frequency/Yeast_freq.csv',
              '../Result/Frequency/Mouse_freq.csv','../Result/Frequency/Human_freq.csv']

# Define a color map for amino acids
amino_acid_colors = {
    'Ala': 'red',
    'Arg': 'blue',
    'Asn': 'green',
    'Asp': 'purple',
    'Cys': 'orange',
    'Gln': 'brown',
    'Glu': 'pink',
    'Gly': 'cyan',
    'His': 'yellow',
    'Ile': 'magenta',
    'Leu': 'lime',
    'Lys': 'olive',
    'Phe': 'navy',
    'Pro': 'gold',
    'Ser': 'coral',
    'Thr': 'indigo',
    'Tyr': 'peru',
    'Val': 'salmon',
    'STOP': 'grey'
}

# Loop through each file
with PdfPages('../Visualization/B/Synonymous_codon_bar.pdf') as pdf:
    for file_path in file_list:
        Species = re.search(r'/([^/]*)_', file_path)

        # Read the CSV file
        data = pd.read_csv(file_path, index_col='codon')
        data = data.drop(['protein1', 'codon_counts', 'protein_counts', 'protein1_counts', 'codon_frequency'], axis=1)
        data = data[~data['protein3'].isin(['Trp', 'Met'])]

        # Plot the data with different colors for each amino acid
        plt.rcParams['font.size'] = 14
        fig, axes = plt.subplots(nrows=4, ncols=5, figsize=(20, 15))
        i = 0
        j = 0
        for amino_acid, color in amino_acid_colors.items():
            amino_acid_data = data[data['protein3'] == amino_acid]
            axes[i, j].bar(amino_acid_data.index, amino_acid_data['protein1_frequency'], color=color)

            # Customize the plot
            axes[i, j].set_xlabel('Codons')
            axes[i, j].set_ylabel('Frequency')
            #axes.set_title(f'Amino Acid Preference - {file_path}')
            axes[i, j].set_title(f'{Species.group(1)}-{amino_acid}')

            # Update indices
            if i < 3:
                i += 1
            elif i == 3 and j < 4:
                i = 0
                j += 1
        # Add an empty subplot for the 21st amino acid
        if j == 4 and i == 3:
            axes[i, j].bar(amino_acid_data.index, amino_acid_data['protein1_frequency'], color=color)

        plt.tight_layout()
        pdf.savefig()

        plt.close()
