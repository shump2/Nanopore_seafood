import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

parser = argparse.ArgumentParser(description='Generate figures from input data file')
parser.add_argument('-i', '--input', help='input CSV file', required=True)
parser.add_argument('-o', '--output', help='output figure file prefix', required=True)
args = parser.parse_args()

# Read the input data
data = pd.read_csv(args.input)

# Get the unique sample IDs and sort them
sample_ids = sorted(data['sample_id'].unique(), key=lambda x: int(x[7:]))

# Set up the plot layout
n_samples = len(sample_ids)
ncols = 2
nrows = n_samples // ncols + (n_samples % ncols > 0)

# Get the max and min for x and y axes
x_max, x_min = data['read_length'].max(), data['read_length'].min()
y_max, y_min = data['base_quality_score'].max(), data['base_quality_score'].min()

# Create the subplots
fig, axes = plt.subplots(nrows, ncols, figsize=(12, 6 * nrows))
axes = axes.flatten()

# Generate the plots
for i, sample_id in enumerate(sample_ids):
    sample_data = data[data['sample_id'] == sample_id]
    ax = axes[i]
    sns.scatterplot(
        x='read_length',
        y='base_quality_score',
        data=sample_data,
        hue='read_accuracy',
        palette='coolwarm',
        ax=ax
    )
    ax.set_title(f'Sample ID: {sample_id}')
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])

    # Add Pearson's correlation
    r, p = pearsonr(sample_data['read_length'], sample_data['base_quality_score'])
    ax.annotate(f"r = {r:.2f}", xy=(0.8, 0.2), xycoords='axes fraction', fontsize=10)
    ax.annotate(f"p = {p:.2e}", xy=(0.8, 0.1), xycoords='axes fraction', fontsize=10)

# Remove any unused subplots
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

# Adjust layout and save the figure
#fig.suptitle('Figure : Read Length vs Base Quality Score')
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(f'{args.output}_read_length_vs_base_quality_score.png', dpi=300)

