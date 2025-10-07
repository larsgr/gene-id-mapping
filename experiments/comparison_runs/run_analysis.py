from pathlib import Path
import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

DATA_PATH = Path('ens_lift_vs_native.tsv')
FIG_DIR = Path('figures')
FIG_DIR.mkdir(parents=True, exist_ok=True)
MPL_CACHE = Path('mpl-cache')
MPL_CACHE.mkdir(parents=True, exist_ok=True)

import os
os.environ['MPLCONFIGDIR'] = str(MPL_CACHE.resolve())

sns.set_theme(style='whitegrid')
plt.rcParams['figure.figsize'] = (8, 5)

columns = ['feature', 'annA', 'geneA', 'txA', 'annB', 'geneB', 'txB', 'stats']
df = pd.read_csv(DATA_PATH, sep='\t', names=columns, header=0)

genes = df[df['feature'] == 'gene'].copy()

def stats_to_dict(stat_string: str) -> dict:
    items = {}
    for chunk in str(stat_string).split(';'):
        if not chunk or '=' not in chunk:
            continue
        key, value = chunk.split('=', 1)
        items[key] = value
    return items

stats_df = genes['stats'].apply(stats_to_dict).apply(pd.Series)
merged = pd.concat([genes.reset_index(drop=True), stats_df], axis=1)
for column in ['best_jaccard_exon', 'best_jaccard_cds_phase']:
    merged[column] = pd.to_numeric(merged[column], errors='coerce')

merged['stableA'] = merged['geneA'].str.split(':').str[-1]
merged['stableB'] = merged['geneB'].str.split(':').str[-1]
merged['stable_same'] = merged['stableA'] == merged['stableB']

class_order = ['Green', 'Yellow', 'Red', 'NotMapped']
class_counts = merged['class'].value_counts().reindex(class_order).fillna(0).astype(int)
class_counts.to_csv(FIG_DIR / 'class_counts.csv', header=['count'])

ax = sns.barplot(x=class_counts.index, y=class_counts.values, palette='deep')
ax.set_title('Gene-level classification counts')
ax.set_xlabel('Class')
ax.set_ylabel('Gene pairs')
for container in ax.containers:
    ax.bar_label(container, fmt='%d', fontsize=10)
plt.tight_layout()
plt.savefig(FIG_DIR / 'class_counts.png', dpi=150)
plt.close()

stable_summary = merged.groupby(['class', 'stable_same']).size().unstack(fill_value=0)
stable_summary.to_csv(FIG_DIR / 'stable_id_comparison.csv')

stable_pct = (stable_summary.T / stable_summary.sum(axis=1)).T.reset_index().melt(
    id_vars='class', var_name='stable_same', value_name='percentage'
)
ax = sns.barplot(data=stable_pct, x='class', y='percentage', hue='stable_same', palette='Set2')
ax.set_title('Stable Ensembl ID agreement by class')
ax.set_ylabel('Percentage (%)')
ax.set_xlabel('Class')
ax.legend(title='Same stable ID')
plt.tight_layout()
plt.savefig(FIG_DIR / 'stable_id_comparison.png', dpi=150)
plt.close()

for column, title, filename in [
    ('best_jaccard_cds_phase', 'Best CDS-phase Jaccard', 'best_jaccard_cds_hist.png'),
    ('best_jaccard_exon', 'Best exon Jaccard', 'best_jaccard_exon_hist.png'),
]:
    ax = sns.histplot(merged[column].dropna(), bins=40, color='#66aa55')
    ax.set_title(title)
    ax.set_xlabel('Jaccard index')
    ax.set_ylabel('Gene pairs')
    plt.tight_layout()
    plt.savefig(FIG_DIR / filename, dpi=150)
    plt.close()

summary = {
    'total_gene_pairs': int(len(merged)),
    'class_counts': class_counts.to_dict(),
    'stable_same': int(merged['stable_same'].sum()),
    'stable_different': int((~merged['stable_same']).sum()),
}
(Path('summary.json')).write_text(json.dumps(summary, indent=2))
print('Analysis complete:', summary['total_gene_pairs'], 'gene pairs')
