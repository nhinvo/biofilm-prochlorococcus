import pandas as pd 
from pathlib import Path 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    HTseq_dir = "data/2.parsed_depth/"

    dfs = [] 
    for fpath in Path(HTseq_dir).glob("9301*"):
        df = pd.read_table(fpath)

        # add exp name, treatment, timepoint
        df['sname'] = fpath.stem.replace('9301-', '')
        df['treatment'] = df['sname'].str.split('-').str[0].str[:-1]
        df['timepoint'] = df['sname'].str.split('-').str[0].str[-1]
        df['rep'] = df['sname'].str.split('-').str[1].str[-1].astype(int)

        dfs.append(df)

    df = pd.concat(dfs)
    df = df.sort_values(by=['treatment', 'timepoint', 'rep'])
    outlier_df = df[(df['mean_bin_depth_zscore'] > 3) | (df['mean_bin_depth_zscore'] < -3)].copy()
    outlier_count = outlier_df['bin_id'].value_counts()
    outlier_df = outlier_df[~outlier_df['bin_id'].isin(outlier_count[outlier_count < 50])]
    print(outlier_count)
    print(outlier_df)

    # plotting 
    for col in ['mean_bin_depth', 'mean_bin_depth_zscore']:
        plt.figure(figsize=(20, 8))

        # Create a boxplot
        sns.violinplot(x='sname', y=col, data=df, palette='Pastel1')
        # if col == 'mean_bin_depth_zscore':
        #     palette = sns.color_palette("Set1", n_colors=outlier_df['bin_id'].nunique()) 
        #     sns.stripplot(
        #         x='sname', y=col, data=outlier_df, 
        #         color='red', size=6, jitter=True, marker='o', 
        #         hue='bin_id', palette=palette, dodge=True
        #     ) 

        #     plt.legend(title='bin_id', bbox_to_anchor=(1.05, 0.5), loc='center left')

        plt.xticks(rotation=90)

        # Show plot
        title = col.replace('_', ' ').title()
        plt.title(f'Violin plot of {title}')
        plt.savefig(f"data/violinplot_{col}.png", dpi=300, bbox_inches="tight")
        plt.close()


main()