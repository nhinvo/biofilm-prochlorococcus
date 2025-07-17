"""
"""
import pandas as pd 
from pathlib import Path 
import matplotlib.pyplot as plt

def main():
    Path('data').mkdir(exist_ok=True)

    df = pd.read_table('../CalculateAbsolute/data/AbsoluteGE/CladeAbosluteGenomeEquivalent.tsv')

    # number of plot rows based on number of depths sequenced 
    fig, axes = plt.subplots(nrows=df['depth'].nunique(), ncols=1, figsize=(5, 4 * df['depth'].nunique()))
    axes = axes.flatten()  # flatten for cycling through 

    df=df[df['Filter size (um)']>0.2].copy()

    # create unique index (x-axis label)
    df['Filter size (um) str'] = df['Filter size (um)'].astype(str)
    df['bar_label'] =  df['sample_name'] + "_" + df['Filter size (um) str']

    # sort and groupby depth
    df = df.sort_values(by='depth', ascending=True)

    groups = df.groupby('depth', sort=False)  # group while keeping previous sort
    for num, (depth, gdf) in enumerate(groups): 
        ax = axes[num]  # obtain plotting axes 

        # pivot table for bar plotting 
        gdf = gdf.pivot_table(
            index='bar_label', 
            columns='classification', 
            values='Absolute Genome Equivalent', 
            fill_value=0
        ).reset_index()

        # sort bar by filter size
        gdf['filter_size'] = gdf['bar_label'].str.split('_').str[-1]
        gdf['size_int'] = gdf['filter_size'].str[:-2].astype(float)
        gdf = gdf.sort_values(by='size_int', ascending=True)

        # obtain bar labels 
        bar_labels= gdf['bar_label']
        x = range(len(bar_labels)) 

        # plot HL and LL bars, centering around the same x-tick 
        ax.bar(x, gdf['HL'], label='HL', width=0.4, align='center', color='#afd48a')
        ax.bar([i+0.4 for i in x], gdf['LL'], label='LL', width=0.4, align='center', color='#8a9bcf')

        # ticks and labels
        ax.set_xticks([i + 0.2 for i in x])
        ax.set_xticklabels(bar_labels, rotation=90)

        ax.set_title(f'{depth}m')
        ax.set_ylabel('Absolute Genome Equivalent')
        ax.set_xlabel('')

        ax.legend()

    # Adjust layout
    plt.tight_layout()
    plt.savefig(f'data/No_0.2_AbsoluteGenomeEquivalentBar.png', dpi=300)  # Save as PNG with 300 DPI
    plt.close()


main()