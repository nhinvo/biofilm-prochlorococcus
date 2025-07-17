"""
Purpose: to plot mapping coverage. 

Based on Konnor's genome depth scripts: 
    - https://github.com/konnorve/Variant-finding-workflow/blob/main/workflow/scripts/make_dup_del_plot.py

Nhi Vo - 06/10/24
"""
from pathlib import Path 
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt

def depth_histogram(fpath, df):
    """
    Plot mapping coverage. 
    """
    fname = fpath.name.split('.')[0].split('_ref')[0]

    plt.hist(df.mean_bin_depth, bins=20)
    plt.title('Histogram of Mean Bin Depth (bins=20)')
    plt.xlabel('Mean Bin Depth')
    plt.ylabel('Count')
    plt.savefig(f"data/3.depth_plots/depth_histogram/{fname}.png")
    plt.close()

def depth_genome(fpath, df):
    """
    """
    fname = fpath.name.split('.')[0].split('_ref')[0]
    Path('data/3.depth_plots/depth_genome').mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(12, 6))
    plt.scatter(df.bin_midpoint, df.mean_bin_depth, c='b', marker='.')
    plt.title('Mean Bin Depth Across Genome')
    plt.xlabel('Position')
    plt.ylabel('Mean Bin Depth')
    plt.savefig(f"data/3.depth_plots/depth_genome/{fname}.png")
    plt.close()

def main():
    parsed_depth_dir = Path("data/2.parsed_depth")

    # output dir for plots 
    Path('data/3.depth_plots/depth_histogram').mkdir(parents=True, exist_ok=True)
    Path('data/3.depth_plots/depth_genome').mkdir(parents=True, exist_ok=True)
    Path('data/3.depth_plots/treatment_zscore').mkdir(parents=True, exist_ok=True)

    ### 1. Plot individual scatter/histogram plots ###     
    # for fpath in parsed_depth_dir.glob("93*.tsv"):
    #     # import binned depth
    #     df = pd.read_table(fpath)

    #     # plot histogram of mean bin depth 
    #     depth_histogram(fpath, df)

    #     # plot scatter of mean bin depth 
    #     depth_genome(fpath, df)

    ### 2. Plot aggreagte treatment line plot ###
    df = pd.read_table('data/2.parsed_depth/aggregate_cov.tsv')
    treatment_groups = df.groupby(['treatment'])
    for index, tdf in treatment_groups:
        treatment = index[0]

        # get treatment min & max (for xlim)
        global_min = tdf['bin_id'].min()-2
        global_max = tdf['bin_id'].max()+2

        # split by rep 
        rep_groups = tdf.groupby(['rep'])
        num_reps = len(rep_groups)

        fig, axs = plt.subplots(num_reps, 1, figsize=(15, 5 * num_reps), sharex=True, sharey=True)


        for i, (index, rdf) in enumerate(rep_groups):
            rep = index[0]

            axs[i].plot(rdf['bin_id'], rdf['mean_bin_depth_zscore'], label=rep, alpha=0.9)
            axs[i].set_ylabel('z-score normalized depth')
            axs[i].set_xlabel('bin_id')
            axs[i].set_xlim(global_min, global_max)
            axs[i].legend()

        fig.suptitle(f'z-score normalized depth for {treatment}.', fontsize=16, y=1.02)
        plt.savefig(f"data/3.depth_plots/treatment_zscore/{treatment}.png")
        plt.close()

if __name__ == "__main__":
    main()