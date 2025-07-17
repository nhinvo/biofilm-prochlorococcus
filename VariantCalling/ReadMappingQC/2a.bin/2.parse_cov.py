"""
To obtain binned depth and z-score. 

Based on Konnor's genome depth scripts: 
    - https://github.com/konnorve/Variant-finding-workflow/blob/main/workflow/scripts/gen_genome_depth_variant_data.py
"""
from pathlib import Path 
import pandas as pd 
import numpy as np

def process_cov_file(fpath, outdir):
    """
    Import output from samtools depth (cov across all positions).
    """
    fname = fpath.name.split('.')[0].split('_ref')[0]

    df = pd.read_table(fpath, names=['contig', 'pos', 'depth'])
    depth_arr = df.depth.to_numpy()

    bin_size = 1000

    genome_size = len(depth_arr)

    bins = np.linspace(0, genome_size, int(genome_size/bin_size), dtype=int)

    bin_info = []
    # for each bin (1000bp windows)
    for i in range(len(bins)-1):
        bin_start = bins[i]
        bin_end = bins[i+1]
        bin_size = bin_end - bin_start
        bin_midpoint = (bin_end + bin_start) / 2
        mean_bin_depth = np.mean(depth_arr[bin_start:bin_end])

        bin_info.append([bin_start, bin_end, bin_size, bin_midpoint, mean_bin_depth])

    df = pd.DataFrame(bin_info, columns=["bin_start", "bin_end", "bin_size", "bin_midpoint", "mean_bin_depth"])
    
    # assign unique ID to bin 
    df['bin_id'] = range(1, (len(df) + 1))

    # calculate z-score
    df['mean_bin_depth_zscore'] = (df['mean_bin_depth'] - df['mean_bin_depth'].mean()) / df['mean_bin_depth'].std()
    
    df.to_csv(f'{outdir}/{fname}.tsv', sep='\t', index=False)


def main():
    # input path 
    samtools_depth_dir = Path("../data/1.samtools_depth_output")

    ### 1. Obtain binned depth: raw and zscore normalized ### 
    outdir = 'data/2.parsed_depth'
    Path(outdir).mkdir(parents=True, exist_ok=True)

    for fpath in samtools_depth_dir.glob("9301*.cov"):
        # process_cov_file(fpath, outdir)
        pass

    ### 2. Aggregate all 9301 results by treatment for plotting ### 
    dfs = []
    for fpath in Path(outdir).glob('9301*tsv'):
        sname = fpath.stem
        treatment = sname[:-3]
        rep = sname[-1]

        df = pd.read_table(fpath)[['bin_id', 'mean_bin_depth_zscore', 'mean_bin_depth']]
        df['treatment'] = treatment
        df['rep'] = rep

        dfs.append(df)

    df = pd.concat(dfs)
    df.to_csv('data/aggregate_cov.tsv', sep='\t', index=False)

    ### 3. obtain pos where std > 2
    dfs = []
    for fpath in Path(outdir).glob('9301*tsv'):
        sname = fpath.stem

        df = pd.read_table(fpath)[['bin_id', 'mean_bin_depth']]
        df['sname'] = sname

        mean = df['mean_bin_depth'].mean()
        std = df['mean_bin_depth'].std()

        # obtain positions with coverage of >2 standard deviation 
        df['> 2std bin_depth'] = (df['mean_bin_depth'] > mean + 2 * std) | (df['mean_bin_depth'] < mean - 2 * std)

        # filter for outliers 
        df = df[df['> 2std bin_depth'] == True]

        dfs.append(df)

    df = pd.concat(dfs)

    # count number of times the bin shows up 
    bin_count_df = df['bin_id'].value_counts().reset_index()
    bin_count_df.columns = ['position', 'sample_count']

    # 
    sample_count = df['sname'].nunique()
    # filter for pos where >2std is not in all samples (sample count=53)
    print(bin_count_df[bin_count_df['sample_count'] < sample_count])
    print(len(bin_count_df[bin_count_df['sample_count'] < sample_count]))

    # count sample_count
    sample_count_df = bin_count_df['sample_count'].value_counts().reset_index()
    sample_count_df.columns = ['sample_count', 'count']

    sample_count_df.to_excel('data/std_count.xlsx', index=False)



main()