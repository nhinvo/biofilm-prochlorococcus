"""
Purpose: 

"""
from pathlib import Path 
import pandas as pd 
import numpy as np

def main():
    samtools_depth_dir = "../data/1.samtools_depth_output/"

    dfs = []
    for fpath in Path(samtools_depth_dir).glob('9301*cov'):
        df = pd.read_table(fpath, names=['contig', 'position', 'count'])
        
        sname = fpath.name.split('.')[0].split('_ref')[0]
        treatment = sname[:-3]
        rep = sname[-1]
        df['treatment'] = treatment
        df['rep'] = rep        

        mean = df['count'].mean()
        std = df['count'].std()

        # obtain positions with coverage of >2 standard deviation 
        df['> 2std'] = (df['count'] > mean + 2 * std) | (df['count'] < mean - 2 * std)

        # filter for outliers 
        df = df[df['> 2std'] == True]

        dfs.append(df)

    df = pd.concat(dfs)

    # count number of times position is present in df 
    pos_count_df = df['position'].value_counts().reset_index()
    pos_count_df.columns = ['position', 'sample_count']

    # filter for pos where >2std is not in all samples (sample count=53)
    print(pos_count_df[pos_count_df['sample_count'] < 53])
    print(len(pos_count_df[pos_count_df['sample_count'] < 53]))

    # count sample_count
    sample_count_df = pos_count_df['sample_count'].value_counts().reset_index()
    sample_count_df.columns = ['sample_count', 'count']

    sample_count_df.to_excel('data/std_count.xlsx', index=False)


main()