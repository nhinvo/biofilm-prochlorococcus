"""
Purpose: reformat results tsv files for ease of plotting. 
"""
import pandas as pd 
import numpy as np
from pathlib import Path 

def parse_normalized(df):
    """
    Returns df with cols: [sample_name classification  genome_equivalents  percentage]
        -  classification: HL, LL, or unclassified 
            - Note: AMZs are considered unclassified 
        - percentage: percent of that classification's genome_equivalents out of the sum 

    """
    # obtain Pro rows only 
    df = df[df['genus'] == 'Prochlorococcus'].copy()

    # handling AMZ rows: 
    df.loc[df['clade'].str.contains('AMZ'), 'clade'] = 'unclassified'  # turn them into unclassified 

    # obtain LL/HL from clade info 
    df['classification'] = df['clade'].apply(
        lambda x: 'HL' if 'HL' in str(x) 
        else ('LL' if 'LL' in str(x) else 'unclassified')
    )

    df = df.groupby(['sample_name', 'classification'])['genome_equivalents'].sum().reset_index()

    # for each sample, obtain percentage 
    df['total_genome_equivalents'] = df.groupby('sample_name')['genome_equivalents'].transform('sum')
    df['percentage'] = (df['genome_equivalents'] / df['total_genome_equivalents']) * 100

    df = df.drop(columns=['total_genome_equivalents'])

    return df

def main():
    # make output dir
    Path('data').mkdir(parents=True, exist_ok=True)

    # import normalized df
    df = pd.read_table('data/AllNormalizedCount.tsv')
    samples_df = df.copy()[['sample_name', 'depth', 'filter size (um)', 'cell_state', 'binned_depth']]
    samples_df = samples_df.drop_duplicates()

    # parse data 
    df = parse_normalized(df)

    # map samples metadata back
    df = pd.merge(df, samples_df, on='sample_name', how='left')
    print(len(df))
    df.to_csv('data/AllNormalizedCounts_PlotFormat.tsv', sep='\t', index=False)
    df.to_excel('data/AllNormalizedCounts_PlotFormat.xlsx', index=False)

main()