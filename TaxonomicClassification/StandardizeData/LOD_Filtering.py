"""
Purpose: to filter for samples that pass LOD. 
    - outputs: 
        - data/LOD_table.tsv
        - data/AllSummaryReadCount.tsv
        - data/AllNormalizedCount.tsv

Filtering parameters (per sample) to keep misclassification rate at <10% (0.1):
    - Pro/Syn ratio of > 0.23
    - Pro Abundance (Pro counts / total classified counts) > 0.0008 > 0.08%
    - Minimum Pro Read Count > 50,000

Note: classified read count col included in case.  
"""
import pandas as pd 
import numpy as np
from pathlib import Path 

# LOD filtering params 
PROSYN_RATIO = 0.24
PRO_ABUNDANCE = 0.0028  # >0.28%
MINIMUM_PRO_COUNT = 50000

def import_snakemake_output(dir_path):
    """
    Import summary_read_count and normalized_counts tables for all Snakemake runs. 
    """
    # import summary (raw read count) df
    summary_df = pd.read_table(f'{dir_path}/summary_read_count.tsv')

    # import normalized df 
    normalized_df = pd.read_table(f'{dir_path}/normalized_counts.tsv')

    return summary_df, normalized_df

def perform_LOD(df, sample_df):
    """
    Returns a list of samples that passed the LOD filtering. 
    LOD filtering criteria: 
        - Pro/Syn ratio of > PROSYN_RATIO AND
            - Calculated by: Classified Pro / Classified Syn 
        - Pro Abundance of > PRO_ABUNDANCE
            - Calculated by: Classified Pro / Total Classified Reads 
        - Pro Count > MINIMUM_PRO_COUNT
    
    Save a table of all the criteria. 
    """
    # filter for read count rows only (remove percent rows)
    df = df[df['summary_type'] == 'reads'].copy()  

    ## 1a. Obtain Pro/Syn Ratio ### 
    df['Pro/Syn Ratio'] = df['Prochlorococcus'] / df['Synechococcus']
    df['Pro/Syn Ratio Condition'] = (df['Pro/Syn Ratio'] > PROSYN_RATIO)

    ## 1b. Obtain Pro/Syn Abundance ## 
    df['Total Classified Read'] = df['Prochlorococcus'] + df['Synechococcus'] + df['other_genus']
    df['Pro Abundance'] = df['Prochlorococcus'] / df['Total Classified Read']
    df['Pro Abundance Condition'] = df['Pro Abundance'] > PRO_ABUNDANCE

    ## 1c. Pro Count ## 
    df['Pro Count Condition'] = df['Prochlorococcus'] >= MINIMUM_PRO_COUNT

    ### 2. Filter Col ### 
    ## 2a. merge to obtain sample metadata ## 
    df = pd.merge(df, sample_df, on=['sample_name'], how='inner')
    df = df[[
        'sample_name', 
        'cell_state', 
        'depth',
        'binned_depth', 
        'Total Classified Read', 
        'Prochlorococcus', 
        'Pro/Syn Ratio', 
        'Pro Abundance', 
        'Pro/Syn Ratio Condition', 
        'Pro Abundance Condition', 
        'Pro Count Condition', 
    ]]

    ### 3. Final LOD ### 
    df['Sample Pass'] = df[[
        'Pro/Syn Ratio Condition', 
        'Pro Abundance Condition', 
        # 'Pro Count Condition',
    ]].all(axis=1)  # both cols have to be True for sample to pass (True)

    # save table of whether samples passed/failed
    df.to_csv('data/LOD_table.tsv', sep='\t', index=False)
    df.to_excel('data/LOD_table.xlsx', index=False)

    # obtain list of samples that passed filter
    samples = df[df['Sample Pass'] == True]['sample_name'].tolist()

    return samples

def LOD_filter(df, LOD_samples, sample_df, fname):
    """
    Filter to remove samples that did not pass LOD filterting and save df. 
    """
    # filter table for samples that passed LOD
    df = df[df['sample_name'].isin(LOD_samples)]

    # merge with samples metadata 
    df = pd.merge(df, sample_df, on='sample_name', how='inner')

    # save table
    df.to_csv(f'data/{fname}.tsv', sep='\t', index=False)


def main():
    Path('data').mkdir(exist_ok=True)

    # import result tables from Snakemake runs
    summary_df, normalized_df = import_snakemake_output('../Modified-ProSynTax-Workflow/results')

    # import samples metadata 
    sample_df = pd.read_table('data/StandardizedSamples.tsv')

    # obtain list of samples that passed LOD filtering 
    LOD_samples = perform_LOD(summary_df, sample_df)

    # filter samples and save df 
    LOD_filter(summary_df, LOD_samples, sample_df, 'AllSummaryReadCount')
    LOD_filter(normalized_df, LOD_samples, sample_df, 'AllNormalizedCount')

main()