"""
Purpose: to make a final samples.tsv file for ease of use downstream. 
    - outputs: 
        - data/StandardizedSamples.tsv
        - data/StandardizedSamples.xlsx

output StandardizedSamples.tsv columns: 
    - ['sample', 'depth', 'filter size (um)', 'cell_state']
    - cell_state: 'Free-living' (<5um filters) or 'Particle-bound' (>= 5um filters)
"""
import pandas as pd 
import numpy as np
from pathlib import Path 

CELL_STATE_THRESHOLD = 5

def process_TARA(df):
    """
    Reformat TARA df to match with remaining datasets. 
    """
    # change col name 
    df = df.rename(columns={'filter_size': 'filter size (um)'}) 

    # remove rows without metadata 
    df = df.dropna(subset=['depth', 'filter size (um)'])

    # fix filter_size values in TARA samples (i.e. ['20' '0.22' from ['20-180' '0.22-1.6' '0.22-3'])
    df['filter size (um)'] = df['filter size (um)'].str.split('-').str[0]

    # remove samples with lat-long > 40 and < -40 
    # print(f'Length of TARA df before latitude filtering: {len(df)}')
    df = df[(df['latitude'] < 40) & (df['latitude'] > -40)]
    # print(f'Length of TARA df after latitude filtering: {len(df)}')

    return df 

def import_samples(dir_path):
    """
    Return samples_df of all *_samples.tsv files. 
    Make all samples_df contain the same columns and content. 
    """
    dfs = [] 

    for fpath in Path(dir_path).glob('*_samples.tsv'):
        df = pd.read_table(fpath)
        
        # convert column names into lowercase
        df.columns = df.columns.str.lower()  

        # process TARA df
        if 'TARA' in str(fpath):
            df = process_TARA(df)

        df = df[['sample', 'depth', 'filter size (um)']]

        dfs.append(df)

    df = pd.concat(dfs)

    return df

def process_samples(df):
    """
    Obtain cell state and binned depths. 
    """
    # obtain sample cell state
    df['filter size (um)'] = df['filter size (um)'].astype(float)
    df['cell_state'] = np.where(
        df['filter size (um)'] >= CELL_STATE_THRESHOLD, 
        'Particle-bound', 'Free-living'
    )

    # bin depths 
    df['binned_depth'] = pd.cut(
        df['depth'], bins=[0, 5, 66, 90, 160, 300, 2000, 5000], 
        labels=[5, 45, 75, 150, 300, 1000, 4000], right=True
    )

    # depth filtering
    df = df[df['depth'] <= 1000]

    # change col name
    df = df.rename(columns={'sample': 'sample_name'})

    return df


def main():
    Path('data').mkdir(exist_ok=True)

    # import samples.tsv tables 
    df = import_samples('data/SamplesTsv')

    # standardize df
    df = process_samples(df)

    # save 
    df.to_csv('data/StandardizedSamples.tsv', sep='\t', index=False)
    # df.to_excel('data/StandardizedSamples.xlsx', index=False)

main()