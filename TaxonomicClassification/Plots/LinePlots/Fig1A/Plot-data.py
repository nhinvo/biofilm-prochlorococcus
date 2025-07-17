"""
Purpose: to obtain data plotted. 
"""
import pandas as pd 
from pathlib import Path 

def main():
    df = pd.read_table("../../../StandardizeData/data/AllSummaryReadCount.tsv")

    # filter for read count only (out of percent and reads) & depth 
    df = df[(df['summary_type'] == 'reads')].copy()

    # obtain total read count 
    df['classified_reads'] = df['other_genus'] + df['Prochlorococcus'] + df['Synechococcus'] 

    # obtain pro percent (out of all classified reads )
    df['pro_percent'] = (df['Prochlorococcus'] / df['classified_reads']) * 100 

    groups = df.groupby(['cell_state', 'binned_depth'])
    data = []
    for index, gdf in groups:
        group_data = {
            'cell_state': gdf['cell_state'].iloc[0], 
            'binned_depth': gdf['binned_depth'].iloc[0], 
            'mean': gdf['pro_percent'].mean(), 
            'std': gdf['pro_percent'].std(), 
            'sem': gdf['pro_percent'].sem(), 
            'sample_size': len(gdf), 
        }

        data.append(group_data)

    df = pd.DataFrame(data)
    df.to_excel('Fig1A-Data.xlsx', index=False)


main()