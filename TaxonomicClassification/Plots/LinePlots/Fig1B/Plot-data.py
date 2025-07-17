"""
Purpose: to obtain data plotted. 
"""
import pandas as pd 
from pathlib import Path 

def main():
    plot_data = '../../../StandardizeData/data/AllNormalizedCounts_PlotFormat.tsv'
    df = pd.read_table(plot_data)
    df = df[~df['classification'].str.contains('unclassified')]

    groups = df.groupby(['cell_state', 'binned_depth', 'classification'])
    data = []
    for index, gdf in groups:
        group_data = {
            'cell_state': gdf['cell_state'].iloc[0], 
            'binned_depth': gdf['binned_depth'].iloc[0], 
            'classification': gdf['classification'].iloc[0], 
            'mean': gdf['percentage'].mean(), 
            'std': gdf['percentage'].std(), 
            'sem': gdf['percentage'].sem(), 
            'sample_size': len(gdf), 
        }
        
        data.append(group_data)

    df = pd.DataFrame(data)
    df.to_excel('Fig1B-Data.xlsx', index=False)
    
main()