"""
"""
import pandas as pd 
import numpy as np
from pathlib import Path 

def main():
    Path('data/AbsoluteGE_AverageEfficiency').mkdir(exist_ok=True)

    # 1. Import table of all extraction efficiencies 
    efficiencies = "data/EfficiencyMetadata/FilterSizeAveragedStandard.tsv"
    efficiencies_df = pd.read_table(efficiencies)

    # 2. Obtain table of genome equivalent 
    standardized_counts = "../../StandardizeData/data/AllNormalizedCount.tsv"
    df = pd.read_table(standardized_counts)
    df = df.rename(columns={'filter size (um)': 'Filter size (um)'})  # change to match 
    df = df[
        ~(df['sample_name'].str.contains('ERR')) &
        ~(df['sample_name'].str.contains('SRR'))
    ]  # remove TARA samples 

    # 3. merge efficiency with genome equivalent
    df = pd.merge(df, efficiencies_df, on='Filter size (um)', how='left')

    # 4. obtain absolute genome count from efficiencies
    df['Absolute Genome Equivalent'] = df['genome_equivalents'] / df['Efficiency']
    
    # Save subclade absolute genome quivalent 
    df.to_csv('data/AbsoluteGE_AverageEfficiency/SubcladeAbosluteGenomeEquivalent.tsv', sep='\t', index=False)
    df.to_excel('data/AbsoluteGE_AverageEfficiency/SubcladeAbosluteGenomeEquivalent.xlsx', index=False)

    # # end script here. doesn't make sense to obtain absolute GE of subclade (HL/LL)
    # return False 

    # 5. Obtain clade level genome equivalent (HL/LL)
    df = df[df['genus'] == 'Prochlorococcus'].copy()  # Pro rows only 

    # handling AMZ rows: 
    df.loc[df['clade'].str.contains('AMZ'), 'clade'] = 'unclassified'  # turn them into unclassified 

    # obtain LL/HL from clade info 
    df['classification'] = df['clade'].apply(
        lambda x: 'HL' if 'HL' in str(x) 
        else ('LL' if 'LL' in str(x) else 'unclassified')
    )

    samples_df = df.copy()  # save table before count sum()
    samples_df = samples_df.drop_duplicates(subset='sample_name')

    df = df.groupby(['sample_name', 'classification'])['Absolute Genome Equivalent'].sum().reset_index()
    df = pd.merge(df, samples_df[['sample_name', 'depth', 'Filter size (um)']], on='sample_name', how='left')
    df['depth'] = df['depth'].astype(int)

    # log absolute GE 
    # df['log10 Absolute Genome Equivalent'] = np.log10(df['Absolute Genome Equivalent'])

    # add PB / FL 
    df['cell_state'] = np.where(df['Filter size (um)'] >= 5, 'Particle-bound', 'Free-living')

    # save clade level absolute genome equivalent 
    df.to_csv('data/AbsoluteGE_AverageEfficiency/CladeAbosluteGenomeEquivalent.tsv', sep='\t', index=False)  
    df.to_excel('data/AbsoluteGE_AverageEfficiency/CladeAbosluteGenomeEquivalent.xlsx', index=False)    

main()