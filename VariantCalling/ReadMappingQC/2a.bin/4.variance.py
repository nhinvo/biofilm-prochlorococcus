import pandas as pd 

def main():
    df = pd.read_table('data/2.parsed_depth/aggregate_cov.tsv')

    # groupby treatment and obtain variance
    # treatment_groups = df.groupby(['treatment'])
    # for index, tdf in treatment_groups:
    #     treatment = index[0]

    df = df.groupby(['bin_id'])['mean_bin_depth_zscore'].var().reset_index()
    df = df.rename(columns={'mean_bin_depth_zscore': 'std_mean_bin_depth_zscore'})
    df.to_csv('data/temp_var_bin_only.tsv', sep='\t', index=False)
        

main()