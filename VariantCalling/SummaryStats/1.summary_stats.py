"""
Purpose: to obtain summary statistics for the variants called. 
"""

import pandas as pd 
from pathlib import Path

def import_cmt_pos(cmt_dir):
    """
    Import and combine the candidate mutation tables for 9301 and 9313.
    """
    pos_data = {}
    for fpath in Path(cmt_dir).glob('*.tsv'):
        group = fpath.stem.split('_')[1]
        df = pd.read_table(fpath)
        pos = set(df['position'].values.tolist())

        pos_data[group] =  list(pos)

    return pos_data

def parse_info_col(row):
    """
    Parse INFO row from vcf file. 
    """
    return(item.split('=') for item in row.split(';'))

def import_vcf_df(vcf_dir_path):
    """
    For .variant.vcf.gz files: 
        - import vcf file as pandas df 
        - combine variant files for all samples
        - split up data in rows of INFO column 
    """
    # cycle through vcf files for all samples 
    dfs = []
    for variant_fpath in Path(vcf_dir_path).glob('*strain.variant.vcf.gz'):
        sample_name = variant_fpath.stem.split('_')[0]
        group = variant_fpath.stem.split('ref_')[1].split('_aligned')[0]

        df = pd.read_csv(
            variant_fpath, 
            compression = 'gzip', 
            delimiter = '\t', 
            comment = '#', 
            names = ['CHROM', 'position', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'OTHER']
        )

        df['sample'] = sample_name
        df['group'] = group
        dfs.append(df)

    # combine vcf of all samples into 1 df 
    df = pd.concat(dfs)

    # split data in the INFO column into their own column 
    data = []
    for index, row in df.iterrows():
        # split data in INFO column 
        parsed_info_data = parse_info_col(row['INFO'])

        # create a new row to add parsed INFO data to 
        new_row = row.to_dict()
        new_row.update(parsed_info_data)

        # add new row dict to list to make new df
        data.append(new_row)

    # make df with parsed INFO data 
    df = pd.DataFrame(data)
    df.drop(columns=['INFO'], inplace=True)

    # obtain unique positions in each ref genome group 
    pos_data = {}
    ref_groups = df.groupby(['group'])
    for index, gdf in ref_groups:
        group = index[0]
        pos = set(gdf['position'].values.tolist())

        pos_data[group] =  list(pos)

    # save parsed vcf df 
    df.to_excel('data/parsed_vcf_data.xlsx', index=False)

    return df, pos_data

def match_check(cmt_pos, vcf_pos, vcf_df):
    """
    Check if there are any SNPs in vcf that were filtered out in cmt. 
    If there are, obtain their vcf stats. 
    """
    ref_genomes = ['MIT9301', 'MIT9313']
    dfs = []
    for ref_genome in ref_genomes:
        cmt_list = cmt_pos[ref_genome]
        vcf_list = vcf_pos[ref_genome]

        unique_to_vcf = [item for item in vcf_list if item not in cmt_list]
        df = vcf_df[vcf_df['position'].isin(unique_to_vcf)].copy()
        df['group'] = ref_genome
        dfs.append(df)

    df = pd.concat(dfs)
    df.to_excel('data/vcf_pos_not_in_cmt.xlsx', index=False)
    

def main():
    """
    """
    # obtain positions from cmt files 
    cmt_dir_path = "../Obtain-cmt/output/" 
    cmt_pos = import_cmt_pos(cmt_dir_path)

    # import vcf data 
    vcf_dir_path = "../Modified-WideVariant-Pipeline/1-Mapping/vcf"
    vcf_df, vcf_pos = import_vcf_df(vcf_dir_path)

    # check if vcf and cmt match 
    match_check(cmt_pos, vcf_pos, vcf_df)


if __name__ == "__main__":
    main()