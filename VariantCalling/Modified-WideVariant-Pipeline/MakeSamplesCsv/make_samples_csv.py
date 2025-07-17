"""
Purpose: to make samples tsv file as input to WV snakemake 
pipeline using data from Maya's spreadsheet. 
"""
import pandas as pd
from pathlib import Path

def format_numbers(num):
    return f'{num:02d}'

def main():
    df = pd.read_excel('dnaLogSheet_submitted20240404.xlsx')

    # convert single digit int (e.g. 1,2,3) into double digits (e.g. 01,02,03)
    # this col will be used for mapping to raw seq samples 
    df['Log index'] = df['Log index'].astype(str).str.zfill(2)

    # create 'Path' column 
    input_sample_dir = "../0-input/"
    df['Path'] = input_sample_dir + "D24-1600" + df['Log index'] + "-6710L/" 

    # create 'Sample' column
    df.rename(columns={'Name': 'Sample'}, inplace=True)

    # create 'FileName' column - e.g. 240403Chi_D24-160027_1_sequence.fastq
    df['FileName'] = "240403Chi_D24-1600" + df['Log index']

    # create 'Reference' and 'Group' columns
    df['Reference'] = "MIT" + df['Ref Genome'].astype(str)
    df['Group'] = "MIT" + df['Ref Genome'].astype(str)

    # create Outgroup column
    df['Outgroup'] = 0  # no outgroups 

    # Remove first 3 rows - Christine's samples 
    df = df.iloc[3:]

    # reorder column (WV pipeline imports csv by col order)
    df = df[['Path', 'Sample', 'FileName', 'Reference', 'Group', 'Outgroup']]

    df.to_csv('../samples.csv', index=False, sep=',')
    
if __name__ == "__main__":
    main()
