"""
Purpose: to obtain percentage of Pro in the community (classified reads) 
for all samples from AllSummaryReadCount.tsv file. 

03/05/25
"""
import pandas as pd
from pathlib import Path 

def main():
    df = pd.read_table('data/AllSummaryReadCount.tsv')

    # obtain read counts only 
    df = df[df['summary_type'] == 'reads']

    # obtain total classified reads
    df['total_classified'] = df['Prochlorococcus'] + df['Synechococcus'] + df['other_genus']

    # obtain percent of pro out of total classified 
    df['pro_percent_in_community'] = (df['Prochlorococcus'] / df['total_classified']) * 100 

    # drop columns
    df = df.drop(columns=['summary_type', 'Prochlorococcus', 'Synechococcus', 'other_genus', 'unclassified'])

    df.to_csv('data/AllProInCommunityPercent.tsv', sep='\t', index=False)
    df.to_excel('data/AllProInCommunityPercent.xlsx', index=False)

main()