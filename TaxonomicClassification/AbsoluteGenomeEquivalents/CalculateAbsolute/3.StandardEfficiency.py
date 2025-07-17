"""
Purpose: to obtain extraction eficiency for each sample. 

Note: DNA standard for the HOT cruise samples = 0.1 ng/uL

N. N. Vo 
10/28/24
"""
import os
import pandas as pd 
from pathlib import Path 
from scipy.stats import zscore
from Bio.SeqUtils import molecular_weight
from Bio import SeqIO


# length of illumina read 
READ_LEN = 150 

# file paths 
GENOME_FASTA = "../../Modified-ProSynTax-Workflow/inputs/Thermus_genome/Thermus_thermophilus_ATCC_BAA_163.fasta" 
SAMPLES_TSV = "../../Modified-ProSynTax-Workflow/inputs/MakeHOT346SamplesTsv/HOT346samples.tsv"
STANDARD_LOG = "data/Edited HOT 346 log sheets Oct 28.xlsx"
BLAST_DIR = "data/BlastOutput"


def genome_stats(genome_fasta):
    """
    Returns: 
        - genome_len: length (total base bair count) of genome and plasmid
        - genome_mw: molecular weight of genome and plasmid
    """
    # read in each fasta contig as a biopython record object 
    fasta_records = list(SeqIO.parse(genome_fasta, "fasta"))

    # extract fasta sequences for each contig 
    seqs = [seq.seq for seq in fasta_records]

    # obtain total length of contigs
    genome_len = len(str(seqs[0])) + len(str(seqs[1]))

    # obtain MW 
    contig1_mw = molecular_weight(seqs[0], "DNA")
    contig2_mw = molecular_weight(seqs[1], "DNA")

    # total MW of the 2 contigs
    genome_mw = contig1_mw + contig2_mw

    return genome_len, genome_mw

def process_added(standard_log, samples_tsv, genome_mw):
    """
    Return df with cols: ['sample', 'Standard Added Molecules']

    Obtain number of molecules of DNA standard added: 
        - Obtain Thermus genome molecular weight 
        - Obtain mass of standard added from volume and concentration 
        - Convert mass from ug to g 
        - Obtain moles of standard added: mass / molecular weight 
        - Obtain molecules of standard added: moles * avogadro's number  
    """
    # import log sheet of DNA standard volume added - final df = 2 cols: [Vial #, Standard Added Vol (uL)]
    standard_df = pd.read_excel(standard_log, sheet_name="size fractions")
    standard_df.columns = standard_df.columns.str.strip()
    standard_df = standard_df[standard_df['Exp/play'] == 'exp']  # remove 'play' samples 
    standard_df['Notes'] = standard_df['Notes'].astype(str)
    standard_df['Notes'] = standard_df['Notes'].fillna('')

    # for rows where samples were combined, combine vial number
    standard_df.loc[standard_df['Notes'].str.contains('DNA extraction combined: '), 'Vial #'] = standard_df['Notes'].str.split(' ').str[-1].str.replace('+', '')
    standard_df = standard_df[['Vial #', 'Vol of DNA standard added (uL)']]
    standard_df = standard_df.dropna(subset=['Vol of DNA standard added (uL)'])
    standard_df = standard_df.rename(columns={'Vol of DNA standard added (uL)': 'Standard Added Vol (uL)'})
    standard_df['Vial #'] = standard_df['Vial #'].astype(int)

    # import samples tsv (from pipeline) to map with sample names 
    samples_df = pd.read_table(samples_tsv)[['sample', 'Vial #', 'Depth', 'Filter size (um)', 'Notes', 'date extracted']]
    samples_df['Vial #'] = samples_df['Vial #'].astype(int)

    # merge dfs 
    df = pd.merge(standard_df, samples_df, on='Vial #', how='right')

    # obtain mass of standard added in grams 
    df['Standard Added Concentration (ng/uL)'] = 0.1
    df['Standard Added Mass (ng)'] = df['Standard Added Concentration (ng/uL)'] * df['Standard Added Vol (uL)']
    df['Standard Added Mass (g)'] = df['Standard Added Mass (ng)'] * 1e-9
    
    # obtain moles of standard added 
    df['Standard Added Moles'] = df['Standard Added Mass (g)'] / genome_mw

    # obtain molecules of standard genome added 
    df['Standard Added Molecules'] = df['Standard Added Moles'] * 6.022e23

    return df

def process_sequenced(dir_path, genome_len):
    """
    """
    parsed_blast_outpath = f'{dir_path}/parsed_blast.tsv'

    # import file if already parsed 
    if os.path.exists(parsed_blast_outpath):
        df = pd.read_table(parsed_blast_outpath)
        return df
    
    # if file doesnt exist, make it 
    dfs = []
    count = 0
    for fpath in Path(dir_path).glob('K*tsv'):
        # obtain sample name and read type (i.e. '1' or '2') from file path 
        sname = fpath.stem[:-2]
        read_type = fpath.stem[-1]

        # import blast output df 
        cols = ['read_name', 'genome_contig', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        df = pd.read_table(fpath, names=cols)[['read_name', 'pident', 'length']]

        # sort hits by pident and alignment length (highest first)
        df = df.sort_values(by=['pident', 'length'], ascending=[False, False])

        # obtain top hits 
        df = df.drop_duplicates(subset=['read_name'], keep='first')

        # add sample name 
        df['sample'] = sname

        dfs.append(df)
        count += 1
        # if count == 6: 
        #     break

    df = pd.concat(dfs)

    # sort hits by pident and alignment length (highest first)
    df = df.sort_values(by=['pident', 'length'], ascending=[False, False])

    # obtain top hits (between fwd and rev reads)
    df = df.drop_duplicates(subset=['sample', 'read_name'], keep='first')

    # obtain count of hits per sample 
    df = df.groupby('sample').size().reset_index(name='Standard Recovered Read Count')
    
    # obtain total base pair of mapped reads 
    df['Standard Recovered Base Pair'] = df['Standard Recovered Read Count'] * READ_LEN 

    # obtain genome_equivalent (coverage)
    df['Standard Recovered Molecules'] = df['Standard Recovered Base Pair'] / genome_len

    df.to_csv(parsed_blast_outpath, sep='\t', index=False)

    return df


def main():
    Path('data/EfficiencyMetadata').mkdir(exist_ok=True, parents=True)

    # obtain reference genome data: length, molecular weight 
    genome_len, genome_mw = genome_stats(GENOME_FASTA)

    # obtain standard genome (molecules) added 
    added_df = process_added(STANDARD_LOG, SAMPLES_TSV, genome_mw)

    # obtain standard genome recovered (sequenced/blasted)
    sequenced_df = process_sequenced(BLAST_DIR, genome_len)

    # merge added vs. recovered data and obtain efficiency 
    df = pd.merge(added_df, sequenced_df, on='sample', how='inner')
    df['Efficiency'] = df['Standard Recovered Molecules'] / df['Standard Added Molecules']

    # obtain z-score of efficiencies - for determining outliers  
    df['Efficiency Z-scores'] = zscore(df['Efficiency'])
    df['Absolute Efficiency Z-scores'] = df['Efficiency Z-scores'].abs()

    # Save all efficiencies 
    df.to_csv('data/EfficiencyMetadata/StandardData.tsv', sep='\t', index=False)



    # Note: this part is done after script 4.EfficiencyCheck.py where 
    # the conclusion was to group by filter size 
    # remove the sample with high Z-score 
    df = df[df['Absolute Efficiency Z-scores'] < 7].copy()

    # group by filter size and obtain average efficiency 
    df = df.groupby('Filter size (um)', as_index=False)['Efficiency'].mean()
    df.to_csv('data/EfficiencyMetadata/FilterSizeAveragedStandard.tsv', sep='\t', index=False)

main()