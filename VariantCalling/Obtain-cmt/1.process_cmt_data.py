"""
Extract candidate mutation data from WideVariant snakemake workflow
into a .tsv file. 
"""
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import numpy as np

def get_base_at_position(fasta_file, position):
    # Parse the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Access the sequence from the record
        sequence = record.seq
        # Get the base at the given position (1-based index)
        base = sequence[position - 1]
        return base

def obtain_ref(cmtfile, ref_dir):
    """
    Obtain reference genome as string. 
    For mapping with SNP pos in cmt df later to obtain ref base. 
    """
    group = cmtfile.stem.split('_')[1]
    group_ref_path = f'{ref_dir}/{group}/genome.fasta'

    for record in SeqIO.parse(group_ref_path, "fasta"):
        ref_genome_sequence = record.seq

    return ref_genome_sequence 



def import_cmt(fpath, ref_genome):
    """
    Import unfiltered cmt data from WideVariant pipeline. 
    Add ref base to cmt. 
    """
    data = {}
    # obtain data from npx cmy file 
    f = np.load(fpath)  # load np file
    data_groups = list(f.keys())  # list of names of files 
    for data_group in data_groups:
        data[data_group] = f[data_group]
    f.close()

    # load data into pandas df
    dfs = []
    for index, sname in enumerate(data['sample_names']):
        A,T,C,G,a,t,c,g = data['counts'][index].T
        cov = data['counts'][index].sum(axis=1)
        df = pd.DataFrame(data=[
            data['p'],
            abs(data['quals'][index]),
            cov,
            A,T,C,G,a,t,c,g,
        ]).T
        df.columns = ['position', 'qual','coverage', 'A','T','C','G','a','t','c','g']
        df['sample'] = sname
        dfs.append(df)

    df = pd.concat(dfs)

    # maps reference to position 
    def get_base_at_position(position):
        return ref_genome[position - 1]  # Convert 1-based to 0-based index

    df['ref'] = df['position'].apply(get_base_at_position)

    return df


def main():
    ### Paths ### 
    cmt_dir = Path("../Modified-WideVariant-Pipeline/2-Case/candidate_mutation_table")  # path to snakemake cmt result 
    ref_dir = "../1-run-WV/ref_genome"
    outdir = Path("output")
    outdir.mkdir(parents=True, exist_ok=True)

    ### Process cmt file ### 
    for cmtfile in cmt_dir.glob("*table.npz"):
        fname = cmtfile.stem

        # obtain path to ref genome 
        ref_genome = obtain_ref(cmtfile, ref_dir)

        df = import_cmt(cmtfile, ref_genome)

        df_out = f"{outdir}/{fname}_all_calls_cmt.tsv"
        df.to_csv(df_out, sep='\t', index=False)


if __name__ == "__main__":
    main()