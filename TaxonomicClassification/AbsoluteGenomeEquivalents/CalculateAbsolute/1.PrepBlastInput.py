"""
Purpose: make input table for bowtie2 read mapping against Thermus. 

"""
import pandas as pd 
from pathlib import Path 

def main():
    Path('data').mkdir(exist_ok=True)
    
    lab_trimmed = "../../Modified-ProSynTax-Workflow/scratch/trimmed_reads"

    data = {}
    # cycle through all read files 
    for fpath in Path(lab_trimmed).glob('*'):
        fname = fpath.name.replace('_trimmed.fastq.gz', '')  # sample name + read type 
        read_type = fname.split('_')[-1]  # '1' or '2' 
        sample_name = '_'.join(fname.split('_')[:-1])  # sample name

        # init list to store data for sample if not exist 
        if sample_name not in data.keys():
            data[sample_name] = {}
        
        data[sample_name][f'read_{read_type}'] = fpath

    df = pd.DataFrame.from_dict(data, orient='index').reset_index()
    df.columns = ['sample_name', 'forward_read', 'reverse_read']
    df.to_csv('data/blast_input.tsv', sep='\t', index=False)

main()