#!/usr/bin/env bash
#SBATCH --job-name=blastn
#SBATCH --time 12:00:00                        
#SBATCH -p sched_mit_hill      
#SBATCH --exclusive
#SBATCH -c 16
#SBATCH --mem 64000
#SBATCH -o logs/blastn.%a.%j.out
#SBATCH -e logs/blastn.%a.%j.err
#SBATCH --array=18# 2-4,11,18-22,26,30,37,42-43,53,65,68,78,80-81,83,85,86,89

echo "$(date)"
eval "$(conda shell.bash hook)"

### 1. Obtain Inputs ### 
thermus_ref_genome=../../Modified-ProSynTax-Workflow/inputs/Thermus_genome/Thermus_thermophilus_ATCC_BAA_163.fasta
sample_name=$(sed "${SLURM_ARRAY_TASK_ID}q;d" data/blast_input.tsv | cut -f1 | tr -d '\r')
fwd_read=$(sed "${SLURM_ARRAY_TASK_ID}q;d" data/blast_input.tsv | cut -f2 | tr -d '\r')
rev_read=$(sed "${SLURM_ARRAY_TASK_ID}q;d" data/blast_input.tsv | cut -f3 | tr -d '\r')

### 2. Convert fastq to fasta ###
fasta_output_dir=data/TrimmedFastaOutput
mkdir -p ${fasta_output_dir}

conda activate seqtk 

echo "$(date): Converting .fastq reads into .fasta..."
seqtk seq -a ${fwd_read} > ${fasta_output_dir}/${sample_name}_1.fasta
seqtk seq -a ${rev_read} > ${fasta_output_dir}/${sample_name}_2.fasta
echo "$(date): Completed .fastq to .fasta conversion!"

conda deactivate

### 3. Run Blastn ### 
blast_output_dir=data/BlastOutput
mkdir -p ${blast_output_dir}

conda activate blast-2.16.0

echo "$(date): Beginning sample ${sample_name} forward read blast..."
blastn \
    -num_threads ${SLURM_CPUS_PER_TASK} \
    -query ${fasta_output_dir}/${sample_name}_1.fasta \
    -subject ${thermus_ref_genome} \
    -outfmt 6 -perc_identity 95 \
    -out ${blast_output_dir}/${sample_name}_1.tsv 

echo "$(date): Beginning sample ${sample_name} reverse read blast..."
blastn \
    -num_threads ${SLURM_CPUS_PER_TASK} \
    -query ${fasta_output_dir}/${sample_name}_2.fasta \
    -subject ${thermus_ref_genome} \
    -outfmt 6 -perc_identity 95 \
    -out ${blast_output_dir}/${sample_name}_2.tsv 

echo "$(date): Completed sample: ${sample_name}"
conda deactivate
