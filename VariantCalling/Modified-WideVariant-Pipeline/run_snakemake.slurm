#!/bin/bash
#SBATCH --job-name widevariant.main
#SBATCH -p sched_mit_chisholm
#SBATCH -n 1
#SBATCH --time=4-0
#SBATCH --mem=8GB
#SBATCH -o logs/snakemake.%j.out
#SBATCH -e logs/snakemake.%j.err


# Activate conda environment (may need to change name of env)
eval "$(conda shell.bash hook)"  
conda activate widevariant-snakemake

snakemake --profile .

# Print "Done!!!" at end of main log file
echo Done!!!
