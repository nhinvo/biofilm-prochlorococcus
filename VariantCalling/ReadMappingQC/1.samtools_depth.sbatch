#!/usr/bin/env bash
#SBATCH --job-name=depth
#SBATCH --time 1-0                         
#SBATCH -p sched_mit_chisholm             
#SBATCH -c 1                            
#SBATCH -N 1
#SBATCH --mem 10G
#SBATCH --array=1-58%40
#SBATCH -o logs/cov.%a.%j.out
#SBATCH -e logs/cov.%a.%j.err

# path to dir with all the bam mapping files 
bam_dir=../Modified-WideVariant-Pipeline/1-Mapping/bowtie2

# obtain path to genome (based on slurm task)
bam_fpath=$(realpath ${bam_dir}/*.bam | sed -n ${SLURM_ARRAY_TASK_ID}p)  

# file name without .bam extention 
name=`basename $bam_fpath .bam` 

# path to samtools depth output 
depth_outdir=data/1.samtools_depth_output
mkdir -p ${depth_outdir}

eval "$(conda shell.bash hook)"
conda activate samtools

echo Running samtools depth on ${name}
echo File path: ${bam_fpath}
samtools depth -o ${depth_outdir}/${name}.cov ${bam_fpath}
# samtools idxstats ${bam_fpath} > ${depth_outdir}/${name}.tsv
echo Completed samtools depth on ${name}

conda deactivate 