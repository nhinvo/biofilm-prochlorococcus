#!/bin/bash

sequencing_dir=/orcd/data/chisholm/001/chisholmlab/experiment_repository/2024/240403Chi  
for sample_dir in $sequencing_dir/*/; do 
    sample_name=`basename $sample_dir`
    echo $sample_name
    mkdir -p $sample_name  # make directory to copy input over if does not exist
    cp $sample_dir/*fastq $sample_name
done