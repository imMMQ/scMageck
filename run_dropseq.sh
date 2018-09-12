#!/bin/bash


#$ -j y
#$ -V 
#$ -pe pvm 4
#$ -cwd

#python run_dropseq.py bam/CROP-seq_Jurkat_TCR_stimulated_r1.bam outputbam/TCR_sti_r1
#python run_dropseq.py bam/CROP-seq_Jurkat_TCR_stimulated_r2.bam outputbam/TCR_sti_r2
echo "bam: $1"
echo "output: $2"

python run_dropseq.py $1 $2

