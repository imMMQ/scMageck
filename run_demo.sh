#!/bin/bash



#$ -j y
#$ -V 
#$ -pe pvm 4
#$ -cwd

# step 1: generate reference index
python generate_ref.py metadata/config.yaml 


# step 2: process bam files

python run_dropseq.py -c metadata/config.yaml

#python run_dropseq.py $1 $2

