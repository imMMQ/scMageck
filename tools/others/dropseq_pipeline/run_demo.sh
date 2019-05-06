#!/bin/bash



#$ -j y
#$ -V 
#$ -pe pvm 4
#$ -cwd

# step 0: before running this demo, go to the demo folder, 
# and make sure all the reference genome sequence and annotation 
# is properly downloaded
# refer to prepare_demo.sh to download these files

# step 1: generate reference index
python generate_ref.py metadata/config.yaml 


# step 2: process bam files

python run_dropseq.py -c metadata/config.yaml

#python run_dropseq.py $1 $2

