#!/bin/bash


# unzip the compressed fastq files

tar  xvzf fastq/cropseq_testfq.tar.gz -C fastq/ 
# download reference genome
# here, we just use chr22 as an example
# to download full reference genome, use 

# source for download:
#wget -O ref/chr22.fa.gz ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
tar xvzf ref/chr22.fa.tar.gz -C ref/

# upzip reference gtf file
# here, we just use chr22 as an example
# to download full reference gtf, use

# source for download:
# wget ftp://ftp.ensembl.org/pub/release-77/gtf/homo_sapiens/Homo_sapiens.GRCh38.77.gtf.gz 

tar xvzf ref/Homo_sapiens.GRCh38.77.chr22.gtf.tar.gz -C ref/



