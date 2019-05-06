# scMAGeCK

scMAGeCK is a computational model to interrogate genes associated multiple expression phenotypes from CRISPR screens coupled with single-cell sequencing data.



# Installation

## Dependencies 

scMAGeCK depends on the following softwares or packages, all available via bioconda channel.

* pysam
* R 
* Seurat

# Demos

## Demo 1: running scMAGeCK from known single cell identity

This demo provides a mini example to run both modules of scMAGeCK: RRA and LR. 
In the terminal, type

    bash run_scmageck_rra.sh   

or

    bash run_scmageck_lr.sh

to run both demos.

## Demo 2: collect single cell identity information from bam files

This demo provides a mini example to collect single cell identity from a CROP-seq experiment. The single-cell RNA-seq is performed on 10X Genomics platform.
In the terminal, type

    bash run_count_bam.sh

to run this demo.


# History

* Version 0.1


