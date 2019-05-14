# scMAGeCK

scMAGeCK is a computational model to interrogate genes associated multiple expression phenotypes from CRISPR screens coupled with single-cell sequencing data.

scMAGeCK is based on our previous [MAGeCK](http://mageck.sourceforge.net) and MAGeCK-VISPR model for pooled CRISPR screens.

# Installation

Currently, no installation is required to run scMAGeCK pipeline (see the demo). We will deposit the pipeline into bioconda channel soon.

## Dependencies 

scMAGeCK depends on the following softwares or packages, all available via bioconda channel.

* MAGeCK
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

# Usage


## collect CROP-seq cell identity information from various scRNA-seq platforms 

The pipeline includes a program, "cropseq_count.py", to collect cell identity information from various file types (bam/fastq) and platforms.

The second demo (demo2) contains an example of using this program.

    usage: cropseq_count.py [-h] [-v] [-n OUTPUT_PREFIX] --lib-grna LIB_GRNA
                        --file-type {fastq,bam,paired-fastq}
                        [--no-reverse-complement] [-m MAX_MISMATCH]
                        [--anchor-before ANCHOR_BEFORE]
                        [--anchor-after ANCHOR_AFTER]
                        [--files FILES [FILES ...]]

cropseq-count: counting sgRNAs from CROPseq experiments.

required arguments:

    --lib-grna LIB_GRNA   A gRNA library file containing the list of sgRNA
                        names, their sequences and associated genes, separated
                        by tab.

optional arguments:

    -h, --help            show this help message and exit
    -v, --version         show program's version number and exit
    -n OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        The prefix of the output file(s). Default sample1.
    --file-type {fastq,bam,paired-fastq}
                        The type of the file to search guide RNA sequence.
    --no-reverse-complement
                        Do not perform reverse complement search of sgRNAs.
    -m MAX_MISMATCH, --max-mismatch MAX_MISMATCH
                        Maximum number of mismatches to be considered in sgRNA search. Default 2. Not recommended for values greater than 2. Decrease this value to speed up the search.
    --anchor-before ANCHOR_BEFORE
                        Anchor sequence before the sgRNA. Default GAAACACCG (at the end of U6 promoter).
    --anchor-after ANCHOR_AFTER
                        Anchor sequence after the sgRNA. Default GTTTTAGAG.
    --files FILES [FILES ...]
                        A list of fastq files, SAM/BAM files, or a wildcard of filenames.


## Use RRA to test the association of gene knockout with certain marker expression

The R script "run_individual_genes.R"  provides an easy-to-run interface to test the association of gene knockout with certain marker expression.

The first demo (demo1) contains an example of running this R script.

    usage: Rscript run_individual_genes.R OPTIONS

required arguments: 

    --BARCODE
                       A txt file to include cell identity information, generated from the cell identity collection step. 
    --RDS
                       An R RDS file of [Seurat](https://satijalab.org/seurat/) object containing the scRNA-seq dataset. Note that the dataset has to be normalized and scaled.
    --GENE
                       The specific gene name to be tested for association (e.g., "MKI67").

optional arguments:

    --RRAPATH
                       The path to the RRA program, if RRA cannot be found in the PATH environment variable.
    --LABEL
                       The label of the output file.



## Use linear regression to test the association of gene knockout with all possible genes

The R script, "run_lr.R", provides an easy-to-run interface to use linear regression for the association of gene knockout with the expression of all the genes.

The first demo (demo1) contains an example of running this R script.

    usage: Rscript run_lr.R OPTIONS

required arguments: 

    --BARCODE
                       A txt file to include cell identity information, generated from the cell identity collection step. 
    --RDS
                       An R RDS file of [Seurat](https://satijalab.org/seurat/) object containing the scRNA-seq dataset. Note that the dataset has to be normalized and scaled.
    --NEGCTRL
                       The name of the genes (separated by ",") served as negative controls.

optional arguments:

    --LABEL
                       The label of the output file.
    --PERMUTATION
                       The number of permutations for p value calculation.



# History

* Version 0.1


