# scMAGeCK

scMAGeCK is a computational model to identify genes associated with multiple expression phenotypes from CRISPR screening coupled with single-cell RNA sequencing data (CROP-seq).

scMAGeCK is based on our previous [MAGeCK](http://mageck.sourceforge.net) and MAGeCK-VISPR models for pooled CRISPR screens. The scMAGeCK manuscript can be found at [bioRxiv](https://www.biorxiv.org/content/10.1101/658146v1/).

Questions? Comments? Join the [MAGeCK Google group](https://groups.google.com/forum/?hl=en#!forum/mageck).

# Installation

Currently, no installation is required to run scMAGeCK pipeline (see the demo). We will deposit the pipeline into bioconda channel soon.

## Dependencies 

scMAGeCK depends on the following softwares or packages, all available via [bioconda](https://bioconda.github.io/) channel.

* [MAGeCK](https://mageck.sourceforge.net)
* pysam
* R 
* [Seurat](https://satijalab.org/seurat/)

To install these dependencies via bioconda, install miniconda (Python v3) first (see instructions [here](https://docs.conda.io/en/latest/miniconda.html)), then use the following command:

    conda create -n scmageck -c bioconda mageck pysam r-seurat
    
This will create an environment named "scmageck" that includes all dependencies. Use the following command to activate this environment:

    source activate scmageck

Sometimes Seurat may not be installed successfully using bioconda. If this is the case, follow the instructions of [Seurat installation](https://satijalab.org/seurat/install.html).


# Demos

## Demo 1: running scMAGeCK from known single cell identity

This demo provides a mini example to run both modules of scMAGeCK: RRA and LR. 
In the terminal, type

    bash run_scmageck_rra.sh   

or

    bash run_scmageck_lr.sh

to run both demos.

The explanations of each input/outp files can be found in the Output files section below.

## Demo 2: collect single cell identity information from bam files

This demo provides a mini example to collect single cell identity from a CROP-seq experiment. The single-cell RNA-seq is performed on 10X Genomics platform.
In the terminal, type

    bash run_count_bam.sh

to run this demo.


The explanations of each input/outp files can be found in the Output files section below.

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

    --lib-grna LIB_GRNA   
                        A gRNA library file containing the list of sgRNA names, their sequences and associated genes, separated by tab.
    --files FILES [FILES ...]
                        A list of fastq files, SAM/BAM files, or a wildcard of filenames.

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
                       Genes whose expressions are to be tested. Multiple genes can be provided, separated by ",". For example, "MKI67,TP53"

optional arguments:

    --RRAPATH
                       The path to the RRA program, if RRA cannot be found in the PATH environment variable.
    --LABEL
                       The label of the output file.
    --NEGCTRL
                       The name of the negative control gene. For example, "NonTargetingControlGuideForHuman". Default is NULL (do not use any negative controls).
    --KEEPTMP
                       Keep temporary files.
    --PATHWAY
                       Treat genes in --GENE option as a pathway. In other words, the averaged expression of these genes will be used for testing.



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

# R environment

If you'd like to use scMAGeCK in your R script, scMAGeCK package can be installed by following command:

    R CMD install scmageck_0.99.0.tar.gz

Here are the examples: 

    BARCODE = system.file("extdata", "barcode_rec.txt", package = "scmageck")
    
    RDS = system.file("extdata", "singles_dox_mki67_v3.RDS", package = "scmageck")
    
    RRAPATH = system.file("extdata", "RRA", package = "scmageck")
  
    scmageck_rra(BARCODE = BARCODE, RDS = RDS, GENE = "MKI67", RRAPATH = RRAPATH,  
             LABEL = 'dox_mki67', NEGCTRL = NULL, KEEPTMP = FALSE, PATHWAY = FALSE)

    scmageck_lr(BARCODE = BARCODE, RDS = RDS, LABEL = 'dox_scmageck_lr', 
            NEGCTRL = 'NonTargetingControlGuideForHuman', PERMUTATION = 1000)


# Output files

## scMAGeCK-RRA

scMAGeCK-RRA will output the ranking and p values of each perturbed genes, using the RRA program in MAGeCK. Users familiar with the MAGeCK program may find it similar with the [gene_summary](https://sourceforge.net/p/mageck/wiki/output/#gene_summary_txt) output in MAGeCK.

Here is the example output of scMAGeCK-RRA:

    Row.names       items_in_group.low      lo_value.low    p.low   FDR.low goodsgrna.low   items_in_group.high     lo_value.high   p.high  FDR.high        goodsgrna.high
    TP53    271     0.11832 0.95619 1       48      271     1.014e-83       4.9975e-06      0.00015 184

Explanations of each column are below:

|Column|Content|
|------|-------|
|Row.names| Perturbed gene name|
|items_in_group.low| The number of single-cells with each gene perturbed |
|lo_value.low | The RRA score in negative selection (reducing the marker expression if this gene is perturbed). The RRA score uses a p value from rank order statistics to measure the degree of selection; the smaller score, the stronger the selection is. More information on the calculation of RRA score can be found in our original [MAGeCK](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0554-4) paper. |      
|p.low  | The raw p-value (using permutation) of this gene in negative selection  |  
|FDR.low | The false discovery rate of this gene in negative selection |
|goodsgrna.low | The number of single-cells that passes the threshold and is considered in the RRA score calculation in negative selection|
|items_in_group.high| The same as items_in_group.low: the number of single-cells with each gene perturbed) |
|lo_value.high| The RRA score in positive selection (increasing the marker expression if this gene is perturbed|      
|p.high| The raw p-value (using permutation) of this gene in positive selection  |  
|FDR.high| The false discovery rate of this gene in positive selection |
|goodsgrna.high| The number of single-cells that passes the threshold and is considered in the RRA score calculation in positive selection|


## scMAGeCK-LR

scMAGeCK-LR will generate several files below:

|File|Description|
|----|----------|
|_score.txt|The score (similar with log fold change) of each perturbed gene (rows) on each marker gene (columns)|
|_score.pval.txt|The associated p values of each score|
|LR.RData|An R object to store scores and p values|

The format of score.txt and score.pval.txt is a simple table file with rows corresponding to perturbed genes and columns corresponding to marker genes. For example in the score.txt,

    Perturbedgene        APC     ARID1A    TP53    MKI67
    APC     0.138075836476524       -0.0343441660045313     0.214449590551132       -0.150287676553705


This row records the effects of perturbing APC gene on the expressions of APC, ARID1A, TP53 and MKI67.

This file can also be imported directly into R (using read.table). 


## cropseq_count (cell identity)

The output of cell identity contains a single text file below:


    cell    barcode sgrna   gene    read_count      umi_count
    AAATCAACGGGTGA-1        NF1_sg_118      AGTCAGTACTGAGCACAACA    NF1     1187    30
    AAATCAACGGGTGA-1        CDKN2A_sg_70    TCTTGGTGACCCTCCGGATT    CDKN2A  1       1
    AAATCAACGGGTGA-1        SETD2_sg_157    AGTTCTTCTCGGTGTCCAAA    SETD2   1       1
    AAATCAACGGGTGA-1        ARID1B_sg_14    GGAAGCAACCAGTCTCGATC    ARID1B  1       1


|Column|Description|
|------|-----------|
|cell|Single-cell name (cell barcode)|
|barcode|The sgRNA ID detected|
|sgrna|The sgRNA sequence|
|gene|The target gene|
|read_count|Number of reads|
|umi_count|Number of UMIs|


# History

* Version 0.1


