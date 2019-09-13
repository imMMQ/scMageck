\name{scmageck-package}
\alias{scmageck-package}
\alias{scmageck}
\docType{package}
\title{
\packageTitle{scmageck}
}
\description{
\packageDescription{scmageck}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{scmageck}
\packageIndices{scmageck}
    scMAGeCK is a computational model to identify genes associated with multiple 
    expression phenotypes from CRISPR screening coupled with single-cell RNA sequencing
    data (CROP-seq).scMAGeCK is based on our previous MAGeCK and MAGeCK-VISPR models for 
    pooled CRISPR screens. 
    
    The scMAGeCK manuscript can be found at bioRxiv(https://www.biorxiv.org/content/10.1101/658146v1/).
}
\author{
\packageAuthor{scmageck}

Maintainer: \packageMaintainer{scmageck}
}

\keyword{ package }
\examples{
    ### BARCODE file contains cell identity information, generated from the cell identity collection step
    BARCODE <- system.file("extdata","barcode_rec.txt",package = "scmageck")
    
    ### RDS can be a Seurat object or local RDS file path that contains the scRNA-seq dataset
    RDS <- system.file("extdata","singles_dox_mki67_v3.RDS",package = "scmageck")
    
    target_gene <- "MKI67"
    
    rra_result <- scmageck_rra(BARCODE=BARCODE, RDS=RDS, GENE=target_gene,
             LABEL='dox_mki67', NEGCTRL=NULL, KEEPTMP=FALSE, PATHWAY=FALSE, SAVEPATH=NULL)
    head(rra_result)
    
    lr_result <- scmageck_lr(BARCODE=BARCODE, RDS=RDS, LABEL='dox_scmageck_lr', 
            NEGCTRL = 'NonTargetingControlGuideForHuman', PERMUTATION = 1000, SAVEPATH=NULL)
    lr_score <- lr_result[1]
    lr_score_pval <- lr_result[2]
    head(lr_score_pval)
}
