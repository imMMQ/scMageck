library(scMAGeCK)
library(R.utils)
# set the BARCODE and RDS file path 
BARCODE = system.file("extdata","barcode_rec.txt",package = "scMAGeCK")
## RDS can be a Seurat object or local RDS file path that contains the scRNA-seq dataset
RDS = system.file("extdata","singles_dox_mki67_v3.RDS",package = "scMAGeCK")
## set the gene signature file path 
SIGNATURE = system.file("extdata", "test_symbols.gmt.txt", package = "scMAGeCK")
# Run scmageck_lr function
# By default, the result will be saved to the current working directory. 
lr_result <- scmageck_lr(BARCODE=BARCODE, RDS=RDS, SIGNATURE=SIGNATURE, LABEL='dox_scmageck_lr', 
                         NEGCTRL = 'NonTargetingControlGuideForHuman', PERMUTATION = 500)
lr_signature <- as.data.frame(lr_result)

#head(lr_score)