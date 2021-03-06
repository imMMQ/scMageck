library(scMAGeCK)

# set the BARCODE and RDS file path 
BARCODE = system.file("extdata","barcode_rec.txt",package = "scMAGeCK")
## RDS can be a Seurat object or local RDS file path that contains the scRNA-seq dataset
RDS = system.file("extdata","singles_dox_mki67_v3.RDS",package = "scMAGeCK")

# Run scmageck_lr function
# By default, the result will be saved to the current working directory. 
lr_result <- scmageck_lr(BARCODE=BARCODE, RDS=RDS, LABEL='dox_scmageck_lr', 
                         NEGCTRL = 'NonTargetingControlGuideForHuman', PERMUTATION = 1000)
lr_score <- lr_result[1]
lr_score_pval <- lr_result[2]

#head(lr_score)
