library(scMAGeCK)

# set the BARCODE and RDS file path 
BARCODE = system.file("extdata","barcode_rec.txt",package = "scMAGeCK")
## RDS can be a Seurat object or local RDS file path that contains the scRNA-seq dataset
RDS = system.file("extdata","singles_dox_mki67_v3.RDS",package = "scMAGeCK")

# Run scmageck_rra function
# By default, the result will be saved to the current working directory. 
rra_result <- scmageck_rra(BARCODE=BARCODE, RDS=RDS, GENE="MKI67",   
                           LABEL='dox_mki67', NEGCTRL="NonTargetingControlGuideForHuman", KEEPTMP=FALSE, PATHWAY=FALSE)
#head(rra_result)
