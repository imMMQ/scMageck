library(scmageck)
library(R.utils)
# set the BARCODE and RDS file path 
BARCODE = system.file("extdata","barcode_rec.txt",package = "scmageck")
## RDS can be a Seurat object or local RDS file path that contains the scRNA-seq dataset
RDS = system.file("extdata","singles_dox_mki67_v3.RDS",package = "scmageck")
## set the gene signature file path 
SIGNATURE = system.file("extdata", "test_symbols.gmt.txt", package = "scmageck")

# Run scmageck_rra function
# By default, the result will be saved to the current working directory. 
scmageck_rra(BARCODE=BARCODE, RDS=RDS, SIGNATURE = SIGNATURE, NEGCTRL="NonTargetingControlGuideForHuman", KEEPTMP = FALSE)

