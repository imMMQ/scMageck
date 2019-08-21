library(scmageck)

# set the BARCODE and RDS file path 
BARCODE = system.file("extdata","barcode_rec.txt",package = "scmageck")
RDS = system.file("extdata","singles_dox_mki67_v3.RDS",package = "scmageck")

# Run scmageck_rra function
# By default, the result will be saved to the current working directory. 
rra_result <- scmageck_rra(BARCODE=BARCODE, RDS=RDS, GENE="MKI67",   
                           LABEL='dox_mki67', NEGCTRL=NULL, KEEPTMP=FALSE, PATHWAY=FALSE)
#head(rra_result)
