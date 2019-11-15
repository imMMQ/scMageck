#!/bin/bash


Rscript ../../tools/R/run_lr.R --RDS ../demo1_scmageck_R/singles_dox_mki67_v3.RDS --BARCODE ../demo1_scmageck_R/barcode_rec.txt  --SIGNATURE gene_signature.txt --LABEL dox_scmageck_lr --NEGCTRL NonTargetingControlGuideForHuman --PERMUTATION 1000
