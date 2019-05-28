#!/bin/bash


Rscript ../../tools/R/run_lr.R --RDS singles_dox_mki67_v3.RDS  --LIBRARY sgrna_lib.txt  --BARCODE barcode_rec.txt   --LABEL dox_scmageck_lr --NEGCTRL NonTargetingControlGuideForHuman --PERMUTATION 1000
