#!/bin/bash


Rscript ../../tools/R/run_lr.R --RDS singles_dox_mki67.RDS   --BARCODE barcode_rec.txt   --LABEL dox_scmageck_lr --NEGCTRL NonTargetingControlGuideForHuman --PERMUTATION 1000
