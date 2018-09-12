#!/bin/bash

indir=`ls fastq/`
#indir="Jurkat_TCR_stimulated_run1 Jurkat_TCR_stimulated_run2"
#indir="Jurkat_TCR_stimulated_run3 Jurkat_TCR_unstimulated_run1"
#indir="Jurkat_TCR_unstimulated_run2 Jurkat_TCR_unstimulated_run4"
indir="Jurkat_TCR_stimulated_run4 Jurkat_TCR_stimulated_run5 Jurkat_TCR_stimulated_run6 Jurkat_TCR_unstimulated_run5"

tmpdir=/scratch/wl948/tmp

for folder in $indir; do
  # echo $folder
  fq1=`ls fastq/$folder/*_1.fastq`
  fq2=`ls fastq/$folder/*_2.fastq`
  outf=fastq/$folder/output.bam
  sample="CROP-seq_${folder/run/r}"
  # command="picard FastqToSam F1=SRR5128078_1.fastq F2=SRR5128078_2.fastq O=SRR5128078.bam SM=CROP-seq_Jurkat_TCR_stimulated_r1"
  command="picard FastqToSam F1=$fq1 F2=$fq2 O=$outf SM=$sample TMP_DIR=$tmpdir"
  echo $command
  eval $command
done
