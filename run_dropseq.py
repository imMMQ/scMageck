#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import re

class parameters:
  dropseq_root="/home/wl948/datarun/wl948/project/crispr/cropseq/dropseqbin"
  picard_jar="picard"
  star="STAR"
  star_index='/scratch/wl948/hg38_spiked_Tcrlibrary'
  cores=4
  refgenome='/scratch/wl948/hg38_spiked_Tcrlibrary/Homo_sapiens.GRCh38.dna.primary_assembly.spiked.fa'
  refflat='/scratch/wl948/hg38_spiked_Tcrlibrary/Homo_sapiens.GRCh38.dna.primary_assembly.spiked.refFlat'
  output_dir="outputbam"
  tmp_dir="/scratch/wl948/tmp"
  cell_barcode_bases="1-12"
  umi_barcode_bases="13-20"
  min_base_quality=10
  min_bases_below_quality=1
  trim_sequence="AAGCAGTGGTATCAACGCAGAGTGAATGGG"
  trim_sequence_length=5
  polya_size=6
  min_genes_per_cell=[500,100,10]
  #  - 100
  #  # - 10
  repair_barcodes=True
  number_seq_error_barcodes_check=10000000
  bead_primer_sequence="AAGCAGTGGTATCAACGCAGAGTAC"
  distance_to_bead_primer_seq=0
  max_number_barcode_bases_to_repair=4

def systemcall(cmd):
    print(cmd)
    os.system(cmd)
    


def step1(input_file,output_prefix):
    # Stage 1: pre-alignment tag and trim
    # Tag with cell barcode
    print("## Step 1: Tagging BAM file with cell barcode")
    cmd = os.path.join(parameters.dropseq_root, "TagBamWithReadSequenceExtended")
    cmd += " TMP_DIR=" + parameters.tmp_dir
    cmd += " BASE_RANGE={}".format(parameters.cell_barcode_bases)
    cmd += " BASE_QUALITY={}".format(parameters.min_base_quality)
    cmd += " BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XC NUM_BASES_BELOW_QUALITY={}".format(parameters.min_bases_below_quality)
    cmd += " INPUT=" + input_file
    #cmd += " SUMMARY=" + os.path.join(output_dir, "unaligned_tagged_Cellular.bam_summary.txt")
    #cmd += " OUTPUT=" + os.path.join(output_dir, "unaligned_tagged_Cell.bam")
    cmd += " SUMMARY=" + output_prefix + "_unaligned_tagged_Cellular.bam_summary.txt"
    cmd += " OUTPUT=" + output_prefix+ "_unaligned_tagged_Cell.bam"
    systemcall(cmd)
    # pipe.clean_add(os.path.join(output_dir, "unaligned_tagged_Cell.bam"), manual=True)

def step1_5(output_prefix):
    # Tag with molecule barcode
    print("## Step 1.5: Tagging BAM file with molecule barcode (UMI)")
    cmd = os.path.join(parameters.dropseq_root, "TagBamWithReadSequenceExtended")
    cmd += " TMP_DIR=" +parameters.tmp_dir 
    cmd += " SUMMARY=" + output_prefix + "_unaligned_tagged_Molecular.bam_summary.txt"
    cmd += " BASE_RANGE={}".format(parameters.umi_barcode_bases)
    cmd += " BASE_QUALITY={}".format(parameters.min_base_quality)
    cmd += " BARCODED_READ=1 DISCARD_READ=true TAG_NAME=XM NUM_BASES_BELOW_QUALITY={}".format(parameters.min_bases_below_quality)
    cmd += " INPUT=" + output_prefix+ "_unaligned_tagged_Cell.bam"
    cmd += " OUTPUT=" + output_prefix+ "_unaligned_tagged_CellMolecular.bam"
    #pipe.run(cmd, os.path.join(output_dir, "unaligned_tagged_CellMolecular.bam"))
    systemcall(cmd)


def step2(output_prefix):
    # Filter bam
    print("## Step 2: Filtering BAM file")
    cmd = os.path.join(parameters.dropseq_root, "FilterBAM")
    cmd += " TAG_REJECT=XQ"
    cmd += " INPUT=" + output_prefix+ "_unaligned_tagged_CellMolecular.bam"
    cmd += " OUTPUT=" + output_prefix+ "_unaligned_tagged_filtered.bam"
    systemcall(cmd)
    # pipe.clean_add(os.path.join(output_dir, "unaligned_tagged_filtered.bam"), manual=True)

def step3(output_prefix):
    # Trim starting sequence
    print("## Step 3: Triming starting sequence")
    cmd = os.path.join(parameters.dropseq_root, "TrimStartingSequence")
    cmd += " SEQUENCE={}".format(parameters.trim_sequence)
    cmd += " MISMATCHES=0 NUM_BASES={}".format(parameters.trim_sequence_length)
    cmd += " OUTPUT_SUMMARY=" + output_prefix+ "_adapter_trimming_report.txt"
    cmd += " INPUT=" + output_prefix+ "_unaligned_tagged_filtered.bam"
    cmd += " OUTPUT=" + output_prefix+"_unaligned_tagged_trimmed_smart.bam"
    #pipe.run(cmd, os.path.join(output_dir, "unaligned_tagged_trimmed_smart.bam"))
    systemcall(cmd)
    #pipe.clean_add(os.path.join(output_dir, "unaligned_tagged_trimmed_smart.bam"), manual=True)


def step4(output_prefix):
    # Trim polyA tail
    print("## Step 4: Trimming polyA tail")
    cmd = os.path.join(parameters.dropseq_root, "PolyATrimmer")
    cmd += " MISMATCHES=0 NUM_BASES={}".format(parameters.polya_size)
    cmd += " OUTPUT_SUMMARY=" + output_prefix+ "_polyA_trimming_report.txt"
    cmd += " INPUT=" + output_prefix+ "_unaligned_tagged_trimmed_smart.bam"
    cmd += " OUTPUT=" + output_prefix+ "_unaligned_mc_tagged_polyA_filtered.bam"
    #pipe.run(cmd, os.path.join(output_dir, "unaligned_mc_tagged_polyA_filtered.bam"))
    #pipe.clean_add(os.path.join(output_dir, "unaligned_mc_tagged_polyA_filtered.bam"), manual=True)
    systemcall(cmd)

def step5(output_prefix):
    # Stage 2: alignment
    # Convert to fastq
    print("## Step 5: Converting to Fastq")
    #cmd = "java -Xmx{}g -jar {} SamToFastq".format(int(args.mem) / 1000, paramters.picard_jar)
    cmd=parameters.picard_jar
    cmd += " SamToFastq"
    cmd += " INPUT=" + output_prefix+ "_unaligned_mc_tagged_polyA_filtered.bam"
    cmd += " FASTQ=" + output_prefix+ "_unaligned_mc_tagged_polyA_filtered.fastq"
    cmd += " TMP_DIR="+parameters.tmp_dir
    #pipe.run(cmd, os.path.join(output_dir, "unaligned_mc_tagged_polyA_filtered.fastq"))
    #pipe.clean_add(os.path.join(output_dir, "unaligned_mc_tagged_polyA_filtered.fastq"), manual=True)
    systemcall(cmd)


def step6(output_prefix):
    # Align reads
    print("## Step 6: Aligning reads with STAR")
    cmd = parameters.star
    cmd += " --genomeDir {}".format(parameters.star_index)
    cmd += " --runThreadN {}".format(parameters.cores)
    cmd += " --outFileNamePrefix " + output_prefix+ "_star."
    cmd += " --readFilesIn " + output_prefix+ "_unaligned_mc_tagged_polyA_filtered.fastq"
    #pipe.run(cmd, os.path.join(output_dir, "star.Aligned.out.sam"))
    #pipe.clean_add(os.path.join(output_dir, "star.Aligned.out.sam"), manual=True)
    systemcall(cmd)

def step7(output_prefix):
    # Stage 3: sort aligned reads (STAR does not necessarily emit reads in the same order as the input)
    print("## Step 7: Sorting aligned BAM file")
    #cmd = "java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx{}g".format(int(args.mem) / 1000)
    #cmd += " -jar {} SortSam".format(pipe.config.tools.piccard_jar)
    cmd=""
    cmd += "{} SortSam".format(parameters.picard_jar)
    cmd += " INPUT=" + output_prefix+ "_star.Aligned.out.sam"
    cmd += " OUTPUT=" + output_prefix+ "_aligned.sorted.bam"
    cmd += " SORT_ORDER=queryname"
    cmd += " TMP_DIR=" + parameters.tmp_dir
    #pipe.run(cmd, os.path.join(output_dir, "aligned.sorted.bam"))
    #pipe.clean_add(os.path.join(output_dir, "aligned.sorted.bam"), manual=True)
    systemcall(cmd)

def step8(output_prefix):
    # Stage 4: merge and tag aligned reads
    # Merge
    print("## Step 8: Merging aligned with unaligned reads")
    # cmd = "java -Djava.io.tmpdir={} -Xmx{}g -jar {} MergeBamAlignment".format(output_dir, int(args.mem) / 1000, pipe.config.tools.piccard_jar)
    cmd = "{} MergeBamAlignment".format(parameters.picard_jar)
    cmd += " REFERENCE_SEQUENCE={}".format( parameters.refgenome)
    cmd += " UNMAPPED_BAM=" + output_prefix+ "_unaligned_mc_tagged_polyA_filtered.bam"
    cmd += " ALIGNED_BAM=" + output_prefix+ "_aligned.sorted.bam"
    cmd += " INCLUDE_SECONDARY_ALIGNMENTS=false"
    cmd += " ALIGNED_READS_ONLY=false"
    cmd += " PAIRED_RUN=false"
    cmd += " OUTPUT=" + output_prefix+ "_merged.bam"
    cmd += " TMP_DIR=" + parameters.tmp_dir
    #pipe.run(cmd, os.path.join(output_dir, "merged.bam"))
    #pipe.clean_add(os.path.join(output_dir, "merged.bam"), manual=True)
    systemcall(cmd)

def step9(output_prefix):
    # Tag reads with exon
    print("## Step 9: Tagging reads with exon")
    cmd = os.path.join(parameters.dropseq_root, "TagReadWithGeneExon")
    cmd += " OUTPUT=" + output_prefix+ "_star_gene_exon_tagged.bam"
    #cmd += " ANNOTATIONS_FILE={}".format(getattr(pipe.config.resources.refflat, sample.genome))
    cmd += " ANNOTATIONS_FILE={}".format(parameters.refflat)
    cmd += " TAG=GE CREATE_INDEX=true"
    cmd += " INPUT=" + output_prefix+"_merged.bam"
    systemcall(cmd)
    #pipe.run(cmd, os.path.join(output_dir, "star_gene_exon_tagged.bam"))

def step10(output_prefix):
    if parameters.repair_barcodes:
        # Detect and fix bead synthesis errors
        print("## Reporting and fixing bead synthesis errors")
        cmd = os.path.join(parameters.dropseq_root, "DetectBeadSynthesisErrors")
        cmd += " INPUT=" + output_prefix+ "_star_gene_exon_tagged.bam"
        cmd += " OUTPUT=" + output_prefix+ "_star_gene_exon_tagged.clean.bam"
        cmd += " OUTPUT_STATS=" + output_prefix+ "_synthesis_statistics.txt"
        cmd += " SUMMARY=" + output_prefix+ "_synthesis_statistics.summary.txt"
        cmd += " NUM_BARCODES={}".format(parameters.number_seq_error_barcodes_check)
        cmd += " PRIMER_SEQUENCE={}".format(parameters.bead_primer_sequence)
        cmd += " EDIT_DISTANCE={}".format(parameters.distance_to_bead_primer_seq)
        cmd += " MAX_NUM_ERRORS={}".format(parameters.max_number_barcode_bases_to_repair)
        cmd += " TMP_DIR=" + parameters.tmp_dir
        #pipe.run(cmd, os.path.join(output_dir, "star_gene_exon_tagged.clean.bam"))
        systemcall(cmd)
        bam_file = output_prefix+ "_star_gene_exon_tagged.clean.bam"
    else:
        bam_file = output_prefix+ "_star_gene_exon_tagged.bam"
    
    
    # Distribution of read quality
    # cell barcode
    print("## Read quality in cell barcodes")
    cmd = os.path.join(parameters.dropseq_root, "GatherReadQualityMetrics")
    cmd += " INPUT=" + bam_file
    cmd += " OUTPUT=" + output_prefix+ "_quality_distribution.cell_barcode.txt"
    cmd += " TAG=XC"
    #pipe.run(cmd, os.path.join(output_dir, "quality_distribution.cell_barcode.txt"))
    systemcall(cmd)
    # UMI
    print("## Read quality in molecule barcodes")
    cmd = os.path.join(parameters.dropseq_root, "GatherReadQualityMetrics")
    cmd += " INPUT=" + bam_file
    cmd += " OUTPUT=" + output_prefix+ "_quality_distribution.mol_barcode.txt"
    cmd += " TAG=XM"
    #pipe.run(cmd, os.path.join(output_dir, "quality_distribution.mol_barcode.txt")) 
    systemcall(cmd)
    
    # Distribution of bases in reads
    # cell barcode
    print("## Distribution of bases in cell barcodes")
    cmd = os.path.join(parameters.dropseq_root, "BaseDistributionAtReadPosition")
    cmd += " INPUT=" + bam_file
    cmd += " OUTPUT=" + output_prefix+ "_base_distribution.cell_barcode.txt"
    cmd += " TAG=XC"
    #pipe.run(cmd, os.path.join(output_dir, "base_distribution.cell_barcode.txt"))   
    systemcall(cmd)
    # UMI
    print("## Distribution of bases in molecule barcodes (UMI)")
    cmd = os.path.join(parameters.dropseq_root, "BaseDistributionAtReadPosition")
    cmd += " INPUT=" + bam_file
    cmd += " OUTPUT=" + output_prefix+ "_base_distribution.mol_barcode.txt"
    cmd += " TAG=XM"
    #pipe.run(cmd, os.path.join(output_dir, "base_distribution.mol_barcode.txt"))
    systemcall(cmd)
    
    # Reads per cell summary
    print("## Reporting summary of reads per cell")
    cmd = os.path.join(parameters.dropseq_root, "BAMTagHistogram")
    cmd += " INPUT=" + bam_file
    cmd += " OUTPUT=" + output_prefix+ "_cell_readcounts.txt"
    cmd += " FILTER_PCR_DUPLICATES=true"
    cmd += " TAG=XC"
    #pipe.run(cmd, os.path.join(output_dir, "cell_readcounts.txt"))
    systemcall(cmd)



def step11(output_prefix):

    if parameters.repair_barcodes:
        bam_file = output_prefix+ "_star_gene_exon_tagged.clean.bam"
    else:
        bam_file = output_prefix+ "_star_gene_exon_tagged.bam"

    for n_genes in parameters.min_genes_per_cell:
        print("## Perform digital gene expression analysis for cells with at least {} genes covered".format(n_genes))
        cmd = os.path.join(parameters.dropseq_root, "DigitalExpression")
        #cmd += " -m {}g".format(int(args.mem) / 1000)
        cmd += " TMP_DIR=" +parameters.tmp_dir 
        cmd += " INPUT=" + bam_file
        cmd += " OUTPUT=" + output_prefix+ "_digital_expression.{}genes.tsv".format(n_genes)
        cmd += " SUMMARY=" + output_prefix+ "_digital_expression.summary.{}genes.tsv".format(n_genes)
        cmd += " MIN_NUM_GENES_PER_CELL={}".format(n_genes)
        #pipe.run(cmd, os.path.join(output_dir, "digital_expression.{}genes.tsv".format(n_genes)), nofail=True)
        systemcall(cmd)

    # Report how often the same UMI is found per cell per gene --> estimate of PCR duplicates
    for n_genes in parameters.min_genes_per_cell:
        print("## Report UMI count per cell per gene for cells with at least {} genes covered".format(n_genes))
        cmd = os.path.join(parameters.dropseq_root, "GatherMolecularBarcodeDistributionByGene")
        #cmd += " -m {}g".format(int(args.mem) / 1000)
        cmd += " TMP_DIR=" +parameters.tmp_dir 
        cmd += " INPUT=" + bam_file
        cmd += " OUTPUT=" + output_prefix+ "cell_umi_barcodes.{}genes.tsv".format(n_genes)
        cmd += " MIN_NUM_GENES_PER_CELL={}".format(n_genes)
        #pipe.run(cmd, os.path.join(output_dir, "cell_umi_barcodes.{}genes.tsv".format(n_genes)))
        systemcall(cmd)


    




if __name__ == '__main__':
  try:
    input_file=sys.argv[1]
    output_prefix=sys.argv[2]
    print('Input:'+input_file)
    print('Output_prefix:'+output_prefix)
    #step1(input_file,output_prefix)
    #step1_5(output_prefix)
    #step2(output_prefix)
    #step3(output_prefix)
    #step4(output_prefix)
    #step5(output_prefix)
    step6(output_prefix)
    step7(output_prefix)
    step8(output_prefix)
    step9(output_prefix)
    step10(output_prefix)
    step11(output_prefix)
  except KeyboardInterrupt:
    sys.stderr.write("Interrupted.\n")
    sys.exit(0)


