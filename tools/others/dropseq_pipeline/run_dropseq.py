#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import re
import yaml
import argparse
import logging

        
def systemcall(command, cmsg=True,run=True):
  if cmsg:
    logging.info('SYSTEM CMD: '+command)
  if run:
    return os.system(command)
  else:
    return 0




class parameters:
  picard_jar="picard"
  star="STAR"
  cores=4
  
  cell_barcode_bases="1-12"
  umi_barcode_bases="13-20"
  min_base_quality=10
  min_bases_below_quality=1
  trim_sequence="AAGCAGTGGTATCAACGCAGAGTGAATGGG"
  trim_sequence_length=5
  polya_size=6
  min_genes_per_cell=[10,50,100,500]
  #  - 100
  #  # - 10
  repair_barcodes=False
  number_seq_error_barcodes_check=10000000
  bead_primer_sequence="AAGCAGTGGTATCAACGCAGAGTAC"
  distance_to_bead_primer_seq=0
  max_number_barcode_bases_to_repair=4
  # to be replaced by config file 
  dropseq_root="/home/wl948/datarun/wl948/project/crispr/cropseq/dropseqbin"
  star_index='/scratch/wl948/hg38_spiked_Tcrlibrary'
  refgenome='/scratch/wl948/hg38_spiked_Tcrlibrary/Homo_sapiens.GRCh38.dna.primary_assembly.spiked.fa'
  refflat='/scratch/wl948/hg38_spiked_Tcrlibrary/Homo_sapiens.GRCh38.dna.primary_assembly.spiked.refFlat'
  output_dir="outputbam"
  tmp_dir="/scratch/wl948/tmp"

def process_args():
    parser=argparse.ArgumentParser(description='Processing crop-seq data')
    parser.add_argument('-v', '--version',action='version',version='0.1')
    parser.add_argument('-c','--config-file',required=True,help='The config yaml file (required)')
    parser.add_argument('-n','--no-run',action='store_true',help='Do not actually run the steps. Just print the command.')
    parser.add_argument('-s','--step-only',choices=[-1,0,1,2,3,4,5,6,7,8,9,10,11,12],nargs='+',type=int,help='Only execute some of the steps. For example, use "-s 0 1 2" if you just want to run steps 0-2. By default, run all the steps.')
    
    args=parser.parse_args()
    config=read_config(args.config_file)
    # set up logging

    logmode="a"
    logging.basicConfig(level=10,
      format='%(levelname)-5s @ %(asctime)s: %(message)s ',
      datefmt='%a, %d %b %Y %H:%M:%S',
      # stream=sys.stderr,
      filename=os.path.join(config["paths"]["output_dir"],'run_dropseq.log'),
      filemode=logmode
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(levelname)-5s @ %(asctime)s: %(message)s ','%a, %d %b %Y %H:%M:%S')
    #formatter.formatTime('%a, %d %b %Y %H:%M:%S')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)
    logging.info('## ')
    logging.info('## ')
    logging.info('## ')
    logging.info('## Starting the pipeline.')
    logging.info('Parameters: '+' '.join(sys.argv))
  

  

    #
    if args.step_only==None:
      args.step_only=[0,1,2,3,4,5,6,7,8,9,10,11,12]
      logging.info('Running all the steps.')
    else:
      logging.info('Only run step(s): '+','.join([str(x) for x in args.step_only]))
    
    return (args,config)

  
    

def fastq2bam(config,args):
  fq1=config["fastq"]["fq1"]
  fq2=config["fastq"]["fq2"]
  #outf=fastq/$folder/output.bam
  outf=os.path.join(parameters.output_dir,"in.bam")
  sample=config["fastq"]["label"]
  # command="picard FastqToSam F1=SRR5128078_1.fastq F2=SRR5128078_2.fastq O=SRR5128078.bam SM=CROP-seq_Jurkat_TCR_stimulated_r1"
  cmd=parameters.picard_jar
  cmd+=" FastqToSam"
  cmd+=" F1="+fq1+" F2="+fq2
  cmd+=" O="+outf
  cmd+=" SM="+sample
  if "tmp_dir" in config["paths"]:
    cmd+=" TMP_DIR="+config["paths"]["tmp_dir"]
  return systemcall(cmd,run=not args.no_run)


def step1(input_file,output_prefix,args,config):
    # Stage 0: pre-alignment tag and trim
    # Tag with cell barcode
    logging.info("## Step 1: Tagging BAM file with cell barcode")
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
    return systemcall(cmd,run=not args.no_run)
    # pipe.clean_add(os.path.join(output_dir, "unaligned_tagged_Cell.bam"), manual=True)

def step2(output_prefix,args,config):
    # Tag with molecule barcode
    logging.info("## Step 2: Tagging BAM file with molecule barcode (UMI)")
    cmd = os.path.join(parameters.dropseq_root, "TagBamWithReadSequenceExtended")
    cmd += " TMP_DIR=" +parameters.tmp_dir 
    cmd += " SUMMARY=" + output_prefix + "_unaligned_tagged_Molecular.bam_summary.txt"
    cmd += " BASE_RANGE={}".format(parameters.umi_barcode_bases)
    cmd += " BASE_QUALITY={}".format(parameters.min_base_quality)
    cmd += " BARCODED_READ=1 DISCARD_READ=true TAG_NAME=XM NUM_BASES_BELOW_QUALITY={}".format(parameters.min_bases_below_quality)
    cmd += " INPUT=" + output_prefix+ "_unaligned_tagged_Cell.bam"
    cmd += " OUTPUT=" + output_prefix+ "_unaligned_tagged_CellMolecular.bam"
    #pipe.run(cmd, os.path.join(output_dir, "unaligned_tagged_CellMolecular.bam"))
    return systemcall(cmd,run=not args.no_run)


def step3(output_prefix,args,config):
    # Filter bam
    logging.info("## Step 3: Filtering BAM file")
    cmd = os.path.join(parameters.dropseq_root, "FilterBAM")
    cmd += " TAG_REJECT=XQ"
    cmd += " INPUT=" + output_prefix+ "_unaligned_tagged_CellMolecular.bam"
    cmd += " OUTPUT=" + output_prefix+ "_unaligned_tagged_filtered.bam"
    return systemcall(cmd,run=not args.no_run)
    # pipe.clean_add(os.path.join(output_dir, "unaligned_tagged_filtered.bam"), manual=True)

def step4(output_prefix,args,config):
    # Trim starting sequence
    logging.info("## Step 4: Triming starting sequence")
    cmd = os.path.join(parameters.dropseq_root, "TrimStartingSequence")
    cmd += " SEQUENCE={}".format(parameters.trim_sequence)
    cmd += " MISMATCHES=0 NUM_BASES={}".format(parameters.trim_sequence_length)
    cmd += " OUTPUT_SUMMARY=" + output_prefix+ "_adapter_trimming_report.txt"
    cmd += " INPUT=" + output_prefix+ "_unaligned_tagged_filtered.bam"
    cmd += " OUTPUT=" + output_prefix+"_unaligned_tagged_trimmed_smart.bam"
    #pipe.run(cmd, os.path.join(output_dir, "unaligned_tagged_trimmed_smart.bam"))
    return systemcall(cmd,run=not args.no_run)
    #pipe.clean_add(os.path.join(output_dir, "unaligned_tagged_trimmed_smart.bam"), manual=True)


def step5(output_prefix,args,config):
    # Trim polyA tail
    logging.info("## Step 5: Trimming polyA tail")
    cmd = os.path.join(parameters.dropseq_root, "PolyATrimmer")
    cmd += " MISMATCHES=0 NUM_BASES={}".format(parameters.polya_size)
    cmd += " OUTPUT_SUMMARY=" + output_prefix+ "_polyA_trimming_report.txt"
    cmd += " INPUT=" + output_prefix+ "_unaligned_tagged_trimmed_smart.bam"
    cmd += " OUTPUT=" + output_prefix+ "_unaligned_mc_tagged_polyA_filtered.bam"
    #pipe.run(cmd, os.path.join(output_dir, "unaligned_mc_tagged_polyA_filtered.bam"))
    #pipe.clean_add(os.path.join(output_dir, "unaligned_mc_tagged_polyA_filtered.bam"), manual=True)
    return systemcall(cmd,run=not args.no_run)

def step6(output_prefix,args,config):
    # Stage 2: alignment
    # Convert to fastq
    logging.info("## Step 6: Converting to Fastq")
    #cmd = "java -Xmx{}g -jar {} SamToFastq".format(int(args.mem) / 1000, paramters.picard_jar)
    cmd=parameters.picard_jar
    cmd += " SamToFastq"
    cmd += " INPUT=" + output_prefix+ "_unaligned_mc_tagged_polyA_filtered.bam"
    cmd += " FASTQ=" + output_prefix+ "_unaligned_mc_tagged_polyA_filtered.fastq"
    cmd += " TMP_DIR="+parameters.tmp_dir
    #pipe.run(cmd, os.path.join(output_dir, "unaligned_mc_tagged_polyA_filtered.fastq"))
    #pipe.clean_add(os.path.join(output_dir, "unaligned_mc_tagged_polyA_filtered.fastq"), manual=True)
    return systemcall(cmd,run=not args.no_run)


def step7(output_prefix,args,config):
    # Align reads
    logging.info("## Step 7: Aligning reads with STAR")
    cmd = parameters.star
    cmd += " --genomeDir {}".format(parameters.star_index)
    cmd += " --runThreadN {}".format(parameters.cores)
    cmd += " --outFileNamePrefix " + output_prefix+ "_star."
    cmd += " --readFilesIn " + output_prefix+ "_unaligned_mc_tagged_polyA_filtered.fastq"
    if "tmp_dir" in config["paths"]:
        tmpdir=os.path.join(config["paths"]["tmp_dir"],'STARtmp')
        if os.path.exists(tmpdir):
            print('Error: STAR tmp already exists. Please manually remove it:'+tmpdir)
            sys.exit(-1)
        cmd += " --outTmpDir {}".format(tmpdir)
    #pipe.run(cmd, os.path.join(output_dir, "star.Aligned.out.sam"))
    #pipe.clean_add(os.path.join(output_dir, "star.Aligned.out.sam"), manual=True)
    return systemcall(cmd,run=not args.no_run)

def step8(output_prefix,args,config):
    # Stage 3: sort aligned reads (STAR does not necessarily emit reads in the same order as the input)
    logging.info("## Step 8: Sorting aligned BAM file")
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
    return systemcall(cmd,run=not args.no_run)

def step9(output_prefix,args,config):
    # Stage 4: merge and tag aligned reads
    # Merge
    logging.info("## Step 9: Merging aligned with unaligned reads")
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
    return systemcall(cmd,run=not args.no_run)

def step10(output_prefix,args,config):
    # Tag reads with exon
    logging.info("## Step 10: Tagging reads with exon")
    cmd = os.path.join(parameters.dropseq_root, "TagReadWithGeneExon")
    cmd += " OUTPUT=" + output_prefix+ "_star_gene_exon_tagged.bam"
    #cmd += " ANNOTATIONS_FILE={}".format(getattr(pipe.config.resources.refflat, sample.genome))
    cmd += " ANNOTATIONS_FILE={}".format(parameters.refflat)
    cmd += " TAG=GE CREATE_INDEX=true"
    cmd += " INPUT=" + output_prefix+"_merged.bam"
    return systemcall(cmd,run=not args.no_run)
    #pipe.run(cmd, os.path.join(output_dir, "star_gene_exon_tagged.bam"))

def step11(output_prefix,args,config):
    logging.info("## Step 11")
    if parameters.repair_barcodes:
        # Detect and fix bead synthesis errors
        logging.info("## Reporting and fixing bead synthesis errors")
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
        if systemcall(cmd,run=not args.no_run)!=0:
          return -1
        bam_file = output_prefix+ "_star_gene_exon_tagged.clean.bam"
    else:
        bam_file = output_prefix+ "_star_gene_exon_tagged.bam"
    
    
    # Distribution of read quality
    # cell barcode
    logging.info("## Read quality in cell barcodes")
    cmd = os.path.join(parameters.dropseq_root, "GatherReadQualityMetrics")
    cmd += " INPUT=" + bam_file
    cmd += " OUTPUT=" + output_prefix+ "_quality_distribution.cell_barcode.txt"
    cmd += " TAG=XC"
    #pipe.run(cmd, os.path.join(output_dir, "quality_distribution.cell_barcode.txt"))
    if systemcall(cmd,run=not args.no_run)!=0:
      return -1
    # UMI
    logging.info("## Read quality in molecule barcodes")
    cmd = os.path.join(parameters.dropseq_root, "GatherReadQualityMetrics")
    cmd += " INPUT=" + bam_file
    cmd += " OUTPUT=" + output_prefix+ "_quality_distribution.mol_barcode.txt"
    cmd += " TAG=XM"
    #pipe.run(cmd, os.path.join(output_dir, "quality_distribution.mol_barcode.txt")) 
    if systemcall(cmd,run=not args.no_run)!=0:
      return -1
    
    # Distribution of bases in reads
    # cell barcode
    logging.info("## Distribution of bases in cell barcodes")
    cmd = os.path.join(parameters.dropseq_root, "BaseDistributionAtReadPosition")
    cmd += " INPUT=" + bam_file
    cmd += " OUTPUT=" + output_prefix+ "_base_distribution.cell_barcode.txt"
    cmd += " TAG=XC"
    #pipe.run(cmd, os.path.join(output_dir, "base_distribution.cell_barcode.txt"))   
    if systemcall(cmd,run=not args.no_run)!=0:
      return -1
    # UMI
    logging.info("## Distribution of bases in molecule barcodes (UMI)")
    cmd = os.path.join(parameters.dropseq_root, "BaseDistributionAtReadPosition")
    cmd += " INPUT=" + bam_file
    cmd += " OUTPUT=" + output_prefix+ "_base_distribution.mol_barcode.txt"
    cmd += " TAG=XM"
    #pipe.run(cmd, os.path.join(output_dir, "base_distribution.mol_barcode.txt"))
    if systemcall(cmd,run=not args.no_run) !=0:
      return -1
    
    # Reads per cell summary
    logging.info("## Reporting summary of reads per cell")
    cmd = os.path.join(parameters.dropseq_root, "BAMTagHistogram")
    cmd += " INPUT=" + bam_file
    cmd += " OUTPUT=" + output_prefix+ "_cell_readcounts.txt"
    cmd += " FILTER_PCR_DUPLICATES=true"
    cmd += " TAG=XC"
    #pipe.run(cmd, os.path.join(output_dir, "cell_readcounts.txt"))
    return systemcall(cmd,run=not args.no_run)



def step12(output_prefix,args,config):

    logging.info("## Step 12")
    if parameters.repair_barcodes:
        bam_file = output_prefix+ "_star_gene_exon_tagged.clean.bam"
    else:
        bam_file = output_prefix+ "_star_gene_exon_tagged.bam"

    for n_genes in parameters.min_genes_per_cell:
        logging.info("## Perform digital gene expression analysis for cells with at least {} genes covered".format(n_genes))
        cmd = os.path.join(parameters.dropseq_root, "DigitalExpression")
        #cmd += " -m {}g".format(int(args.mem) / 1000)
        cmd += " TMP_DIR=" +parameters.tmp_dir 
        cmd += " INPUT=" + bam_file
        cmd += " OUTPUT=" + output_prefix+ "_digital_expression.{}genes.tsv".format(n_genes)
        cmd += " SUMMARY=" + output_prefix+ "_digital_expression.summary.{}genes.tsv".format(n_genes)
        cmd += " MIN_NUM_GENES_PER_CELL={}".format(n_genes)
        #pipe.run(cmd, os.path.join(output_dir, "digital_expression.{}genes.tsv".format(n_genes)), nofail=True)
        if systemcall(cmd,run=not args.no_run)!=0:
          return -1

    # Report how often the same UMI is found per cell per gene --> estimate of PCR duplicates
    for n_genes in parameters.min_genes_per_cell:
        logging.info("## Report UMI count per cell per gene for cells with at least {} genes covered".format(n_genes))
        cmd = os.path.join(parameters.dropseq_root, "GatherMolecularBarcodeDistributionByGene")
        #cmd += " -m {}g".format(int(args.mem) / 1000)
        cmd += " TMP_DIR=" +parameters.tmp_dir 
        cmd += " INPUT=" + bam_file
        cmd += " OUTPUT=" + output_prefix+ "cell_umi_barcodes.{}genes.tsv".format(n_genes)
        cmd += " MIN_NUM_GENES_PER_CELL={}".format(n_genes)
        #pipe.run(cmd, os.path.join(output_dir, "cell_umi_barcodes.{}genes.tsv".format(n_genes)))
        if systemcall(cmd,run=not args.no_run)!=0:
          return -1
    return 0


    

def read_config(configfile):
    config=yaml.load(open(configfile))
    parameters.dropseq_root=""
    parameters.star_index=os.path.join(config["paths"]["output_dir"],"spiked_genomes")
    parameters.refgenome=os.path.join(config["paths"]["output_dir"],"spiked_genomes","spiked.fa")
    parameters.refflat=os.path.join(config["paths"]["output_dir"],"spiked_genomes","spiked.refFlat")
    parameters.output_dir=os.path.join(config["paths"]["output_dir"],"outputbam")
    if not os.path.exists(parameters.output_dir):
        os.makedirs(parameters.output_dir)
    parameters.tmp_dir=config["paths"]["tmp_dir"]
    return config



if __name__ == '__main__':
  try:
    (args,config)=process_args()
    input_file=os.path.join(parameters.output_dir,"in.bam")
    output_prefix=os.path.join(parameters.output_dir,"out")
    # convert fastq to bam
    if 0 in args.step_only:
      if fastq2bam(config,args)!=0:
        logging.error('Error converting fastq to bam (step 0).')
        sys.exit(-1)
    #input_file=sys.argv[1]
    #output_prefix=sys.argv[2]
    logging.info('Input:'+input_file)
    logging.info('Output_prefix:'+output_prefix)
    if 1 in args.step_only:
      if step1(input_file,output_prefix,args,config)!=0:
        logging.error('Error in step 1.')
        sys.exit(-1)
    if 2 in args.step_only:
      if step2(output_prefix,args,config)!=0:
        logging.error('Error in step 2.')
        sys.exit(-1)
    if 3 in args.step_only:
      if step3(output_prefix,args,config)!=0:
        logging.error('Error in step 3.')
        sys.exit(-1)
    if 4 in args.step_only:
      if step4(output_prefix,args,config)!=0:
        logging.error('Error in step 4.')
        sys.exit(-1)
    if 5 in args.step_only:
      if step5(output_prefix,args,config)!=0:
        logging.error('Error in step 5.')
        sys.exit(-1)
    if 6 in args.step_only:
      if step6(output_prefix,args,config)!=0:
        logging.error('Error in step 6.')
        sys.exit(-1)
    if 7 in args.step_only:
      if step7(output_prefix,args,config)!=0:
        logging.error('Error in step 7.')
        sys.exit(-1)
    if 8 in args.step_only:
      if step8(output_prefix,args,config)!=0:
        logging.error('Error in step 8.')
        sys.exit(-1)
    if 9 in args.step_only:
      if step9(output_prefix,args,config)!=0:
        logging.error('Error in step 9.')
        sys.exit(-1)
    if 10 in args.step_only:
      if step10(output_prefix,args,config)!=0:
        logging.error('Error in step 10.')
        sys.exit(-1)
    if 11 in args.step_only:
      if step11(output_prefix,args,config)!=0:
        logging.error('Error in step 11.')
        sys.exit(-1)
    if 12 in args.step_only:
      if step12(output_prefix,args,config)!=0:
        logging.error('Error in step 12.')
        sys.exit(-1)
  except KeyboardInterrupt:
    sys.stderr.write("Interrupted.\n")
    sys.exit(0)


