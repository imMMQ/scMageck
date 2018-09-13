#!/usr/bin/env python

#from looper.models import Project
import os
import sys
#import pandas as pd
import yaml
import re

def systemcall(command, cmsg=True):
  if cmsg:
    print('SYSTEM CMD: '+command)
    os.system(command)


def write_annotation(sgannotation, config, output_fasta, output_gtf, cas9=True):
    # create fasta and gtf entries for each gRNA
    fasta_entries = []
    gtf_entries = []
    sglist=[]
    
    nline=0
    for line in open(sgannotation):
        nline+=1
        field=line.strip().split()
        if len(field)==1:
            field=line.strip(',').split()
        if len(field)<3:
            print('Error: cannot parse annotation file at line '+str(nline))
        if len(re.findall('[^ATGC]',field[1].upper()))>0:
            print('Warning: non sequence field found at line '+str(nline))
        else:
            sglist+=[field]

    for sgentry in sglist:
        oligo_name = sgentry[0]
        guide_sequence =sgentry[1] 

        # add fasta entry
        sequence = config['crop-seq']['u6'] + guide_sequence + config['crop-seq']['rest']
        header = fasta_header_template.format(chrom=oligo_name, length=len(sequence))
        fasta_entries.append(header)
        fasta_entries.append(sequence)

        # add gtf entry
        gtf_entries.append(gtf_template.format(chrom=oligo_name, id=oligo_name, length=len(sequence)))

    if cas9:
        # add Cas9 vector
        seq_name = "Cas9_blast"
        sequence = "".join([
            config['crop-seq']['cas9'],
            config['crop-seq']['nls'],
            config['crop-seq']['flag'],
            config['crop-seq']['p2a'],
            config['crop-seq']['blast'],
            config['crop-seq']['space'],
            config['crop-seq']['virus_ltr']])
        header = fasta_header_template.format(chrom=seq_name, length=len(sequence))
        # add cas9 entry
        fasta_entries.append(header)
        fasta_entries.append(sequence)
        gtf_entries.append(gtf_template.format(chrom=seq_name, id=seq_name, length=len(sequence)))

    # write to file
    with open(output_fasta, "w") as fasta_handle:
        fasta_handle.writelines("\n".join(fasta_entries))
    with open(output_gtf, "w") as gtf_handle:
        gtf_handle.writelines(gtf_entries)


fasta_header_template = ">{chrom}_chrom dna:chromosome chromosome:GRCh38:{chrom}_chrom:1:{length}:1 REF"

gtf_template = """{chrom}_chrom\thavana\tgene\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA";
{chrom}_chrom\thavana\ttranscript\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; transcript_id "{id}_transcript"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; transcript_name "{id}_transcript"; transcript_source "havana";
{chrom}_chrom\thavana\texon\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; transcript_id "{id}_transcript"; exon_number "1"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; transcript_name "{id}_transcript"; transcript_source "havana"; exon_id "{id}_exon";
"""

# initialize project
#prj = Project(os.path.join("metadata", "config.yaml"))

# read in guide annotation
#annotation = pd.read_csv(os.path.join("metadata", "guide_annotation.csv"))

def read_config(configfile):
    config=yaml.load(open(configfile))
    return config
  

def generate_annotations(config,sgrna_annotation):
        prefix=config["paths"]["output_dir"]
        output_dir = os.path.join(prefix, "spiked_genomes")
    
        #spiked_dir = os.path.join(output_dir, library) 
        spiked_dir=output_dir
    
        for _dir in [prefix, spiked_dir]:
            if not os.path.exists(_dir):
                os.makedirs(_dir)
    
        output_fasta = os.path.join(spiked_dir, "gRNA_spikes.fa")
        output_gtf = os.path.join(spiked_dir, "gRNA_spikes.gtf")
        # write gRNA library annotation
        write_annotation(sgrna_annotation, config, output_fasta, output_gtf)
        build_star_index(config,output_fasta,output_gtf,spiked_dir)
    
def get_hg38_ref():
    # Get fasta genome
    os.system(
        "wget -O {} ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        .format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")))
    os.system("gzip -d {}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")))
    # Get gtf annotation
    os.system(
        "wget -O {} ftp://ftp.ensembl.org/pub/release-77/gtf/homo_sapiens/Homo_sapiens.GRCh38.77.gtf.gz"
        .format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.gtf.gz")))
    os.system("gzip -d {}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.gtf.gz")))


def build_star_index(config,output_fasta,output_gtf,spiked_dir):
    reference_fa=config['genomes']['fa']
    #os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.fa")
    #os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.fa")
    final_fa=os.path.join(spiked_dir,'spiked.fa')
    
    reference_gtf=config['genomes']['gtf']
    final_gtf=os.path.join(spiked_dir,'spiked.gtf')
    #os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.gtf")
    #os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.spiked.gtf")
    
    # Add extra chromosomes (CROP-seq constructs) to genome
    systemcall("cat {} {} > {}".format(reference_fa, output_fasta, final_fa))
    systemcall("cat {} {} > {}".format(reference_gtf, output_gtf, final_gtf ))

    # Build STAR index (contruct + spiked with gRNAs)
    # cmd = "srun --mem 80000 -p develop /cm/shared/apps/star/2.4.2a/STAR"
    cmd = "STAR"
    cmd += " --runThreadN 8"
    cmd += " --runMode genomeGenerate"
    cmd += " --genomeDir {}".format(spiked_dir)
    cmd += " --genomeFastaFiles {}".format(final_fa)
    cmd += " --sjdbGTFfile {}".format(final_gtf)
    cmd += " --sjdbOverhang 74"
    if "tmp_dir" in config["paths"]:
        cmd += " --outTmpDir {}".format(config["paths"]["tmp_dir"])
    systemcall(cmd)

    # Create sequence dictionaries (for piccard)
    # cmd = "srun --mem 80000 -p develop java -Xmx8g -jar /cm/shared/apps/picard-tools/1.140/picard.jar"
    picard_jar="picard" # can be installed via bioconda
    output_dict=os.path.join(spiked_dir, "spiked.dict")
    cmd = picard_jar
    cmd += " CreateSequenceDictionary"
    cmd += " REFERENCE={}".format(final_fa)
    cmd += " OUTPUT={}".format(output_dict)
    #cmd += " GENOME_ASSEMBLY={}".format(genome)
    #cmd += " SPECIES=human"
    systemcall(cmd)

    # Create reflat files
    #cmd = "srun --mem 80000 -p develop java -Xmx80g -jar ~/Drop-seq_tools-1.12/jar/dropseq.jar ConvertToRefFlat"
    #dropseq_jar=config['program']['dropseqjar']
    #cmd = "java -Xmx80g -jar "+dropseq_jar+" ConvertToRefFlat"
    #cmd = "java -Xmx80g -jar ~/Drop-seq_tools-1.12/jar/dropseq.jar ConvertToRefFlat"
    #dict:os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.dict")
    #refflat:os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.refFlat")
    output_refflat=os.path.join(spiked_dir, "spiked.refFlat")
    cmd = "ConvertToRefFlat "
    cmd += " SEQUENCE_DICTIONARY={}".format(output_dict)
    cmd += " ANNOTATIONS_FILE= {}".format(final_gtf)
    cmd += " OUTPUT={}".format(output_refflat)
    systemcall(cmd)

    # Remove vanilla genome
    #os.remove(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.fa"))
    #os.remove(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.gtf"))

def main():
  config=read_config(sys.argv[1])
  sgrnafile=config["metadata"]["sgrna_table"]
  generate_annotations(config,sgrnafile)

if __name__ == '__main__':
  try:
    main()
  except  KeyboardInterrupt:
    sys.stderr.write("User interrupt me! ;-) Bye!\n")
    sys.exit(0)
