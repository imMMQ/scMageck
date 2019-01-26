#!/usr/bin/env python3
'''
count occurences in cropseq
Usage:
python custom_count.py filtered_fastq/bc0014.fastq
or
python custom_count.py "filtered_fastq/bc001*"
'''

import sys
import re
import glob
import os
import itertools
import argparse
import logging

def crop_parseargs():
  """
  Parsing arguments
  """
  parser = argparse.ArgumentParser(description='cropseq-count: counting sgRNAs from CROPseq experiments.')
  
  parser.add_argument('-v', '--version',action='version',version='0.1')
  
  parser.add_argument('-n','--output-prefix',default='sample1',help='The prefix of the output file(s). Default sample1.')

  parser.add_argument('-l','--lib-grna',required=True,help='A gRNA library file containing the list of sgRNA names, their sequences and associated genes, separated by tab.')
  parser.add_argument('--no-reverse-complement',action='store_true',help='Do not perform reverse complement search of sgRNAs.')
  parser.add_argument('-m','--max-mismatch',type=int,default=2,help='Maximum number of mismatches to be considered in sgRNA search. Default 2. Not recommended for values greater than 2. Decrease this value to speed up the search.')
  parser.add_argument('--anchor-before',default='GAAACACCG',help='Anchor sequence before the sgRNA. Default GAAACACCG (at the end of U6 promoter).')
  parser.add_argument('--anchor-after',default='GTTTTAGAG',help='Anchor sequence after the sgRNA. Default GTTTTAGAG.')
  parser.add_argument('--files',nargs='+',help='A list of fastq files, SAM/BAM files, or a wildcard of filenames.')
  


  # post-processing
  args=parser.parse_args()
    
  logging.basicConfig(level=10,
    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
    datefmt='%a, %d %b %Y %H:%M:%S',
    # stream=sys.stderr,
    filename=args.output_prefix+'.log',
    filemode='w'
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
  
  # add paramters
  logging.info('Parameters: '+' '.join(sys.argv))
  
  return args


# two versions of rev comp, depending on different versions of python
'''
Reverse complement
'''
if sys.version_info >= (3,1):
  trans_table=str.maketrans("ACGT","TGCA")
  def count_revcomp(x):
    return x.translate(trans_table)[::-1]
else:
  trans_table=string.maketrans("ACGT","TGCA")
  def count_revcomp(x):
    return x.translate(trans_table)[::-1]
#u6_before='GAAACACCG'
#u6_after='GTTTTAGAG'
# get grna_list
#grna_list_file='/groups/ligrp/weili/JinZhang/180912_RawDataForSingleCell/180912/metadata/grna_list.txt'

def gen_mismatches(sequence,num_mismatches):
  """
  generate mismatches of a certain sequence
  borrow from https://github.com/shendurelab/single-cell-ko-screens/blob/master/get_barcodes.py
  """
  letters = 'ACGT'
  mismatches=[]
  for locs in itertools.combinations(range(len(sequence)), num_mismatches):
      sequence_list = [[char] for char in sequence]
      for loc in locs:
          orig_char = sequence[loc]
          sequence_list[loc] = [l for l in letters if l != orig_char]

      for poss in itertools.product(*sequence_list):
          mismatches.append(''.join(poss))

  return mismatches

def gen_mismatch_dict(grna_table,num_mismatches):
  """
  Generate a dictionary of grna and its mismatches
  """

  grdict={}
  gmismatchdict={}
  g_orig_seq_cnt={}
  g_mis_seq_cnt={}
  for gr in grna_table:
    orig_seq=gr[1].upper()
    if orig_seq in g_orig_seq_cnt:
      logging.warning('Identical sequences '+orig_seq+ ' for '+gr[0]+' and '+g_orig_seq_cnt[orig_seq]+'. Skip '+gr[0]+'.')
      continue
    g_orig_seq_cnt[orig_seq]=gr[0]
    #grdict[gr[0]]=u6_before+gr[1]+u6_after
    grdict[gr[0]]=[orig_seq]
    nms=gen_mismatches(orig_seq,num_mismatches)
    gmismatchdict[gr[0]]=nms
    for sq in nms+[orig_seq]:
      g_mis_seq_cnt[sq]=(lambda s: g_mis_seq_cnt[s]+1 if s in g_mis_seq_cnt else 1)(sq)
  # remove those that overlap with each other
  n_skipped=0
  n_totalcnt=0
  for gr in grna_table:
    if gr[0] not in grdict:
      continue
    orig_seq=gr[1].upper()
    nms=gen_mismatches(orig_seq,num_mismatches)
    for sq in nms:
      if g_mis_seq_cnt[sq]>1:
        n_skipped=1
      else:
        gmismatchdict[gr[0]]+=[sq]
        n_totalcnt+=1
  logging.info('Total sgRNAs:'+str(len(grdict)))
  logging.info('Total sgRNAs with mismatches:'+str(n_totalcnt))
  logging.info(str(n_skipped)+ ' mismatched sequences are skipped due to overlap with other sequences.')
  
  
  return (grdict,gmismatchdict)

def process_input_files(args):
 
  in_fqfilelist=[]
  for sf in args.files:
    if os.path.isfile(sf):
      in_fqfilelist+=[sf]
    else:
      for sf2 in glob.glob(sf):
        if os.path.isfile(sf2):
          in_fqfilelist+=[sf2]
        else:
          logging.error(sf+' does not exist.')
          sys.exit(-1)
  return in_fqfilelist

def process_library_file(args): 
  grna_list_file=args.lib_grna
  grna=[line.strip().split() for line in open(grna_list_file)]
  grna=grna[1:]
  logging.info('Total gRNA:'+str(len(grna)))
  return grna

def search_sequence(line,gr_count,grdict,gmismatchdict,args):
  """
  Perform sequence search between anchor_before and anchor_after
  Parameters
    line
      a sequence to be search
    gr_count
      a dictionary to be updated if a hit was found
    grdict
    gmismatchdict
      a dictionary of {sgrna_name:[sequence]}
    args
      arguments
  Returns
    True/False whether a hit was found
	  
  """
  #line=count_revcomp(line)
  # u6 before or u6 after must be present
  u6_before=args.anchor_before
  u6_after=args.anchor_after
  if u6_before not in line and u6_after not in line:
    return 0
  # first, count perfect matches
  found_rec=False
  for (gn,gseqlist) in grdict.items():
    for gseq in gseqlist:
      if u6_before+gseq in line or gseq+u6_after in line:
        gr_count[gn]+=1
        found_rec=True
        break
  # if perfect matches are not found, go to the next time-consuming step for mismatches 
  if found_rec == False:
    for (gn,gseqlist) in gmismatchdict.items():
      for gseq in gseqlist:
        if u6_before+gseq in line or gseq+u6_after in line:
          gr_count[gn]+=1
          found_rec=True
          break
  return found_rec
 
def process_fastq_files(in_fqfile,grdict,gmismatchdict,args):
    nl=0
    gr_count={s:0 for s in grdict.keys()}
    for line in open(in_fqfile):
      nl+=1
      if nl%4!=2:
        continue
      hashit=search_sequence(line,gr_count,grdict,gmismatchdict,args)
      if hashit==False and args.no_reverse_complement == False:
        line=count_revcomp(line)
        search_sequence(line,gr_count,grdict,gmismatchdict,args)
    
    nct=0
    for (sg,ct) in gr_count.items():
      if ct>0:
        print(sg+'\t'+str(ct))
        nct+=ct
    logging.info(str(nct)+' hits found..')
    return gr_count
 
def output_to_file(args,gr_count_dict,grdict):
  """
  Write count output to file
  """  
  out_file=args.output_prefix+'.output.txt'
  outff=open(out_file,'w')
  grdict_k=[x for x in grdict.keys()]
  print('\t'.join(['sample']+grdict_k),file=outff)
  for in_fqfile in gr_count_dict.keys():
    gr_count=gr_count_dict[in_fqfile]
    print('\t'.join([in_fqfile]+[str(gr_count[sg]) for sg in grdict_k]),file=outff)
  
  outff.close()


  
if __name__ == '__main__':
  args=crop_parseargs()
  in_fqfilelist=process_input_files(args)

  logging.info('Total number of files:'+str(len(in_fqfilelist)))
  grna=process_library_file(args) 



  (grdict,gmismatchdict)=gen_mismatch_dict(grna,args.max_mismatch)



  gr_count_dict={}

  # process file
  for in_fqfile in in_fqfilelist:
    logging.info('Processing file '+in_fqfile)
    gr_count=process_fastq_files(in_fqfile,grdict,gmismatchdict,args)
    gr_count_dict[in_fqfile]=gr_count
   
  # output to file
  output_to_file(args,gr_count_dict,grdict)
  
    
    
          




