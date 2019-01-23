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

u6_before='GAAACACCG'
u6_after='GTTTTAGAG'
# get grna_list
grna_list_file='/groups/ligrp/weili/JinZhang/180912_RawDataForSingleCell/180912/metadata/grna_list.txt'

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
  Generate a dictionary of mismatches
  """

  grdict={}
  gmismatchdict={}
  g_orig_seq_cnt={}
  g_mis_seq_cnt={}
  for gr in grna_table:
    orig_seq=gr[1].upper()
    if orig_seq in g_orig_seq_cnt:
      print('Identical sequences '+orig_seq+ ' for '+gr[0]+' and '+g_orig_seq_cnt[orig_seq]+'. Skip '+gr[0]+'.')
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
  print(str(n_skipped)+ ' mismatched sequences are skipped due to overlap with other sequences.')
  
  
  return (grdict,gmismatchdict)
  
  


in_fqfilelist=sys.argv[1]
if os.path.isfile(in_fqfilelist):
  in_fqfilelist=[in_fqfilelist]
else:
  in_fqfilelist=glob.glob(in_fqfilelist)

if len(sys.argv)>2:
  out_file=sys.argv[2]
  if os.path.isfile(out_file) and out_file.endswith('txt')==False:
    print('Error: file '+out_file+' already exists.')
    sys.exit(-1)
else:
  out_file='output.txt'
#out_countfile=sys.argv[2]

print('Total number of files:'+str(len(in_fqfilelist)))

#sys.exit(0)


grna=[line.strip().split() for line in open(grna_list_file)]
grna=grna[1:]


(grdict,gmismatchdict)=gen_mismatch_dict(grna,2)

print('Total gRNA:'+str(len(grna)))


gr_count_dict={}

for in_fqfile in in_fqfilelist:
  print('Processing file '+in_fqfile)
  nl=0
  gr_count_dict[in_fqfile]={s:0 for s in grdict.keys()}
  gr_count=gr_count_dict[in_fqfile]
  for line in open(in_fqfile):
    nl+=1
    if nl%4!=2:
      continue
    # u6 before or u6 after must be present
    if u6_before not in line and u6_after not in line:
      continue
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
            break
  
  for (sg,ct) in gr_count.items():
    if ct>0:
      print(sg+'\t'+str(ct))
  
# output to file

outff=open(out_file,'w')
grdict_k=[x for x in grdict.keys()]
print('\t'.join(['sample']+grdict_k),file=outff)
for in_fqfile in gr_count_dict.keys():
  gr_count=gr_count_dict[in_fqfile]
  print('\t'.join([in_fqfile]+[str(gr_count[sg]) for sg in grdict_k]),file=outff)

outff.close()


  
  
        




