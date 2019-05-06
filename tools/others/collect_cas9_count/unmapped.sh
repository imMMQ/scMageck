
# samtools view -b -f 4 ../possorted_genome_bam.bam > unmapped/unmapped.bam

# bedtools bamtofastq -i unmapped/unmapped.bam -fq unmapped/unmapped.fq

# bowtie2 -p 4 --no-unal -x /lustre/groups/ligrp/weili/Shendure/nmeth18/pipeline/out/cas9_only/bowtie2/cas9_seq -U unmapped/unmapped.fq -S unmapped/unmapped_to_cas9.sam   

cd unmapped
samtools view unmapped.bam | ./collect_cas9_count.py > cas9count.txt 
