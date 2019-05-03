# command line: run_RRA.R RDS label start_ind end_ind
#
library(R.utils)

# parsing the command 
args <- commandArgs(trailingOnly=TRUE,asValues=TRUE)
for(ag in args){
  print(ag)
}

# determine the location of virtual_facs_function.R 
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

print(script.basename)




#source('gseafunc.R')
function_script=file.path(script.basename,'virtual_facs_functions.R')
print(paste('func script:',function_script))
source(function_script)

quit()

RRAPATH='/home/wei_li/.conda/envs/py3k/bin/RRA'




targetobj=readRDS(args[1])
data_label=args[2]

bc_dox=read.table(paste('../pairko/run_in_cluster/wl_get_barcodes_output_',data_label,'.txt',sep=''),header=T,as.is = T)


bc_gene=read.table('test_sgrna.txt',header=T,as.is = T)
rownames(bc_gene)=bc_gene[,2]
sum(bc_dox$barcode%in%bc_gene$sequence)
dim(bc_dox)

bc_dox[,c('oligo','gene')]=bc_gene[bc_dox$barcode,c('oligo','gene')]
bc_dox[,1]=sub('-\\d','',bc_dox[,1])

guide_count=table(bc_dox$cell)
ncnt=table(table(bc_dox$cell))
ncnt
barplot(ncnt)

# only leave cells with unique guides

dupsq=bc_dox[duplicated(bc_dox$cell),1]
bc_dox_uq=bc_dox[!bc_dox[,1]%in%dupsq,]
rownames(bc_dox_uq)=bc_dox_uq[,1]

dim(bc_dox)
dim(bc_dox_uq)



target_gene_list=c('BAK1','BAD','BAX', 'BCL2','BCL2L1', 'CASP9', 'CASP2', 'CASP8','CASP7','CASP3')
target_gene_list=unique(bc_dox_uq[bc_dox_uq$gene!='NonTargetingControlGuideForHuman','gene'])
target_gene_list=target_gene_list[!is.na(target_gene_list)]

#target_gene='CASP8'
for(target_gene in target_gene_list){
  if(!target_gene%in%rownames(targetobj@scale.data)){
    next
  }
  texp=targetobj@scale.data[target_gene,]
  texp=sort(texp)
  texp_withg=texp[names(texp)%in%rownames(bc_dox_uq) & !is.na(bc_dox_uq[names(texp),'oligo'])]

  other_table=get_rank_tables_from_rra(texp_withg,bc_dox_uq,tmpprefix=paste('sample_',runif(1,1,10000),sep=''),rrapath = RRAPATH)

  write.table(other_table,file=paste('individual_table/',data_label,'_',target_gene,'_RRA.txt',sep=''),sep='\t',quote=F,row.names=F)

}



