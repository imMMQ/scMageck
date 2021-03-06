# command line: Rscript run_individual_genes.R OPTIONS
# Required arguments 
# --BARCODE
#	A txt file to include cell identity information, generated from the cell identity collection step.
# --RDS
#	An R RDS file of [Seurat](https://satijalab.org/seurat/) object containing the scRNA-seq dataset. Note that the dataset has to be normalized and scaled.
# --GENE
#	Genes whose expressions are to be tested. Multiple genes can be provided, separated by ",". For example, "MKI67,TP53"
#
# Optional arguments
# --RRAPATH
#	The path to the RRA program, if RRA cannot be found in the PATH environment variable.
# --LABEL
#	The label of the output files.
# --NEGCTRL 
#	The name of the negative control gene. For example, "NonTargetingControlGuideForHuman". Default is NULL (do not use any negative controls).
# --KEEPTMP
#	Keep temporary files.
# --PATHWAY
#	Treat genes in --GENE option as a pathway. In other words, the averaged expression of these genes will be used for testing.
#
library(R.utils)
library(Seurat)

# parsing the command  ####
args <- commandArgs(trailingOnly=TRUE,asValues=TRUE)
print('Arguments:')
for(i in 1:length(args)){
  print(paste(names(args)[i],':',args[i]))
}

# determine the location of virtual_facs_function.R  ####
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

print(script.basename)

#source('gseafunc.R')
function_script=file.path(script.basename,'virtual_facs_functions.R')
#print(paste('func script:',function_script))
source(function_script)

# determine the path of RRA

#RRAPATH='/home/wei_li/.conda/envs/py3k/bin/RRA'
RRAPATH=ifelse(is.null(args[['RRAPATH']]),'RRA',args[['RRAPATH']])
print('Checking RRA...')
#print(RRAPATH)
if(system(RRAPATH,ignore.stdout = TRUE, ignore.stderr = TRUE)!=0){
  print(paste('Error: cannot find RRA in ',RRAPATH,'. Please specify the path of RRA (in MAGeCK)'))
  quit()
}


#quit()


#data_label=ifelse(is.null(args['LABEL']),'sample1',args['LABEL'])
if('LABEL'%in%names(args)){data_label=args[['LABEL']]}else{data_label='sample1'}

# read cell assignment and libray file ####

#bc_dox=read.table(paste('../pairko/run_in_cluster/wl_get_barcodes_output_',data_label,'.txt',sep=''),header=T,as.is = T)
#bc_gene=read.table('test_sgrna.txt',header=T,as.is = T)
bc_dox=read.table(args[['BARCODE']],header=T,as.is=T)

keep_tmp=ifelse('KEEPTMP'%in%names(args),T,F)
print(paste('keep_tmp:',keep_tmp))
ispathway=ifelse('PATHWAY'%in%names(args),T,F)
print(paste('ispathway:',ispathway))


if('NEGCTRL'%in%names(args)){negctrl_gene=args[['NEGCTRL']]}else{negctrl_gene=NULL}

#bc_gene=read.table(args[['LIBRARY']],header=T,as.is=T)
#rownames(bc_gene)=bc_gene[,2]
#n_doxinlib=sum(bc_dox$barcode%in%bc_gene$sequence)
#print(paste('Number of barcodes in library:',n_doxinlib))
#dim(bc_dox)

#bc_dox[,c('oligo','gene')]=bc_gene[bc_dox$barcode,c('oligo','gene')]
bc_dox[,1]=sub('-\\d$','',bc_dox[,1])

guide_count=table(bc_dox$cell)
ncnt=table(table(bc_dox$cell))
#print('Cell count:')
#print(ncnt)
#barplot(ncnt)

# only leave cells with unique guides ####

dupsq=bc_dox[duplicated(bc_dox$cell),1]
bc_dox_uq=bc_dox[!bc_dox[,1]%in%dupsq,]
rownames(bc_dox_uq)=bc_dox_uq[,1]

print(paste('Total barcode records:',nrow(bc_dox)))
print(paste('Unique barcode records:',nrow(bc_dox_uq)))

# test target genes ####

target_gene_list=strsplit(args[['GENE']],',')[[1]]
print(paste('Target gene:',paste(target_gene_list,collapse=';')))


# read Seurat RDS file ####

print(paste("Reading RDS file:",args[['RDS']]))
targetobj=readRDS(args[['RDS']])

# run RRA ####

if('scale.data'%in%names(attributes(targetobj))){
  scalef=targetobj@scale.data # for version 2
}else{
  scalef=GetAssayData(object = targetobj, slot = "scale.data")
}

if(ispathway==TRUE){
  for(target_gene in target_gene_list){
    if(!target_gene%in%rownames(scalef)){
      print(paste('Error: gene ',target_gene,' not in expression list.'))
      quit()
    }
  }
  texp=colMeans(scalef[target_gene_list,])
  texp=sort(texp)
  texp_withg=texp[names(texp)%in%rownames(bc_dox_uq) & !is.na(bc_dox_uq[names(texp),'barcode'])]
  
  other_table=get_rank_tables_from_rra(texp_withg,bc_dox_uq,tmpprefix=paste('sample_',runif(1,1,10000),sep=''),rrapath = RRAPATH,keeptmp=keep_tmp,negctrlgenelist=negctrl_gene)
  
  write.table(other_table,file=paste(data_label,'_',paste(target_gene_list,collapse='_'),'_RRA.txt',sep=''),sep='\t',quote=F,row.names=F)
  
}else{
  # treat genes separately

  for(target_gene in target_gene_list){
    if(!target_gene%in%rownames(scalef)){
      print(paste('Warning: gene ',target_gene,' not in expression list.'))
      next
    }else{
      print(paste('Testing gene ',target_gene,'...'))
    }
    texp=scalef[target_gene,]
    texp=sort(texp)
    texp_withg=texp[names(texp)%in%rownames(bc_dox_uq) & !is.na(bc_dox_uq[names(texp),'barcode'])]
  
    other_table=get_rank_tables_from_rra(texp_withg,bc_dox_uq,tmpprefix=paste('sample_',runif(1,1,10000),sep=''),rrapath = RRAPATH,keeptmp=keep_tmp,negctrlgenelist=negctrl_gene)
  
    write.table(other_table,file=paste(data_label,'_',target_gene,'_RRA.txt',sep=''),sep='\t',quote=F,row.names=F)
  
  }

}



