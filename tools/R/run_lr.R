# command line: Rscript run_individual_genes.R OPTIONS
# Required options
# --BARCODE
# --RDS
# --NEGCTRL
#
# Optional options
# --LABEL
# --PERMUTATION
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
source(function_script)
function_script=file.path(script.basename,'multiple_guides_function.R')
source(function_script)

#quit()

# optional parameters ####

data_label=ifelse(is.null(args['LABEL']),'sample1',args['LABEL'])
n_permutation=ifelse(is.null(args['PERMUTATION']),10000,as.integer(args['PERMUTATION']))

# read cell assignment and libray file ####

#bc_dox=read.table(paste('../pairko/run_in_cluster/wl_get_barcodes_output_',data_label,'.txt',sep=''),header=T,as.is = T)
#bc_gene=read.table('test_sgrna.txt',header=T,as.is = T)
bc_dox=read.table(args[['BARCODE']],header=T,as.is=T)
#bc_gene=read.table(args[['LIBRARY']],header=T,as.is=T)

#rownames(bc_gene)=bc_gene[,2]
#n_doxinlib=sum(bc_dox$barcode%in%bc_gene$sequence)
#print(paste('Number of barcodes in library:',n_doxinlib))
#dim(bc_dox)

#bc_dox[,c('oligo','gene')]=bc_gene[bc_dox$barcode,c('oligo','gene')]
bc_dox[,1]=sub('-\\d','',bc_dox[,1])

guide_count=table(bc_dox$cell)
ncnt=table(table(bc_dox$cell))
#print('Cell count:')
#print(ncnt)
#barplot(ncnt)


print(paste('Total barcode records:',nrow(bc_dox)))

# load neg control guides ####

#ngctrlgenelist=read.table(args[['NEGCTRL']],header=F,as.is=T)[,1]
ngctrlgenelist=strsplit(args[['NEGCTRL']],',')[[1]]
print(paste('Neg Ctrl guide:',paste(ngctrlgenelist,collapse=';')))

# read Seurat RDS file ####

print(paste("Reading RDS file:",args[['RDS']]))
targetobj=readRDS(args[['RDS']])

# convert to ind_matrix ####

ind_matrix<-frame2indmatrix(bc_dox,targetobj)

print(paste('Index matrix dimension:',nrow(ind_matrix),',',ncol(ind_matrix)))


#save(list=ls(),file='tmp.RData')
# try to perform matrix regresson on single genes ####

mat_for_single_reg=single_gene_matrix_regression(targetobj,ngctrlgene=ngctrlgenelist,indmatrix=ind_matrix,high_gene_frac=-1.00)
Xmat=mat_for_single_reg[[1]]
#
# Xmat[,which(colnames(Xmat)%in%ngctrlgenelist)[1]]=1 # already integrated into function

Ymat=mat_for_single_reg[[2]]
Amat_pm_lst=getsolvedmatrix_with_permutation_cell_label(Xmat,Ymat,npermutation = n_permutation)
Amat=Amat_pm_lst[[1]]
Amat_pval=Amat_pm_lst[[2]]

save(Amat,Amat_pval,Xmat,Ymat,ind_matrix,ngctrlgenelist,bc_dox,file=paste(data_label,'_LR.RData',sep=''))
write.table(Amat,file=paste(data_label,'_score.txt',sep=''),sep='\t',quote=F,row.names=T)
write.table(Amat_pval,file=paste(data_label,'_score_pval.txt',sep=''),sep='\t',quote=F,row.names=T)



