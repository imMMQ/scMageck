# command line: Rscript run_lr.R OPTIONS
# Required arguments
# --BARCODE
# --RDS
# --NEGCTRL
#
# Optional arguments
# --LABEL
# --PERMUTATION
# --SIGNATURE
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

#data_label=ifelse(is.null(args['LABEL']),'sample1',args['LABEL'])
if('LABEL'%in%names(args)){data_label=args[['LABEL']]}else{data_label='sample1'}
#n_permutation=ifelse(is.null(args['PERMUTATION']),10000,as.integer(args['PERMUTATION']))
if('PERMUTATION'%in%names(args)){n_permutation=as.integer(args[['PERMUTATION']])}else{n_permutation=10000}
# To identify whether use the gmt file
if('SIGNATURE'%in%names(args)){data_signature=args[['SIGNATURE']]}
run_signature=ifelse('SIGNATURE'%in%names(args),T,F)
print(paste('run_signature:',run_signature))

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
bc_dox[,1]=sub('-\\d$','',bc_dox[,1])

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

# Optional function
# Get the results based on gmt file
if(run_signature==TRUE){
  gmt <- read.delim(args[['SIGNATURE']], header = FALSE)
  gmt <- t(as.matrix(gmt))
  colnames(gmt) <- gmt[1,]
  gmt <- gmt[-1:-2,]
  print(paste('Total signature records:',ncol(gmt)))
  sig_mat <- getsigmat(Ymat, gmt_file = gmt)
  if(ncol(sig_mat) > 0) {
    Amat_sig_lst=getsolvedmatrix_with_permutation_cell_label(Xmat,sig_mat,npermutation = 1000)
    sig_score=Amat_sig_lst[[1]]
    sig_pval=Amat_sig_lst[[2]]
    sig_re <- getsigresult(signature_score = sig_score, signature_pval = sig_pval)
    write.table(data.frame(sig_re),file=paste(data_label,'_signature.txt',sep=''),sep='\t',quote=F,row.names=F)
  }
}

save(Amat,Amat_pval,Xmat,Ymat,ind_matrix,ngctrlgenelist,bc_dox,file=paste(data_label,'_LR.RData',sep=''))
write.table(data.frame(Perturbedgene=rownames(Amat),Amat),file=paste(data_label,'_score.txt',sep=''),sep='\t',quote=F,row.names=F)
write.table(data.frame(Perturbedgene=rownames(Amat),Amat_pval),file=paste(data_label,'_score_pval.txt',sep=''),sep='\t',quote=F,row.names=F)



