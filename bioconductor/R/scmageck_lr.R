scmageck_lr <-
function(BARCODE,RDS,NEGCTRL,LABEL=NULL,PERMUTATION=NULL,SAVEPATH='./'){
  if(!is.null(LABEL)){
    data_label=LABEL}
  else{data_label='sample1'}
  
  if(!is.null(PERMUTATION)){
    n_permutation=as.integer(PERMUTATION)}
  else{n_permutation=10000}
  
  # read cell assignment and libray file ####
  bc_dox=read.table(BARCODE,header=T,as.is=T)
  bc_dox[,1]=sub('-\\d','',bc_dox[,1])
  
  guide_count=table(bc_dox$cell)
  ncnt=table(table(bc_dox$cell))
  print(paste('Total barcode records:',nrow(bc_dox)))
  
  # load neg control guides ####
  ngctrlgenelist=strsplit(NEGCTRL,',')[[1]]
  print(paste('Neg Ctrl guide:',paste(ngctrlgenelist,collapse=';')))
  
  # read Seurat RDS file ####
  print(paste("Reading RDS file:",RDS))
  targetobj=readRDS(RDS)
  
  # convert to ind_matrix ####
  ind_matrix<-frame2indmatrix(bc_dox,targetobj)
  print(paste('Index matrix dimension:',nrow(ind_matrix),',',ncol(ind_matrix)))
  
  # try to perform matrix regresson on single genes ####
  mat_for_single_reg=single_gene_matrix_regression(targetobj,ngctrlgene=ngctrlgenelist,indmatrix=ind_matrix,high_gene_frac=-1.00)
  Xmat=mat_for_single_reg[[1]]

  # Xmat[,which(colnames(Xmat)%in%ngctrlgenelist)[1]]=1 # already integrated into function
  Ymat=mat_for_single_reg[[2]]
  Amat_pm_lst=getsolvedmatrix_with_permutation_cell_label(Xmat,Ymat,npermutation = n_permutation)
  Amat=Amat_pm_lst[[1]]
  Amat_pval=Amat_pm_lst[[2]]
  #save(Amat,Amat_pval,Xmat,Ymat,ind_matrix,ngctrlgenelist,bc_dox,file=paste(data_label,'_LR.RData',sep=''))
  if(!is.null(SAVEPATH)){
    write.table(data.frame(Perturbedgene=rownames(Amat),Amat),file=paste(SAVEPATH,data_label,'_score.txt',sep=''),sep='\t',quote=F,row.names=F)
    write.table(data.frame(Perturbedgene=rownames(Amat),Amat_pval),file=paste(SAVEPATH,data_label,'_score_pval.txt',sep=''),sep='\t',quote=F,row.names=F)
  }
  return(list(data.frame(Perturbedgene=rownames(Amat),Amat), data.frame(Perturbedgene=rownames(Amat),Amat_pval)))
}
