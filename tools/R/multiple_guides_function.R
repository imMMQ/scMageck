# function definitions #####
# version: 02-22-2019
# should sync from the version in macbook:
# /Users/weili/Dropbox/work/cropseq/Shendure/nmeth18/multiple_guides_function.R


getscaledata<-function(targetobj,scaled=TRUE){
# if scaled=FALSE, return raw.data
  if('scale.data'%in%names(attributes(targetobj))){
    if(scaled){
      scalef=targetobj@scale.data # for version 2
    }else{
      scalef=targetobj@raw.data # for version 2
    }
  }else{
    if(scaled){
      scalef=GetAssayData(object = targetobj, slot = "scale.data")
    }else{
      scalef=GetAssayData(object = targetobj, slot = "counts")
    }
  }
  return (scalef)
}


# perform matrix decomposition
# get A for Y=XA
# A= (X^TX)^-1X^TY
# Y:  (cells * expressed genes)
# X: design matrix, (cells * KO genes), can also be (cells * double KO genes)
# A: (KO genes * expressed genes)
getsolvedmatrix<-function(Xm,Ym,lambda=0.01){
  #Amat=solve(Xmat,Ymat) # solve AX=B, or Xmat * A =Ymat
  TMmat_g= (t(Xm) %*% Xm) + lambda * diag(ncol(Xm))
  
  Amat_g= solve(TMmat_g) %*% t(Xm) %*% Ym
  return (Amat_g)
}

getsolvedmatrix_with_permutation_cell_label<-function(Xm,Ym,lambda=0.01,npermutation=1000){
  Amat_ret=getsolvedmatrix(Xm,Ym,lambda=lambda)
  Amat_ret_higher=matrix(rep(0,ncol(Amat_ret)*nrow(Amat_ret)),nrow = nrow(Amat_ret))
  rownames(Amat_ret_higher)=rownames(Amat_ret)
  colnames(Amat_ret_higher)=colnames(Amat_ret)
  # permute N times
  # randomly shuffle cell labels
  for(npm in 1:npermutation){
    if(npm%%100==0){
      print(paste('Permutation:',npm,'/',npermutation,'...'))
    }
    cells_shu=sample(rownames(Ym),nrow(Ym))
    Xm_s=Xm[cells_shu,]
    Ym_s=Ym # [cells_shu,]
    rownames(Ym_s)=cells_shu
    Amat_random=getsolvedmatrix(Xm_s,Ym_s,lambda=lambda)
    
    Amat_ret_higher=Amat_ret_higher+ (abs(Amat_random)>abs(Amat_ret))*1.0
    #browser()
  }
  Amat_ret_higher=Amat_ret_higher/npermutation
  return (list(Amat_ret,Amat_ret_higher))
}


# construct a matrix of Y=XA, Y= (cells*expressed genes), X=(cells* KO genes), A=(KO genes * expressed genes)

single_gene_matrix_regression<-function(targetobj,ngctrlgene=c('NonTargetingControlGuideForHuman'),indmatrix=NULL,high_gene_frac=0.01){
  # return X matrix and Y matrix for regression
  # note that all the ngctrlgene are merged into one column, "NegCtrl"
  # if indmatrix is provided, the Xmat will be constructed from indmatrix
  rawf=getscaledata(targetobj,scaled=F)
  select_genes=rownames(rawf)[ which(rowSums(as.matrix(rawf)!=0)>=ncol(rawf)*high_gene_frac)]
  print(paste('Selected genes:',length(select_genes)))
  # browser()

  scalef=getscaledata(targetobj)

  if(is.null(indmatrix)){
    select_cells=rownames(targetobj@meta.data)[which(!is.na(targetobj@meta.data$geneID))]
  }else{
    select_cells=rownames(indmatrix)
    select_cells=select_cells[select_cells %in% colnames(scalef)]
  }
  YmatT=scalef[select_genes,select_cells]
  
  Ymat=as.matrix(t(YmatT)) # (cells * expressed genes)
  if(is.null(indmatrix)){
    tgf=targetobj@meta.data[select_cells,'geneID']
    tgf[tgf%in%ngctrlgene]='NegCtrl'
    tgphenotype=as.factor(tgf)
    Xmat=matrix(rep(0,length(select_cells)*length(unique(tgphenotype))),nrow=length(select_cells))
    rownames(Xmat)=select_cells
    colnames(Xmat)=levels(tgphenotype)
    Xmat[as.matrix(cbind(1:nrow(Xmat),as.numeric(tgphenotype)))]=1
    Xmat[,'NegCtrl']=1 # set up base line
  }else{
    tgf=colnames(indmatrix)
    tgf[tgf%in%ngctrlgene]='NegCtrl'
    tgphenotype=as.factor(tgf)
    
    Xmat=matrix(rep(0,length(select_cells)*length(unique(tgphenotype))),nrow=length(select_cells))
    rownames(Xmat)=select_cells
    colnames(Xmat)=levels(tgphenotype)
    for(cnl in colnames(indmatrix)){
      cellns=which(indmatrix[,cnl]==TRUE) #make sure indmatrix 
      if(cnl%in%ngctrlgene){
        Xmat[cellns,'NegCtrl']=1
      }else{
        Xmat[cellns,cnl]=1
      }
    }
    Xmat[,'NegCtrl']=1
    
  }# end if
  
  return (list(Xmat,Ymat))
}


plot_single_genes<-function(targetob,gene1,targetgene,haslog=TRUE,plotfigure=TRUE,ngctrlgene=c('NonTargetingControlGuideForHuman'),indmatrix=NULL){
  # plot single gene ko effect
  # if indmatrix =NULL, use targetobj@metadata to identify cell groups
  # else, use indmatrix where cells that do not bear gene1 target as controls

  if(is.null(indmatrix)){
    cell_ctrl=rownames(targetob@meta.data)[which(targetob@meta.data$geneID%in% ngctrlgene)]
    cell_gene1=rownames(targetob@meta.data)[which(targetob@meta.data$geneID==gene1)]
    if(length(cell_gene1)==0){
      gene1=sub('-','_',gene1)
      cell_gene1=rownames(targetob@meta.data)[which(targetob@meta.data$geneID==gene1)]
    }
  }else{
     cell_ctrl=rownames(indmatrix)[indmatrix[,gene1]==FALSE]
     cell_gene1=rownames(indmatrix)[indmatrix[,gene1]==TRUE]
  }
  #cell_gene2=rownames(targetob@meta.data)[which(targetob@meta.data$geneID==gene2)]
  #mg_geneid=paste(gene1,gene2,sep='_')
  #cell_merged=select_pair_genes[[mg_geneid]]
  scalef=getscaledata(targetobj) 
  t_exp=scalef[targetgene,c(cell_ctrl,cell_gene1)]
  if(min(t_exp)<0){
    t_exp=t_exp-(min(t_exp))+0.1
  }
  t_type=c(rep('NegCtrl',length(cell_ctrl)),rep(gene1,length(cell_gene1)))
  t_type=factor(t_type,levels = c('NegCtrl',gene1))
  ds=data.frame(Expression=t_exp,
                Type=t_type, Cells=c(cell_ctrl,cell_gene1))
  
  if(plotfigure){
    p<-ggplot(ds, aes(x=Type, y=Expression)) + 
      geom_violin()+
      geom_point(position=position_jitter(w=0.1,h=0)) +
      ggtitle(paste(targetgene,'expression'))
    if(haslog){
      p=p+  scale_y_log10()
    }
    #ggtitle(paste(ensemblID,geneID))
    print(p)
  }
  
  return (ds)
}

plot_gi_genes<-function(targetobj,gene1,gene2,targetgene,select_pair_genes,haslog=TRUE,plotfigure=TRUE,ngctrlgene=c('NonTargetingControlGuideForHuman')){
  if(gene1>gene2){
    tx=gene2
    gene2=gene1
    gene1=tx
  }
  cell_ctrl=rownames(targetobj@meta.data)[which(targetobj@meta.data$geneID%in% ngctrlgene)]
  cell_gene1=rownames(targetobj@meta.data)[which(targetobj@meta.data$geneID==gene1)]
  cell_gene2=rownames(targetobj@meta.data)[which(targetobj@meta.data$geneID==gene2)]
  mg_geneid=paste(gene1,gene2,sep='_')
  cell_merged=select_pair_genes[[mg_geneid]]
  cell_merged=cell_merged[cell_merged%in%rownames(targetobj@meta.data)]
  
  #browser()

  scalef=getscaledata(targetobj)
  t_exp=scalef[targetgene,c(cell_ctrl,cell_gene1,cell_gene2,cell_merged)]
  if(min(t_exp)<0){
    t_exp=t_exp-(min(t_exp))+0.1
  }
  t_type=c(rep('NegCtrl',length(cell_ctrl)),rep(gene1,length(cell_gene1)),rep(gene2,length(cell_gene2)),rep(mg_geneid,length(cell_merged)))
  t_type=factor(t_type,levels = c('NegCtrl',gene1,gene2,mg_geneid))
  ds=data.frame(Expression=t_exp,
                Type=t_type, Cells=c(cell_ctrl,cell_gene1,cell_gene2,cell_merged))
  
  if(plotfigure){
    p<-ggplot(ds, aes(x=Type, y=Expression)) + 
      geom_violin()+
      geom_point(position=position_jitter(w=0.1,h=0)) +
      ggtitle(paste(targetgene,'expression'))
    if(haslog){
      p=p+  scale_y_log10()
    }
    #ggtitle(paste(ensemblID,geneID))
    print(p)
  }
  
  return (ds)
}


prepare_matrix_for_pair_reg<-function(targetobj,bc_dox,Xmat,Ymat,Amat,cell_cutoff=4,ngctrlgene=c('NonTargetingControlGuideForHuman')){
  # this function returls (X, Y) for paired KO regression
  # Xmat, Ymat, Amat: these are regression for single genes
  # the gene and cell column in bc_dox is used to determine the target of each cell
  bc_dox_nonuq=bc_dox[!is.na(bc_dox$gene),]
  dupsq=bc_dox_nonuq[duplicated(bc_dox_nonuq$cell),1]
  bc_dox_nonuq=bc_dox_nonuq[bc_dox_nonuq[,1]%in%dupsq,]
  print(paste('non_uq genes:',nrow(bc_dox_nonuq)))
  
  # get the targeting genes per cell
  cell_list=list()
  for(i in 1:nrow(bc_dox_nonuq)){
    cellname=bc_dox_nonuq[i,1]
    targetgene=bc_dox_nonuq[i,'gene']
    if(cellname %in% names(cell_list)){
      cell_list[[cellname]]=append(cell_list[[cellname]],targetgene)
    }else{
      cell_list[[cellname]]=c(targetgene)
    }
  }
  print(paste('finish compiling',length(cell_list),'cells'))
  #browser()
  # filter pairs
  library(utils)
  pair_genes=list()
  
  for(si in 1:length(cell_list)){
    gls=unique(cell_list[[si]])
    gls_cell_name=names(cell_list)[si]
    if(length(gls)<2){
      next
    }
    if(length(gls)>2){
      #print(paste(length(gls),'combinations...'))
    }
    zi_cb=combn(sort(gls),2)
    for(zi in ncol(zi_cb)){
      zi_c=zi_cb[,zi]
      zi_c_id=paste(zi_c,collapse = '_')
      if(zi_c_id %in% names(pair_genes)){
        pair_genes[[zi_c_id]]=append(pair_genes[[zi_c_id]],gls_cell_name)
      }else{
        pair_genes[[zi_c_id]]=c(gls_cell_name)
      }
    }
  }
  print(paste('finished calculating ',length(pair_genes),'pairs'))

 
  scalef=getscaledata(targetobj)
 
  pair_genes_length=unlist(lapply(pair_genes,length))
  
  select_pair_genes=names(pair_genes_length)[pair_genes_length>=cell_cutoff]
  
  # construct Y and X matrix with only double KO
  select_pair_genes=pair_genes[(select_pair_genes)]
  
  cells_combine=unique(unlist(select_pair_genes))
  cells_combine=cells_combine[cells_combine%in% colnames(scalef)]
  
  print(paste('gene pairs left:',length(cells_combine)))
  
  # construct a matrix of Y=XA, Y= (cells*expressed genes), X=(cells* KO genes), A=(KO genes * expressed genes)
  #select_genes=rownames(targetobj@raw.data)[ which(rowSums(targetobj@raw.data!=0)>ncol(targetobj@raw.data)/100)]
  select_genes=colnames(Ymat)
  select_cells=rownames(Ymat)
  YmatT_db=scalef[select_genes,cells_combine]
  
  Ymat_db=as.matrix(t(YmatT_db)) # (cells * expressed genes)
  #tgphenotype=targetobj@meta.data[select_cells,'geneID']
  #tgf=targetobj@meta.data[select_cells,'geneID']
  #tgf[tgf%in%ngctrlgene]='NegCtrl'
  #tgphenotype=as.factor(tgf)
  tgphenotype=as.factor(colnames(Xmat))
  # need to calculate residue
  Xmat_db_res=matrix(rep(0,length(cells_combine)*length(unique(tgphenotype))),nrow=length(cells_combine))
  
  rownames(Xmat_db_res)=cells_combine
  colnames(Xmat_db_res)=levels(tgphenotype)
  #cells_combine[as.matrix(cbind(1:nrow(cells_combine),as.numeric(tgphenotype)))]=1
  #browser()
  for(i in 1:nrow(bc_dox_nonuq)){
    t_cellname=bc_dox_nonuq[i,1]
    t_gene_target=bc_dox_nonuq[i,'gene']
    if(t_cellname%in%cells_combine & t_gene_target %in% colnames(Xmat_db_res)){
      Xmat_db_res[t_cellname,t_gene_target]=1
    }
  }
  # residule of Xmat * A =Ymat
  Ymat_db_residule=Ymat_db- Xmat_db_res %*% Amat
  
  print('creating matrix for paired gene...')
  # now, create X matrix
  
  Xmat_db=matrix(rep(0,length(cells_combine)*length(names(select_pair_genes))),nrow=length(cells_combine))
  rownames(Xmat_db)=cells_combine
  colnames(Xmat_db)=names(select_pair_genes)
  
  for ( i in 1:length(select_pair_genes)){
    t_pair_name=names(select_pair_genes)[i]
    for(cn in select_pair_genes[[i]]){
      if(cn %in% cells_combine){
        Xmat_db[cn,t_pair_name]=1
      }
    }
  }
  
  return (list(Xmat_db,Ymat_db_residule,select_pair_genes))
}


prepare_matrix_for_pair_reg_indmat<-function(targetobj,ind_matrix,Xmat,Ymat,Amat,cell_cutoff=4,ngctrlgene=c('NonTargetingControlGuideForHuman')){
  # this function returns (X, Y) for paired KO regression
  # Xmat, Ymat, Amat: these are regression for single genes
  # the gene and cell column in ind_matrix is used to determine the target of each cell

  scalef=getscaledata(targetobj)
  # new code
  s_pair_x=c()
  s_pair_y=c()
  cells_s=rep(0,nrow(ind_matrix));names(cells_s)=rownames(ind_matrix)
  for(i in 1:(ncol(ind_matrix)-1)){
    if(i%%10==1){
      print(paste(i,'...'))
    }
    for(j in (i+1):ncol(ind_matrix)){
      nct=sum(ind_matrix[,i] & ind_matrix[,j])
      if(nct>=cell_cutoff){
        s_pair_x=append(s_pair_x,i)
        s_pair_y=append(s_pair_y,j)
        cells_s=cells_s+ (ind_matrix[,i] & ind_matrix[,j])
      }
    }
  }
  cells_combine=names(cells_s)[which(cells_s>0)]
  cells_combine=cells_combine[cells_combine%in% colnames(scalef)]
  
  
  # construct a matrix of Y=XA, Y= (cells*expressed genes), X=(cells* KO genes), A=(KO genes * expressed genes)
  #select_genes=rownames(targetobj@raw.data)[ which(rowSums(targetobj@raw.data!=0)>ncol(targetobj@raw.data)/100)]
  select_genes=colnames(Ymat)
  select_cells=rownames(Ymat)
  YmatT_db=scalef[select_genes,cells_combine]
  
  Ymat_db=as.matrix(t(YmatT_db)) # (cells * expressed genes)
  #tgphenotype=targetobj@meta.data[select_cells,'geneID']
  #tgf=targetobj@meta.data[select_cells,'geneID']
  #tgf[tgf%in%ngctrlgene]='NegCtrl'
  #tgphenotype=as.factor(tgf)
  tgphenotype=as.factor(colnames(Xmat))
  # need to calculate residue
  #Xmat_db_res=matrix(rep(0,length(cells_combine)*length(unique(tgphenotype))),nrow=length(cells_combine))
  #rownames(Xmat_db_res)=cells_combine
  #colnames(Xmat_db_res)=levels(tgphenotype)
  Xmat_db_res=Xmat[cells_combine,]
  
  #cells_combine[as.matrix(cbind(1:nrow(cells_combine),as.numeric(tgphenotype)))]=1
  #browser()

  # residule of Xmat * A =Ymat
  Ymat_db_residule=Ymat_db- Xmat_db_res %*% Amat
  
  print('creating matrix for paired gene...')
  # now, create X matrix
  
  Xmat_db=matrix(rep(0,length(cells_combine)*length(s_pair_x)),nrow=length(cells_combine))
  rownames(Xmat_db)=cells_combine
  colnames(Xmat_db)=paste(colnames(ind_matrix)[s_pair_x],'_',colnames(ind_matrix)[s_pair_y],sep='')
  
  print(paste('selected gene pairs:',length(s_pair_x)))
  select_pair_genes=list()
  for ( i in 1:length(s_pair_x)){
      nct=rownames(ind_matrix)[which(ind_matrix[,s_pair_x[i]] & ind_matrix[,s_pair_y[i]])]
      cn=nct[nct%in%cells_combine]
      Xmat_db[cn,i]=1
      namex=colnames(ind_matrix)[s_pair_x[i]]
      namey=colnames(ind_matrix)[s_pair_y[i]]
      genepairname=paste(namex,'_',namey,sep='')
      select_pair_genes[[genepairname]]=cn
  }
  
  return (list(Xmat_db,Ymat_db_residule,select_pair_genes))
}


frame2indmatrix<-function(bc_d,targetobj){

  scalef=getscaledata(targetobj)
  colnames(scalef) = sub('-\\d$','',colnames(scalef))
  rnm=unique(bc_d$cell)
  cnm=unique(bc_d$gene)
  rnm=rnm[!is.na(rnm)]
  rnm=rnm[rnm%in%colnames(scalef)]
  cnm=cnm[!is.na(cnm)]
  ind_matrix=matrix(rep(FALSE,length(rnm)*length(cnm)),nrow=length(rnm))
  rownames(ind_matrix)=rnm
  colnames(ind_matrix)=cnm
  for(si in 1:nrow(bc_d)){
    t_r=bc_d[si,'cell']
    t_c=bc_d[si,'gene']
    if((t_r%in%rnm) & (t_c%in% cnm)){
      ind_matrix[t_r,t_c]=TRUE
    }
  }
  return (ind_matrix)
}
