# virtura_facs_functions

getscaledata<-function(targetobj){
  if('scale.data'%in%names(attributes(targetobj))){
    scalef=targetobj@scale.data # for version 2
  }else{
    scalef=GetAssayData(object = targetobj, slot = "scale.data")
  }
  return (scalef)
}

get_rank_tables<-function(genes_to_rank,negctrlgenelist='NonTargetingControlGuideForHuman'){
  score_test=seq(-1.0,1.0,length.out = length(genes_to_rank))
  score_tggene=paste(genes_to_rank,1:length(genes_to_rank),sep='_')
  names(score_test)=score_tggene
  
  candidate_gene=unique(genes_to_rank)
  candidate_gene=candidate_gene[!is.na(candidate_gene)]
  gsea_score=rep(1,length(candidate_gene))
  p_vsall=gsea_score
  p_vsnegctrl=gsea_score
  for(tg_i in 1:length(candidate_gene) ){
    
    tggene=candidate_gene[tg_i]
    # print(paste(tggene,'...'))
    ghit=score_tggene[genes_to_rank==tggene]
    
    gs_score=getoriginalgseascore(score_test,ghit)
    z=max(gs_score)
    gsea_score[tg_i]=z
    
    z=ks.test(score_test[!score_tggene%in%ghit],score_test[ghit])
    p_vsall[tg_i]=z$p.value
    
    
    
    z=ks.test(score_test[genes_to_rank%in%negctrlgenelist],score_test[ghit])
    p_vsnegctrl[tg_i]=z$p.value
    
  }
  
  report_f=data.frame(gene=candidate_gene,gsea=gsea_score,ks_p=p_vsall,ks_negctrl=p_vsnegctrl)
  
  report_f[,'ks_p_adj']=p.adjust(report_f$ks_p,method = 'fdr')
  
  report_f[,'ks_negctrl_adj']=p.adjust(report_f$ks_negctrl,method = 'fdr')
  return (report_f)
}


get_rank_tables_from_rra<-function(rankexp,bc_dox_u,rrapath=NULL,pcutoff=0.3,tmpprefix=paste('sample_',runif(1,1,10000),sep=''),negctrlgenelist='NonTargetingControlGuideForHuman',more_rra='',negsel=T,possel=T){
  rankexp=rankexp[names(rankexp)%in%rownames(bc_dox_u) & !is.na(bc_dox_u[names(rankexp),'barcode'])]
  if(length(rankexp)<3){
    print('Error: cannot find enough cells.')
    return (NULL)
  }
  rankexp=sort(rankexp)
  
  texp_guide_ass=bc_dox_u[names(rankexp),'barcode']
  texp_gene_ass=bc_dox_u[names(rankexp),'gene']
  texp_guide_ass1=paste(texp_guide_ass,1:length(texp_guide_ass),sep='_r')
  
  rra_oframe=data.frame(guide=texp_guide_ass1,gene=texp_gene_ass,
                        list=rep('list',length(texp_guide_ass)),value=rankexp,
                        prob=rep(1,length(texp_guide_ass)),chosen=rep(1,length(texp_guide_ass)))
  low_file=paste(tmpprefix,'_rra_low.txt',sep='')
  write.table(rra_oframe,file=low_file,row.names = F,quote=F,sep='\t')
  
  
  rra_oframe_h=rra_oframe
  rra_oframe_h[,'value']=-1*rra_oframe_h[,'value']
  rra_oframe_h=rra_oframe_h[order(rra_oframe_h[,'value']),]
  high_file=paste(tmpprefix,'_rra_high.txt',sep='')
  write.table(rra_oframe_h,file=high_file,row.names = F,quote=F,sep='\t')
  ngguidefile=paste(tmpprefix,'_negctrl.txt',sep='')
  
  if(!is.null(negctrlgenelist)){
    ngguidelist=texp_guide_ass1[texp_gene_ass%in%negctrlgenelist]
    write.table(ngguidelist,file=ngguidefile,sep='\t',row.names = F,col.names = F,quote = F)
    ngguidecommand=paste('--control',ngguidefile)
  }else{
    ngguidecommand=''
  }
  if(is.null(rrapath)){
    rracommand='RRA'
  }else{
    rracommand=rrapath
  }
  
  rra_low_out=paste(tmpprefix,'_rra_low.out',sep='')
  rra_c=paste(rracommand,'-i', low_file,
              '-o',rra_low_out,
              ngguidecommand,
              '-p',pcutoff,
              '--max-sgrnapergene-permutation 10000 ',more_rra)
  if(negsel){
  print(rra_c)
  system(rra_c,ignore.stdout = TRUE,ignore.stderr = TRUE)
  }
  
  rra_high_out=paste(tmpprefix,'_rra_high.out',sep='')
  rra_c=paste(rracommand,'-i', high_file,
              '-o',rra_high_out,
              ngguidecommand,
              '-p',pcutoff,
              '--max-sgrnapergene-permutation 10000 ',more_rra)
  if(possel){
  print(rra_c)
  system(rra_c,ignore.stdout = TRUE,ignore.stderr = TRUE)
  }
  
  # merge both
  if(negsel){
    frame_l=read.table(rra_low_out,header = T,as.is = T,row.names = 1,na.strings = '')
  }
  if(possel){
    frame_h=read.table(rra_high_out,header = T,as.is = T,row.names = 1,na.strings = '')
  }
  
  if(negsel & !possel){
    system(paste('rm',low_file,rra_low_out))
    if(!is.null(negctrlgenelist)){
      system(paste('rm',ngguidefile))
    }
    return (frame_l)
  }
  if(!negsel & possel){
    system(paste('rm',high_file,rra_high_out))
    if(!is.null(negctrlgenelist)){
      system(paste('rm',ngguidefile))
    }
    return (frame_h)
  }
  report_f=merge(frame_l,frame_h,by=0,suffixes=c('.low','.high'))
  
  system(paste('rm',low_file,high_file,rra_low_out,rra_high_out))
  if(!is.null(negctrlgenelist)){
    system(paste('rm',ngguidefile))
  }
  return (report_f)
}


test_sc_in_gsea<-function(targetobj, gs_test=gs_c2_exps, select_gs_names=NULL,padjcutoff=0.01,gsea_cutoff=0.3,negctrlgenelist='NonTargetingControlGuideForHuman'){
  #select_gs_names=grep('^KEGG|^BIOCARTA|^REACTOME|^PID',gs_c2$names,invert = T)
  #padjcutoff=0.01
  if(is.null(select_gs_names)){
    select_gs_names=1:length(gs_test$names)
  }
  r_cs=list()
  r_cs$rpt=list()
  r_cs$dir=c()
  r_cs$pathway=c()

  scalef=getscaledata(targetobj)
  
  for(gs_pw in select_gs_names){
    
    search_ind=which(select_gs_names==gs_pw)
    if(search_ind%%10==1){
      print(paste(search_ind,'/',length(select_gs_names)))
    }
    
    gs_target=gs_test$geneset[gs_pw]
    gs_name=gs_test$names[gs_pw]
    print(gs_name)
    #target_gene='MKI67'
    texp_mat=scalef[rownames(scalef)%in%gs_target[[1]],]
    texp=colSums(texp_mat)
    texp=sort(texp)
    
    #texp_sel_cells=which(names(texp)%in%rownames(bc_dox_uq) & !is.na(bc_dox_uq[names(texp),'oligo']))
    texp_sel_cells=which(!is.na(targetobj@meta.data[names(texp),'target.gene']))
    texp_withg=texp[texp_sel_cells]
    #texp_guide_ass=bc_dox_uq[names(texp_withg),'oligo']
    #texp_gene_ass=bc_dox_uq[names(texp_withg),'gene']
    texp_guide_ass=targetobj@meta.data[names(texp_withg),'target.gene']
    texp_gene_ass=targetobj@meta.data[names(texp_withg),'geneID']
    
    
    # negative
    
    genes_to_rank=(texp_gene_ass)
    report_neg=get_rank_tables(genes_to_rank,negctrlgenelist=negctrlgenelist)
    
    
    # positive
    
    genes_to_rank=rev(texp_gene_ass)
    report_pos=get_rank_tables(genes_to_rank,negctrlgenelist=negctrlgenelist)
    
    nsel=which(report_neg$gsea>gsea_cutoff & 
                 #report_neg$gene!='TP53'& 
                 (  report_neg$ks_p_adj<padjcutoff | report_neg$ks_negctrl_adj<padjcutoff))
    if(length(nsel)>0){
      
      print(paste(gs_name,'negative'))
      print(report_neg[nsel,])
      r_cs$rpt=rbind(r_cs$rpt,report_neg[nsel,])
      r_cs$dir=c(r_cs$dir,rep('negative',length(nsel)))
      r_cs$pathway=c(r_cs$pathway,rep(gs_name,length(nsel)))
    }
    
    nsel=which(report_pos$gsea>gsea_cutoff& 
                 #report_pos$gene!='TP53' & 
                 (report_pos$ks_p_adj<padjcutoff| report_pos$ks_negctrl_adj<padjcutoff))
    if(length(nsel)>0){
      
      print(paste(gs_name,'positive'))
      print(report_pos[nsel,])
      
      
      r_cs$rpt=rbind(r_cs$rpt,report_pos[nsel,])
      r_cs$dir=c(r_cs$dir,rep('positive',length(nsel)))
      
      r_cs$pathway=c(r_cs$pathway,rep(gs_name,length(nsel)))
    }
    
  }
  
  report_fall=r_cs$rpt
  report_fall[,'direction']=r_cs$dir
  report_fall[,'pathway']=r_cs$pathway
  report_fall[,'ks_p_adj.2']=p.adjust(report_fall$ks_p,method='fdr')
  report_fall[,'ks_negctrl_adj.2']=p.adjust(report_fall$ks_negctrl,method='fdr')
  return (report_fall)
}


test_sc_in_rra<-function(targetobj,  gs_test=gs_c2_exps, select_gs_names=NULL,padjcutoff=0.25,gsea_cutoff=0.3,negctrlgenelist='NonTargetingControlGuideForHuman',rra_path=NULL,tmpprefix='sample1'){
  #select_gs_names=grep('^KEGG|^BIOCARTA|^REACTOME|^PID',gs_c2$names,invert = T)
  #padjcutoff=0.01
  if(is.null(select_gs_names)){
    select_gs_names=1:length(gs_test$names)
  }
  r_cs=list()
  r_cs$rpt=list()
  r_cs$dir=c()
  r_cs$pathway=c()


  scalef=getscaledata(targetobj)

 
  bc_dx_uq=colnames(scalef)
  bc_guide_ass=targetobj@meta.data[bc_dx_uq,'target.gene']
  bc_gene_ass=targetobj@meta.data[bc_dx_uq,'geneID']
  
  bc_dx_uq=bc_dx_uq[!is.na(bc_gene_ass)]
  bc_guide_ass=bc_guide_ass[!is.na(bc_gene_ass)]
  bc_gene_ass=bc_gene_ass[!is.na(bc_gene_ass)]
  
  bc_dx_ttb=data.frame(cell=bc_dx_uq,barcode=bc_guide_ass,gene=bc_gene_ass)
  rownames(bc_dx_ttb)=bc_dx_uq
  
  for(gs_pw in select_gs_names){
    
    search_ind=which(select_gs_names==gs_pw)
    if(search_ind%%10==1){
      print(paste(search_ind,'/',length(select_gs_names)))
    }
    
    gs_target=gs_test$geneset[gs_pw]
    gs_name=gs_test$names[gs_pw]
    print(gs_name)
    #target_gene='MKI67'
    texp_mat=scalef[rownames(scalef)%in%gs_target[[1]],]
    texp=colSums(texp_mat)
    texp=sort(texp)
    
    rept_frm=get_rank_tables_from_rra(texp,bc_dx_ttb,rrapath = rra_path,negctrlgenelist = negctrlgenelist,tmpprefix=tmpprefix )
    
    
    nsel=which(#report_neg$gene!='TP53'& 
      (  rept_frm$FDR.low<padjcutoff))
    if(length(nsel)>0){
      
      print(paste(gs_name,'negative'))
      print(rept_frm[nsel,])
      r_cs$rpt=rbind(r_cs$rpt,rept_frm[nsel,])
      r_cs$dir=c(r_cs$dir,rep('negative',length(nsel)))
      r_cs$pathway=c(r_cs$pathway,rep(gs_name,length(nsel)))
    }
    
    nsel=which(#report_pos$gene!='TP53' & 
      (rept_frm$FDR.high<padjcutoff))
    if(length(nsel)>0){
      
      print(paste(gs_name,'positive'))
      print(rept_frm[nsel,])
      
      
      r_cs$rpt=rbind(r_cs$rpt,rept_frm[nsel,])
      r_cs$dir=c(r_cs$dir,rep('positive',length(nsel)))
      
      r_cs$pathway=c(r_cs$pathway,rep(gs_name,length(nsel)))
    }
    
  }
  
  report_fall=r_cs$rpt
  report_fall[,'direction']=r_cs$dir
  report_fall[,'pathway']=r_cs$pathway
  #report_fall[,'ks_p_adj.2']=p.adjust(report_fall$ks_p,method='fdr')
  #report_fall[,'ks_negctrl_adj.2']=p.adjust(report_fall$ks_negctrl,method='fdr')
  return (report_fall)
}

get_sc_signature<-function(targetobj, gs_test=gs_c2_exps, padjcutoff=0.01,gsea_cutoff=0.3,negctrlgenelist='NonTargetingControlGuideForHuman'){
  #select_gs_names=grep('^KEGG|^BIOCARTA|^REACTOME|^PID',gs_c2$names,invert = T)
  #padjcutoff=0.01
  select_gs_names=gs_test$names

  scalef=getscaledata(targetobj)

  r_cs=matrix(rep(0,length(select_gs_names)*ncol(scalef)),nrow=length(select_gs_names))
  rownames(r_cs)=select_gs_names
  colnames(r_cs)=colnames(scalef)
  for(gs_pw in 1:length(select_gs_names)){
    
    #search_ind=which(select_gs_names==gs_pw)
    if(gs_pw%%10==1){
      print(paste(gs_pw,'/',length(select_gs_names)))
    }
    
    gs_target=gs_test$geneset[gs_pw]
    gs_name=gs_test$names[gs_pw]
    #print(gs_name)
    #target_gene='MKI67'
    texp_mat=scalef[rownames(scalef)%in%gs_target[[1]],]
    texp=colMeans(texp_mat)
    r_cs[gs_pw,]=texp
    #texp=sort(texp)
  }
  return (r_cs)
  
}

