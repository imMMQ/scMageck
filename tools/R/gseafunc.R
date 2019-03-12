
# gene sets and functions

# read a geneset from GCT file
readgmtfile<-function(filename){
  geneset=strsplit(as.character(read.table(filename,sep=' ')[,1]), '\t')
  geneset_names=unlist(lapply(geneset,function(x){return (x[1])}))
  geneset_sizes=unlist(lapply(geneset,length))-2
  gs=list();
  gs$geneset=geneset;
  gs$names=geneset_names;
  gs$sizes=geneset_sizes;
  return (gs);
}

gs_all=readgmtfile('msigdb.v4.0.symbols.gmt')
gs_c2=readgmtfile('c2.all.v4.0.symbols.gmt')
# calculate the normalized t value
# input parameter:
#   ntval: a ranked list of t-value, with its name as gene name
#   geneset: the list of the query gene set
#   method: the method to compute p-value. 'normal': use normal distribution; others: use chi-sq distribution
getgstval<-function(ntval,geneset,method='normal'){
    gset2=geneset[!is.na(geneset) & geneset%in% names(ntval)];
    n=length(gset2)
    if(n==0){
        print('Error: the length of effective gene set is zero.');
        return (0);
    }
    if(method=='normal'){
      nt=sum(ntval[gset2])/sqrt(n);
      return (nt);
    }else{
      mv=mean(ntval[gset2]);
      nt=(sum((ntval[gset2]-mv)^2)-(n-1))/(2*(n-1));
      return (nt);
    }
}



# for a given vector of t-value (with the gene name as its name), return GSEA results
# parameters:
# genetval: the t-value for each gene, with its name as gene name
# ranktransform: whether to convert t-values into a vector of norm distributions
# standarize: whether to standarize the gene expression
# genesetname: the name of the gene set. 'all': all gene set; 'c2': c2 gene set. others: provided in givengeneset
getgseaframe<-function(genetval,ranktransform=F,standardize=T,genesetname='c2',givengeneset=NULL,sorted=F){

  if(genesetname=='all'){
    gsobj=gs_all;
  }else if(genesetname=='c2'){
    gsobj=gs_c2;
  }else{
    gsobj=givengeneset;
    if(is.null(givengeneset)){
      print('Error: gene set must be provided!')
      return (0);
    }
  }
  
  gsnames=gsobj$names;
  gssz=gsobj$sizes;
  genesets=gsobj$geneset;
  
  # account for adjustments of the correlation
  if("adj" %in% names(gsobj)){
    gsadj=gsobj$adj;
  }else{
    gsadj=rep(0,length(gssz));
  }
  genetval_norm=genetval[!is.na(genetval)];
  if(ranktransform==T){
      ux=genetval;
      ux[order(ux)]=sort(rnorm(length(ux)))
      genetval_norm=ux;
  }
  if(standardize==T){
      genetval_norm=(genetval_norm-mean(genetval_norm,na.rm=T))/sqrt(var(genetval_norm,na.rm=T))
  }
  genesetsnt=unlist(lapply(genesets,function(x){ return (getgstval(genetval_norm,x))}))
  usetval=genesetsnt;
  if(standardize==T){
      # why standarize again?
      # usetval=(usetval-mean(usetval))/sqrt(var(usetval))
  }

  adjsd=sqrt(1+gsadj);
  genesetsframe=data.frame(sig=gsnames,size=gssz,
  normt=usetval,
  pval=pnorm(abs(usetval),lower.tail=F),
  logp=-1*pnorm(abs(usetval),lower.tail=F,log.p=T),
  adjpval=2*pnorm(abs(usetval),sd=adjsd,lower.tail=F),
  adjlogp=-1*pnorm(abs(usetval),sd=adjsd,lower.tail=F,log.p=T),
  adjp.up=pnorm((usetval),sd=adjsd,lower.tail=F),
  adjp.down=pnorm((usetval),sd=adjsd,lower.tail=T),
  adjweight=gsadj)
  
  rownames(genesetsframe)=genesetsframe[,'sig']
  genesetsframe[,'adj.p']=p.adjust(genesetsframe[,'adjpval'])
  
  if(sorted){
    genesetsframe=genesetsframe[order(-genesetsframe[,'logp']),]
  }
  return (genesetsframe);
}


# for a given gene rank, perform GSEA test
getgsearankframe<-function(generank,genesetname='c2',givengeneset=NULL,sorted=F){
  z_tval=1:length(generank);
  names(z_tval)=toupper(as.character(generank));
  z_tval=sort(z_tval);
  ztb_gsea=getgseaframe(z_tval,ranktransform=T,genesetname=genesetname,givengeneset=givengeneset,sorted=sorted)
  return (ztb_gsea);
}


# perform GSEA test on a matrix of t-values
# parameters:
#   tval: return the normalized t value instead of p-value
tmat_test<-function(genesetobj,tmat,rankT=F,std=T,correct=T,tval=F){
  ntk=apply(tmat,2,function(query){
    tf=getgseaframe(query,ranktransform=rankT,standardize=std,genesetname='',givengeneset=genesetobj,sorted=F);
    
    retnn= (tf[,'pval'])
    if(tval){
      retnn=(tf[,'normt'])
    }
    names(retnn)=tf[,'sig'];
    return (retnn)
  })
  
  if(correct & tval==F){
    gsmat0=p.adjust(ntk)
    gsmat0=matrix(gsmat0,ncol=ncol(ntk))
    rownames(gsmat0)=rownames(ntk)
    colnames(gsmat0)=colnames(ntk)
    ntk=gsmat0
  }
  return (ntk);
}


# perform a overlap test based on the gene set and gene query
# parameter:
#   geneset: a list of genes
#   genequery: a list of gene sets
#   ntotal: the total number of genes in the query
#   method: hyper (hypergeometric), chi (chi-squared) or fisher (fisher's exact) test
#   lowert: whether to perform under-represented p-values (not valid for chi-squared test)
#   noprint: do not print debug information
overlaptest<-function(geneset,genequery,ntotal,method='hyper',lowert=F,logp=F,noprint=F){
  # for fisher's exact test
  a=length(intersect(geneset,genequery))
  b=length(genequery)-a;
  c=length(geneset)-a;
  d=ntotal-b-c+a;
  rbmat=rbind(c(a,b),c(c,d))
  # for hypergenomic test
  x=a;
  k=a+b; # total query
  m=a+c; # total geneset
  n=b+d; # total genes not in the geneset 
  if(method=='chi'){
    zp=chisq.test(rbmat)
    t=zp$p.value;
    if(!noprint){
      print(paste('Chi-squared test: a=',a,',b=',b,',c=',c,',d=',d,',p=',t,',lower=',lowert))
    }
    if(logp){
      t=log(t);
    }
  }else if(method=='fisher'){
    if(lowert){
      zp=fisher.test(rbmat,alternative='less');
    }else{
      zp=fisher.test(rbmat,alternative='greater');
    }
    t=zp$p.value;
    if(!noprint){
      print(paste('Fishers exact test: a=',a,',b=',b,',c=',c,',d=',d,',p=',t,',lower=',lowert))
    }
    if(logp){
      t=log(t)
    }
  }else{
    t= (phyper(x,m,n,k,lower.tail=lowert,log.p=logp))
    if(!noprint){
      print(paste('Hypergeometric test: x=',x,',k=',k,',m=',m,',n=',n,',p=',t,',lower=',lowert))
    }
  }
  return (t);
}



# for a given vector genes, return the overlap test results
# parameters:
#   querygene: the query geneset
#   ntotal: the total # of genes
#   genesetname: the name of the gene set. 'all': all gene set; 'c2': c2 gene set. others: provided in givengeneset
#   testmethod: method for testing
#   lowertest: whether to perform under-represented test
#   logtransform: whether p-value is log-transformed?
#   sorted: whether to sort the frame according to p-value
getoverlaptestframe<-function(querygene,ntotal,genesetname='c2',givengeneset=NULL,
                              testmethod='hyper',lowertest=F,logtransform=F,sorted=F){
  
  if(genesetname=='all'){
    gsnames=gs_all$names;
    gssz=gs_all$sizes;
    genesets=gs_all$geneset;
  }else if(genesetname=='c2'){
    gsnames=gs_c2$names;
    gssz=gs_c2$sizes;
    genesets=gs_c2$geneset;
  }else{
    if(is.null(givengeneset)){
      print('Error: gene set must be provided!')
      return (0);
    }
    gsnames=givengeneset$names;
    gssz=givengeneset$sizes;
    genesets=givengeneset$geneset;
  }
  
  genepval=unlist(lapply(genesets,function(x){ 
      return (overlaptest(x,querygene,ntotal,method=testmethod,lowert=lowertest,logp=logtransform))
    }))
  
  genesetsframe=data.frame(sig=gsnames,size=gssz,
                           pval=genepval)
  rownames(genesetsframe)=genesetsframe[,'sig']
  genesetsframe[,'adj.p']=p.adjust(genesetsframe[,'pval'])
  
  if(sorted){
    genesetsframe=genesetsframe[order(-genesetsframe[,'adj.p']),]
  }
  return (genesetsframe);
}



# perform non-GSEA test between a given gene set and a binary matrix
# can be test on both over-represented and under-represented
binmat_test<-function(genesetobj,binmat,overrepresented=T,meth='hyper',correct=T){
  
  ntk=(apply(binmat,2,function(query){
    qgene=rownames(binmat)[query!=0];
    retfm= (getoverlaptestframe(qgene,nrow(binmat),genesetname='',givengeneset=genesetobj,
                                lowert=!overrepresented,testmethod=meth));
    
    retnn=retfm[,'pval']
    names(retnn)=retfm[,'sig']
    return (retnn)
  }))
  
  if(correct){
    gsmat0=p.adjust(ntk)
    gsmat0=matrix(gsmat0,ncol=ncol(ntk))
    rownames(gsmat0)=rownames(ntk)
    colnames(gsmat0)=colnames(ntk)
    ntk=gsmat0
  }
  return (ntk)
}


###########
# get gsea score
getoriginalgseascore<-function(corvalue,geneset){
  # binary vector
  tgframeu_vec=(names(corvalue)%in%geneset)*1
  
  
  nh=sum(tgframeu_vec)
  n=length(tgframeu_vec)
  # old method?
  #incvalue=sqrt((n-nh)/nh)
  #decvalue=sqrt(nh/(n-nh))
  # new GSEA 
  incvalue=1/nh;
  decvalue=1/(n-nh);
  gsvec=rep(0,n)
  
  #corvec=abs(tgframeu[,2])/max(abs(tgframeu[,2]))
  corvec=abs(corvalue)
  sumcorins=sum(corvec[which(tgframeu_vec==1)])
  for(i in 1:n){
    incvalue=corvec[i]/sumcorins;
    if(i==1){
      if(tgframeu_vec[1]==1){
        gsvec[1]=incvalue;
      }else{
        gsvec[1]=-1*decvalue;
      }
    }else{
      if(tgframeu_vec[i]==1){
        gsvec[i]=gsvec[i-1]+incvalue;
      }else{
        gsvec[i]=gsvec[i-1]-decvalue;
      }
    }
  }
  return (gsvec);
}
# shuffle GSEA score
shufflegseascore<-function(corvalue,geneset,n=1000){
  corvname=names(corvalue);
  nval=rep(0,n);
  for(i in 1:n){
    corvalue=sample(corvalue);
    names(corvalue)=corvname;
    gz=getoriginalgseascore(corvalue,geneset);
    nval[i]=max(gz);
    if(i%%10==1){
      print(paste(i,'...'))
    }
  }
  return (nval)
}


# plot GSEA function
# parameters:
# corvalue: a numeric vector of correlations. The names must be gene names. Must be between [-1,1]
# geneset: a vector of gene names
# greyscalebar: whether the bar is plotted in greyscale?
plotrankedgsea<-function(corvalue,geneset, greyscalebar=F, steps=100, ...){
  
  # binary vector
  tgframeu_vec=(names(corvalue)%in%geneset)*1
  
  gsvec=getoriginalgseascore(corvalue,geneset)
  
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  
  nf <- layout(matrix(c(1,0,2,0),2,2,byrow=TRUE), c(7,0), c(4,1), TRUE)
  #layout.show(nf)
  
  par(mar=c(1,5,3,1))
  plot(1:length(gsvec),gsvec, xlim = c(1,length(gsvec)), type='l',lwd=4,col='green', ylab='Enrichment Score (ES)',
       xaxt='n',cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.5, ...) # 
  
  # steps=100;
  pvgroup=floor((1:length(tgframeu_vec))/steps)+1
  pvrs=tapply(tgframeu_vec,pvgroup,sum)
  if(greyscalebar==F){
    pvrs=(pvrs>0)*1
  }
  par(mar=c(1,5,1,1))
  image(as.matrix(pvrs),col=rev(gray.colors(max(pvrs)+1,start=0,end=1)),axes=F)
  box()
  
  
  par(def.par)#- reset to default
  
  return (max(gsvec))
}
