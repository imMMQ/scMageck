# Get the gene signature expression matrix from Ymat
getsigmat <- function(Ymat, gmt_file){
  sig_mat <- data.frame(row.names = rownames(Ymat))
  sig <- colnames(gmt_file)
  for (gene_sig in sig) {
    genes <- as.character(gmt_file[,gene_sig])
    if(any(genes%in%colnames(Ymat))){
      genes <- genes[genes%in%colnames(Ymat)]
      Ymat_sig <- as.data.frame(Ymat[,genes])
      Ymat_sig$m <- rowMeans(Ymat_sig)
      sig_mat <- cbind(sig_mat, as.data.frame(Ymat_sig$m)) # identify whether the genome is mouse or human
    }else{
      genes <- capitalize(tolower(genes))
      if(any(genes%in%colnames(Ymat))){
        genes <- genes[genes%in%colnames(Ymat)]
        Ymat_sig <- as.data.frame(Ymat[,genes])
        Ymat_sig$m <- rowMeans(Ymat_sig)
        sig_mat <- cbind(sig_mat, as.data.frame(Ymat_sig$m))
      }else{
        print("Didn't find any signatures in this dataset")
        break
      }
    }
  }
  if(ncol(sig_mat) > 0) {
    colnames(sig_mat) <- sig
    sig_mat <- as.matrix(sig_mat)
  }
  return(sig_mat)
}
