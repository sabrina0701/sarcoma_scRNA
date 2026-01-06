# ******************************************************************************************** 
# Identifying heterogeneity programs using nonnegative matrix factorization - human tumors
# (i.e. continuous + discrete heterogeneity) 

library(Seurat)


source("/home/d_qian/scripts/nmf_programs.R")
setwd("/mnt/beegfs01/scratch/d_qian/SingleCell")

tumor<-readRDS("/mnt/beegfs01/scratch/d_qian/SingleCell/combined_datasets/tumor.rds")
Dong<-subset(tumor,study=="Dong_NB")


# perform NMF with ranks ranging from 4 to 9       
w_basis_tumor <- list() # nmf gene scores
h_coef_tumor <- list() # nmf cell scores


expr_tumor<-SplitObject(Dong,split.by="sample_ID")

for(i in names(expr_tumor)) {
  print(i)
  w <- NULL
  h <- NULL
  sparse_matrix <- GetAssayData(expr_tumor[[i]], assay = "RNA", layer = "counts")
  sparse_matrix<-as.matrix(sparse_matrix)
  #nonzero <- sparse_matrix  > 0
  #keep_genes<-names(sort(Matrix::rowSums(nonzero),decreasing= T)[1:10000])
  #sparse_matrix <- sparse_matrix[keep_genes, ]
  sparse_matrix <- t(t(sparse_matrix)/colSums(sparse_matrix))*1000000
  for(j in 4:9)  {
    n <- nmf_programs( sparse_matrix, is.log=F, rank=j)
    colnames(n$w_basis) <- paste0(i, "_rank4_9_nruns10.RDS.", j, ".", 1:j)
    colnames(n$h_coef) <- paste0(i, "_rank4_9_nruns10.RDS.", j, ".", 1:j)
    w <- cbind(w, n$w_basis)
    h <- cbind(h, n$h_coef)
  }
  
  w_basis_tumor[[i]] <- w
  h_coef_tumor[[i]] <- h
}

# save output
saveRDS(w_basis_tumor, "nmf/Zhou/nmf_w_basis_tumor.RDS")
saveRDS(h_coef_tumor, "nmf/Zhou/nmf_h_coef_tumor.RDS")



