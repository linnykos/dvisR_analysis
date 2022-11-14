rm(list=ls())
library(Seurat)
library(npregfast)
library(alphahull)
library(dbscan)
library(R.utils)

load("../../../../../out/dvisR_analysis/kevin_tmp/Writeup1/Writeup1_bm_saver.RData")
# load("../../../out/dvisR_analysis/kevin_tmp/Writeup1/Writeup1_bm_saver.RData")
source("criterion.R")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

mat <- t(saver_res$estimate)
quantile(apply(mat, 2, max))
length(which(mat >= 10))/prod(dim(mat))
length(which(mat >= 100))/prod(dim(mat))
mat <- pmin(mat, 100)
quantile(mat)
quantile(apply(mat, 2, max))
table(bm$celltype.l2)

mat <- mat[which(bm$celltype.l2 == "CD4 Naive"),]
quantile(apply(mat, 2, max))

.nonzero_col <- function(mat, col_idx, bool_value){
  stopifnot(inherits(mat, "dgCMatrix"), col_idx %% 1 == 0,
            col_idx > 0, col_idx <= ncol(mat))
  
  val1 <- mat@p[col_idx]
  val2 <- mat@p[col_idx+1]
  
  if(val1 == val2) return(numeric(0))
  if(bool_value){
    # return the value
    mat@x[(val1+1):val2]
  } else {
    # return the row index
    mat@i[(val1+1):val2]+1
  }
}

tmp <- Matrix::t(bm[["RNA"]]@counts[Seurat::VariableFeatures(bm),])
nonzero_percentage_vec <- sapply(1:ncol(tmp), function(j){
  length(.nonzero_col(tmp, col_idx = j, bool_value = F))/nrow(tmp)
})
names(nonzero_percentage_vec) <- colnames(tmp)
round(100*quantile(nonzero_percentage_vec))
length(which(nonzero_percentage_vec >= 0.2))
mat <- mat[,names(nonzero_percentage_vec)[which(nonzero_percentage_vec >= 0.2)]]
dim(mat)
mat <- mat[,which(matrixStats::colVars(mat) >= 1e-4)]
mat <- scale(mat)
rm(list = c("saver_res", "bm"))
gc(T)

#######

criterion_all <- function(dat){
  dbscan_val <- criterion_dbscan_number(dat)
  kl_val <- criterion_KLdiv_diagonalGaussian(dat)
  ks_val <- criterion_kolmogrov(dat)
  nlm_val <- tryCatch({criterion_linear_vs_nonparametric_fit(dat)}, error = function(e){NA})
  
  pearson_val <- correlation_pearson(dat)
  spearman_val <- correlation_spearman(dat)
  kendall_val <- correlation_kendall(dat)
  xi_val <- correlation_xi(dat)
  beta_val <- correlation_beta(dat)
  energy_val <- correlation_energy(dat)
  
  vec <- c( beta = beta_val,
            dbscan = dbscan_val,
            energy = energy_val,
            kendall = kendall_val,
            kl = kl_val,
            ks = ks_val,
            nlm = nlm_val,
            pearson = pearson_val,
            spearman = spearman_val,
            xi = xi_val)
  vec
}

# dat <- mat[,c(1:2)]; tmp <- criterion_all(dat)

#######

set.seed(10)
num_pairs <- 5000
p <- ncol(mat)
sampling_pairs <- cbind(sample(1:p, size = num_pairs, replace = T), 
                        sample(1:p, size = num_pairs, replace = T))
sampling_pairs <- sampling_pairs[apply(sampling_pairs, 1, function(x){abs(diff(x))})!=0,]
sampling_pairs <- t(apply(sampling_pairs, 1, function(x){sort(x)}))
mapping_val <- sampling_pairs[,1] + sampling_pairs[,2]*(p+1)
rm_idx <- which(duplicated(mapping_val))
if(length(rm_idx) > 0) sampling_pairs <- sampling_pairs[-rm_idx,]
mapping_val <- sampling_pairs[,1] + sampling_pairs[,2]*(p+1)
sampling_pairs <- sampling_pairs[order(mapping_val, decreasing = F),]
dim(sampling_pairs)

num_pairs <- nrow(sampling_pairs)
criterion_mat <- matrix(NA, nrow = num_pairs, ncol = 10)
for(i in 1:num_pairs){
  print(i)
  criterion_mat[i,] <- criterion_all(mat[,sampling_pairs[i,]])
  
  if(i %% 100 == 0) {
    print("Saving")
    save(mat, criterion_mat, sampling_pairs,
         date_of_run, session_info,
         file = "../../../../../out/dvisR_analysis/kevin_tmp/Writeup1/Writeup1_initial_scan.RData")
  }
}

save(mat, criterion_mat, sampling_pairs,
     date_of_run, session_info,
     file = "../../../../../out/dvisR_analysis/kevin_tmp/Writeup1/Writeup1_initial_scan.RData")


