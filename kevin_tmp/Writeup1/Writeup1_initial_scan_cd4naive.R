rm(list=ls())
library(dbscan)
library(energy)
library(Hmisc)
library(mclust)
library(minerva)
library(npregfast)
library(pcaPP)
library(R.utils)
library(Seurat)
library(TauStar)
library(XICOR)

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

celltype <- "CD4 Naive"
mat <- mat[which(bm$celltype.l2 == celltype),]
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

tmp <- Matrix::t(bm[["RNA"]]@counts[Seurat::VariableFeatures(bm),which(bm$celltype.l2 == celltype)])
nonzero_percentage_vec <- sapply(1:ncol(tmp), function(j){
  length(.nonzero_col(tmp, col_idx = j, bool_value = F))/nrow(tmp)
})
names(nonzero_percentage_vec) <- colnames(tmp)
round(100*quantile(nonzero_percentage_vec))
length(which(nonzero_percentage_vec >= 0.1))
mat <- mat[,names(nonzero_percentage_vec)[which(nonzero_percentage_vec >= 0.1)]]
mat <- mat[,which(matrixStats::colVars(mat) >= 1e-4)]
mat <- scale(mat)
dim(mat)
rm(list = c("saver_res", "bm"))
gc(T)

#######

criterion_all <- function(dat){
  dbscan_val <- criterion_dbscan_number(dat)
  em_val <- correlation_em(dat)
  kl_val <- criterion_KLdiv_diagonalGaussian(dat)
  ks_val <- criterion_kolmogrov(dat)
  nlm_val <- criterion_linear_vs_nonparametric_fit(dat)
  
  beta_val <- correlation_beta(dat)
  energy_val <- correlation_energy(dat)
  hoeffd_val <- correlation_hoeffd(dat)
  hsic_val <- correlation_hsic(dat)
  kendall_val <- correlation_kendall(dat)
  mic_val <- correlation_mic(dat)
  nmi_val <- correlation_nmi(dat)
  pearson_val <- correlation_pearson(dat)
  spearman_val <- correlation_spearman(dat)
  taustar_val <- correlation_taustar(dat)
  xi_val <- correlation_xi(dat)
  
  vec <- c(beta = beta_val,
           dbscan = dbscan_val,
           em = em_val,
           energy = energy_val,
           hsic = hsic_val,
           hoeff = hoeffd_val,
           kendall = kendall_val,
           kl = kl_val,
           ks = ks_val,
           mic = mic_val,
           nlm = nlm_val,
           nmi = nmi_val,
           pearson = pearson_val,
           spearman = spearman_val,
           taustar = taustar_val,
           xi = xi_val)
  vec
}

# dat <- mat[,c(1:2)]; tmp <- criterion_all(dat)

#######

set.seed(10)
num_pairs <- 10000
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
criterion_mat <- matrix(NA, nrow = num_pairs, ncol = 16)
for(i in 1:num_pairs){
  print(i)
  tmp <- criterion_all(mat[,sampling_pairs[i,]])
  if(i == 1) colnames(criterion_mat) <- names(tmp)
  criterion_mat[i,] <- tmp
  
  if(i %% 100 == 0) {
    print("Saving")
    save(mat, criterion_mat, sampling_pairs,
         date_of_run, session_info,
         file = paste0("../../../../../out/dvisR_analysis/kevin_tmp/Writeup1/Writeup1_initial_scan_", celltype, ".RData"))
  }
}

save(mat, criterion_mat, sampling_pairs,
     date_of_run, session_info,
     file = paste0("../../../../../out/dvisR_analysis/kevin_tmp/Writeup1/Writeup1_initial_scan_", celltype, ".RData"))


