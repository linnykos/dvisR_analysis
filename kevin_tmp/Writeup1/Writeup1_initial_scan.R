rm(list=ls())
library(Seurat)
library(np)
library(npregfast)
library(alphahull)
library(dbscan)
library(R.utils)

load("../../../../../out/dvisR_analysis/kevin_tmp/Writeup1/Writeup1_bm_saver.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

mat <- t(saver_res$estimate)
mat <- pmin(mat, 10)
quantile(mat)
quantile(apply(mat, 2, max))
table(bm$celltype.l2)

mat <- mat[which(bm$celltype.l2 == "CD4 Naive"),]
quantile(apply(mat, 2, max))
rm(list = c("saver_res", "bm"))
gc(T)

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

#######

criterion_dbscan_number <- function(dat){
  res <- dbscan::dbscan(dat, eps = 0.1, minPts = 10, borderPoints = F)
  max(res$cluster)
}

criterion_kolmogrov <- function(dat){
  dat[,1] <- jitter(dat[,1])
  dat[,2] <- jitter(dat[,2])
  as.numeric(stats::ks.test(dat[,1], "pnorm")$statistic +
               stats::ks.test(dat[,2], "pnorm")$statistic)
}

criterion_KLdiv_diagonalGaussian <- function(dat){
  cov.full <- stats::cov(dat)
  cov.diag <- diag(c(stats::var(dat[,1]), stats::var(dat[,2])))
  
  res <- .5 * (.trace(base::solve(cov.diag) %*% cov.full) +
                 log(base::det(cov.diag) / base::det(cov.full)))
  
  ifelse(res > 1e4, 1e4, res)
}

.trace <- function(mat){sum(diag(mat))}

.areaahull <- function (x, timeout = 5) {
  area <- R.utils::withTimeout(try(alphahull::areaahulleval(x), silent = TRUE),
                               timeout = timeout)
  if (!is.numeric(area)) {
    warning("Problem in area computation (Returns NA)")
    area <- NA
  }
  if (is.numeric(area) & area < 0) {
    warning("Problem in area computation (Returns NA)")
    area <- NA
  }
  return(area)
}

criterion_alpha_area <- function(dat){
  dat[,1] <- jitter(dat[,1])
  dat[,2] <- jitter(dat[,2])
  tot <- prod(apply(apply(dat, 2, range), 2, diff))
  a <- .areaahull(alphahull::ahull(dat, alpha = 0.5))
  a/tot
}

criterion_linear_vs_nonparametric_fit <- function(dat){
  dat[,1] <- jitter(dat[,1])
  dat[,2] <- jitter(dat[,2])
  
  tmp <- data.frame(x = dat[,1], y = dat[,2])
  nlm_res1 <- npregfast::frfast(y ~ x, data = tmp)
  nlm_vec1 <- nlm_res1$p[,1,1]
  lm_res1 <- stats::lm(y ~ x, data = tmp)
  lm_vec1 <- stats::predict(lm_res1, newdata = data.frame(x = nlm_res1$x))
  
  nlm_res2 <- npregfast::frfast(x ~ y, data = tmp)
  nlm_vec2 <- nlm_res2$p[,1,1]
  lm_res2 <- stats::lm(x ~ y, data = tmp)
  lm_vec2 <- stats::predict(lm_res2, newdata = data.frame(y = nlm_res2$x))
  
  mean((lm_vec1 - nlm_vec1)^2) + mean((lm_vec2 - nlm_vec2)^2)
}

criterion_all <- function(dat){
  cat("Dbscan")
  dbscan_val <- criterion_dbscan_number(dat)
  cat("Ks")
  ks_val <- criterion_kolmogrov(dat)
  cat("KL")
  kl_val <- criterion_KLdiv_diagonalGaussian(dat)
  cat("Alpha")
  alpha_val <- criterion_alpha_area(dat)
  cat("NLM")
  nlm_val <- tryCatch({criterion_linear_vs_nonparametric_fit(dat)}, error = function(e){NA})
  
  vec <- c(alpha = alpha_val,
           dbscan = dbscan_val,
           kl = kl_val,
           ks = ks_val,
           nlm = nlm_val)
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
criterion_mat <- matrix(NA, nrow = num_pairs, ncol = 5)
for(i in 1:num_pairs){
  print(i)
  criterion_mat[i,] <- criterion_all(mat[,sampling_pairs[i,]])
  
  if(i %% 100 == 0) {
    print("Saving")
    save(mat, criterion_mat,
         date_of_run, session_info,
         file = "../../../../../out/dvisR_analysis/kevin_tmp/Writeup1/Writeup1_initial_scan.RData")
  }
}

save(mat, criterion_mat,
     date_of_run, session_info,
     file = "../../../../../out/dvisR_analysis/kevin_tmp/Writeup1/Writeup1_initial_scan.RData")




