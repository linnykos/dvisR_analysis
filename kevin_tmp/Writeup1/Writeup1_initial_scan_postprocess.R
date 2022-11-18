rm(list=ls())
load( "../../../out/dvisR_analysis/kevin_tmp/Writeup1/Writeup1_initial_scan.RData")

summary(round(criterion_mat,2))
n <- nrow(criterion_mat)
rank_mat <- sapply(1:ncol(criterion_mat), function(i){
  rank(criterion_mat[,i])/n
})
colnames(rank_mat) <- colnames(criterion_mat)
image(abs(cor(rank_mat)), zlim = c(0,1))
par(mar = c(0,0,0,0)); image(abs(cor(rank_mat[,c("beta", "energy", "hsic",
                                                 "hoeff", "kendall", "mic",
                                                 "nmi", "pearson", "spearman", 
                                                 "taustar", "xi")])), asp = T, zlim = c(0,1))
par(mar = c(0,0,0,0)); image(abs(cor(rank_mat[,c("beta", "kendall", "pearson", "spearman")])), asp = T, zlim = c(0,1))
par(mar = c(0,0,0,0)); image(abs(cor(rank_mat[,c("energy", "hsic", "hoeff",
                                                 "mic", "nmi", "taustar", "xi")])), asp = T, zlim = c(0,1))


for(i in 1:ncol(criterion_mat)){
  idx_vec <- sampling_pairs[order(criterion_mat[,i], decreasing = T)[1:6],]
  par(mfrow = c(2,3), mar = c(4, 4, 0.5, 0.5))
  for(kk in 1:nrow(idx_vec)){
    idx <- idx_vec[kk,]
    tmp <- mat[,idx[1:2]]
    print(colnames(tmp))
    plot_density(tmp, n = 100,
                 xlab = colnames(tmp)[1], ylab = colnames(tmp)[2],
                 xlim = quantile(tmp[,1], probs = c(0.01,0.99), na.rm = T),
                 ylim = quantile(tmp[,2], probs = c(0.01,0.99), na.rm = T),
                 main = paste0(colnames(criterion_mat)[i], ": Max"),
                 point_color = grDevices::rgb(0.5, 0.5, 0.5, 0.2))
  }
  
  idx_vec <- sampling_pairs[order(criterion_mat[,i], decreasing = F)[1:6],]
  par(mfrow = c(2,3), mar = c(4, 4, 0.5, 0.5))
  for(kk in 1:nrow(idx_vec)){
    idx <- idx_vec[kk,]
    tmp <- mat[,idx[1:2]]
    print(colnames(tmp))
    plot_density(tmp, n = 100,
                 xlab = colnames(tmp)[1], ylab = colnames(tmp)[2],
                 xlim = quantile(tmp[,1], probs = c(0.01,0.99), na.rm = T),
                 ylim = quantile(tmp[,2], probs = c(0.01,0.99), na.rm = T),
                 main = paste0(colnames(criterion_mat)[i], ": Min"))
  }
}

######################################

n <- nrow(criterion_mat)
ordering_1 <- rank(abs(criterion_mat[,"kendall"]))/n
ordering_2 <- rank(abs(criterion_mat[,"xi"]))/n
plot(ordering_1, ordering_2, asp = T, col = rgb(0.5,0.5,0.5,0.2), pch = 16)

which_pairs <- which(ordering_1 <= 0.3 & ordering_2 >= 0.7)
length(which_pairs)
which_pair <- which_pairs[1]
idx <- sampling_pairs[which_pair,]
tmp <- mat[,idx[1:2]]
print(colnames(tmp))
plot_density(tmp, n = 500,
             xlab = colnames(tmp)[1], ylab = colnames(tmp)[2],
             xlim = quantile(tmp[,1], probs = c(0.01,0.99), na.rm = T),
             ylim = quantile(tmp[,2], probs = c(0.01,0.99), na.rm = T),
             point_color = grDevices::rgb(0.5, 0.5, 0.5, 0.1),
             main = paste0("Kendall: ", round(criterion_mat[which_pair,"kendall"],2), 
                           ", Quantile = ", round(ordering_1[which_pair],2), 
                           "\nXi: ", round(criterion_mat[which_pair,"xi"],2), 
                           ", Quantile = ", round(ordering_2[which_pair],2)))

########################

criterion_mat2 <- criterion_mat[,c("beta", "energy", "hsic",
                                  "hoeff", "kendall", "mic",
                                  "nmi", "pearson", "spearman", 
                                  "taustar", "xi")]
criterion_mat2 <- criterion_mat2 >= 0.5
summary(criterion_mat2)

which_pairs <- which(!criterion_mat2[,"xi"] & criterion_mat2[,"kendall"])
ordering_1 <- rank(abs(criterion_mat[,"kendall"]))/n
ordering_2 <- rank(abs(criterion_mat[,"xi"]))/n
length(which_pairs)
which_pair <- which_pairs[1]
idx <- sampling_pairs[which_pair,]
tmp <- mat[,idx[1:2]]
print(colnames(tmp))
plot_density(tmp, n = 500,
             xlab = colnames(tmp)[1], ylab = colnames(tmp)[2],
             xlim = quantile(tmp[,1], probs = c(0.01,0.99), na.rm = T),
             ylim = quantile(tmp[,2], probs = c(0.01,0.99), na.rm = T),
             point_color = grDevices::rgb(0.5, 0.5, 0.5, 0.1),
             main = paste0("Pearson: ", round(criterion_mat[which_pair,"pearson"],2), 
                           ", Quantile = ", round(ordering_1[which_pair],2), 
                           "\nXi: ", round(criterion_mat[which_pair,"xi"],2), 
                           ", Quantile = ", round(ordering_2[which_pair],2)))

