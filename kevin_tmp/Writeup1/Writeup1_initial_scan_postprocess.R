rm(list=ls())
load( "../../../out/dvisR_analysis/kevin_tmp/Writeup1/Writeup1_initial_scan.RData")

colnames(criterion_mat) <- c("beta", "dbscan",
                             "energy", "kendall",
                             "kl", "ks",
                             "nlm", "pearson",
                             "spearman", "xi")
summary(criterion_mat)
image(cor(criterion_mat))
image(cor(criterion_mat[,c("beta", "energy", "kendall", "pearson", "spearman", "xi")]))

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
ordering_1 <- order(abs(criterion_mat[,"pearson"]), decreasing = T)/n
ordering_2 <- order(criterion_mat[,"xi"], decreasing = T)/n
plot(ordering_1, ordering_2, asp = T, col = rgb(0.5,0.5,0.5,0.2), pch = 16)

which_pairs <- which(ordering_1 <= 0.05, ordering_2 >= 0.95)
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

