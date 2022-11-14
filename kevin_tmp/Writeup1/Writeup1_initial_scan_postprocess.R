rm(list=ls())
load( "../../../out/dvisR_analysis/kevin_tmp/Writeup1/Writeup1_initial_scan.RData")

# idx <- sampling_pairs[which.max(criterion_mat[,1]),]
# plot(mat[,idx[1]], mat[,idx[2]], pch = 16, col = rgb(0.5,0.5,0.5,0.5),
#      xlab = colnames(mat)[idx[1]], ylab = colnames(mat)[idx[2]],
#      xlim = quantile(mat[,idx[1]], probs = c(0.01,0.99), na.rm = T),
#      ylim = quantile(mat[,idx[2]], probs = c(0.01,0.99), na.rm = T))
# 
# idx <- sampling_pairs[which.min(criterion_mat[,1]),]
# plot(mat[,idx[1]], mat[,idx[2]], pch = 16, col = rgb(0.5,0.5,0.5,0.5),
#      xlab = colnames(mat)[idx[1]], ylab = colnames(mat)[idx[2]],
#      xlim = quantile(mat[,idx[1]], probs = c(0.01,0.99), na.rm = T),
#      ylim = quantile(mat[,idx[2]], probs = c(0.01,0.99), na.rm = T))

#####

idx <- sampling_pairs[which.max(criterion_mat[,2]),]
plot(mat[,idx[1]], mat[,idx[2]], pch = 16, col = rgb(0.5,0.5,0.5,0.3),
     xlab = colnames(mat)[idx[1]], ylab = colnames(mat)[idx[2]],
     xlim = quantile(mat[,idx[1]], probs = c(0.01,0.99), na.rm = T),
     ylim = quantile(mat[,idx[2]], probs = c(0.01,0.99), na.rm = T),
     main = "DBSCAN Max")

idx <- sampling_pairs[which.min(criterion_mat[,2]),]
plot(mat[,idx[1]], mat[,idx[2]], pch = 16, col = rgb(0.5,0.5,0.5,0.3),
     xlab = colnames(mat)[idx[1]], ylab = colnames(mat)[idx[2]],
     xlim = quantile(mat[,idx[1]], probs = c(0.01,0.99), na.rm = T),
     ylim = quantile(mat[,idx[2]], probs = c(0.01,0.99), na.rm = T),
     main = "DBSCAN Min")

#####

idx <- sampling_pairs[which.max(criterion_mat[,3]),]
plot(mat[,idx[1]], mat[,idx[2]], pch = 16, col = rgb(0.5,0.5,0.5,0.3),
     xlab = colnames(mat)[idx[1]], ylab = colnames(mat)[idx[2]],
     xlim = quantile(mat[,idx[1]], probs = c(0.01,0.99), na.rm = T),
     ylim = quantile(mat[,idx[2]], probs = c(0.01,0.99), na.rm = T),
     main = "KL Max")

idx <- sampling_pairs[which.min(criterion_mat[,3]),]
plot(mat[,idx[1]], mat[,idx[2]], pch = 16, col = rgb(0.5,0.5,0.5,0.5),
     xlab = colnames(mat)[idx[1]], ylab = colnames(mat)[idx[2]],
     xlim = quantile(mat[,idx[1]], probs = c(0.01,0.99), na.rm = T),
     ylim = quantile(mat[,idx[2]], probs = c(0.01,0.99), na.rm = T),
     main = "KL Min")

#####

idx <- sampling_pairs[which.max(criterion_mat[,4]),]
plot(mat[,idx[1]], mat[,idx[2]], pch = 16, col = rgb(0.5,0.5,0.5,0.3),
     xlab = colnames(mat)[idx[1]], ylab = colnames(mat)[idx[2]],
     xlim = quantile(mat[,idx[1]], probs = c(0.01,0.99), na.rm = T),
     ylim = quantile(mat[,idx[2]], probs = c(0.01,0.99), na.rm = T),
     main = "KS Max")

idx <- sampling_pairs[which.min(criterion_mat[,4]),]
plot(mat[,idx[1]], mat[,idx[2]], pch = 16, col = rgb(0.5,0.5,0.5,0.3),
     xlab = colnames(mat)[idx[1]], ylab = colnames(mat)[idx[2]],
     xlim = quantile(mat[,idx[1]], probs = c(0.01,0.99), na.rm = T),
     ylim = quantile(mat[,idx[2]], probs = c(0.01,0.99), na.rm = T),
     main = "KS Min")

#####

idx <- sampling_pairs[which.max(criterion_mat[,5]),]
plot(mat[,idx[1]], mat[,idx[2]], pch = 16, col = rgb(0.5,0.5,0.5,0.3),
     xlab = colnames(mat)[idx[1]], ylab = colnames(mat)[idx[2]],
     xlim = quantile(mat[,idx[1]], probs = c(0.01,0.99), na.rm = T),
     ylim = quantile(mat[,idx[2]], probs = c(0.01,0.99), na.rm = T),
     main = "NLM Max")

idx <- sampling_pairs[which.min(criterion_mat[,5]),]
plot(mat[,idx[1]], mat[,idx[2]], pch = 16, col = rgb(0.5,0.5,0.5,0.3),
     xlab = colnames(mat)[idx[1]], ylab = colnames(mat)[idx[2]],
     xlim = quantile(mat[,idx[1]], probs = c(0.01,0.99), na.rm = T),
     ylim = quantile(mat[,idx[2]], probs = c(0.01,0.99), na.rm = T),
     main = "NLM Min")

