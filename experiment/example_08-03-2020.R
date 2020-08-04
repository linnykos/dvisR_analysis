rm(list=ls())
set.seed(10)
mat <- matrix(rcauchy(100), 10, 10)
mat <- abs(mat)
clockwise90 = function(a) { t(a[nrow(a):1,]) } 

# "naive" way to visualize matrices
image(clockwise90(abs(mat)), asp = T)
image(clockwise90(abs(mat)), asp = T, 
      col = hcl.colors(12, "YlOrRd", rev = TRUE),
      breaks = seq(min(mat), max(mat), length.out = 13))
# the above two plots are exactly the same

plot(1:12, col = hcl.colors(12, "YlOrRd", rev = TRUE), pch = 16, cex = 5)


# "naive" way to visualize matrices
plot(sort(mat), pch = 16)
quantile(mat, probs = seq(0, 1, length.out = 13))
image(clockwise90(abs(mat)), asp = T, 
      col = hcl.colors(12, "YlOrRd", rev = TRUE),
      breaks = quantile(mat, probs = seq(0, 1, length.out = 13)))

# if we had to change the values in mat directly for the sake of this visualization,
#  we can do something like the following:
mat2 <- mat
break_vals <- quantile(mat, probs = seq(0, 1, length.out = 13))
for(i in 1:12){
 idx <- intersect(which(mat >= break_vals[i]), which(mat <= break_vals[i+1]))
 mat2[idx] <- i
}

image(clockwise90(mat2),  asp = T, 
      col = hcl.colors(12, "YlOrRd", rev = TRUE))

##########################################

set.seed(10)
# 3 cell-types, 10 cells per cell-type, and 5 genes total
dat_list <- lapply(1:3, function(x){
 set.seed(x)
 MASS::mvrnorm(10, rep(0,5), diag(5))
})

cor_list <- lapply(1:3, function(x){
 stats::cor(dat_list[[x]])
})

# this is analogous to what you've been currently doing:
plot(dat_list[[3]][,2], dat_list[[3]][,4], asp = T, pch = 16)

# this is more literal on what to do:
# you also include all the other cells
dat_all <- do.call(rbind, dat_list)
plot(NA, xlim = range(dat_all[,2]), ylim = range(dat_all[,4]))
for(i in 1:3){
 points(dat_list[[i]][,2], dat_list[[i]][,4], asp = T, pch = 16, col = i)
}

#########

# how do you find the pairs of genes that give the "largest contrast"
array_all <- array(NA, dim = c(5, 5, 3))
for(i in 1:3){
 array_all[,,i] <- cor_list[[i]]
}

# the following diff_mat quantities for each pair of genes, what the difference between the largest
## and smallest pearson corrleation is among the 3 cell-types
diff_mat <- apply(array_all, MARGIN = c(1,2), FUN = function(x){diff(range(x))})
