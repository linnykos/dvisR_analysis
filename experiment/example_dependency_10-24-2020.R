set.seed(10)
n <- 100
x1 <- stats::rnorm(n); y1 <- x1^2 + stats::rnorm(n, sd = 0.1)
x2 <- stats::rnorm(n); y2 <- x2^2 + stats::rnorm(n, sd = 0.1)
dat <- cbind(x1, y1, x2, y2)

par(mfrow = c(1,3))
plot(dat[,1], dat[,2], asp = T, pch = 16)
plot(dat[,3], dat[,4], asp = T, pch = 16)
plot(dat[,1], dat[,3], asp = T, pch = 16)

# category 1 of ways to induce dependency
# Method 1.A
## Instead of strictly ordering, we "randomize" the ordering slightly
# example of strict ordering
dat2 <- dat
dat2[,c(1,2)] <- dat2[order(dat2[,1]), c(1,2)]
dat2[,c(3,4)] <- dat2[order(dat2[,3]), c(3,4)]

par(mfrow = c(1,3))
plot(dat2[,1], dat2[,2], asp = T, pch = 16)
plot(dat2[,3], dat2[,4], asp = T, pch = 16)
plot(dat2[,1], dat2[,3], asp = T, pch = 16)

# instead, we could do a "randomized" ordering somewhat
# Method 1.A.I
shuffling_function <- function(vec, window = 5){
  order_vec <- order(vec)
  
  for(i in 1:(length(order_vec)-window)){
    order_vec[i:(i+window)] <- sample(order_vec[i:(i+window)])
  }
  
  order_vec
}

dat2 <- dat
new_order1 <- shuffling_function(dat2[,1], window = 10)
new_order3 <- shuffling_function(dat2[,3], window = 10)

dat2[,c(1,2)] <- dat[new_order1, c(1,2)]
dat2[,c(3,4)] <- dat[new_order3, c(3,4)]

par(mfrow = c(1,3))
plot(dat2[,1], dat2[,2], asp = T, pch = 16)
plot(dat2[,3], dat2[,4], asp = T, pch = 16)
plot(dat2[,1], dat2[,3], asp = T, pch = 16)

# Method 1.A.II
quadratic_ordering <- function(vec){
  order_vec <- order(vec)
  
  n <- length(order_vec)
  odd_idx <- seq(1, n, by = 2)
  even_idx <- seq(2, n, by = 2)
  
  new_order_vec <- c(order_vec[odd_idx], rev(order_vec[even_idx]))
  new_order_vec
}

dat2 <- dat
new_order3 <- quadratic_ordering(dat2[,3])

dat2[,c(1,2)] <- dat[order(dat[,1]), c(1,2)]
dat2[,c(3,4)] <- dat[new_order3, c(3,4)]

par(mfrow = c(1,3))
plot(dat2[,1], dat2[,2], asp = T, pch = 16)
plot(dat2[,3], dat2[,4], asp = T, pch = 16)
plot(dat2[,1], dat2[,3], asp = T, pch = 16)

# Method 1.B
## All we're going to do is introduce the notion of "cell-types"/"clusters"
set.seed(10)
cluster_label <- sample(1:6, 100, replace = T)

dat2 <- dat
for(k in 1:max(cluster_label)){
  idx <- which(cluster_label == k)
  
  dat2[idx,c(1,2)] <- dat2[idx[order(dat2[idx,1])], c(1,2)]
  dat2[idx,c(3,4)] <- dat2[idx[order(dat2[idx,3])], c(3,4)]
}

par(mfrow = c(1,3))
plot(dat2[,1], dat2[,2], asp = T, pch = 16)
plot(dat2[,3], dat2[,4], asp = T, pch = 16)
plot(dat2[,1], dat2[,3], asp = T, pch = 16, col = cluster_label)


###########################################
###########################################

rm(list=ls())
set.seed(10)
p <- 4
mean_vec <- rep(0, 4)
off_val <- 0.3
cov_mat <- matrix(c(1, 0.9, off_val, off_val,
                    0.9, 1, off_val, off_val, 
                    off_val, off_val, 1, 0.9,
                    off_val, off_val, 0.9, 1), 4, 4)
n <- 500
# (In the notation of https://jmlr.csail.mit.edu/papers/volume10/liu09a/liu09a.pdf, section 3)
## dat represent the n-samples of "Z"
dat <- MASS::mvrnorm(n = n, mu = mean_vec, Sigma = cov_mat)

par(mfrow = c(1,3))
plot(dat[,1], dat[,2], asp = T, pch = 16)
plot(dat[,3], dat[,4], asp = T, pch = 16)
plot(dat[,1], dat[,3], asp = T, pch = 16)

# (In the notation of https://jmlr.csail.mit.edu/papers/volume10/liu09a/liu09a.pdf, section 3)
## dat2 represent the n-samples of "X"
dat2 <- dat
dat2[,1] <- sign(dat2[,1])*abs(dat2[,1])^0.9
dat2[,2] <- sign(dat2[,2])*abs(dat2[,2])^0.8
dat2[,3] <- sign(dat2[,3])*abs(dat2[,3])^0.5
dat2[,4] <- sign(dat2[,4])*abs(dat2[,4])^0.4

par(mfrow = c(1,3))
plot(dat2[,1], dat2[,2], pch = 16)
plot(dat2[,3], dat2[,4], pch = 16)
plot(dat2[,1], dat2[,3], pch = 16)


