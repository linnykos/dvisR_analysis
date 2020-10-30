rm(list=ls())
set.seed(10)
p <- 4
mean_vec <- rep(0, 4)
off_val <- 0
cov_mat <- matrix(c(1, 0.9, off_val, off_val,
                    0.9, 1, off_val, off_val, 
                    off_val, off_val, 1, 0.9,
                    off_val, off_val, 0.9, 1), 4, 4)
n <- 500
# (In the notation of https://jmlr.csail.mit.edu/papers/volume10/liu09a/liu09a.pdf, section 3)
## dat represent the n-samples of "Z"
dat <- MASS::mvrnorm(n = n, mu = mean_vec, Sigma = cov_mat)

dat2 <- dat
dat2[,1] <- sign(dat2[,1])*abs(dat2[,1])^0.9 # these are very specific functions, can we generate these functions more generally?
dat2[,2] <- sign(dat2[,2])*abs(dat2[,2])^0.8
dat2[,3] <- sign(dat2[,3])*abs(dat2[,3])^0.5
dat2[,4] <- sign(dat2[,4])*abs(dat2[,4])^0.4

par(mfrow = c(1,3))
plot(dat2[,1], dat2[,2], pch = 16)
plot(dat2[,3], dat2[,4], pch = 16)
plot(dat2[,1], dat2[,3], pch = 16)

####################################

# the overall strategy is:
## we need a monotonic transformation, so the most flexible monotonic function we can think of is: the CDF
## so: when we need a general monotonic function -- let's just estimate the CDF of some other dataset's variable
##   and use that CDF as the CDF we use for our data-generating process

.generate_clusters_L <- function(n){
  t(sapply(1:n, function(i){
    k <- sample(1:4, 1, prob = rep(1/4, 4))
    
    if (k == 1){
      x <- stats::rnorm(1, mean = -1, sd = 0.1)
      y <- x + stats::rnorm(1, mean = 1, sd = 1)
    } else if (k == 2){
      x <- stats::rnorm(1, mean = 0.5, sd = 0.5)
      y <- stats::rnorm(1, mean = -2, sd = 0.1)
    } else if (k == 3){
      x <- stats::rnorm(1, mean = -1, sd = 0.1)
      y <- x + stats::rnorm(1, mean = 0.5, sd = 0.5)
    } else {
      x <- stats::rnorm(1, mean = -.5, sd = 0.2)
      y <- stats::rnorm(1, mean = -2, sd = 0.1)
    }
    
    c(x,y)
  }))
}

.generate_clusters2 <- function(n){
  t(sapply(1:n, function(i){
    k <- sample(1:5, 1, prob = rep(1/5, 5))
    if (k == 1) {
      x <- stats::runif(1, min = -0.5, max = 2)
      y <- -0.5 * x + stats::rnorm(1, mean = 1.25, sd = 0.1)
    } else if (k == 2) {
      x <- stats::rnorm(1, mean = -0.25, sd = 0.1)
      y <- stats::rnorm(1, mean = 0, sd = 0.1)
    } else if (k == 3) {
      x <- stats::rnorm(1, mean = 1, sd = 0.75)
      y <- stats::rnorm(1, mean = 0, sd = 0.1)
    } else if (k == 4) {
      x <- stats::runif(1, min = -0.5, max = 0.25)
      y <- x + stats::rnorm(1, mean = 1.25, sd = 0.25)
    } else {
      x <- stats::rnorm(1, mean = 0.25, sd = 0.25)
      y <- stats::rnorm(1, mean = 0.5, sd = 0.1)
    }
    c(x,y)
  }))
}

set.seed(10)
dat_test1 <- .generate_clusters_L(n)
dat_test2 <- .generate_clusters2(n)

# now that we have some datasets 
tmp_dataset <- cbind(dat_test1, dat_test2)

dat3 <- dat
for(j in 1:4){
  tmp <- dat3[,j]
  tmp2 <- rep(NA, n)
  
  # now we do the matching
  for(i in 1:n){
    quant_val <- length(which(tmp <= tmp[i]))/n
    map_to_value <- stats::quantile(tmp_dataset[,j], probs = quant_val)
    tmp2[i] <- map_to_value
  }
  
  dat3[,j] <- tmp2
}

# now-generate nonparanormal (i.e., the copula) dataset
par(mfrow = c(1,3))
plot(dat3[,1], dat3[,2], pch = 16)
plot(dat3[,3], dat3[,4], pch = 16)
plot(dat3[,1], dat3[,3], pch = 16)

# for sake comparison, let's look at tmp_dataset
par(mfrow = c(1,3))
plot(tmp_dataset[,1], tmp_dataset[,2], pch = 16)
plot(tmp_dataset[,3], tmp_dataset[,4], pch = 16)
plot(tmp_dataset[,1], tmp_dataset[,3], pch = 16)


