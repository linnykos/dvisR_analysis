rm(list=ls())

generate_clusters_L <- function(n){
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

variable_transformation <- function(vec1, vec2, mean_val = 0, sd_val = 1){
  val <- stats::pnorm(vec1, mean = mean_val, sd = sd_val)
  
  stats::quantile(vec2, probs = val)
}

##################################

set.seed(10)
n <- 1000
dat_test1 <- generate_clusters_L(n)
dat_test2 <- generate_clusters_L(n)

dat1 <- cbind(dat_test1, dat_test2)
den_list <- lapply(1:4, function(x){stats::density(dat1[,x])})

mean_vec <- colMeans(dat1); cov_mat <- matrix(0.8, 4, 4); diag(cov_mat) <- 1
dat2 <- MASS::mvrnorm(n, mean_vec, cov_mat)

# nonparanormal transform
dat2_new <- do.call(cbind, lapply(1:4, function(x){
  variable_transformation(dat2[,x], dat1[,x], mean_val = mean_vec[x],
                          sd_val = sqrt(diag(cov_mat)[x]))
}))

# this is the proxy dataset we make to grab monotone functions from
par(mfrow = c(1,3))
plot(dat1[,1], dat1[,2], pch = 16, main = "Original data", asp = T, col = rgb(0.5,0.5,0.5,0.5))
plot(dat1[,3], dat1[,4], pch = 16, asp = T, col = rgb(0.5,0.5,0.5,0.5))
plot(dat1[,1], dat1[,3], pch = 16, asp = T, col = rgb(0.5,0.5,0.5,0.5))

##############

# this is the gaussian data we need to make nonparanormals
par(mfrow = c(1,3))
plot(dat2[,1], dat2[,2], pch = 16, main = "Gaussian data", asp = T, col = rgb(0.5,0.5,0.5,0.5))
plot(dat2[,3], dat2[,4], pch = 16, asp = T, col = rgb(0.5,0.5,0.5,0.5))
plot(dat2[,1], dat2[,3], pch = 16, asp = T, col = rgb(0.5,0.5,0.5,0.5))

##############

# this is the nonparanormal 
par(mfrow = c(1,3))
plot(dat2_new[,1], dat2_new[,2], pch = 16, main = "Final data", asp = T, col = rgb(0.5,0.5,0.5,0.5))
plot(dat2_new[,3], dat2_new[,4], pch = 16, asp = T, col = rgb(0.5,0.5,0.5,0.5))
plot(dat2_new[,1], dat2_new[,3], pch = 16, asp = T, col = rgb(0.5,0.5,0.5,0.5))


########

# this is the property we ensure by nonparanormals -- the marginal densities 
## between the nonparanormal and the proxy dataset are the same
for(i in 1:4){
  par(mfrow = c(1,2))
  hist(dat1[,i], col = "gray", main = i, breaks = 50)
  hist(dat2_new[,i], col = "gray", main = i, breaks = 50)
}


