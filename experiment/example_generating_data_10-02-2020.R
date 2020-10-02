rm(list=ls())

.generate_lollipop <- function(n){
 t(sapply(1:n, function(i){
  # generate "cluster" assignment
  k <- sample(c(1,2), 1)
  
  # generate "ball"
  if(k == 1){
   x <- stats::rnorm(1, mean = 0, sd = 1)
   y <- stats::rnorm(1, mean = 0, sd = 1)
  } else {
   # k == 2, generate "stick"
   x <- stats::runif(1, min = -1, max = 5)
   y <- 1*x + stats::rnorm(1, mean = 0, sd = 0.5)
  }
  
  c(x,y)
 }))
}

.generate_v <- function(n){
 t(sapply(1:n, function(i){
  # generate "cluster" assignment
  k <- sample(c(1,2), 1, prob = c(0.7, 0.3))
  
  # generate "ball"
  if(k == 1){
   x <- stats::runif(1, min = 0, max = 5)
   y <- 2*x + stats::rnorm(1, mean = 0, sd = 0.5)
  } else {
   # k == 2, generate "stick"
   x <-  stats::runif(1, min = 0, max = 5)
   y <- 0.5*x + stats::rnorm(1, mean = 0, sd = 0.1)
  }
  
  c(x,y)
 }))
}

##########################

generate_data <- function(n, p){
 stopifnot(p %% 2 == 0)
 
 # generate p/2 2-dimensional datasets
 dat_list <- lapply(1:(p/2), function(round){
  type <- sample(c(1,2), 1)
  
  if(type == 1){
   .generate_lollipop(n)
  } else {
   .generate_v(n)
  }
 })
 
 # concatenate column-wise
 do.call(cbind, dat_list)
}

set.seed(10)
zz <- generate_data(500, 10)

# 5 out of the 45 pairs of dependency
# would work for c(1,2), c(3,4), c(5,6), c(7,8), c(9,10)
plot(zz[,1], zz[,2], asp = T) # of lollipop
plot(zz[,3], zz[,4], asp = T) # of v

# any other pair looks like independent variables
plot(zz[,1], zz[,3], asp = T)
