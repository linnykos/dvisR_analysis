rm(list=ls())

vine_func <- function(x, y, mode = "pkg"){
  if(mode == "pkg"){
    VineCopula::BetaMatrix(cbind(x,y))[1,2]
  } else {
    # following Blomqvist, N. (1950).  On a measure of dependence between two
    # random variables. The Annals of Mathematical Statistics, 21(4), 593-600.
    mx <- median(x); my <- median(y)
    n1 <- sum(x < mx & y > my) +  sum(x > mx & y < my)
    n2 <- sum(x < mx & y < my) +  sum(x > mx & y > my)
    (n1 - n2)/(n1 + n2)
  }
 
}

### 
n <- 50
set.seed(10)
xtmp <- runif(n); ytmp <- runif(n)
dat <- cbind(c(xtmp, xtmp, xtmp+2, xtmp+2), c(ytmp, ytmp+2, ytmp, ytmp+2))
plot(dat[,1], dat[,2], pch = 16, asp = T)

set.seed(10)
vine_func(dat[,1], dat[,2], mode = "pkg")
vine_func(dat[,1], dat[,2], mode = "custom")


vine_func(dat[,1]+5, dat[,2]+5, mode = "pkg")
vine_func(dat[,1]/3, dat[,2]/3, mode = "pkg")

