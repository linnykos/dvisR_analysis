
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

correlation_pearson <- function(dat){
  stats::cor(dat[,1], dat[,2], method = "pearson")
}

correlation_spearman <- function(dat){
  stats::cor(dat[,1], dat[,2], method = "spearman")
}

correlation_kendall <- function(dat){
  pcaPP::cor.fk(dat)[1,2]
}

correlation_xi <- function(dat){
  XICOR::calculateXI(dat[,1], dat[,2])
}

correlation_beta <- function(dat){
  median_x <- stats::median(dat[,1])
  median_y <- stats::median(dat[,2])
  tab <- table(dat[,1] <= median_x, dat[,2] <= median_y)
  ((tab[1,1]+tab[2,2]) - (tab[1,2]+tab[2,1]))/sum(tab)
}

correlation_energy <- function(dat){
  as.numeric(energy::dcor2d(dat[,1], dat[,2]))
}