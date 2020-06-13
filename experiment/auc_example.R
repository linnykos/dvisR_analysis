rm(list=ls())

# setup
trials <- 1000

## generate the data
## note: you can put whatever type of way you generate data here
generate_data <- function(seed = 0, n = 100, sd = 2, independent = T){
 set.seed(seed)
 x <- stats::runif(n)
 if(independent){
  y <- stats::rnorm(n, sd = sd)
 } else {
  y <- x + stats::rnorm(n, sd = sd)
 }

 list(x = x, y = y)
}

## compute the pearson correlation's significance (i.e., returns a p-value)
## note: this function doesn't necessary have to return a p-value. it can be the magnitude of the pearson
##  correlation, for example. As long as the function returns something between [0,1], it will be easy
##  to compute the ROC (and AUC) for this statistic
compute_pearson_significance <- function(x, y){
 res <- stats::cor.test(x, y, method = "pearson")
 res$p.value
}

##########################

# generate a lot data
## note: for half the datasets, we make the two variable independent. 
##   for the other half, we make the two variables dependent
## here, 0 = dependent. 1 = independent
label_vec <- c(rep(0, round(trials/2)), rep(1, trials - round(trials/2)))
dat_list <- lapply(1:trials, function(trial){
 if(label_vec[trial] == 1){
  generate_data(seed = trial, n = 100, sd = 2, independent = T)
 } else {
  generate_data(seed = trial, n = 100, sd = 2, independent = F)
 }
})

## visualize one of the datasets just so we know what it looks like
## add red lines for reference
graphics::plot(NA, asp = T, xlab = "x", ylab = "y", xlim = range(dat_list[[1]]$x),
               ylim = range(dat_list[[1]]$y))
graphics::lines(rep(0,2), c(-1e4, 1e4), col = "red", lty = 2) 
graphics::lines(c(-1e4, 1e4), rep(0,2), col = "red", lty = 2) 
graphics::points(dat_list[[1]]$x, dat_list[[1]]$y, pch = 16)

## compute the pearson significance for each dataset in dat_list
pval_vec <- sapply(dat_list, function(dat){
 compute_pearson_significance(dat$x, dat$y)
})

## visualize the p-values
graphics::par(mfrow = c(1,2))
## first, the histogram
graphics::hist(pval_vec, col = "gray", breaks = 20, xlab = "p-value", main = "Histogram of p-values")
## next, the QQ-plot
graphics::plot(NA, asp = T,
               xlim = c(0,1), ylim = c(0,1), xlab = "Observed quantiles", ylab = "Theoretical quantiles",
               main = "QQ plot of Pearson p-values")
graphics::lines(c(0,1), c(0,1), col = "red", lty = 2)
graphics::points(sort(pval_vec), seq(0, 1, length.out = length(pval_vec)), pch = 16)

#########################

# How to understand ROC and AUC

## first make the ROC curve. Note there are MANY ways to make an ROC curve. In general, you just want
##   to ask -- in my setting, what makes sense to put on x-axis and y-axis such that when I plot
##   my results (one point per dataset), if my method had NO power (i.e., the same as randomly guessing
##   whether or not this dataset has x and y that are dependent), then the method would have points
##   that lie exactly along the diagonal.
## In our case (and often for our simulations), we want to perform the following thought-experiment:
##   If I were to set a threshold (say, threshold = 0.3) and call any p-value less than threshold as 
##   "the dataset corresponding to this p-value displays dependent data" and any p-value larger than the
##   threshold as "the dataset corresponding to this p-value displays independent data", how often
##   will I be correct or incorrect? This is where we calculate the true positive rate and false positive
##   rate (TPR and FPR). Then, if I varied this threshold from 0 to 1, my values of TPR and FPR will
##   "trace out a curve." This curve is our ROC curve.
## This variability in how to design an ROC curve is partially why the Wiki page of ROC curve is so
##   complicated. https://en.wikipedia.org/wiki/Receiver_operating_characteristic. Recall, if we were doing
## classification (ex: labeling whether or not a patient has cancer), this setup with be slightly different.

## making the ROC curve (again, in our case, it's exactly the same as the QQ-plot)
## first, compute the TPR and FPR across many values of the threshold
## here, we will define a "true" condition as "dataset is independent", 
##   and the "positive" condition as "data is predicted to be independent"
threshold_vec <- seq(0, 1, length.out = 500) # the length of this vector doesn't have to be equal to trials
tpr_vec <- sapply(threshold_vec, function(threshold){
 idx <- which(pval_vec >= threshold)
 length(intersect(idx, which(label_vec == 1)))/length(which(label_vec == 1))
})
fpr_vec <- sapply(threshold_vec, function(threshold){
 idx <- which(pval_vec >= threshold)
 length(intersect(idx, which(label_vec == 0)))/length(which(label_vec == 0))
})

graphics::plot(NA, asp = T,
               xlim = c(0,1), ylim = c(0,1), xlab = "FPR", ylab = "TPR ",
               main = "ROC curve for Pearson test of significance")
graphics::lines(c(0,1), c(0,1), col = "red", lty = 2)
graphics::points(fpr_vec, tpr_vec,  pch = 16)


## the AUC (area-under-the-curve) then is simply "what is the area between the curve and the diagonal line?"
## there are many ways to compute this. Fancier methods use a "trapezoid-integration method" to approximately
##   integrate this area. 
## Coding this integration is a pain, so we'll just the AUC package to do this for us

library(AUC)
roc_res <- AUC::roc(pval_vec, labels = as.factor(label_vec))
## we can plot the roc_res to get a similar to plot to above. it's slightly different since AUC::roc uses
##   a slightly different set of thresholds (see roc_res$cutoffs)  
graphics::plot(roc_res, asp = T)

AUC::auc(roc_res) # our AUC

#####################

# final notes:
## so clearly, when you play around with this, there are a lot of things you can vary.
## First: what's the "signal-to-noise" (SNR) ratio? For example, if you add the SNR
##   (decrease the value of sd or equivalently, increase n), the setting becomes easier
## Second: there are different ways to generate dependent data and independent data. Within this AUC
##   example, we need to generate many datasets (given by "trials"), some of which are dependent, and some
##   of which are independent. Changing literally how these are generated can change the AUC
## Regardless, the AUC is a good benchmarking way so you can compare across different methods. Here,
##   we only do Pearson's significance as the method. But it could be using the magnitude of Pearson correlation
##   (i.e., no hypothesis testing, so no p-value) or Kendall's tau p-value, etc. The main point is --
##   AUC measures the aggregate performance of a method by seeing how it performs across many datasets
##   and comparing its performance against a method that randomly guesses.

