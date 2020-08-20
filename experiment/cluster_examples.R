# 1) smallest p-value in two-sample testing
# 2) clustering objective compared to random assignment
# 3) classification error

dat <- read.csv("../../data/Zeisel_preprocessed.csv", row.names = 1)
cell_type <- read.table("../../data/Zeisel_cell_info.txt", sep = "\t", header = 1)

# grab an example scatterplot
idx1 <- grep("Sema5a", colnames(dat))
idx2 <- grep("Lamp5", colnames(dat))
cell_type <- as.numeric(as.factor(cell_type$level1class))
plot(dat[,idx1], dat[,idx2], 
     col = cell_type, asp = T, pch =16)
plot(1:7, col = 1:7, pch = 16, cex = 5)

tmp <- which(cell_type %in% c(3,5))
cell_type2 <- cell_type[tmp]
dat2 <- dat[tmp,c(idx1, idx2)]
plot(dat2[,1], dat2[,2], asp = T, col = cell_type2, pch = 16)

# now that we have our 2 classes, we can try our different methods

## method 1 : smallest p-value in two-sample testing
## you could: use any type of two-sample hypothesis test here
## the simplest one could be a two-sample t-test
x <- dat2[which(cell_type2 == 3),]
y <- dat2[which(cell_type2 == 5),]
ttest_res <- stats::t.test(x, y, alternative = "two.sided")
ttest_res
ttest_res$p.value
# if the p-values is small, then the clusters are "different"
# of course, we want to pick "reasonably good tests" to use, as the t-test here
#  is only testing the mean

## method 2 : classification error
## we can imagine using a logistic regression
## you could: use any type of classification algorithm here
dat_df <- data.frame(x1 = dat2[,1], x2 = dat2[,2], label = as.factor(cell_type2))
log_res <- stats::glm(label ~ x1 + x2 + 1, data = dat_df, family = "binomial")
predictions <- (stats::predict(log_res, newdata = dat_df, type = "response") > 0.5)
table(cell_type2, predictions) # can compute the classification error from here

## method 3 : clustering objective compared to random assignment
## what we're going to do is use the k-means objective function
## you could: use any metric that is "small" when the cluster labels "conforms"/"agrees" with the data
.l2norm <- function(x){sqrt(sum(x^2))}

compute_ss <- function(dat, cluster_label){
 stopifnot(nrow(dat) == length(cluster_label))
 
 uniq_val <- unique(cluster_label)
 
 mean_list <- lapply(uniq_val, function(x){
  idx <- which(cluster_label == x)
  mean_vec <- colMeans(dat[idx,])
 })
 
 sum(sapply(1:nrow(dat), function(i){
  .l2norm(dat[i,] - mean_list[[which(uniq_val == cluster_label[i])]])
 }))
}

compute_ss(dat2, cell_type2)
set.seed(10)
new_label <- sample(1:2, size = nrow(dat2), replace = T)
compute_ss(dat2, new_label)

par(mfrow = c(1,2))
plot(dat2[,1], dat2[,2], asp = T, col = as.numeric(as.factor(cell_type2)), pch = 16)
plot(dat2[,1], dat2[,2], asp = T, col = new_label, pch = 16)

