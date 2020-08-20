set.seed(10)
rm(list=ls())

## http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
## https://jschoeley.github.io/ced18-tidyr/visualization.html
## https://stackoverflow.com/questions/48424682/how-do-i-limit-the-range-of-the-viridis-colour-scale

cor_mat <- stats::cor(MASS::mvrnorm(n = 4, mu = rep(0,10), Sigma = diag(10)))
cor_mat <- abs(cor_mat)
cor_mat <- cor_mat^2+5

melted_matrix <- reshape::melt(cor_mat, na.rm = T)
colnames(melted_matrix) = c("Var1", "Var2", "value")

len <- 3
break_vec <- as.numeric(quantile(melted_matrix[,3], probs = seq(0, 1,length.out = len)))
break_vec[1] <- break_vec[1]-2
break_vec[length(break_vec)] <- break_vec[length(break_vec)]+2

ggplot2::ggplot(data = melted_matrix, aes(x = Var2, y = Var1, fill = value)) + 
 geom_tile() + 
 # scale_fill_gradientn(colors = scales::viridis_pal()(len), breaks = break_vec,
 #                      labels = round(break_vec, 2))
 ggplot2::scale_fill_viridis_b(
  #rescaler = function(x, from = NULL){sapply(x,function(y){(which.min(abs(y - break_vec))-1)/(len-1)})},
  breaks = break_vec,
  # begin = 0, end = 1,
  # values = seq(0, 1, length.out = len-1),
  option = "D", direction = -1
  )

###################################
# just playing around below

# library(ggplot2)
# dsamp <- diamonds[sample(nrow(diamonds), 1000), ]
# ggplot(dsamp, aes(carat, price)) +
#  geom_point(aes(colour = clarity))
# 
# x <- LETTERS[1:20]
# y <- paste0("var", seq(1,20))
# data <- expand.grid(X=x, Y=y)
# data$Z <- runif(400, 0, 5)
# 
# # Heatmap 
# ggplot(data, aes(X, Y, fill= Z)) + 
#  geom_tile()


# mydata <- mtcars[, c(1,3,4,5,6,7)]
# cormat <- round(cor(mydata),2)
# library(reshape2)
# melted_cormat <- melt(cormat)
# library(ggplot2)
# ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
#  geom_tile()
