set.seed(10)
rm(list=ls())

## http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
## https://jschoeley.github.io/ced18-tidyr/visualization.html
## https://stackoverflow.com/questions/48424682/how-do-i-limit-the-range-of-the-viridis-colour-scale
## https://stackoverflow.com/questions/10981324/ggplot2-heatmap-with-colors-for-ranged-values

# ## playing around
# len <- 3
# col_vec <- scales::viridis_pal()(len)
# plot(1:len, col = col_vec, pch = 16, cex = 5)

###################

cor_mat <- stats::cor(MASS::mvrnorm(n = 4, mu = rep(0,10), Sigma = diag(10)))
cor_mat <- abs(cor_mat)
# cor_mat <- cor_mat^2+5
cor_mat[lower.tri(cor_mat, diag = T)] <- NA

melted_matrix <- reshape::melt(cor_mat, na.rm = T)
colnames(melted_matrix) = c("Var1", "Var2", "value")
melted_matrix <- melted_matrix[!is.na(melted_matrix$value),]

len <- 11
break_vec <- round(as.numeric(quantile(cor_mat, probs = seq(0, 1,length.out = len), na.rm = T)), 2)
break_vec[1] <- break_vec[1]-2
break_vec[length(break_vec)] <- break_vec[length(break_vec)]+2

melted_matrix$value <- cut(melted_matrix$value, breaks = break_vec)

ggplot2::ggplot(data = melted_matrix, aes(x = Var2, y = Var1, fill = value)) + 
 geom_tile() + 
 ggplot2::scale_fill_manual(breaks = sort(unique(melted_matrix$value)), 
                            values = rev(scales::viridis_pal(begin = 0, end = 1)(len-1)))
# 
#  ggplot2::scale_fill_viridis_b(
#   #rescaler = function(x, from = NULL){sapply(x,function(y){(which.min(abs(y - break_vec))-1)/(len-1)})},
#   breaks = break_vec,
#   # begin = 0, end = 1,
#   # values = seq(0, 1, length.out = len-1),
#   option = "D", direction = -1
#   )

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
