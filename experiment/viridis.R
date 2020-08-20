set.seed(10)
rm(list=ls())

## http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
## https://jschoeley.github.io/ced18-tidyr/visualization.html
## https://stackoverflow.com/questions/48424682/how-do-i-limit-the-range-of-the-viridis-colour-scale
## https://stackoverflow.com/questions/10981324/ggplot2-heatmap-with-colors-for-ranged-values

# step 1: generate correlation matrix
cor_mat <- stats::cor(MASS::mvrnorm(n = 4, mu = rep(0,10), Sigma = diag(10)))
cor_mat <- abs(cor_mat)
# cor_mat <- cor_mat^2+5 ## <- [uncomment this line to convince yourself this code does what we think it does]
cor_mat[lower.tri(cor_mat, diag = T)] <- NA

# step 2: melt the matrix
melted_matrix <- reshape::melt(cor_mat, na.rm = T)
colnames(melted_matrix) = c("Var1", "Var2", "value")
melted_matrix <- melted_matrix[!is.na(melted_matrix$value),]

# step 3: generate cutoffs
len <- 3 ## <- [you can change this when you use this yourself, it's set to 3 just for ease of demonstration]
break_vec <- round(as.numeric(quantile(cor_mat, probs = seq(0, 1,length.out = len), na.rm = T)), 2)
break_vec[1] <- break_vec[1]-2
break_vec[length(break_vec)] <- break_vec[length(break_vec)]+2

# step 4: cut the values off into discrete bins
melted_matrix$value <- cut(melted_matrix$value, breaks = break_vec)

# step 5: using scale_fill_manual, construct the heatmap
ggplot2::ggplot(data = melted_matrix, aes(x = Var2, y = Var1, fill = value)) + 
 ggplot2::geom_tile() + 
 ggplot2::scale_fill_manual(breaks = sort(unique(melted_matrix$value)), 
                            values = rev(scales::viridis_pal(begin = 0, end = 1)(len-1)))
