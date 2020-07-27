# This file is designed to draw a heatmap for each correlation matrix.

# Import Part #

## Import the stored RData
load(file = "Correlation_data.RData")

## Import ggplot library to draw heatmap.
library(ggplot2) # Visualization tool
library(reshape2) # Matrix melting tool

## Store the matrices in a list and make a vector to label
matrice <- list(cor_pearson_mat, cor_spearman_mat, faster_kendall_mat,
                hoeff_dist, dist_cor_mat)
label <- c("Pearson", "Spearman", "Kendall", "Hoeff", "Dist")

# Drawing heatmap Part #

for (i in 1:length(matrice)){
  cor_mat <- matrice[[i]]
  ## melt the matrix to draw a heatmap.
  melted_matrix <- melt(cor_mat, na.rm = T)
  
  ## ggplot to draw a heatmap with the viridis palette
  cor_heatmap <- ggplot(data = melted_matrix, 
                        aes(x = Var2, y = Var1, fill = value)) + 
    geom_tile() +
    scale_fill_viridis_b(option = "D", direction = -1) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = "transparent")) +
    labs(title = paste0("Heatmap of ", label[i], " Correlation")) 
  
  ## Store the heatmap as a file
  jpeg(filename = paste0("Heatmap_", label[i], "Corr.jpeg"),
       width = 1080,
       height = 1080)
  print(cor_heatmap)
  dev.off()
}

# Drawing heatmaps of absolute values Part #

for (i in 1:length(matrice)){
  cor_mat <- abs(matrice[[i]])
  ## melt the matrix to draw a heatmap.
  melted_matrix <- melt(cor_mat, na.rm = T)
  
  ## ggplot to draw a heatmap with the viridis palette
  cor_heatmap <- ggplot(data = melted_matrix, 
                        aes(x = Var2, y = Var1, fill = value)) + 
    geom_tile() +
    scale_fill_viridis_b(option = "D", direction = -1) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = "transparent")) +
    labs(title = paste0("Abs_Heatmap of ", label[i], " Correlation")) 
  
  ## Store the heatmap as a file
  jpeg(filename = paste0("Abs_Heatmap_", label[i], "Corr.jpeg"),
       width = 1080,
       height = 1080)
  print(cor_heatmap)
  dev.off()
}




