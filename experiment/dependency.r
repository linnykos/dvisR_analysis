install.packages("energy", repos = "http://cran.us.r-project.org")

library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)
library(SAVER)
library(reshape2) # melt function
library(ggplot2) # ggplot function
library(pcaPP) # Fast Kendall function
library(energy) # Distance Correlation
library(Hmisc) # Hoeffding's D measure
library(zebu) # Normalized Mutual Information
library(XICOR) # Chatterjee's Coefficient
library(VineCopula) # Blomqvist's Beta

# bm <- LoadData(ds = "bmcite")
# orig_count <- Matrix::t(bm@assays$RNA@counts)
# 
# load("full_saver.rds")
# saver_mat <- Matrix::t(saver_dat$estimate)
# 
# CD4_ind <- which(as.factor(bm@meta.data$celltype.l2) == "CD4 Naive")
# 
# sub_dat <- saver_mat[CD4_ind, ]
# 
# first_filter <- apply(sub_dat, 2, function(x) {sd(x)})
# first_filter_ind <- which(first_filter > 0)
# sub_dat <- sub_dat[, first_filter_ind]
# 
# print("filter is done")
# 
# dat_hclust <- hclust(dist(t(as.matrix(sub_dat))))
# dat_index <- dat_hclust$order
# 
# sub_dat <- sub_dat[, dat_index]
# save_date <- date()
# save(sub_dat, dat_index, save_date, file = "CD4_dat.rds")

load("CD4_dat.rds")
make_cormat <- function(dat_mat){
  matrix_dat <- matrix(nrow = ncol(dat_mat), ncol = ncol(dat_mat))
  cor_mat_list <- list()
  
  basic_cor <- c("pearson", "spearman")
  # find each of the correlation matrices with Pearson, Spearman, Kendall Correlation Coefficients
#  for (i in 1:2){
#    print(i)
#    cor_mat <- stats::cor(dat_mat, method = basic_cor[i])
#    cor_mat[upper.tri(cor_mat, diag = T)] <- NA
#    cor_mat_list[[i]] <- cor_mat
#    save(cor_mat_list, file = "CD4_temp_cor.rds")
#  }
  
  # functions that take matrix or data.frame as input
  no_loop_function <- c(pcaPP::cor.fk, Hmisc::hoeffd, 
                        VineCopula::BetaMatrix)
  for (i in 3:5){
    if (i == 3){
      next
    }
    print(i)
    fun <- no_loop_function[[i-2]]
    cor_mat <- fun(dat_mat)
    if (i == 4){ # Hoeffding's D
      cor_mat <- cor_mat$D
    }
    cor_mat[upper.tri(cor_mat, diag = T)] <- NA
    cor_mat_list[[i]] <- cor_mat
    save(cor_mat_list, file = "CD4_temp_cor1_2.rds")
  }

  return (cor_mat_list)
  # functions that take two variables as input to calculate correlations.
  need_loop <- c(zebu::lassie, energy::dcor2d, XICOR::calculateXI)
  
  for (i in 6:8){
    break
    print(i)
    fun <- need_loop[[i-5]]
    
    cor_mat <- matrix(nrow = ncol(dat_mat),
                      ncol = ncol(dat_mat))
    
    for (j in 2:ncol(dat_mat)){
      for (k in 1:(j-1)){
        if (i == 6){
          cor_mat[j, k] <- fun(cbind(dat_mat[, j], dat_mat[, k]), continuous=c(1,2), breaks = 6, measure = "npmi")$global
          
        } else {
          cor_mat[j, k] <- fun(as.numeric(dat_mat[, j]),
                               as.numeric(dat_mat[, k]))
        }
      }
    }
    
    cor_mat[upper.tri(cor_mat, diag = T)] <- NA
    cor_mat_list[[i]] <- cor_mat
  #  save(cor_mat_list, file = "CD4_temp_cor.rds")
  }
  return(cor_mat_list)
}

draw_heatmap <- function(cor_mat){
  len <- 5
  melted_cormat <- melt(cor_mat)
  melted_cormat <- melted_cormat[!is.na(melted_cormat$value),]
  break_vec <- round(as.numeric(quantile(melted_cormat$value,
                                         probs = seq(0, 1, length.out = len),
                                         na.rm = T)),
                     4)
  break_vec[1] <- break_vec[1]-1
  break_vec[len] <- break_vec[len]+1
  melted_cormat$value <- cut(melted_cormat$value, breaks = break_vec)
  heatmap_color <- unique(melted_cormat$value)
  
  heatmap <- ggplot(data = melted_cormat, aes(x = Var2, y = Var1, fill = value))+
    geom_tile(colour = "Black") +
    ggplot2::scale_fill_manual(breaks = sort(heatmap_color), 
                               values = rev(scales::viridis_pal(begin = 0, end = 1)
                                            (length(heatmap_color)))) +
    theme_bw() + # make the background white
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.ticks = element_blank(),
          # erase tick marks and labels
          axis.text.x = element_blank(), axis.text.y = element_blank())
  
  return (heatmap)
}

make_cor_heatmap <- function(dat_mat, cor_mat_list){
  fun_lable <- c("Pearson's Correlation", "Spearman's Correlation", "Kendall's Correlation",
                 "Hoeffding's D", "Blomqvist's Beta", "NMI",
                 "Distance Correlation", "XI Correlation")
  
  cor_heatmap_list <- list()
  cor_abs_heatmap_list <- list()
  
  # make correlation matrices
  #cor_mat_list <- make_cormat(dat_mat)
  
  for (i in 1:8){
    print(i)
    if(i == 5) {
      cor_heatmap_list[[i]] <- NULL
      cor_abs_heatmap_list[[i]] <- NULL
      next 
    }
    cor_mat <- abs(cor_mat_list[[i]])
    
    # get heatmaps
    cor_heatmap <- draw_heatmap(cor_mat)
    
    # ggplot labels
    ggplot_labs <- labs(title = paste("Heatmap of", fun_lable[i]),
                        x = "",
                        y = "",
                        fill = "Coefficient") # change the title and legend label
    
    cor_heatmap_list[[i]] <- cor_heatmap + ggplot_labs
    
    if (i %in% c(1,2,3,4,6)){
      cor_abs_mat <- abs(cor_mat_list[[i]])
      cor_abs_heatmap <- draw_heatmap(cor_abs_mat)
      ggplot_abs_labs <- labs(title = paste("Abs Heatmap of", fun_lable[i]),
                              x = "", # change the title and legend label
                              y = "", 
                              fill = "Coefficient") 
      cor_abs_heatmap_list[[i]] <- cor_abs_heatmap + ggplot_abs_labs
    } else {
      ggplot_abs_labs <- labs(title = paste("Abs Heatmap of", fun_lable[i]),
                              subtitle = "Equivalent to Non-Abs Heatmap",
                              x = "", # change the title and legend label
                              y = "", 
                              fill = "Coefficient") 
      cor_abs_heatmap_list[[i]] <- cor_heatmap + ggplot_abs_labs
    }
  }
  
  ans <- list(cor_heatmap_list, cor_abs_heatmap_list)
  
  return (ans)
}

cormat_list <- make_cormat(sub_dat)
save(cormat_list, sub_dat, dat_index, store_date,
     file = "CD4_seurat_corr2_3.RData")
