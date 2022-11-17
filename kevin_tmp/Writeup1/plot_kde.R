#' Plotting density
#'
#' @param mat matrix
#' @param n integer for \code{stats::kde}
#' @param max_num integer
#' @param density_colors color vector
#' @param point_color color
#' @param contour_color color
#' @param grid_color color
#' @param pch integer
#' @param contour_lwd numeric
#' @param grid_lty integer
#' @param grid_lwd numeric
#' @param ... additional graphical parameters
#'
#' @return nothing
#' @export
plot_density <- function(mat,
                         xlim,
                         ylim,
                         n = 100,
                         max_num = 20000,
                         density_colors = grDevices::heat.colors(100, alpha = 0.8),
                         point_color = grDevices::rgb(0.5, 0.5, 0.5, 0.8),
                         contour_color = grDevices::rgb(1,1,1,0.5),
                         grid_color = "gray",
                         xaxt = "n",
                         yaxt = "n",
                         axes = F,
                         xlab = "",
                         ylab = "",
                         pch = 16,
                         contour_lwd = 2,
                         grid_lty = 2,
                         grid_lwd = 0.5,
                         ...){
  stopifnot(ncol(mat) == 2)
  
  mat2 <- rbind(mat,
                c(xlim[1], ylim[1]),
                c(xlim[1], ylim[2]),
                c(xlim[2], ylim[1]),
                c(xlim[2], ylim[2]))
  f1 <- MASS::kde2d(mat2[,1], mat2[,2], n = n)
  graphics::image(f1,
                  col = density_colors,
                  xlim = xlim,
                  ylim = ylim,
                  xlab = xlab,
                  ylab = ylab,
                  ...)
  
  if(nrow(mat) > max_num) {
    idx <- sample(1:nrow(mat), size = max_num)
  } else {
    idx <- 1:nrow(mat)
  }
  graphics::points(mat[idx,1], mat[idx,2],
                   col = point_color,
                   pch = pch)
  
  graphics::contour(f1, add = T,
                    drawlabels = F,
                    col = contour_color,
                    lwd = contour_lwd)
  
  invisible()
}