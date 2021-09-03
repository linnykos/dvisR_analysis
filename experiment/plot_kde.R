#' Plotting density
#'
#' @param mat matrix
#' @param grid_spacing numeric vector
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
plot_density <- function(mat, grid_spacing,
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
  f1 <- MASS::kde2d(mat2[,1], mat2[,2], n = 100)
  graphics::image(f1,
                  col = density_colors,
                  xlim = xlim,
                  ylim = ylim,
                  ...)

  for(i in grid_spacing){
    graphics::lines(rep(i,2), c(-2,2)*max(abs(grid_spacing)),
                    lty = grid_lty,
                    lwd = grid_lwd,
                    col = grid_color)
    graphics::lines(c(-2,2)*max(abs(grid_spacing)), rep(i,2),
                    lty = grid_lty,
                    lwd = grid_lwd,
                    col = grid_color)
  }

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

##############################

# set.seed(10)
# mat <- do.call(rbind, lapply(c(0, 4), function(x){
#   MASS::mvrnorm(n = 100, mu = c(x,0), Sigma = diag(2))
# }))
# plot_density(mat, grid_spacing = seq(-3,7,by=2),
#              xlim = c(-3,7), ylim = c(-3,7),
#              xlab = "Example x-label",
#              ylab = "Custom y-label",
#              main = "Title here")
