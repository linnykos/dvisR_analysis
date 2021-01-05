#############################################################################

# Import Ziesel dataset
dat <- read.csv("Zeisel_preprocessed.csv", row.names = 1)
cell_type <- read.table("Zeisel_cell_info.txt", sep = "\t", header = 1)

# Get the labels for each cell
cluster_labels <- as.numeric(as.factor(cell_type$level1class))

cell_labels <- unique(as.factor(cell_type$level1class))
subset_data <- list()
subset_heatmap <- list()
load("meaningful_indices.RData")

sub_full <- dat[, indices]

for (i in 1:length(cell_labels)){
  label <- cell_labels[[i]]
  extracted <- dat[which(cell_type$level1class == label), indices]
  subset_data[[i]] <- extracted
  
  hclustered <- hclust(dist(t(extracted)))$order
  reordered <- extracted[, hclustered]
  subset_heatmap[[i]] <- reordered
}

col <- ncol(subset_data[[1]])
cell_num <- length(cell_labels)

#############################################################################

### manually generate data with 7 different shapes

# 1. lollipop shape
.generate_lollipop <- function(n){
  t(sapply(1:n, function(i){
    # generate "cluster" assignment
    k <- sample(c(1,2), 1)
    
    # generate "ball"
    if(k == 1){
      x <- stats::rnorm(1, mean = 0, sd = 1)
      y <- stats::rnorm(1, mean = 0, sd = 1)
    } else {
      # k == 2, generate "stick"
      x <- stats::rnorm(1, mean = 3, sd = 0.8)
      y <- x + stats::rnorm(1, mean = 1, sd = 0.5)
    }
    
    c(x,y)
  }))
}

# 2. V shape
.generate_v <- function(n){
  t(sapply(1:n, function(i){
    # generate "cluster" assignment
    k <- sample(c(1,2,3), 1, prob = c(0.45, 0.45, 0.1))
    
    # generate "stick"
    if(k == 1){
      x <- stats::rnorm(1, mean = 2.5, sd = 1)
      y <- 2 * x + stats::rnorm(1, mean = 0, sd = 0.5)
    } else if (k == 2) {
      # k == 2, generate "stick"
      x <- stats::rnorm(1, mean = 2.5, sd = 1)
      y <- 0.5 * x + stats::rnorm(1, mean = 0, sd = 0.1)
    } else {
      x <- stats::rnorm(1, mean = 0.5, sd = 0.2)
      y <- stats::rnorm(1, mean = 0.5, sd = 0.3)
    }
    
    c(x,y)
  }))
}

# 3. Local correlated shape
.generate_local <- function(n){
  t(sapply(1:n, function(i){
    k <- sample(1:3, 1, prob = c(0.4, 0.4, 0.2))
    
    # generate "left_ball"
    if (k == 1){
      x <- stats::rnorm(1, mean = -3, sd = 1)
      y <- stats::rnorm(1, mean = 2, sd = 1)
    } else if (k == 2) {
      # generate "right_ball"
      x <- stats::rnorm(1, mean = 3, sd = 1)
      y <- stats::rnorm(1, mean = 5, sd = 1)
    } else{
      # generate "local correlated data"
      x <- stats::rnorm(1, mean = 0, sd = 0.5)
      y <- x + stats::rnorm(1, mean = 3.5, sd = 0.1) 
    }
    c(x,y)
  }))
}

# 4. Quadratic Shape
.generate_quadratic <- function(n){
  t(sapply(1:n, function(i){
    x <- stats::rnorm(1, mean = 0, sd = 1)
    y <- x^2 + stats::rnorm(1, mean = 0, sd = 0.5) 
    c(x,y)
  }))
}

# 5. L Shape
.generate_clusters_L <- function(n){
  t(sapply(1:n, function(i){
    k <- sample(1:4, 1, prob = rep(1/4, 4))
    
    if (k == 1){
      x <- stats::rnorm(1, mean = -1, sd = 0.1)
      y <- x + stats::rnorm(1, mean = 1, sd = 1)
    } else if (k == 2){
      x <- stats::rnorm(1, mean = 0.5, sd = 0.5)
      y <- stats::rnorm(1, mean = -2, sd = 0.1)
    } else if (k == 3){
      x <- stats::rnorm(1, mean = -1, sd = 0.1)
      y <- x + stats::rnorm(1, mean = 0.5, sd = 0.5)
    } else {
      x <- stats::rnorm(1, mean = -.5, sd = 0.2)
      y <- stats::rnorm(1, mean = -2, sd = 0.1)
    }
    
    c(x,y)
  }))
}

# 6. Three clusters with a line
.generate_clusters1 <- function(n){
  t(sapply(1:n, function(i){
    k <- sample(1:4, 1, prob = seq(0.25, 1, by = 0.25))
    
    if (k == 1){
      x <- stats::rnorm(1, mean = -1.5, sd = 0.5)
      y <- -x + stats::rnorm(1, mean = 0.2, sd = 0.2) 
    } else if (k == 2){
      x <- stats::rnorm(1, mean = -3, sd = 0.5)
      y <- stats::rnorm(1, mean = -2, sd = 0.5)
    } else if (k == 3){
      x <- stats::rnorm(1, mean = 0, sd = 0.5)
      y <- x + -1 * stats::rnorm(1, mean = 1, sd = 0.5)
    } else{
      x <- stats::rnorm(1, mean = 2, sd = 0.5)
      y <- stats::rnorm(1, mean = 0.5, sd = 0.5)
    }
    
    c(x,y)
  }))
}

# 7. Triangle Shape with a cluster
.generate_clusters2 <- function(n){
  t(sapply(1:n, function(i){
    k <- sample(1:5, 1, prob = rep(1/5, 5))
    if (k == 1) {
      x <- stats::rnorm(1, mean = 0.75, sd = 0.5)
      y <- -0.5 * x + stats::rnorm(1, mean = 1.25, sd = 0.1)
    } else if (k == 2) {
      x <- stats::rnorm(1, mean = -0.25, sd = 0.1)
      y <- stats::rnorm(1, mean = 0, sd = 0.1)
    } else if (k == 3) {
      x <- stats::rnorm(1, mean = 1, sd = 0.75)
      y <- stats::rnorm(1, mean = 0, sd = 0.1)
    } else if (k == 4) {
      x <- stats::rnorm(1, mean = -0.25, sd = 0.1)
      y <- x + stats::rnorm(1, mean = 1.25, sd = 0.25)
    } else {
      x <- stats::rnorm(1, mean = 0.25, sd = 0.25)
      y <- stats::rnorm(1, mean = 0.5, sd = 0.1)
    }
    c(x,y)
  }))
}

# 8. Linear Shape (pos, neg)

.generate_pos_linear <- function(n){
  t(sapply(1:n, function(i){
    # generate "cluster" assignment
    x <- stats::rnorm(1, mean = 2.5, sd = 1)
    y <- 2 * x + stats::rnorm(1, mean = 0, sd = 1)
    c(x,y)
  }))
}

.generate_neg_linear <- function(n){
  t(sapply(1:n, function(i){
    # generate "cluster" assignment
    x <- stats::rnorm(1, mean = 2.5, sd = 1)
    y <- - 2 * x + stats::rnorm(1, mean = 0, sd = 1)
    c(x,y)
  }))
}

# 9. Independent Shape 

.generate_independent <- function(n){
  t(sapply(1:n, function(i){
    # generate "cluster" assignment
    x <- stats::rnorm(1, mean = 1, sd = 0.1)
    y <- 2 * x + stats::rnorm(1, mean = 0, sd = 1)
    c(x,y)
  }))
}
