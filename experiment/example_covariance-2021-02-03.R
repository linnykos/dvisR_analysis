### FOR THE COVARIANCE MATRIX
# we can think of our covariance matrix as 
# a 2x2 block matrix. So let \Sigma \in \mathbb{R}^{9x9}
# matrix. So we can pictorally describe this as:
# [A   , B ]
# [B^T , C ], where A is the 4x4 matrix from dataset 1, 
# C is the 5x5 matris from dataset 2, and B is the new matrix
# we need to construct (but we can set it to be all 0's as a default)

# Of course, B (a 4x5 submatrix) doesn't have to be 0 though. So there are various
# things we can play with: 
# - We can play with B, A or C 
# - Most likely, we'll like to change only B (because we "like"
#    what A and C currently look like)
# - Some tricky things to remind ourselves though: Changing B
#    might result in Sigma no longer being positive semidefinite (PSD)
#    (So realistically, we just need to be mindful of this, and guess-and-check)
# - However, it's a little easier to ensure Sigma is still PSD
#    when we change A or C

# There are general strategies on how to design a covariance matrix
# https://online.stat.psu.edu/stat502/lesson/10/10.3
# https://support.sas.com/resources/papers/proceedings/proceedings/sugi30/198-30.pdf
# https://www.spratings.com/documents/20184/86990/BlockFinalsp/876fb460-00d2-4baa-8136-3403d84c8ee3
# https://www.stat.berkeley.edu/~bickel/BL2008-banding.pdf

# by and large, there are a few ways to think about constructing covariance matrices
# - 1) diagonal matrix
# - 2) autoregressive (AR) matrix (this is a special case of Toeplitz matrices): matrix with "decaying off-diagonal terms"
p <- 4; rho <- .5
mat <- matrix(0, p, p)
for(i in 1:p){
  for(j in 1:p){
    mat[i,j] <- rho^(abs(i-j)) #note that we can change the power function slight (this gives us the spatial power structure)
  }
}
mat
# however, in this case, all the diagonal entries are exactly 1. So what do we do?
# Answer: We can multiply this matrix on the left and right by a diagonal matrices
diag_vec <- c(5, 10, 5, 10)
mat2 <- diag(sqrt(diag_vec)) %*% mat %*% diag(sqrt(diag_vec))
round(mat2, 2)
eigen(mat2)$values # this transform perserved PSD
# So if you wanted to use this as the B matrix, one idea
# is to acknowledge that the B matrix is an off-diagonal block
# so perhaps we want the largest values of B to be closest to the
# lower-left corner
#  - 3) Banded covariance matrix (all this means is that it's a "tri-diagonal matrix")
#    (in general, very similar to the AR matrix in the sense that you have different "bands" of off-diagonal terms)
#  - 4) Block covariance matrix (the only difficult part about this
#   is that you need to chose the values of the blocks appropriately to ensure the matrix is still PSD)

# concluding thoughts: 
## - you might need to change A and C slightly 
#      (by making terms originally 0 something slightly not 0)
#      (The reason: The more non-zero A and C are, the more freedom you have in setting B for Sigma to still be PSD)
## - you can then have more freedom in choosing B. I think setting B to be a block matrix would be 
##   reasonable [so for example: The B matrix could look as follows]
mat <- matrix(0, nrow = 4, ncol = 5)
mat[1:2,1:2] <- 0.1
mat[3:4,3:5] <- 0.1
mat
## I think for now, you can keep the "off-diagonal" blocks to be 0 still
## - as per usual, then you can go in and change any term you want
## manually (as you've done before -- just check that the overall matrix is still PSD).

### MEAN VECTORS AND MONOTONIC TRANSFORMATION
# Recall: Our end goal is to design a meaningful
# mean-vector and covariance-matrix for each of the
# 3 cell types. And all we're trying to do right now
# is to "pick" the "best looking parts" of both
# datasets and slapping them together to get
# our desired mean vectors/covariance matrices.
# (Of course, we'll need to choose the mean vector as well)
# this is akin to:
p <- 10
idx <- c(1,2,5,7,9)
corpula_mean_first <- c(10, 0, 10, 0, 10,
                        0, 10, 0, 10, 0) 
corpula_mean_first[idx] # these will be the 5 values in the new mean vector we want
# so copula function could be different among the 3 cell-types
# but when you're "merging" the data-generation-process
# of the two datasets, you can apply different copula functions to different
# variables (So realistically, you'll use the same copula functions for each 
# of the respective variables that you've already been using)
# so the way you write corpula_third would need to be a little bit more complicated now:
# for example, it might look like this:
corpula_third <- sapply(1:ncol(corpula_third), function(j){
  if(j %in% c(1:4)){
    0.4 * sign(x) * abs(x)^1.4
  } else { 
    0.6 * sign(x) * abs(x)^1.4
  }
})

