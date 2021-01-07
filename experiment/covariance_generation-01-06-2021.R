rm(list=ls())

set.seed(10)
n <- 500
mean_vec_1 <- c(20,0)
cov_mat_1 <- matrix(c(100, 0,
                      0, 50), 2, 2)
####
mean_vec_2 <- c(0,20)
cov_mat_2 <- matrix(c(50, 0,
                      0, 100), 2, 2)

dat_1 <- MASS::mvrnorm(n, mu = mean_vec_1, Sigma = cov_mat_1)
dat_1 <- apply(dat_1, 2, function(x){x^2}) # copula-part (monotonic transformation)
dat_2 <- MASS::mvrnorm(n, mu = mean_vec_2, Sigma = cov_mat_2)
dat_2 <- apply(dat_2, 2, function(x){x^2})

dat <- rbind(dat_1, dat_2)
# dat_3 <- MASS::mvrnorm(n, mu = c(0,0), Sigma = diag(2))
# dat <- rbind(dat_1, dat_2, dat_3)

# par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(dat[,1], dat[,2], asp = T,
     pch = 16, col = rgb(0.5, 0.5, 0.5, 0.1))
# lines(rep(20,2), c(-100,100), col = "red", lty = 2)
# lines(rep(-20,2), c(-100,100), col = "red", lty = 2)

# Thing to remember #1: the mean and covariance are closely
## related to make this "L-shape" work.
## Remember you can play around with all 4 values, but
## changing some of the values might require you to
## also change the mean to get the "L-shape" to work.

# Thing to remember #2: The main appeal of using copulas
## is to get "non-gaussian-looking" data, of which
## this might be useful in getting "skewed" distributions
## so that most of the mass in the "L-shape" are closer
## to the origin. Then, to achieve this, one could consider
## a monotonic function such as f(x) = x^2

# Thing to remember #3: This is currently a 2-dimensional
## example, but if we want to consider the 8-dimensional
## example, we can just think of each 2x2 submatrix
## at a time (within the 8x8 covariance matrix)

# Thing to remember #4: Remember your covariance matrices
## need to be PSD (positive semi-definite). 
## You can always explicitly check the eigenvalues of cov_mat 
## to make sure it's PSD.
## Example:
cov_mat <- matrix(c(1, 5,
                    5, 1), 2, 2)
eigen(cov_mat)$values # this means this is a indefinite matrix (some positive, some negative eigenvalues)
MASS::mvrnorm(n, mu = c(0,0), Sigma = cov_mat)

# https://en.wikipedia.org/wiki/Sylvester%27s_criterion
### what is this: This is a convenient "trick" to quickly see if a matrix is PSD or not
### key ideas behind the trick: We have a simple formula for the determinant of a 2x2 matrix
### For a matrix mat(c(a,b,c,d),2,2), the determinant is ad-bc
### Why is this even relevant: We know there's a close connection b/w determinants and eigenvalues
### The relation: the product of the eigenvalues is the determinant
### So why do we care about this for PSD matrices: 
### Because if the determinant is negative for a 2x2 matrix, that MUST mean
### you had one positive and one negative eigenvalue (and hence, not PSD)
## Note: For 2x2 matrices, ad-bc > 0 means that the product of diagonal entries 
## must be larger than the product of the off-diagonal <- the most important line.
#### Recall though, we're working with an 8x8 matrix, so this is where 
## we need Slyvester's criterion to help us along the way.
#### There two ways to think about this:
## 1) [the literal way, on what the criterion says -- the wiki page] 
## A matrix is PSD if and only if (i.e is an equivalent definition essentially)
## if all the upper-left square matrices have a positive determinant
## Remark 1: This is how we rule of "two negative eigenvalues resulting in a positive determinant"
## Remark 2: It's kinda hard to use in practical tho, because we don't have nice simple formulas
## for determinant for 3x3, 4x4, etc. matrices
## 2) [the practical way]
## A /necessary/ (but not sufficient) condition for a matrix to be PSD
## (in simple words: If this trick says it's a PSD matrix, it might not actually be,
### but all PSD matrices will satisfy this trick)
## For an 8x8 matrix, pick 2 diagonal entries and check all the 
## determinant for their corresponding 2x2 submatrix. If ALL possible 2x2
## matrices chosen in this way have a positive determinant, then we can
## probably intuit that the 8x8 matrix is PSD.
## This works: because the determinant is unchanged when we swap rows&columns accordingly
## Hence, we can make ANY two variables the top-two variables in our 8x8 matrix
## and then apply Slyvester's criterion (meaning we'll be checking the determinant of 
## that top-left 2x2 matrix)

## Takeaway: When we're thinking about whether or not an 8x8 matrix is PSD
## it's "good enough" for us (at least, to get some intuition) to consider
## all possible 2x2 submatrices along the diagonal of the 8x8 matrix.
## [example: take the 2x2 matrix formed by: (1,1), (1,5), (5,1), (5,5)]

## ^^^ this is really just a trick to quickly see if a matrix is PSD (but note,
## is a "necessary but not sufficient", so it's not 100% reliable, but it's good
## enough for us to get some intuition on when a matrix might be PSD)

## Takeaway #2: For the most part, if your covariance matrix is "diagonally dominant",
## then you should fine. [the largest entry in each row, is the diagonal entry]

