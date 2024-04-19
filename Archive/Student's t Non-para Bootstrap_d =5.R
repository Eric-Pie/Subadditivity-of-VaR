### 0 Setup ####################################################################
#install.packages("simsalapar")
#install.packages("copula")
library(simsalapar)
library(xtable)
require(copula)
n.obs <- 1000 #% sample size
alpha <- 1-10^(seq(log10(1-0.8), log10(1-0.999), length.out = 100))
N.sim <- 5000  
d <- 5
tau <- 0.5
doExtras <- simsalapar:::doExtras()

### 1 Generate Data ############################################################
set.seed(123)
## generate a sample from a t-distribution with degree of freedom 4, Kendall's tau 0.5, dimension 5 and t copula.
cop_t <- ellipCopula("t", param=iTau(ellipCopula("t"), tau=tau), dim=d)
original_U <- rCopula(n.obs, cop_t)
qmargin_t4 <- function(p){
  qt(p, df=4)
}
X <- apply(original_U, 2, qmargin_t4)

### 2 Student's t Non-para Bootstrap ###########################################
## Student's t Non-para Bootstrap
t_boot <- function(X, N.sim){
  n <- nrow(X)
  d <- ncol(X)
  t_boot <- list(NA, N.sim)
  for(i in 1:N.sim){
    ## each row in X is a sample with dimiension 5, we resample n rows from X with replacement to generate a new sample space.
    X_boot <- X[sample(1:n.obs, n.obs, replace=TRUE), ]
    ## store the X_boot into t_boot
    t_boot[[i]] <- X_boot
  }
  t_boot
}
Bootstrap_samples <- t_boot(X, N.sim)

### 3 For each b bootstrap sample, calculate Y_b and CI ########################

## First, we calculate s_b
## create the ecdf for each dimension of b th bootstrap sample
#ecdf_b_Xj <- function(X_b){
#  n <- nrow(X_b)
#  d <- ncol(X_b)
#  ecdf_b <- list(NA, d)
#  for(j in 1:d){
#    ecdf_b[[j]] <- ecdf(X_b[, j])
#  }
#  ecdf_b
#}

#List_ecdf_b <- list(data = NA, N.sim)

#for (i in 1:N.sim){
#  List_ecdf_b[[i]] <- ecdf_b_Xj(Bootstrap_samples[[i]])
#}

## function to calculate s_b
s_b <- function(b,boots_samples,alpha){
  Matrix_sorted <- matrix(NA, nrow = n.obs, ncol = d)
  for (j in 1:d){
    data_array <- Bootstrap_samples[[b]][,j]
    data_array_sorted <- sort(data_array)
    Matrix_sorted[,j] <- data_array_sorted
    }
  ## sum up the row of the sorted matrix
  s_b <- sum(Matrix_sorted[ceiling(alpha*n.obs),])
  return(s_b)
}

Matrix_s_b <- matrix(data = NA, nrow = length(alpha), ncol = N.sim)
for (i in 1:length(alpha)){
  for (b in 1:N.sim){
    Matrix_s_b[i,b] <- s_b(b, Bootstrap_samples, alpha[i])
  }
}

## add dimnames for Matrix_s_b
dimnames(Matrix_s_b) <- list(alpha, 1:N.sim)

## Second, we calculate Y_b
## calculate all S_bi
S_bi_list <- list(NA, N.sim)
for (b in 1:N.sim){
  S_bi_list[[b]] <- rowSums(Bootstrap_samples[[b]])
}

## calculate Y_b and CI
Tilde_alpha <- 0.05 ## significance level
Matrix_Y_b <- matrix(data = NA, nrow = length(alpha), ncol = N.sim)
Matrix_bounds <- matrix(data = NA, nrow = length(alpha), ncol = 3)
for (i in 1:length(alpha)){
  for (b in 1:N.sim){
    Matrix_Y_b[i,b] <- ecdf(S_bi_list[[b]])(Matrix_s_b[i,b])
  }
  array_Y_b <- Matrix_Y_b[i,]
  array_Y_b_sorted <- sort(array_Y_b)
  L_n <- array_Y_b_sorted[ceiling(Tilde_alpha*N.sim)]
  U_n <- array_Y_b_sorted[ceiling((1-Tilde_alpha)*N.sim)]
  Matrix_bounds[i,] <- c(L_n, U_n, alpha[i])
}

### 4 Output ###################################################################
## draw a plot
## draw a plot for d = 5,
matplot((1-10^(seq(log10(1-0.8), log10(1-0.999), length.out = 100))), Matrix_bounds, type = "l", lty = 1, col = c("blue", "red","brown"), xlab = expression(alpha), ylab = substitute(Bootstrap~CI~of~F[S](d*{F^{-1}}(alpha))~"for"~l~alpha~"in [0.8,1]"~(d~"="~5), list(l = length(alpha))),xlim = c(0.8,1),ylim = c(0.8,1))

## add the legend
legend("bottomright", legend = c('Lower bound (if '~alpha~'< Lower bound, VaR is subadditive)','Upper bound (if '~alpha~'> Upper bound, VaR is superadditive) ',expression(alpha)), col = c("blue", "red","brown"), lty = 1, cex = 0.5,bty = "n")
