### 0 Setup ####################################################################
#install.packages("simsalapar")
#install.packages("copula")
library(simsalapar)
library(xtable)
require(copula)
n.obs <- 1000 #% sample size
alpha <- 1-10^(seq(log10(1-0.01), log10(1-0.999), length.out = 200))
N.sim <- 5000  
tau <- 0.5
doExtras <- simsalapar:::doExtras()

### 1 Generate Data ############################################################
## d = 5
set.seed(123)
## generate a sample from a t-distribution with degree of freedom 4
cop_t <- ellipCopula("t", param=iTau(ellipCopula("t"), tau=tau), dim=5)
original_U <- rCopula(n.obs, cop_t)
qmargin_t4 <- function(p){
  qt(p, df=4)
}
X_5 <- apply(original_U, 2, qmargin_t4)

## d=10
set.seed(124)
## generate a sample from a t-distribution with degree of freedom 4
cop_t <- ellipCopula("t", param=iTau(ellipCopula("t"), tau=tau), dim=10)
original_U <- rCopula(n.obs, cop_t)
qmargin_t4 <- function(p){
  qt(p, df=4)
}
X_10 <- apply(original_U, 2, qmargin_t4)

## d=20
set.seed(125)
## generate a sample from a t-distribution with degree of freedom 4
cop_t <- ellipCopula("t", param=iTau(ellipCopula("t"), tau=tau), dim=20)
original_U <- rCopula(n.obs, cop_t)
qmargin_t4 <- function(p){
  qt(p, df=4)
}
X_20 <- apply(original_U, 2, qmargin_t4)

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
Bootstrap_samples_5 <- t_boot(X_5, N.sim)
Bootstrap_samples_10 <- t_boot(X_10, N.sim)
Bootstrap_samples_20 <- t_boot(X_20, N.sim)

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
## d = 5
s_b_5 <- function(b,boots_samples,alpha){
  Matrix_sorted <- matrix(NA, nrow = n.obs, ncol = 5)
  for (j in 1:5){
    data_array <- boots_samples[[b]][,j]
    data_array_sorted <- sort(data_array)
    Matrix_sorted[,j] <- data_array_sorted
  }
  ## sum up the row of the sorted matrix
  s_b <- sum(Matrix_sorted[ceiling(alpha*n.obs),])
  return(s_b)
}

Matrix_s_b_5 <- matrix(data = NA, nrow = length(alpha), ncol = N.sim)
for (i in 1:length(alpha)){
  for (b in 1:N.sim){
    Matrix_s_b_5[i,b] <- s_b_5(b, Bootstrap_samples_5, alpha[i])
  }
}

## d= 10
s_b_10 <- function(b,boots_samples,alpha){
  Matrix_sorted <- matrix(NA, nrow = n.obs, ncol = 10)
  for (j in 1:10){
    data_array <- boots_samples[[b]][,j]
    data_array_sorted <- sort(data_array)
    Matrix_sorted[,j] <- data_array_sorted
  }
  ## sum up the row of the sorted matrix
  s_b <- sum(Matrix_sorted[ceiling(alpha*n.obs),])
  return(s_b)
}

Matrix_s_b_10 <- matrix(data = NA, nrow = length(alpha), ncol = N.sim)
for (i in 1:length(alpha)){
  for (b in 1:N.sim){
    Matrix_s_b_10[i,b] <- s_b_10(b, Bootstrap_samples_10, alpha[i])
  }
}

## d= 20
s_b_20 <- function(b,boots_samples,alpha){
  Matrix_sorted <- matrix(NA, nrow = n.obs, ncol = 20)
  for (j in 1:20){
    data_array <- boots_samples[[b]][,j]
    data_array_sorted <- sort(data_array)
    Matrix_sorted[,j] <- data_array_sorted
  }
  ## sum up the row of the sorted matrix
  s_b <- sum(Matrix_sorted[ceiling(alpha*n.obs),])
  return(s_b)
}

Matrix_s_b_20 <- matrix(data = NA, nrow = length(alpha), ncol = N.sim)
for (i in 1:length(alpha)){
  for (b in 1:N.sim){
    Matrix_s_b_20[i,b] <- s_b_20(b, Bootstrap_samples_20, alpha[i])
  }
}


## add dimnames for Matrix_s_b
dimnames(Matrix_s_b_5) <- list(alpha, 1:N.sim)
dimnames(Matrix_s_b_10) <- list(alpha, 1:N.sim)
dimnames(Matrix_s_b_20) <- list(alpha, 1:N.sim)

## Second, we calculate Y_b
## calculate all S_bi
## d = 5
S_bi_list_5 <- list(NA, N.sim)
for (b in 1:N.sim){
  S_bi_list_5[[b]] <- rowSums(Bootstrap_samples_5[[b]])
}

## d = 10
S_bi_list_10 <- list(NA, N.sim)
for (b in 1:N.sim){
  S_bi_list_10[[b]] <- rowSums(Bootstrap_samples_10[[b]])
}

## d =20
S_bi_list_20 <- list(NA, N.sim)
for (b in 1:N.sim){
  S_bi_list_20[[b]] <- rowSums(Bootstrap_samples_20[[b]])
}


## define the t_statistic function
t_statistic <- function(n,dd,alpha)
{
  sqrt(n)*(dd - alpha)/sqrt(dd*(1-dd))
}

## the function defined to see if we reject null hypothesis
reject <- function(T_stat,tilde_alpha)
{
  if (T_stat < qnorm(tilde_alpha)){
    return(1)   ## 1 represents superadditive
  }
  else if (T_stat >= qnorm(tilde_alpha)){
    return(0)   ## 0 represents subadditive
  }
}

## calculate Y_b and do hypothesis test
Tilde_alpha <- 0.05 ## significance level
## d =5
Matrix_Y_b_5 <- matrix(data = NA, nrow = length(alpha), ncol = N.sim)
for (i in 1:length(alpha)){
  for (b in 1:N.sim){
    Matrix_Y_b_5[i,b] <- ecdf(S_bi_list_5[[b]])(Matrix_s_b_5[i,b])
  }
}

## d = 10
Matrix_Y_b_10 <- matrix(data = NA, nrow = length(alpha), ncol = N.sim)
for (i in 1:length(alpha)){
  for (b in 1:N.sim){
    Matrix_Y_b_10[i,b] <- ecdf(S_bi_list_10[[b]])(Matrix_s_b_10[i,b])
  }
}

## d =20
Matrix_Y_b_20 <- matrix(data = NA, nrow = length(alpha), ncol = N.sim)
for (i in 1:length(alpha)){
  for (b in 1:N.sim){
    Matrix_Y_b_20[i,b] <- ecdf(S_bi_list_20[[b]])(Matrix_s_b_20[i,b])
  }
}

## calculate all t-statistics and store in the list
DD_5 <- matrix(0,nrow = length(alpha),ncol = N.sim)
DD_10 <- matrix(0,nrow = length(alpha),ncol = N.sim)
DD_20 <- matrix(0,nrow = length(alpha),ncol = N.sim)

## calculate t-statistics for d = 5
for (i in 1:N.sim){
  for (j in 1:length(alpha)){
    DD_5[j,i] <- t_statistic(n.obs,Matrix_Y_b_5[j,i],alpha[j])
  }
}

## calculate t-statistics for d = 10
for (i in 1:N.sim){
  for (j in 1:length(alpha)){
    DD_10[j,i] <- t_statistic(n.obs,Matrix_Y_b_10[j,i],alpha[j])
  }
}

## calculate t-statistics for d = 20
for (i in 1:N.sim){
  for (j in 1:length(alpha)){
    DD_20[j,i] <- t_statistic(n.obs,Matrix_Y_b_20[j,i],alpha[j])
  }
}

## calculate the rejection rate for d = 5
reject_rate_5 <- matrix(0,nrow = length(alpha),ncol = 1)
for (i in 1:length(alpha)){
  reject_rate_5[i,1] <- mean(sapply(DD_5[i,],reject,Tilde_alpha))
}

## calculate the rejection rate for d = 10
reject_rate_10 <- matrix(0,nrow = length(alpha),ncol = 1)
for (i in 1:length(alpha)){
  reject_rate_10[i,1] <- mean(sapply(DD_10[i,],reject,Tilde_alpha))
}

## calculate the rejection rate for d = 20
reject_rate_20 <- matrix(0,nrow = length(alpha),ncol = 2)
for (i in 1:length(alpha)){
  reject_rate_20[i,1] <- mean(sapply(DD_20[i,],reject,Tilde_alpha))
  reject_rate_20[i,2] <- Tilde_alpha
}

# combine the rejection rate for d =5, d = 10 and d = 20 into one matrix
reject_rate <- cbind(reject_rate_5,reject_rate_10,reject_rate_20)



## transform the plot into a pdf file
pdf("t Hypothesis_bootstrap_graph_n.obs=1000.pdf")
matplot((1-10^(seq(log10(1-0.001), log10(1-0.999), length.out = 200))), reject_rate, type = "l", lty = 1, col = c("blue", "red","green","brown"), xlab = expression(alpha), ylab = substitute(Probability~of~Rejecting~Null~Hypothesis~'for'~d~'='~'5,10,20'~(n.obs~"="~1000)),xlim = c(0,1),ylim = c(0,1))
legend("right", legend = c('rejection probability for d = 5','rejection probability for d = 10','rejection probability for d = 20','significance level = 0.05'), col = c("blue", "red","green","brown"), lty = 1, cex = 1,bty = "n")
abline(v = 0.5, lty = 2)
axis(1, at = 0.5, labels = 0.5)
points(0.5,0.05,pch = 19,col = "black",cex = 1)
dev.off()