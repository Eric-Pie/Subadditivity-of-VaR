### 0 Setup ####################################################################
#install.packages("simsalapar")
#install.packages("copula")
#install.packages("xtable")
library(simsalapar)
library(xtable)
#require(copula)
n.obs <- 2000 #% sample size
## generate a sequence of alpha having more value close to 0
alpha <- 1-10^(seq(log10(1-0.001), log10(1-0.1), length.out = 200))  #% alpha
N.sim <- 2000 #% number of simulations
doExtras <- simsalapar:::doExtras()


### 1 Generate Data ############################################################
## d = 4
set.seed(139)
## generate a sample from a bernoulli distribution
original_U <- matrix(runif(n.obs * 4), nrow = n.obs, ncol = 4)
qber <- function(p) {
  qbinom(p, size = 1, prob = 0.5)
}
X_4 <- apply(original_U, 2, qber)

## d = 5
set.seed(140)
## generate a sample from a bernoulli distribution
original_U <- matrix(runif(n.obs * 5), nrow = n.obs, ncol = 5)
qber <- function(p) {
  qbinom(p, size = 1, prob = 0.5)
}
X_5 <- apply(original_U, 2, qber)

## d = 6
set.seed(143)
## generate a sample from a bernoulli distribution
original_U <- matrix(runif(n.obs * 6), nrow = n.obs, ncol = 6)
qber <- function(p) {
  qbinom(p, size = 1, prob = 0.5)
}
X_6 <- apply(original_U, 2, qber)

## d = 7
set.seed(144)
## generate a sample from a bernoulli distribution
original_U <- matrix(runif(n.obs * 7), nrow = n.obs, ncol = 7)
qber <- function(p) {
  qbinom(p, size = 1, prob = 0.5)
}
X_7 <- apply(original_U, 2, qber)


### 2 Bernoulli Non-para Bootstrap ###########################################
## Bernoulli Non-para Bootstrap
ber_boot <- function(X, N.sim){
  n <- nrow(X)
  d <- ncol(X)
  ber_boot <- list(NA, N.sim)
  for(i in 1:N.sim){
    ## each row in X is a sample with dimiension d, we resample n rows from X with replacement to generate a new sample space.
    X_boot <- X[sample(1:n.obs, n.obs, replace=TRUE), ]
    ## store the X_boot into t_boot
    ber_boot[[i]] <- X_boot
  }
  ber_boot
}
Bootstrap_samples_4 <- ber_boot(X_4, N.sim)
Bootstrap_samples_5 <- ber_boot(X_5, N.sim)
Bootstrap_samples_6 <- ber_boot(X_6, N.sim)
Bootstrap_samples_7 <- ber_boot(X_7, N.sim)

### 3 For each b bootstrap sample, calculate Y_b and CI ########################

## d = 4
s_b_4 <- function(b,boots_samples,alpha){
  Matrix_sorted <- matrix(NA, nrow = n.obs, ncol = 4)
  for (j in 1:4){
    data_array <- boots_samples[[b]][,j]
    data_array_sorted <- sort(data_array)
    Matrix_sorted[,j] <- data_array_sorted
  }
  ## sum up the row of the sorted matrix
  s_b <- sum(Matrix_sorted[ceiling(alpha*n.obs),])
  return(s_b)
}

Matrix_s_b_4 <- matrix(data = NA, nrow = length(alpha), ncol = N.sim)
for (i in 1:length(alpha)){
  for (b in 1:N.sim){
    Matrix_s_b_4[i,b] <- s_b_4(b, Bootstrap_samples_4, alpha[i])
  }
}


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

## d= 6
s_b_6 <- function(b,boots_samples,alpha){
  Matrix_sorted <- matrix(NA, nrow = n.obs, ncol = 6)
  for (j in 1:6){
    data_array <- boots_samples[[b]][,j]
    data_array_sorted <- sort(data_array)
    Matrix_sorted[,j] <- data_array_sorted
  }
  ## sum up the row of the sorted matrix
  s_b <- sum(Matrix_sorted[ceiling(alpha*n.obs),])
  return(s_b)
}

Matrix_s_b_6 <- matrix(data = NA, nrow = length(alpha), ncol = N.sim)
for (i in 1:length(alpha)){
  for (b in 1:N.sim){
    Matrix_s_b_6[i,b] <- s_b_6(b, Bootstrap_samples_6, alpha[i])
  }
}
## d = 7
s_b_7 <- function(b,boots_samples,alpha){
  Matrix_sorted <- matrix(NA, nrow = n.obs, ncol = 7)
  for (j in 1:7){
    data_array <- boots_samples[[b]][,j]
    data_array_sorted <- sort(data_array)
    Matrix_sorted[,j] <- data_array_sorted
  }
  ## sum up the row of the sorted matrix
  s_b <- sum(Matrix_sorted[ceiling(alpha*n.obs),])
  return(s_b)
}

Matrix_s_b_7 <- matrix(data = NA, nrow = length(alpha), ncol = N.sim)
for (i in 1:length(alpha)){
  for (b in 1:N.sim){
    Matrix_s_b_7[i,b] <- s_b_7(b, Bootstrap_samples_7, alpha[i])
  }
}

## add dimnames for Matrix_s_b
dimnames(Matrix_s_b_4) <- list(alpha, 1:N.sim)
dimnames(Matrix_s_b_5) <- list(alpha, 1:N.sim)
dimnames(Matrix_s_b_6) <- list(alpha, 1:N.sim)
dimnames(Matrix_s_b_7) <- list(alpha, 1:N.sim)

## Second, we calculate Y_b
## calculate all S_bi
## d = 4
S_bi_list_4 <- list(NA, N.sim)
for (b in 1:N.sim){
  S_bi_list_4[[b]] <- rowSums(Bootstrap_samples_4[[b]])
}

## d = 5
S_bi_list_5 <- list(NA, N.sim)
for (b in 1:N.sim){
  S_bi_list_5[[b]] <- rowSums(Bootstrap_samples_5[[b]])
}

## d = 6
S_bi_list_6 <- list(NA, N.sim)
for (b in 1:N.sim){
  S_bi_list_6[[b]] <- rowSums(Bootstrap_samples_6[[b]])
}

## d = 7
S_bi_list_7 <- list(NA, N.sim)
for (b in 1:N.sim){
  S_bi_list_7[[b]] <- rowSums(Bootstrap_samples_7[[b]])
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

## d =4
Matrix_Y_b_4 <- matrix(data = NA, nrow = length(alpha), ncol = N.sim)
for (i in 1:length(alpha)){
  for (b in 1:N.sim){
    Matrix_Y_b_4[i,b] <- ecdf(S_bi_list_4[[b]])(Matrix_s_b_4[i,b])
  }
}


## d =5
Matrix_Y_b_5 <- matrix(data = NA, nrow = length(alpha), ncol = N.sim)
for (i in 1:length(alpha)){
  for (b in 1:N.sim){
    Matrix_Y_b_5[i,b] <- ecdf(S_bi_list_5[[b]])(Matrix_s_b_5[i,b])
  }
}

## d = 6
Matrix_Y_b_6 <- matrix(data = NA, nrow = length(alpha), ncol = N.sim)
for (i in 1:length(alpha)){
  for (b in 1:N.sim){
    Matrix_Y_b_6[i,b] <- ecdf(S_bi_list_6[[b]])(Matrix_s_b_6[i,b])
  }
}

##d = 7
Matrix_Y_b_7 <- matrix(data = NA, nrow = length(alpha), ncol = N.sim)
for (i in 1:length(alpha)){
  for (b in 1:N.sim){
    Matrix_Y_b_7[i,b] <- ecdf(S_bi_list_7[[b]])(Matrix_s_b_7[i,b])
  }
}

## calculate all t-statistics and store in the list
DD_4 <- matrix(0,nrow = length(alpha),ncol = N.sim)
DD_5 <- matrix(0,nrow = length(alpha),ncol = N.sim)
DD_6 <- matrix(0,nrow = length(alpha),ncol = N.sim)
DD_7 <- matrix(0,nrow = length(alpha),ncol = N.sim)

## calculate t-statistics for d = 4
for (i in 1:N.sim){
  for (j in 1:length(alpha)){
    DD_4[j,i] <- t_statistic(n.obs,Matrix_Y_b_4[j,i],alpha[j])
  }
}

## calculate t-statistics for d = 5
for (i in 1:N.sim){
  for (j in 1:length(alpha)){
    DD_5[j,i] <- t_statistic(n.obs,Matrix_Y_b_5[j,i],alpha[j])
  }
}

## d = 6
for (i in 1:N.sim){
  for (j in 1:length(alpha)){
    DD_6[j,i] <- t_statistic(n.obs,Matrix_Y_b_6[j,i],alpha[j])
  }
}

## d = 7
for (i in 1:N.sim){
  for (j in 1:length(alpha)){
    DD_7[j,i] <- t_statistic(n.obs,Matrix_Y_b_7[j,i],alpha[j])
  }
}

## calculate the rejection rate for d = 4
reject_rate_4 <- matrix(0,nrow = length(alpha),ncol = 1)
for (i in 1:length(alpha)){
  reject_rate_4[i,1] <- mean(sapply(DD_4[i,],reject,Tilde_alpha))
}

## calculate the rejection rate for d = 5
reject_rate_5 <- matrix(0,nrow = length(alpha),ncol = 1)
for (i in 1:length(alpha)){
  reject_rate_5[i,1] <- mean(sapply(DD_5[i,],reject,Tilde_alpha))
}

## calculate the rejection rate for d = 6
reject_rate_6 <- matrix(0,nrow = length(alpha),ncol = 1)
for (i in 1:length(alpha)){
  reject_rate_6[i,1] <- mean(sapply(DD_6[i,],reject,Tilde_alpha))
}

## calculate the rejection rate for d = 7
reject_rate_7 <- matrix(0,nrow = length(alpha),ncol = 2)
for (i in 1:length(alpha)){
  reject_rate_7[i,1] <- mean(sapply(DD_7[i,],reject,Tilde_alpha))
  reject_rate_7[i,2] <- Tilde_alpha
}

## combine the rejection rate for d =4,5,6,7 into one matrix
reject_rate <- cbind(reject_rate_4,reject_rate_5,reject_rate_6,reject_rate_7)


## plot the rejection rate
matplot(log10(1-10^(seq(log10(1-0.001), log10(1-0.1), length.out = 200))), reject_rate, type = "l", lty = 1, col = c("black","blue","orange","purple","brown"), xlab = expression(log~alpha), ylab = substitute(Probability~of~Rejecting~Null~Hypothesis~'for'~d~'='~'4,5,6,7'~(n.obs~"="~2000)),xlim = c(-3,-1),ylim = c(0,1))

## add the legend
legend("left", legend = c('rejection probability for d = 4','rejection probability for d = 5','rejection probability for d = 6','rejection probability for d = 7','significance level = 0.05'), col = c("black","blue","orange","purple","brown"), lty = 1, cex = 0.6,bty = "n")

## add points on the cross point of the two lines with ptch = 19 and col = black
points(-1.205,0.05,pch = 19,col = "black",cex = 0.8)
points(-1.58,0.05,pch = 19,col = "black",cex = 0.8)
points(-1.845,0.05,pch = 19,col = "black",cex = 0.8)
points(-2.05,0.05,pch = 19,col = "black",cex = 0.8)


## vertical dashed line
abline(v = -1.235,lty = 2)
abline(v = -1.55,lty = 2)
abline(v = -1.845,lty = 2)
abline(v = -2.15,lty = 2)


## add the x value of points at the top of the plot
axis(3, at = -1.235, labels = -1.235)
axis(3, at = -1.55, labels = -1.55)
axis(3, at = -1.845, labels = -1.845)
axis(3, at = -2.15, labels = -2.15)





