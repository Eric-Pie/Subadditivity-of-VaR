### 0 Setup ####################################################################
#install.packages("simsalapar")
#install.packages("copula")
#install.packages("xtable")
library(simsalapar)
library(xtable)
#require(copula)
n.obs <- 1000 #% sample size
alpha <- 1-10^(seq(log10(1-0.5), log10(1-0.5), length.out = 100))  #% alpha
doExtras <- simsalapar:::doExtras()


## list of variables
varList <- varlist(
  ## sample size
  n = list(value = n.obs),
  ## margins
  qmargin = list(type="frozen", expr = quote(F^{-1}), #% F[j] is changed to F^{-1}.
                 value = c(Ber_0.5 = function(p) qbinom(p,size = 1, prob = 0.5))), ## bernoulli quantile function with probability 0.5
  ## VaR confidence levels
  alpha = list(type="inner", value = alpha))

#getEl(varList)

## function defined to compute F_{X_1+..+X_d}(d*F_1^-(\alpha))
doOne <- function(n, d, qmargin, alpha)
{
  U <- matrix(runif(n * d), nrow = n, ncol = d)
  ## compute F_{X_1+..+X_d}(d*F_1^-(\alpha)) for all confidence levels alpha
  ## => VaR_alpha superadditive <=> F_{X_1+..+X_d}(d*F_1^-(\alpha)) - alpha < 0
  t(sapply(qmargin, function(FUN) ecdf(rowSums(FUN(U)))(d*FUN(alpha)))) 
  ## note: t() is important here, since, otherwise, the order of the variables
  ## ----  would not be correct (=> check should reveal this) 
}

## apply doOne to all samples in varList
nonGr <- get.nonGrids(varList)$nonGrids
dd_5 <- doOne(n= min(nonGr$n, 1000), d=5,
              qmargin=nonGr$qmargin, alpha=nonGr$alpha)
dd_10 <- doOne(n= min(nonGr$n, 1000), d=10,
               qmargin=nonGr$qmargin, alpha=nonGr$alpha)
dd_20 <- doOne(n= min(nonGr$n, 1000), d=20,
                qmargin=nonGr$qmargin, alpha=nonGr$alpha)
#stopifnot(dim(dd) == with(nonGr, c(length(qmargin), length(alpha))))

DD <- matrix(data = NA, nrow = length(alpha), ncol = 3)
DD[,1] <- t(dd_5)
DD[,2] <- t(dd_10)
DD[,3] <- t(dd_20)

## the function defined to calculate the Lower bound and Upper bound of CI of F_s at significance level Tilde_alpha 
L_n <- function(dd,n,tilde_alpha){
  dd - qnorm(1-tilde_alpha)*sqrt(dd*(1-dd))/sqrt(n)
}
U_n<- function(dd,n,tilde_alpha){
  dd + qnorm(1-tilde_alpha)*sqrt(dd*(1-dd))/sqrt(n)
}

data_array_5 <- c()## build up a array to store L_n and U_n, d = 5
data_array_10<- c()## build up a array to store L_n and U_n, d = 10
data_array_20<- c()## build up a array to store L_n and U_n, d = 20
Tilde_alpha <- 0.05 ## significance level of CI
## compute all results  
for (j in 1:length(alpha)){
  L_n_dd <- L_n(dd = DD[j,1],n=n.obs,tilde_alpha=Tilde_alpha)
  U_n_dd <- U_n(dd = DD[j,1],n=n.obs,tilde_alpha=Tilde_alpha)
  
  data_array_5 <- append(data_array_5, L_n_dd)
  data_array_5 <- append(data_array_5, U_n_dd)
  data_array_5 <- append(data_array_5, alpha[j])
}

for (j in 1:length(alpha)){
  L_n_dd <- L_n(dd = DD[j,2],n=n.obs,tilde_alpha=Tilde_alpha)
  U_n_dd <- U_n(dd = DD[j,2],n=n.obs,tilde_alpha=Tilde_alpha)
  
  data_array_10 <- append(data_array_10, L_n_dd)
  data_array_10 <- append(data_array_10, U_n_dd)
  data_array_210 <- append(data_array_10, alpha[j])
}

for (j in 1:length(alpha)){
  L_n_dd <- L_n(dd = DD[j,3],n=n.obs,tilde_alpha=Tilde_alpha)
  U_n_dd <- U_n(dd = DD[j,3],n=n.obs,tilde_alpha=Tilde_alpha)
  
  data_array_20 <- append(data_array_20, L_n_dd)
  data_array_20 <- append(data_array_20, U_n_dd)
  data_array_20 <- append(data_array_20, alpha[j])
}


## build up a matrix to store the data
Matrix_5 <- matrix(data = data_array_5, nrow = length(alpha), ncol = 3, byrow = TRUE)
Matrix_10 <- matrix(data = data_array_10, nrow = length(alpha), ncol = 3, byrow = TRUE)
Matrix_20 <- matrix(data = data_array_20, nrow = length(alpha), ncol = 3, byrow = TRUE)


## add the dim names
dimnames(Matrix_5) <- list(alpha,c('Lower bound','Upper bound',expression(alpha)))
dimnames(Matrix_10) <- list(alpha,c('Lower bound','Upper bound',expression(alpha)))
dimnames(Matrix_20) <- list(alpha,c('Lower bound','Upper bound',expression(alpha)))

## draw a plot for d = 5,
matplot(seq(0.6,0.7,length = 100), Matrix_5, type = "l", lty = 1, col = c("blue", "red","brown"), xlab = expression(alpha), ylab = substitute(CI~of~F[S](d*{F^{-1}}(alpha))~"for"~l~alpha~"in [0.8,1]"~(d~"="~5), list(l = length(alpha))),xlim = c(0.6,0.7),ylim = c(0,1))

## add the legend
legend("right", legend = c('Lower bound (if '~alpha~'< Lower bound, VaR is subadditive)','Upper bound (if '~alpha~'> Upper bound, VaR is superadditive) ',expression(alpha)), col = c("blue", "red","brown"), lty = 1, cex = 0.5,bty = "n")

## draw a plot for d = 10,
matplot(seq(0.6,0.7,length = 100), Matrix_20, type = "l", lty = 1, col = c("blue", "red","brown"), xlab = expression(alpha), ylab = substitute(CI~of~F[S](d*{F^{-1}}(alpha))~"for"~l~alpha~"in [0.8,1]"~(d~"="~20), list(l = length(alpha))),xlim = c(0.6,0.7),ylim = c(0,1))

## add the legend
legend("right", legend = c('Lower bound (if '~alpha~'< Lower bound, VaR is subadditive)','Upper bound (if '~alpha~'> Upper bound, VaR is superadditive) ',expression(alpha)), col = c("blue", "red","brown"), lty = 1, cex = 0.5,bty = "n")
## draw a plot for d = 20,
matplot(seq(0.6,0.7,length = 100), Matrix_100, type = "l", lty = 1, col = c("blue", "red","brown"), xlab = expression(alpha), ylab = substitute(CI~of~F[S](d*{F^{-1}}(alpha))~"for"~l~alpha~"in [0.8,1]"~(d~"="~100), list(l = length(alpha))),xlim = c(0.6,0.7),ylim = c(0,1))

## add the legend
legend("right", legend = c('Lower bound (if '~alpha~'< Lower bound, VaR is subadditive)','Upper bound (if '~alpha~'> Upper bound, VaR is superadditive) ',expression(alpha)), col = c("blue", "red","brown"), lty = 1, cex = 0.5,bty = "n")