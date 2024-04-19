### 0 Setup ####################################################################
#install.packages("simsalapar")
#install.packages("copula")
library(simsalapar)
library(xtable)
require(copula)
n.obs <- 1000 #% sample size
alpha <- 1-10^(seq(log10(1-0.8), log10(1-0.999), length.out = 100))
N.sim <- 5000  
doExtras <- simsalapar:::doExtras()

## list of variables
varList <- varlist(
  ## simulation number
  n.sim = list(type = "N", expr = quote(N[sim]), value = N.sim), 
  ## sample size
  n = list(value = n.obs),
  ## dimensions, and weights (vector) for each d
  d = list(type="grid", value = c(5,20,100)),
  ## copula family names
  family = list(type="frozen", expr = quote(C),
                value = "t"), # t = t_4
  ## dependencies by Kendall's tau
  tau = list(type="frozen", value = 0.5),
  ## margins
  qmargin = list(type="frozen", expr = quote(F[j]),
                 value = c(t4   = function(p) qt(p, df=4))), #function only requires one para. p (cum. prob.), since we have known the distribution is t with df = 4.
  ## VaR confidence levels
  alpha = list(type="frozen", value = alpha))

#getEl(varList)

## function defined to compute F_{X_1+..+X_d}(d*F_1^-(\alpha))
doOne <- function(n, d, family, tau, qmargin, alpha)
{
  #stopifnot(require(copula))
  cop <- switch(family,
                "t" =
                  ellipCopula("t", param=iTau(ellipCopula("t"), tau=tau), dim=d),
                stop("unsupported 'family'"))
  U <- rCopula(n, copula=cop)
  
  t(sapply(qmargin, function(FUN) ecdf(rowSums(FUN(U)))(d*FUN(alpha)))) 
}

## apply doOne to all samples
DD <- doMclapply(vList = varList,doOne = doOne,check = FALSE)

Tilde_alpha <- 0.05 ## significance level 0.05

## build up a matrix to store the data and Probability
Matrix_DD <- matrix(data = "NA", nrow = length(alpha), ncol = 3)
Matrix_Prob_sub <- matrix(data = "NA", nrow = length(alpha), ncol = 3)
Matrix_Prob_super <- matrix(data = "NA", nrow = length(alpha), ncol = 3)

## build up a matrix to store all L_n, U_n  and alpha with different alpha and d
Matrix_Bounds <- matrix(data = NA, nrow = length(alpha),ncol = 9)

## calculate the lower bond and the upper bound of one-tailed CI
for (i in 1:3){
  for (j in 1:length(alpha)){
    data_array <- c()
    for (m in 1:N.sim){
      data_array <- append(data_array,DD[i,][m][[1]]$value[j])
    }
    sorted_arr <- sort(data_array)
    L_n <- sorted_arr[ceiling(N.sim*Tilde_alpha)]
    U_n <- sorted_arr[ceiling(N.sim*(1-Tilde_alpha))]
    ## store the lower bounds and upper bounds
    Matrix_Bounds[j,3*i-2] <- L_n
    Matrix_Bounds[j,3*i-1] <- U_n
    Matrix_Bounds[j,3*i] <- alpha[j]
    
    if (alpha[j] < L_n){
      Matrix_DD[j,][i] <- 1  ## 1 represents subadditive
      num_sub <- 0
      for (x in 1:length(sorted_arr)){
        if (sorted_arr[x] >= alpha[j]){
          num_sub <- num_sub + 1
        }
      }
      Matrix_Prob_sub[j,][i] <- round(num_sub/N.sim, digits = 3)
    }
    else if (alpha[j] > U_n){
      Matrix_DD[j,][i] <- 2  ## 2 represents superadditive
      num_super <- 0
      for (x in 1:length(sorted_arr)){
        if (sorted_arr[x] < alpha[j]){
          num_super <- num_super + 1
        }
      }
      Matrix_Prob_super[j,][i] <- round(num_super/N.sim, digits = 3)
      Matrix_Prob_super[j,][i]
    }
    else {
      Matrix_DD[j,][i] <- 0  ## 0 represents both not
    }
  }
}

## add the dim names
dimnames(Matrix_DD) <- list(alpha,c(5,20,100))
dimnames(Matrix_Prob_sub) <- list(alpha,c(5,20,100))
dimnames(Matrix_Prob_super) <- list(alpha,c(5,20,100))

## split the Matrix_Bounds into 3 matrix by d
Matrix_5 <- matrix(data = Matrix_Bounds[,1:3], nrow = length(alpha), ncol = 3)
Matrix_20 <- matrix(data = Matrix_Bounds[,1:3], nrow = length(alpha), ncol = 3)
Matrix_100 <- matrix(data = Matrix_Bounds[,1:3], nrow = length(alpha), ncol = 3)

## transform matrix to table in Latex
df_DD <- data.frame(Matrix_DD)
df_Prob_sub <- data.frame(Matrix_Prob_sub)
df_Prob_super <- data.frame(Matrix_Prob_super)

Table_DD <- xtable(df_DD)
Table_Prob_sub <- xtable(df_Prob_sub)
Table_Prob_super <- xtable(df_Prob_super)

toLatex.xtable(Table_DD)
toLatex.xtable(Table_Prob_sub)
toLatex.xtable(Table_Prob_super)

## add the dim names
dimnames(Matrix_5) <- list(alpha,c('Lower bound','Upper bound',expression(alpha)))
dimnames(Matrix_20) <- list(alpha,c('Lower bound','Upper bound',expression(alpha)))
dimnames(Matrix_100) <- list(alpha,c('Lower bound','Upper bound',expression(alpha)))

## draw a plot for d = 5,
matplot((1-10^(seq(log10(1-0.8), log10(1-0.999), length.out = 100))), Matrix_5, type = "l", lty = 1, col = c("blue", "red","brown"), xlab = expression(alpha), ylab = substitute(Simulation~CI~of~F[S](d*{F^{-1}}(alpha))~"for"~l~alpha~"in [0.8,1]"~(d~"="~5), list(l = length(alpha))),xlim = c(0.8,1),ylim = c(0.8,1))

## add the legend
legend("bottomright", legend = c('Lower bound (if '~alpha~'< Lower bound, VaR is subadditive)','Upper bound (if '~alpha~'> Upper bound, VaR is superadditive) ',expression(alpha)), col = c("blue", "red","brown"), lty = 1, cex = 0.5,bty = "n")

## draw a plot for d = 20,
matplot((1-10^(seq(log10(1-0.8), log10(1-0.999), length.out = 100))), Matrix_5, type = "l", lty = 1, col = c("blue", "red","brown"), xlab = expression(alpha), ylab = substitute(Simulation~CI~of~F[S](d*{F^{-1}}(alpha))~"for"~l~alpha~"in [0.8,1]"~(d~"="~20), list(l = length(alpha))),xlim = c(0.8,1),ylim = c(0.8,1))

## add the legend
legend("bottomright", legend = c('Lower bound (if '~alpha~'< Lower bound, VaR is subadditive)','Upper bound (if '~alpha~'> Upper bound, VaR is superadditive) ',expression(alpha)), col = c("blue", "red","brown"), lty = 1, cex = 0.5,bty = "n")

## draw a plot for d = 100,
matplot((1-10^(seq(log10(1-0.8), log10(1-0.999), length.out = 100))), Matrix_5, type = "l", lty = 1, col = c("blue", "red","brown"), xlab = expression(alpha), ylab = substitute(Simulation~CI~of~F[S](d*{F^{-1}}(alpha))~"for"~l~alpha~"in [0.8,1]"~(d~"="~100), list(l = length(alpha))),xlim = c(0.8,1),ylim = c(0.8,1))

## add the legend
legend("bottomright", legend = c('Lower bound (if '~alpha~'< Lower bound, VaR is subadditive)','Upper bound (if '~alpha~'> Upper bound, VaR is superadditive) ',expression(alpha)), col = c("blue", "red","brown"), lty = 1, cex = 0.5,bty = "n")


## now here is a problem. When alpha is close to 0.999, let's say 1, the VaR is not is both type, it is weird.