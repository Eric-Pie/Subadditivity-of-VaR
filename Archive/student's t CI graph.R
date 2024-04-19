### 0 Setup ####################################################################
#install.packages("simsalapar")
#install.packages("copula")
library(simsalapar)
library(xtable)
require(copula)
n.obs <- 1000 #% sample size
alpha <- 1-10^(seq(log10(1-0.501), log10(1-0.999), length.out = 100))  #% alpha
N.sim <- 5000 #% number of simulations
doExtras <- simsalapar:::doExtras()

## list of variables
varList <- varlist(
  ## simulation number
  n.sim = list(type = "N", expr = quote(N[sim]), value = N.sim),
  ## sample size
  n = list(value = n.obs),
  ## dimensions, and weights (vector) for each d
  d = list(type="grid", value = c(5,10,20)),
  ## copula family names
  family = list(type="frozen", expr = quote(C),
                value = c("t")), # t = t_4
  ## dependencies by Kendall's tau
  tau = list(type="frozen", value = c(0.5)),
  ## margins
  qmargin = list(type="frozen", expr = quote(F^{-1}), #% F[j] is changed to F^{-1}.
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
  
  ## compute F_{X_1+..+X_d}(d*F_1^-(\alpha)) for all confidence levels alpha
  ## => VaR_alpha superadditive <=> F_{X_1+..+X_d}(d*F_1^-(\alpha)) - alpha < 0
  t(sapply(qmargin, function(FUN) ecdf(rowSums(FUN(U)))(d*FUN(alpha)))) 
  ## note: t() is important here, since, otherwise, the order of the variables
  ## ----  would not be correct (=> check should reveal this) 
}

## apply doOne to all samples in varList
DD <- doMclapply(vList = varList,doOne = doOne,check = FALSE)

## the function defined to calculate the Lower bound and Upper bound of CI of F_s at significance level Tilde_alpha 
L_n <- function(dd,n,tilde_alpha){
  dd - qnorm(1-tilde_alpha)*sqrt(dd*(1-dd))/sqrt(n)
}
U_n<- function(dd,n,tilde_alpha){
  dd + qnorm(1-tilde_alpha)*sqrt(dd*(1-dd))/sqrt(n)
}

Tilde_alpha <- 0.05 ## significance level 0.05

## calculate all lower bounds and store in the list
DD_5_L <- matrix(0,nrow = length(alpha),ncol = N.sim)
DD_10_L <- matrix(0,nrow = length(alpha),ncol = N.sim)
DD_20_L <- matrix(0,nrow = length(alpha),ncol = N.sim)

## calculate all upper bounds and store in the list
DD_5_U <- matrix(0,nrow = length(alpha),ncol = N.sim)
DD_10_U <- matrix(0,nrow = length(alpha),ncol = N.sim)
DD_20_U <- matrix(0,nrow = length(alpha),ncol = N.sim)

## compute all bounds for d = 5
for (i in 1:N.sim){
  for (j in 1:length(alpha)){
    L_n_dd <- L_n(dd = DD[1,][i][[1]]$value[j],n=n.obs,tilde_alpha=Tilde_alpha)
    U_n_dd <- U_n(dd = DD[1,][i][[1]]$value[j],n=n.obs,tilde_alpha=Tilde_alpha)
    
    DD_5_L[j,i] <- L_n_dd
    DD_5_U[j,i] <- U_n_dd
  }
}

## compute all bounds for d = 10
for (i in 1:N.sim){
  for (j in 1:length(alpha)){
    L_n_dd <- L_n(dd = DD[2,][i][[1]]$value[j],n=n.obs,tilde_alpha=Tilde_alpha)
    U_n_dd <- U_n(dd = DD[2,][i][[1]]$value[j],n=n.obs,tilde_alpha=Tilde_alpha)
    
    DD_10_L[j,i] <- L_n_dd
    DD_10_U[j,i] <- U_n_dd
  }
}

## compute all bounds for d = 20
for (i in 1:N.sim){
  for (j in 1:length(alpha)){
    L_n_dd <- L_n(dd = DD[3,][i][[1]]$value[j],n=n.obs,tilde_alpha=Tilde_alpha)
    U_n_dd <- U_n(dd = DD[3,][i][[1]]$value[j],n=n.obs,tilde_alpha=Tilde_alpha)
    
    DD_20_L[j,i] <- L_n_dd
    DD_20_U[j,i] <- U_n_dd
  }
}

#to see if alpha < L_n or alpha > U_n.
judgement_L <- function(L_n,alpha){
  if (alpha <= L_n){
    return(1)  ## 1 represents subadditive
  }
  else {
    return(0)  ## 0 represents both not 
  }
}

judgement_U <- function(U_n,alpha){
  if (alpha >= U_n){
    return(1)  ## 1 represents superadditive
  }
  else {
    return(0)  ## 0 represents both not 
  }
}

## calculate the probability in CI for d = 5
Prob_5 <- matrix(0,nrow = length(alpha),ncol = 3)
for (i in 1:length(alpha)){
  Prob_5[i,1] <- mean(sapply(DD_5_L[i,],judgement_L,alpha = alpha[i]))
  Prob_5[i,2] <- mean(sapply(DD_5_U[i,],judgement_U,alpha = alpha[i]))
  Prob_5[i,3] <- 1-Tilde_alpha
}

## calculate the probability in CI for d = 10
Prob_10 <- matrix(0,nrow = length(alpha),ncol = 3)
for (i in 1:length(alpha)){
  Prob_10[i,1] <- mean(sapply(DD_10_L[i,],judgement_L,alpha = alpha[i]))
  Prob_10[i,2] <- mean(sapply(DD_10_U[i,],judgement_U,alpha = alpha[i]))
  Prob_10[i,3] <- 1-Tilde_alpha
}

## calculate the probability in CI for d = 20
Prob_20 <- matrix(0,nrow = length(alpha),ncol = 3)
for (i in 1:length(alpha)){
  Prob_20[i,1] <- mean(sapply(DD_20_L[i,],judgement_L,alpha = alpha[i]))
  Prob_20[i,2] <- mean(sapply(DD_20_U[i,],judgement_U,alpha = alpha[i]))
  Prob_20[i,3] <- 1-Tilde_alpha
}


## draw a plot for d = 5,
matplot((1-10^(seq(log10(1-0.501), log10(1-0.999), length.out = 100))), Prob_5, type = "l", lty = 1, col = c("blue", "red","brown"), xlab = expression(alpha), ylab = substitute("Probability~of~alpha~not~in~CI~for"~l~alpha~"in [0.5,1]"~(d~"="~5), list(l = length(alpha))),xlim = c(0.5,1),ylim = c(0,1))

## add the legend
legend("bottomright", legend = c('Probability of alpha is below lower bound','Probability of alpha is larger than upper bound','1-significance level = 0.95'), col = c("blue", "red","brown"), lty = 1, cex = 0.5,bty = "n")

## draw a plot for d = 10,
matplot((1-10^(seq(log10(1-0.501), log10(1-0.999), length.out = 100))), Prob_10, type = "l", lty = 1, col = c("blue", "red","brown"), xlab = expression(alpha), ylab = substitute("Probability~of~alpha~not~in~CI~for"~l~alpha~"in [0.5,1]"~(d~"="~10), list(l = length(alpha))),xlim = c(0.5,1),ylim = c(0,1))
  
## add the legend
legend("bottomright", legend = c('Probability of alpha is below lower bound','Probability of alpha is larger than upper bound','1-significance level = 0.95'), col = c("blue", "red","brown"), lty = 1, cex = 0.5,bty = "n")

## draw a plot for d = 20,
matplot((1-10^(seq(log10(1-0.501), log10(1-0.999), length.out = 100))), Prob_20, type = "l", lty = 1, col = c("blue", "red","brown"), xlab = expression(alpha), ylab = substitute("Probability~of~alpha~not~in~CI~for"~l~alpha~"in [0.5,1]"~(d~"="~20), list(l = length(alpha))),xlim = c(0.5,1),ylim = c(0,1))

## add the legend
legend("bottomright", legend = c('Probability of alpha is below lower bound','Probability of alpha is larger than upper bound','1-significance level = 0.95'), col = c("blue", "red","brown"), lty = 1, cex = 0.5,bty = "n")
