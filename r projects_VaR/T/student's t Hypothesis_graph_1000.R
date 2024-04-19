### 0 Setup ####################################################################
#install.packages("simsalapar")
#install.packages("copula")
#install.packages("xtable")
library(simsalapar)
library(xtable)
require(copula)
n.obs <- 1000 #% sample size
alpha <- 1-10^(seq(log10(1-0.001), log10(1-0.999), length.out = 200))  #% alpha
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

Tilde_alpha <- 0.05 ## significance level 0.05
#print(qnorm(Tilde_alpha))

## calculate all t-statistics and store in the list
DD_5 <- matrix(0,nrow = length(alpha),ncol = N.sim)
DD_10 <- matrix(0,nrow = length(alpha),ncol = N.sim)
DD_20 <- matrix(0,nrow = length(alpha),ncol = N.sim)

## calculate t-statistics for d = 5
for (i in 1:N.sim){
  for (j in 1:length(alpha)){
    DD_5[j,i] <- t_statistic(n.obs,DD[1,][i][[1]]$value[j],alpha[j])
  }
}

## calculate t-statistics for d = 10
for (i in 1:N.sim){
  for (j in 1:length(alpha)){
    DD_10[j,i] <- t_statistic(n.obs,DD[2,][i][[1]]$value[j],alpha[j])
  }
}

## calculate t-statistics for d = 20
for (i in 1:N.sim){
  for (j in 1:length(alpha)){
    DD_20[j,i] <- t_statistic(n.obs,DD[3,][i][[1]]$value[j],alpha[j])
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

## combine the rejection rate for d =5, d = 10 and d = 20 into one matrix
reject_rate <- cbind(reject_rate_5,reject_rate_10,reject_rate_20)

## transform the plot into a pdf file
pdf("t Hypothesis_graph_n.obs=1000.pdf")
matplot((1-10^(seq(log10(1-0.001), log10(1-0.999), length.out = 200))), reject_rate, type = "l", lty = 1, col = c("blue", "red","green","brown"), xlab = expression(alpha), ylab = substitute(Probability~of~Rejecting~Null~Hypothesis~'for'~d~'='~'5,10,20'~(n.obs~"="~1000)),xlim = c(0,1),ylim = c(0,1))
legend("right", legend = c('rejection probability for d = 5','rejection probability for d = 10','rejection probability for d = 20','significance level = 0.05'), col = c("blue", "red","green","brown"), lty = 1, cex = 1,bty = "n")
abline(v = 0.5, lty = 2)
axis(1, at = 0.5, labels = 0.5)
points(0.5,0.05,pch = 19,col = "black",cex = 1)
dev.off()



