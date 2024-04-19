### 0 Setup ####################################################################
#install.packages("simsalapar")
#install.packages("copula")
#install.packages("xtable")
library(simsalapar)
library(xtable)
#require(copula)
n.obs <- 2000 #% sample size
## generate a sequence of alpha having more value close to 0
alpha <- 1-10^(seq(log10(1-0.001), log10(1-0.1), length.out = 400))  #% alpha
N.sim <- 5000 #% number of simulations
doExtras <- simsalapar:::doExtras()


## list of variables
varList <- varlist(
  ## sample size
  n = list(value = n.obs),
  ## margins
  qmargin = list(type="frozen", expr = quote(F^{-1}), #% F[j] is changed to F^{-1}.
                 value = c(Ber_0.5 = function(p) qbinom(p,size = 1, prob = 0.5))), ## bernoulli quantile function with probability 0.3
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

## N.sim simulation
## build matrices to store the results for d = 4,5,6,7
F_4 <- matrix(data = NA, nrow = length(alpha), ncol = N.sim)
F_5 <- matrix(data = NA, nrow = length(alpha), ncol = N.sim)
F_6 <- matrix(data = NA, nrow = length(alpha), ncol = N.sim)
F_7 <- matrix(data = NA, nrow = length(alpha), ncol = N.sim)




for (i in 1:N.sim){
  dd_4 <- doOne(n= min(nonGr$n, 1000), d=4,
                qmargin=nonGr$qmargin, alpha=nonGr$alpha)
  dd_5 <- doOne(n= min(nonGr$n, 1000), d=5,
                qmargin=nonGr$qmargin, alpha=nonGr$alpha)
  dd_6 <- doOne(n= min(nonGr$n, 1000), d=6,
                qmargin=nonGr$qmargin, alpha=nonGr$alpha)
  dd_7 <- doOne(n= min(nonGr$n, 1000), d=7,
                qmargin=nonGr$qmargin, alpha=nonGr$alpha)


  F_4[,i] <- t(dd_4)
  F_5[,i] <- t(dd_5)
  F_6[,i] <- t(dd_6)
  F_7[,i] <- t(dd_7)


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

Tilde_alpha <- 0.05 ## significance level 0.05
#print(qnorm(Tilde_alpha))

## calculate all t-statistics and store in the list
DD_4 <- matrix(0,nrow = length(alpha),ncol = N.sim)
DD_5 <- matrix(0,nrow = length(alpha),ncol = N.sim)
DD_6 <- matrix(0,nrow = length(alpha),ncol = N.sim)
DD_7 <- matrix(0,nrow = length(alpha),ncol = N.sim)




## calculate t-statistics for d = 4,5,6,7,8
for (i in 1:N.sim){
  for (j in 1:length(alpha)){
    DD_4[j,i] <- t_statistic(n.obs,F_4[j,i],alpha[j])
  }
}

for (i in 1:N.sim){
  for (j in 1:length(alpha)){
    DD_5[j,i] <- t_statistic(n.obs,F_5[j,i],alpha[j])
  }
}

for (i in 1:N.sim){
  for (j in 1:length(alpha)){
    DD_6[j,i] <- t_statistic(n.obs,F_6[j,i],alpha[j])
  }
}

for (i in 1:N.sim){
  for (j in 1:length(alpha)){
    DD_7[j,i] <- t_statistic(n.obs,F_7[j,i],alpha[j])
  }
}






## calculate the rejection rate for d = 4,5,6,7
reject_rate_4 <- matrix(0,nrow = length(alpha),ncol = 1)
for (i in 1:length(alpha)){
  reject_rate_4[i,1] <- mean(sapply(DD_4[i,],reject,Tilde_alpha))
}

reject_rate_5 <- matrix(0,nrow = length(alpha),ncol = 1)
for (i in 1:length(alpha)){
  reject_rate_5[i,1] <- mean(sapply(DD_5[i,],reject,Tilde_alpha))
}

reject_rate_6 <- matrix(0,nrow = length(alpha),ncol = 1)
for (i in 1:length(alpha)){
  reject_rate_6[i,1] <- mean(sapply(DD_6[i,],reject,Tilde_alpha))
}

reject_rate_7 <- matrix(0,nrow = length(alpha),ncol = 2)
for (i in 1:length(alpha)){
  reject_rate_7[i,1] <- mean(sapply(DD_7[i,],reject,Tilde_alpha))
  reject_rate_7[i,2] <- Tilde_alpha
}





## combine the rejection rate for d =4,5,6,7 into one matrix
reject_rate <- cbind(reject_rate_4,reject_rate_5,reject_rate_6,reject_rate_7)


## save the above plot as a pdf
pdf("Bernoulli_different_d_n.obs=2000.pdf")

## plot the rejection rate
matplot(log10(1-10^(seq(log10(1-0.001), log10(1-0.1), length.out = 400))), reject_rate, type = "l", lty = 1, col = c("black","blue","orange","purple","brown"), xlab = expression(log~alpha), ylab = substitute(Probability~of~Rejecting~Null~Hypothesis~'for'~d~'='~'4,5,6,7'~(n.obs~"="~2000)),xlim = c(-3,-1),ylim = c(0,1))

## add the legend
legend("left", legend = c('rejection probability for d = 4','rejection probability for d = 5','rejection probability for d = 6','rejection probability for d = 7','significance level = 0.05'), col = c("black","blue","orange","purple","brown"), lty = 1, cex = 0.6,bty = "n")

## add points on the cross point of the two lines with ptch = 19 and col = black
points(-1.235,0.05,pch = 19,col = "black",cex = 0.8)
points(-1.55,0.05,pch = 19,col = "black",cex = 0.8)
points(-1.845,0.05,pch = 19,col = "black",cex = 0.8)
points(-2.15,0.05,pch = 19,col = "black",cex = 0.8)


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

dev.off()



