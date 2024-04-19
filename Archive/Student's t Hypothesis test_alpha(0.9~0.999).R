### 0 Setup ####################################################################
#install.packages("simsalapar")
#install.packages("copula")
library(simsalapar)
require(copula) # for the following call of doOne()
n.obs <- 1000 ## sample size
alpha <- 1-10^(seq(log10(1-0.9), log10(1-0.999), length.out = 201)) ## alphas from 0.9 tp 0.999
doExtras <- simsalapar:::doExtras()

## list of variables
varList <- varlist(
  ## sample size
  n = list(value = n.obs),
  ## dimensions, and weights (vector) for each d
  d = list(type="frozen", value = 20),
  ## copula family names
  family = list(type="frozen", expr = quote(C),
                value = "t"), # t = t_4
  ## dependencies by Kendall's tau
  tau = list(type="frozen", value = 0.5),
  ## margins
  qmargin = list(type="frozen", expr = quote(F[j]),
                 value = c(t4   = function(p) qt(p, df=4))), #function only requires one para. p (cum. prob.), since we have known the distribution is t with df = 4.
  ## VaR confidence levels
  alpha = list(type="inner", value = c(alpha)))

#getEl(varList)

## Function to Compute F_{X_1+..+X_d}(d*F_1^-(\alpha))
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
}


nonGr <- get.nonGrids(varList)$nonGrids
#print(nonGr)

## apply function doOne, store the results into dd
dd <- doOne(n= min(nonGr$n, 1000), d=4, family="t", tau=0.5,
            qmargin=nonGr$qmargin, alpha=nonGr$alpha)
print(dd)

## function defined to calculate t_statistic
t_statistic <- function(n,dd,alpha)
{
  sqrt(n)*(dd - alpha)/sqrt(dd*(1-dd))
}

## calculate t_statistic
T_statistic <- t_statistic(n = n.obs,dd=dd,alpha=alpha)
print(T_statistic)

## define the function used to judge if we should reject Null Hypothesis or not
reject <- function(T_stats,tilde_alpha)
{
  num_sub <- 0
  num_super <- 0
  for (T_stat in T_stats){
    ## if (T_stat >= qnorm(tilde_alpha)), then we will not reject the Null Hypothesis and VaR is subadditive.
    if (T_stat >= qnorm(tilde_alpha)){
      num_sub <- num_sub + 1
    }
    ## if (T_stat < qnorm(tilde_alpha)), then we will reject the Null Hypothesis and VaR is superadditive.
    else if (T_stat < qnorm(tilde_alpha)){
      num_super <- num_super + 1
  }
  }
  ## show the number of subadditive and superadditive VaR among these alpha
  print("num_super =")
  print(num_super)
  print("num_sub =")
  print(num_sub)
}

## significance level of hypothesis test
Tilde_alpha <- 0.05
## to judge if we reject the Null at a significance level Tilde_alpha
reject(T_stats = T_statistic,tilde_alpha = Tilde_alpha)