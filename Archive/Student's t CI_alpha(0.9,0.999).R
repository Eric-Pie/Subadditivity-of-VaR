### 0 Setup ####################################################################
#install.packages("simsalapar")
#install.packages("copula")
library(simsalapar)
require(copula) # for the following call of doOne()
n.obs <- 10000 #% sample size
alpha <- 1-10^(seq(log10(1-0.9), log10(1-0.999), length.out = 20)) ## alphas from 0.9 tp 0.999
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
  alpha = list(type="inner", value = alpha))

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
dd <- doOne(n= min(nonGr$n, 10000), d=4, family="t", tau=0.5,
            qmargin=nonGr$qmargin, alpha=nonGr$alpha)
print(dd)

## the function defined to calculate the Lower bound and Upper bound of CI of F_s at significance level Tilde_alpha 
L_n <- function(dd,n,tilde_alpha){
  dd - qnorm(1-tilde_alpha)*sqrt(dd*(1-dd))/sqrt(n)
}
U_n<- function(dd,n,tilde_alpha){
  dd + qnorm(1-tilde_alpha)*sqrt(dd*(1-dd))/sqrt(n)
}

## the function defined to see if the VaR is subadditive or superadditive
judgement <- function(L_n,U_n,alpha){
  num_sub <- 0
  num_super <- 0
  for (i in 1:length(alpha)){
    if (alpha[i] < L_n[i]){
      num_sub <- num_sub + 1
    }
    else if (alpha[i] > U_n[i]){
      num_super <- num_super + 1
    }
  }
  print("num_super =")
  print(num_super)
  print("num_sub =")
  print(num_sub)
}

Tilde_alpha <- 0.05 ## significance level of CI
L_n <- L_n(dd =dd, n =n.obs,tilde_alpha =Tilde_alpha)
U_n <- U_n(dd =dd, n =n.obs,tilde_alpha =Tilde_alpha)
print(L_n)
print(U_n)
judgement(L_n = L_n,U_n=U_n,alpha=alpha)