### 0 Setup ####################################################################
#install.packages("simsalapar")
#install.packages("copula")
library(simsalapar)
require(copula)
n.obs <- 10000 #% sample size
alpha <- 0.99
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


nonGr <- get.nonGrids(varList)$nonGrids
#print(nonGr)
dd <- doOne(n= min(nonGr$n, 1000000), d=4, family="t", tau=0.5,
            qmargin=nonGr$qmargin, alpha=nonGr$alpha)
print(dd)

L_n <- function(dd,n,tilde_alpha){
  dd - qnorm(1-tilde_alpha)*sqrt(dd*(1-dd))/sqrt(n)
}
U_n<- function(dd,n,tilde_alpha){
  dd + qnorm(1-tilde_alpha)*sqrt(dd*(1-dd))/sqrt(n)
}

#to see if alpha < L_n or alpha > U_n.
judgement <- function(L_n,U_n,alpha){
  if (alpha < L_n){
    print("with probability 1−α˜, VaR is subadditive.")
  }
  else if (alpha > U_n){
    print("with probability 1−α˜, VaR is superadditive")
  }
}

Tilde_alpha <- 0.05
L_n <- L_n(dd =dd, n =n.obs,tilde_alpha =Tilde_alpha)
U_n <- U_n(dd =dd, n =n.obs,tilde_alpha =Tilde_alpha)
#print(L_n)
#print(U_n)

#To see if n is sufficient large
n_star <- function(tilde_alpha,dd,alpha){
  (qnorm(1-tilde_alpha))^2 *(dd*(1-dd))/(dd-alpha)^2
}
N_star <- n_star(tilde_alpha=Tilde_alpha,dd = dd, alpha = alpha)
print(N_star)


judgement(L_n = L_n,U_n=U_n,alpha=alpha)