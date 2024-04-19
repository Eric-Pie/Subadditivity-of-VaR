### 0 Setup ####################################################################
#install.packages("simsalapar")
#install.packages("copula")
#install.packages("xtable")
library(simsalapar)
library(xtable)
#require(copula)
n.obs <- 1000 #% sample size
alpha <- seq(0.6,0.7,length = 20)  #% alpha
doExtras <- simsalapar:::doExtras()


## list of variables
varList <- varlist(
  ## sample size
  n = list(value = n.obs),
  ## margins
  qmargin = list(type="frozen", expr = quote(F^{-1}), #% F[j] is changed to F^{-1}.
                 value = c(Ber_0.3 = function(p) qbinom(p,size = 1, prob = 0.3))), ## bernoulli quantile function with probability 0.3
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
dd_20 <- doOne(n= min(nonGr$n, 1000), d=20,
               qmargin=nonGr$qmargin, alpha=nonGr$alpha)
dd_100 <- doOne(n= min(nonGr$n, 1000), d=100,
                qmargin=nonGr$qmargin, alpha=nonGr$alpha)
#stopifnot(dim(dd) == with(nonGr, c(length(qmargin), length(alpha))))

DD <- matrix(data = NA, nrow = length(alpha), ncol = 3)
DD[,1] <- t(dd_5)
DD[,2] <- t(dd_20)
DD[,3] <- t(dd_100)

## the function defined to calculate the Lower bound and Upper bound of CI of F_s at significance level Tilde_alpha 
L_n <- function(dd,n,tilde_alpha){
  dd - qnorm(1-tilde_alpha)*sqrt(dd*(1-dd))/sqrt(n)
}
U_n<- function(dd,n,tilde_alpha){
  dd + qnorm(1-tilde_alpha)*sqrt(dd*(1-dd))/sqrt(n)
}

#to see if alpha < L_n or alpha > U_n.
judgement <- function(L_n,U_n,alpha){
  if (alpha < L_n){
    return(1)  ## 1 represents subadditive
  }
  else if (alpha > U_n){
    return(2)  ## 2 represents superadditive
  }
  else {
    return(0)  ## 0 represents both not 
  }
}

Tilde_alpha <- 0.05 ## significance level of CI

data_array <- c()## build up a array to store results

## compute all results
for (i in 1:3){
  for (j in length(alpha)){
    L_n_dd <- L_n(dd = DD[j,i],n=n.obs,tilde_alpha=Tilde_alpha)
    U_n_dd <- U_n(dd = DD[j,i],n=n.obs,tilde_alpha=Tilde_alpha)
    
    data_array <- append(data_array, judgement(L_n=L_n_dd,U_n=U_n_dd,alpha=alpha[j]))
  }
}

## build up a matrix to store the data
Matrix_DD <- matrix(data = data_array, nrow = length(alpha), ncol = 3)

## add the dim names
dimnames(Matrix_DD) <- list(alpha,c(5,20,100))

## transform matrix to table in Latex
df_DD <- data.frame(Matrix_DD)
Table_DD <- xtable(df_DD)
toLatex.xtable(Table_DD)