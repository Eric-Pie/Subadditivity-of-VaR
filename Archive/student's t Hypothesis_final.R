### 0 Setup ####################################################################
#install.packages("simsalapar")
#install.packages("copula")
#install.packages("xtable")
library(simsalapar)
library(xtable)
require(copula)
n.obs <- 1000 #% sample size
alpha <- 1-10^(seq(log10(1-0.9), log10(1-0.999), length.out = 20))  #% alpha
doExtras <- simsalapar:::doExtras()

## list of variables
varList <- varlist(
  ## sample size
  n = list(value = n.obs),
  ## dimensions, and weights (vector) for each d
  d = list(type="grid", value = c(5,20,100)),
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

data_array <- c() ## build up a array to store results

## compute all results 
for (i in 1:3){
  for (j in 1:length(alpha)){
    T_statistic <- t_statistic(n = n.obs,dd=DD[[i]]$value[,j],alpha=alpha[j])
    print(T_statistic)
    data_array <- append(data_array, reject(T_stat = T_statistic,tilde_alpha = Tilde_alpha))
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