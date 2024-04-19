### 0 Setup ####################################################################
#install.packages("simsalapar")
#install.packages("copula")
library(simsalapar)
require(copula)
n.obs <- 1000 #% sample size
alpha <- 0.9  #% alpha
doExtras <- simsalapar:::doExtras()

## list of variables
varList <- varlist(
  ## sample size
  n = list(type = "grid",value = n.obs),
  ## dimensions, and weights (vector) for each d
  d = list(type="grid", value = c(5,20,100)),
  ## copula family names
  family = list(type="grid", expr = quote(C),
                value = c("t")), # t = t_4
  ## dependencies by Kendall's tau
  tau = list(type="grid", value = 0.5),
  ## margins
  qmargin = list(type="grid", expr = quote(F^{-1}), #% F[j] is changed to F^{-1}.
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
  
  ## compute F_{X_1+..+X_d}(d*F_1^-(\alpha)) for all confidence levels alpha
  ## => VaR_alpha superadditive <=> F_{X_1+..+X_d}(d*F_1^-(\alpha)) - alpha < 0
  t(sapply(qmargin, function(FUN) ecdf(rowSums(FUN(U)))(d*FUN(alpha)))) 
  ## note: t() is important here, since, otherwise, the order of the variables
  ## ----  would not be correct (=> check should reveal this) 
}
##doMCLapply
doCheck(doOne, varList)