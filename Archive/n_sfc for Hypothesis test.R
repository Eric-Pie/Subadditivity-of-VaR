### 0 Setup ####################################################################
#install.packages("simsalapar")
#install.packages("copula")
library(simsalapar)
require(copula)
n.obs <- 1000#% sample size
alpha <- 0.9  #% alpha
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
  qmargin = list(type="frozen", expr = quote(F^{-1}), #% F[j] is changed to F^{-1}.
                 value = c(t4   = function(p) qt(p, df=4))), #function only requires one para. p (cum. prob.), since we have known the distribution is t with df = 4.
  ## VaR confidence levels
  alpha = list(type="frozen", value = alpha))

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
## get all S_i
doOne_S_i <- function(n, d, family, tau, qmargin, alpha)
{
  #stopifnot(require(copula))
  cop <- switch(family,
                "t" =
                  ellipCopula("t", param=iTau(ellipCopula("t"), tau=tau), dim=d),
                stop("unsupported 'family'"))
  U <- rCopula(n, copula=cop)
  
  ## compute F_{X_1+..+X_d}(d*F_1^-(\alpha)) for all confidence levels alpha
  ## => VaR_alpha superadditive <=> F_{X_1+..+X_d}(d*F_1^-(\alpha)) - alpha < 0
  t(sapply(qmargin, function(FUN) rowSums(FUN(U)))) 
  ## note: t() is important here, since, otherwise, the order of the variables
  ## ----  would not be correct (=> check should reveal this) 
}

##get s
doOne_s <- function(d, qmargin, alpha)
{
  alphas <- matrix(alpha,nrow = 1, ncol = d)
  t(sapply(qmargin, function(FUN) rowSums(FUN(alphas)))) 
  ## note: t() is important here, since, otherwise, the order of the variables
  ## ----  would not be correct (=> check should reveal this) 
}

nonGr <- get.nonGrids(varList)$nonGrids
#print(nonGr)
dd <- doOne(n= min(nonGr$n, 1000), d=4, family="t", tau=0.5,
            qmargin=nonGr$qmargin, alpha=nonGr$alpha)
dd_S_i <- doOne_S_i(n= min(nonGr$n, 1000), d=4, family="t", tau=0.5,
            qmargin=nonGr$qmargin, alpha=nonGr$alpha)
s <- doOne_s(d=4, qmargin=nonGr$qmargin, alpha=nonGr$alpha)
dd
dd_S_i
s

##calculate the sd. of samples
less_equal <- 0
for (i in 1:length(dd_S_i)){
  if (dd_S_i[i] <= s){
    less_equal <- less_equal+1
  }
}

values <- c(rep(1, times = less_equal), rep(0, times = length(dd_S_i)-less_equal))
samples <- array(values, dim = c(1, length(dd_S_i)))

n_mean <- mean(samples)
n_sd <- sd(samples)
n_sfc <- 1082.22*(n_sd/n_mean)^2
n_sfc
