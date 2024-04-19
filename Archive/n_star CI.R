### 0 Setup ####################################################################
#install.packages("simsalapar")
#install.packages("copula")
n_Times <- 200 #the number of samples for N_star
n_res <- c()
library(simsalapar)
require(copula)
lambda <- 1082.22 ## r =0.05, p =0.9
Tilde_alpha <- 0.05 ## significance level of F_s CI

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

n_star <- function(tilde_alpha,dd,alpha){
  (qnorm(1-tilde_alpha))^2 *(dd*(1-dd))/(dd-alpha)^2
}

## generate n_Times N_star with a fixed n.obs
for (i in 1:n_Times){
  n.obs <- 500 #% sample size
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
  
  
  nonGr <- get.nonGrids(varList)$nonGrids
  dd <- doOne(n= min(nonGr$n, 100000), d=4, family="t", tau=0.5,
              qmargin=nonGr$qmargin, alpha=nonGr$alpha)
  
  N_star <- n_star(tilde_alpha=Tilde_alpha,dd = dd, alpha = alpha)
  n_res <- append(n_res,N_star)
}

print(n_res)
sorted_N_star <- sort(n_res)
upper_N_star <- sorted_N_star[floor(0.95*length(n_res))]
lower_N_star <- sorted_N_star[ceiling(0.05*length(n_res))]
upper_N_star
lower_N_star

