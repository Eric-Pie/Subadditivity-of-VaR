### 0 Setup ####################################################################
#install.packages("simsalapar")
#install.packages("copula")
library(simsalapar)
require(copula)
n.obs <- 1000 #% sample size
alpha <- 0.95 
N.sim <- 100  #5000
doExtras <- simsalapar:::doExtras()

## list of variables
varList <- varlist(
  ## simulation number
  n.sim = list(type = "N", expr = quote(N[sim]), value = N.sim), #? No need?
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
  
  t(sapply(qmargin, function(FUN) ecdf(rowSums(FUN(U)))(d*FUN(alpha)))) 
}


nonGr <- get.nonGrids(varList)$nonGrids
#print(nonGr)

##n.sim times repulicate
dd_list <- rep(0,N.sim)
for (i in 1:N.sim){
  dd <- doOne(n= min(nonGr$n, 1000), d=4, family="t", tau=0.5,
              qmargin=nonGr$qmargin, alpha=nonGr$alpha)
  dd_list[i] <- as.numeric(dd)
}
sorted_dd_list <- sort(dd_list) ## order n.sim Y_b (i.e. dd)

## calculate the lower bond and the upper bound of one-tailed CI
Tilde_alpha <- 0.05
L_n_dd  <- sorted_dd_list[ceiling(Tilde_alpha*N.sim)]
U_n_dd <- sorted_dd_list[ceiling((1-Tilde_alpha)*N.sim)]

## to calculate the approximately probability of subadditive or superadditive
prob_sub <- function(alpha,n.sim){
  num_sub <- 0
  for (i in 1:n.sim){
    if (dd_list[i] >= alpha){
      num_sub <- num_sub + 1
    }
  }
  print(num_sub/n.sim)
}
prob_super <- function(alpha,n.sim){
  num_super <- 0
  for (i in 1:n.sim){
    if (dd_list[i] < alpha){
      num_super <- num_super + 1
    }
  }
  print(num_super/n.sim)
}


## to see if the VaR is subadditive or superadditive
judgement <- function(L_n_dd,U_n_dd,alpha,n.sim){
  if (alpha < L_n_dd){
    print("Approximately with probability 1−α˜, VaR is subadditive.")
    prob_sub(alpha,n.sim)
  }
  else if (alpha > U_n_dd){
    print("Approximately with probability 1−α˜, VaR is superadditive")
    prob_super(alpha,n.sim)
  }
}

judgement(L_n_dd = L_n_dd,U_n_dd = U_n_dd,alpha = alpha,n.sim = N.sim)

