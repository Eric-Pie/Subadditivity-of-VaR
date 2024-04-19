### 0 Setup ####################################################################
#install.packages("simsalapar")
#install.packages("copula")
library(simsalapar)
require(copula)
n.obs <- 10000 #% sample size
alpha <- 1-10^(seq(log10(1-0.9), log10(1-0.999), length.out = 201))
N.sim <- 10 #5000
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

Tilde_alpha <- 0.05 ##significance level
##n.sim times replicate
dd_list <- rep(list(NULL),length(alpha))
sorted_dd_list <- rep(list(NULL),length(alpha))

## calculate the Lower bound and Upper bound with different alpha of CI of F_s at significance level Tilde_alpha 
L_n_dd <- rep(0,length(alpha))
U_n_dd <- rep(0,length(alpha))
for (j in 1:length(alpha)){
  for (i in 1:N.sim){
    dd <- doOne(n= min(nonGr$n, 1000), d=4, family="t", tau=0.5,
              qmargin=nonGr$qmargin, alpha=nonGr$alpha[j])
    dd_list[[j]][i] <- as.numeric(dd)
  }
  sorted_dd_list[[j]] <- sort(dd_list[[j]])## order n.sim Y_b (i.e. dd)
  
  ## calculate the lower bond and the upper bound of one-tailed CI at significance level Tilde_alpha
  L_n_dd[j] <- sorted_dd_list[[j]][ceiling(Tilde_alpha*N.sim)]
  U_n_dd[j] <- sorted_dd_list[[j]][ceiling((1-Tilde_alpha)*N.sim)]
}


## function defined to calculate the approximately probability of subadditive or superadditive
prob_sub <- function(alpha,n.sim,dd_List){
  num_sub <- 0
  for (x in 1:n.sim){
    if (dd_List[x] >= alpha){
      num_sub <- num_sub + 1
    }
  }
  return(num_sub/n.sim)
}
prob_super <- function(alpha,n.sim,dd_List){
  num_super <- 0
  for (x in 1:n.sim){
    if (dd_List[x] < alpha){
      num_super <- num_super + 1
    }
  }
  return(num_super/n.sim)
}

## to see if the VaR is subadditive or superadditive  
judgement <- function(L_n_dd,U_n_dd,alpha,n.sim){
  num_sub <- 0
  num_super <- 0
  prob_sub_list <- rep(0,length(alpha))
  prob_super_list <- rep(0,length(alpha))
  for (m in 1:length(alpha)){
    if (alpha[m] < L_n_dd[m]){
      num_sub <= num_sub + 1 
      prob_sub_list[m] <- prob_sub(alpha[m],n.sim,dd_list[[m]])
    }
    else if (alpha[m] > U_n_dd[m]){
      num_super <- num_super + 1
      prob_super_list[m] <- prob_super(alpha[m],n.sim,dd_list[[m]])
    }
  }
  print("prob_sub_list")
  print(prob_sub_list)
  print("prob_super_list")
  print(prob_super_list)
}

judgement(L_n_dd = L_n_dd,U_n_dd = U_n_dd,alpha = alpha,n.sim = N.sim)

