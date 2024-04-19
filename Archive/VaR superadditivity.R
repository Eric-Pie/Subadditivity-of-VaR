## By Marius Hofert

## Simulation study (demo from R package simsalapar) to determine superadditivity of VaR


### 0 Setup ###############################################################
##install.packages("foreach")
##install.packages("doParallel")
##install.packages("crop")
n.obs <- 100 #% sample size
n.alpha <- 32 #% number of quantile. Dividing the samples into n.alpha equal intervals.
if(n.alpha > n.obs)
  stop("number of quantiles, n.alpha, must be smaller than the sample size, n.obs")
doExtras <- simsalapar:::doExtras()
#%调用simsalapar包中的doExtras()函数，并将其结果赋值给名为doExtras的变量（或对象）doExtras

## load packages
#install.packages("crop")
library(simsalapar)
doPDF <- require(crop) 

## list of variables
varList <-
  varlist(
    ## sample size
    n = list(value = n.obs),
    ## dimensions, and weights (vector) for each d
    d = list(type="grid", value = c(4, 20, 100)),
    ## copula family names
    family = list(type="grid", expr = quote(C),
                  value = c("normal", "t", "Clayton", "Gumbel")), # t = t_4
    ## dependencies by Kendall's tau
    tau = list(type="grid", value = c(0.2, 0.5, 0.8)),
    ## margins
    qmargin = list(type="inner", expr = quote(F[j]),
                   value = c(norm = qnorm,
                             t4   = function(p) qt(p, df=4), #function only requires one para. p (cum. prob.), since we have known the distribution is t with df = 4.
                             Par2 = function(p) (1-p)^(-1/2))), # Pareto(2)
    ## VaR confidence levels
    alpha = list(type="inner", value = 0:n.alpha/n.alpha)) #% alpha is inputted into type="inner" with value...
#% The reason why we define value for "family", "tau" and "qmargin" is because they are our future choices, right?

### 1 Functions ################################################################

##' @title Function to Compute F_{X_1+..+X_d}(d*F_1^-(\alpha))
##' @author Marius Hofert
doOne <- function(n, d, family, tau, qmargin, alpha)
{
  ## checks (and load required packages here for parallel computing later on)
  ## stopifnot(require(copula))
  ## Note: This require() affects run time due to loading copula once on each
  ##       node, but this can easily be identified by *robust* analysis of all
  ##       run times (or is negligible due to much larger run times of the jobs)
  ##       or, as another alternative, can be avoided by using initExpr (see below)
  
  cop <- switch(family,
                "normal" =
                  ellipCopula("normal", param=iTau(ellipCopula("normal"), tau=tau),
                              dim=d),
                "t" =
                  ellipCopula("t", param=iTau(ellipCopula("t"), tau=tau), dim=d),
                "Clayton" =
                  onacopulaL("Clayton", list(th=iTau(archmCopula("clayton"), tau),
                                             seq_len(d))),
                "Gumbel" =
                  onacopulaL("Gumbel", list(th=iTau(archmCopula("gumbel"), tau),
                                            seq_len(d))),
                stop("unsupported 'family'"))
  U <- rCopula(n, copula=cop)
  #% The copula is used to split the dependence and marginal dist. from joint dist. If I use "t" here, does it represent those r.v. are iid.?
  
  ## compute F_{X_1+..+X_d}(d*F_1^-(\alpha)) for all confidence levels alpha
  ## => VaR_alpha superadditive <=> F_{X_1+..+X_d}(d*F_1^-(\alpha)) - alpha < 0
  t(sapply(qmargin, function(FUN) ecdf(rowSums(FUN(U)))(d*FUN(alpha)) - alpha)) #% What is FUN, what does d*FUN(alpha) represent?
  ## note: t() is important here, since, otherwise, the order of the variables
  ## ----  would not be correct (=> check should reveal this) 
}
  #% do we need to mention this in paper?
  #% t() represents transpose


### 2 Main #####################################################################

## check doOne() #% to check what aspects?

## manually
require(copula) # for the following call of doOne()
nonGr <- get.nonGrids(varList)$nonGrids
dd <- doOne(n= min(nonGr$n, 100), d=4, family="Clayton", tau=0.5,
            qmargin=nonGr$qmargin, alpha=nonGr$alpha)
stopifnot(dim(dd) == with(nonGr, c(length(qmargin), length(alpha))))

## with doCheck()
doCheck(doOne, varList)

## if simsalapar no longer *depends* on parallel:
library(parallel) #%
makeCluster <- parallel::makeCluster
## to be CRAN check compatible and as called from tests with explicit doExtras <- FALSE:
nc <- simsalapar:::nCores4test(); nc <- if(doExtras) nc else min(nc, 2) #% number of cores
nc
nc.win <- if(.Platform$OS.type=="windows") 1 else nc # otherwise win-builder fails

## computation
system.time(res <-
              doClusterApply(varList, cluster=makeCluster(nc, type="PSOCK"),
                             sfile = if(n.obs > 1000) "VaR_superadd.rds" else NULL,
                             doOne=doOne, monitor = printInfo[["gfile"]],
                             ## load copula once on each worker, to not affect run time
                             initExpr = require("copula"),
                             timer=mkTimer(gcFirst=TRUE) ))
#% doMCapply()
#% subjob(doOne)
if(doExtras) {
  ## doMclapply() ------------------------------------------------------------
  ## "init expression" for mclapply can happen here locally (see above)
  print(system.time(res2 <- doMclapply(varList, cores=nc.win, doOne=doOne)))
  stopifnot( doRes.equal(res, res2) )
  
  ## doRmpi() ----------------------------------------------------------------
  ## => passing an init expression does not (easily) work
  
  ## doForeach() -------------------------------------------------------------
  print(system.time(res3 <- {
    doForeach(varList, cluster=makeCluster(nc, type="PSOCK"),
              doOne=doOne, extraPkgs = "copula",
              timer=mkTimer(gcFirst=TRUE) )
  }))
  stopifnot( doRes.equal(res, res3) )
}


### 3 Analysis #################################################################

## extract results
val  <- getArray(res) # array of values
err  <- getArray(res, "error") # array of error indicators
warn <- getArray(res, "warning") # array of warning indicators
time <- getArray(res, "time") # array of user times in ms

## warnings, errors
if(any(err > 0))
  ftable(100* err, col.vars="tau") # percentage of errors
if(any(warn > 0))
  ftable(100*warn, col.vars="tau") # percentage of warnings

## run time
ftable(time, row.vars=c("family", "d"), col.vars="tau")

## add 'tau==' (just nicer)
str(val)
dimnames(val)[["tau"]] <- paste0("tau==", dimnames(val)[["tau"]])

## plot of VaR estimates
## plotting to pdf if not interactive graphics (=> R CMD BATCH VaRsuperadd.R)
names(varList[["qmargin"]][["value"]]) # "norm" , "t4", "Par2"
for(m in names(varList[["qmargin"]][["value"]])) {
  if(doPDF) pdf(file = (file <- paste0("VaR_superadd_", m, ".pdf")))
  mayplot(val[qmargin=m,,,,], varList, row.vars="family", col.vars="tau",
          xvar="alpha", ylim=if(n.obs > 1000) "local" else "global",
          panel.first = function(...) abline(h=0, col="gray40"), # gray line
          panel.last = function(x,y, col, ...) {
            ## For demo: write the 'panel.last' string, dependent on (x,y)
            rx <- range(x); dx <- diff(rx)
            ry <- range(y); dy <- diff(ry)
            text(rx[1]+.7*dx, ry[1]+.2*dy, "< 0: superadd.", col=col)
          })
  if(doPDF) dev.off.crop(file)
}