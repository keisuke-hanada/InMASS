library(Rcpp)
library(data.table)


sourceCpp("04_simulation-v1.0/test_4to0/make_models.cpp")

# x2k: normal
make_model1 <- function(n, sigma, strata, nsim) {
  x <- make_models(n=n, sigma=sigma, strata=strata, nsim=nsim, dist_x2="normal", dist_x3="normal")
  return(x)
}

# x2k: chisq
make_model2 <- function(n, sigma, strata, nsim) {
  x <- make_models(n=n, sigma=sigma, strata=strata, nsim=nsim, dist_x2="chi2", dist_x3="normal")
  return(x)
}

# n <- 100
# d01 <- make_model1(n=n, sigma=1, strata=10, nsim=1)
# d02 <- make_model2(n=n, sigma=1, strata=5, nsim=1)
# 
# 
# hist(d01$target_ipd$x3k, freq=F, breaks=100)
# hist(d02$target_ipd$x3k, freq=F, breaks=100)
# c(mean(d01$target_ipd$x3k),
#   mean(d02$target_ipd$x3k))
# c(var(d01$target_ipd$x3k),
#   var(d02$target_ipd$x3k))
# 
# c(mean(d01$target_ipd$x2k),
#   mean(d02$target_ipd$x2k))
# c(var(d01$target_ipd$x2k),
#   var(d02$target_ipd$x2k))
# 
# 
# 
# hist(subset(d01$strata_ipd, nsim==1 & strata==1)$x2k)
# hist(subset(d02$strata_ipd, nsim==1 & strata==1)$x2k)
# 
# boxplot(x2k ~ strata, subset(d01$strata_ipd, nsim==1))
# boxplot(x2k ~ strata, subset(d02$strata_ipd, nsim==1))
# 
# boxplot(x3k ~ strata, subset(d01$strata_ipd, nsim==1))
# boxplot(x3k ~ strata, subset(d02$strata_ipd, nsim==1))
# 
# d11 <- subset(d02$strata_ipd, nsim==1)
# plot(d11$x2k, d11$yik)


