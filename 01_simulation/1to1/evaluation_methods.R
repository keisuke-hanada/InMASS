
meanfunc <- function(res, method){
  res$params <- NULL
  res$x2k_m <- NULL
  beta <- sapply(res, function(x) x$coefficients[2])
  return(beta)
}

varfunc <- function(res, method){
  res$params <- NULL
  res$x2k_m <- NULL
  if(class(res[[1]]) == "iwlm"){
    v2 <- sapply(res, function(x) x$cov[2,2] )
  }else{
    v2 <- sapply(res, function(x) summary(x)$coefficients[2,2]^2)
  }
  return(v2)
}

cpfunc <- function(res, method){
  btrue <- res$params$beta[c(2,4)]
  # theta <- btrue[1] + btrue[2] * res$x2k_m
  theta <- btrue[1]
  res$params <- NULL
  res$x2k_m <- NULL
  ci <- sapply(res, function(x) confint(x)[2,])
  cp <- ci[1,] <= theta & theta <= ci[2,]
  return(cp)
}

powerfunc <- function(res, method, b.max=4, h=0.01){
  btrue <- res$params$beta[c(2,4)]
  nb <- (b.max-btrue[1])/h
  b <- btrue[1] + (0:nb)*h
  theta0 <- 0
  res$params <- NULL
  res$x2k_m <- NULL
  ci <- sapply(res, function(x) confint(x)[2,])
  power <- sapply(1:(nb+1), function(i) mean(b[i] + theta0 <= ci[1,] | ci[2,] <= b[i] + theta0))
  return(list(power=power, x=b))
}




