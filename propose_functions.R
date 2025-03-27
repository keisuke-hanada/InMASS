library(metafor)


iwlm <- function(formula, data, weights) { #Importance weighted linear regression model
  
  Yall <- model.response(model.frame(formula, data))
  Xall <- model.matrix(formula, data)
  Wall <- diag(weights)
  p.all <- ncol(Xall)
  
  InvX <- solve(t(Xall) %*% Wall %*% Xall)
  bh <- InvX %*% t(Xall) %*% Wall %*% Yall
  var.bh <- InvX %*% (t(Xall) %*% Wall %*% diag(weights^2*c((Yall-Xall%*%bh))^2) %*% Wall %*% Xall) %*% InvX
  
  
  fitted <- Xall %*% bh
  residuals <- Yall - fitted
  
  result <- list(
    coefficients = bh,
    cov = var.bh,
    residuals = residuals,
    fitted.values = fitted,
    formula = formula,
    data = data,
    weights = weights
  )
  
  class(result) <- "iwlm"
  return(result)
  
}


print.iwlm <- function(x) {
  cat("Call:\n")
  print(x$formula)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}

summary.iwlm <- function(object) {
  X <- model.matrix(object$formula, object$data)
  y <- model.response(model.frame(object$formula, object$data))
  W <- diag(object$weights)
  n <- nrow(X)
  nw <- sum(W)
  p <- ncol(X)
  residuals <- object$residuals
  coef <- object$coefficients
  var_beta <- object$cov
  

  sigma <- sqrt(sum(W%*%residuals^2) / (nw - p))
  se <- sqrt(diag(var_beta))
  

  t_values <- coef / se
  p_values <- 2 * pt(-abs(t_values), df = nw - p)
  

  p_values_formatted <- ifelse(p_values < 0.0001, "<0.0001", sprintf("%.4f", p_values))
  

  residual_summary <- summary(residuals)
  
  
  r_squared <- round(1 - sum(W%*%residuals^2) / sum(W%*%(y - mean(y))^2), 4)
  adj_r_squared <- round(1 - (1 - (1 - sum(W%*%residuals^2) / sum(W%*%(y - mean(y))^2))) * (nw - 1) / (nw - p), 4)
  

  result <- list(
    call = object$formula,
    residuals = residuals,
    coefficients = data.frame(
      Estimate = coef,
      Std.Error = se,
      t.value = t_values,
      Pr = p_values_formatted,
      row.names = names(coef)
    ),
    sigma = sigma,
    weighted.r.squared = r_squared,
    weighted.adj.r.squared = adj_r_squared,
    df = c(p, nw - p)
  )
  
  class(result) <- "summary.iwlm"
  return(result)
}

print.summary.iwlm <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  
  cat("\nResiduals:\n")
  print(summary(x$residuals))
  
  cat("\nCoefficients:\n")
  print(x$coefficients)
  
  cat("\nResidual standard error:", round(x$sigma, 4),
      "on", x$df[2], "degrees of freedom\n")
  cat("Multiple Weighted R-squared:", round(x$weighted.r.squared, 4), "\n")
  cat("Weighted Adjusted R-squared:", round(x$weighted.adj.r.squared, 4), "\n")
}


confint.iwlm <- function(object, level = 0.95) {
  X <- model.matrix(object$formula, head(object$data))
  coef <- object$coefficients
  var_beta <- object$cov
  
  nw <- sum(object$weights)
  p <- ncol(X)
  se <- sqrt(diag(var_beta))
  

  alpha <- 1 - level
  t_value <- qt(1 - alpha / 2, df = nw - p)
  lower <- coef - t_value * se
  upper <- coef + t_value * se
  

  ci <- data.frame(
    Lower = lower,
    Upper = upper
  )
  
  return(ci)
}



### total of integrate_meta_into_specific_study ###
inmass <- function(formula, formula.ma, data.ma.mean, data.ma.var, data.ipd, strata,
                   ps_meta=1, seed=1234){ 
  
  set.seed(seed)
  
  ################################################################################
  ### meta-analysis
  ################################################################################
  vars <- all.vars(formula.ma)
  yik <- data.ma.mean[[vars[1]]]
  v2.yk <- data.ma.var[[vars[1]]]
  ma <- tryCatch({
    rma(yik, v2.yk, mods = update(formula.ma, NULL ~ .), data = data.ma.mean, method="DL")
  }, error = function(e) {
    list(b = c(NA, NA, NA, NA))
  })
  betah <- ma$b
  p <- length(vars)
  
  
  ################################################################################
  ### reconstructing IPD
  ################################################################################
  dat01ripd <- data.frame()
  arm0 <- unique(data.ma.mean[[vars[2]]])
  
  
  for(st in 1:strata){
    ad02mean <- subset(data.ma.mean, subset=strata==st)
    ad02var <- subset(data.ma.var, subset=strata==st)
    n_strata <- ad02mean$n
    
    ripd03 <- data.frame()
    x_vars <- vars[-c(1:2)]
    for(arm in arm0) {
      ripd01 <- ad02mean[rep(arm+1, each=n_strata[arm+1]),]
      for(var in x_vars) {
        ripd01[[var]] <- rnorm(n_strata[arm+1], mean=ad02mean[[var]][arm+1], sd=sqrt(ad02var[[var]][arm+1]))
      }
      
      
      ripd02 <- model.frame(formula.ma, data=ripd01)
      
      y_hat <- model.matrix(formula.ma, data=ripd02) %*% betah
      sd_y <- sqrt( max(ad02var[[vars[1]]][arm+1] - var(y_hat), 0) )
      ripd02[[1]] <- y_hat + rnorm(n_strata[arm+1], sd=sd_y)
      
      ripd02$strata <- st
      ripd02$vi <- sd_y^2
      
      ripd03 <- rbind(ripd03, ripd02)
    }
    
    dat01ripd <- rbind(dat01ripd, ripd03)
    
  }
  plot(dat01ripd$x, dat01ripd$yi, col=dat01ripd$z+1)
  
  

  data.ipd$vi[data.ipd[vars[2]]==0] <- var(data.ipd[[vars[1]]][data.ipd[vars[2]]==0])
  data.ipd$vi[data.ipd[vars[2]]==1] <- var(data.ipd[[vars[1]]][data.ipd[vars[2]]==1])
  data.ipd$strata <- 0
  data.ipd1 <- data.ipd
  misscol <- setdiff(names(data.ipd), names(dat01ripd))
  for (col in misscol){
    data.ipd1[[col]] <- NA
  }
  data.ipd1 <- data.ipd1[names(dat01ripd)]
  if(ps_meta == 1) {
    dat01all <- rbind(data.ipd1, dat01ripd)
  }else if(ps_meta == 2) {
    dat01all <- rbind(data.ipd1, subset(dat01ripd, dat01ripd[[vars[2]]]==0))
  }else {
    stop("ps_meta need to set 1 or 2. 
         1: treatment and control groups in meta-analysis are used to calculate PS, 
         2: control group in meta-analysis is used to calculate PS.")
  }
  
  
  dat01all$set <- as.numeric(dat01all$strata==0)
  
  n1 <- sum(data.ipd[[vars[2]]])
  n0 <- sum(1-data.ipd[[vars[2]]])
  
  ################################################################################
  ### density ratio
  ################################################################################
  vars_new <- vars[-c(1:2)]
  vars_new2 <- paste("I(", vars_new, "^2)", sep="")
  
  ps_dat <- rbind(data.frame(id=1, subset(dat01all, set==1)), data.frame(id=0, dat01all))
  n_target <- sum(ps_dat$id==1)
  n_source <- length(ps_dat$id==0)
  ps_formula <- reformulate(c("1", vars_new, vars_new2), response = "id")
  ps_model <- glm(ps_formula,
                  family = binomial(link = "logit"),
                  data = ps_dat)
  psval <- predict(ps_model, type="response")[ps_dat$id==0]
  
  dat01all$weight <- psval / (1 - psval) * (n_source/n_target)
  if(ps_meta == 2){
    dat01all$weight[dat01all[[vars[2]]]==1] <- 1
  }
  
  n1h <- sum(dat01all$weight[dat01all[[vars[2]]]==1])
  n0h <- sum(dat01all$weight[dat01all[[vars[2]]]==0])

  
  
  ################################################################################
  ### weighted regression
  ################################################################################
  
  retval <- iwlm(formula, data=dat01all, weights=dat01all$weight)
  confint(retval)
  return(retval)
}




### total of integrate_meta_into_specific_study for control group ###
inmass_control <- function(formula, formula.ma, data.ma.mean, data.ma.var, data.ipd, strata,
                           ps_meta=2){ 
  
  
  ################################################################################
  ### meta-analysis
  ################################################################################
  vars <- all.vars(formula.ma)
  yik <- data.ma.mean[[vars[1]]]
  v2.yk <- data.ma.var[[vars[1]]]
  ma <- tryCatch({
    rma(yik, v2.yk, mods = update(formula.ma, NULL ~ .), data = data.ma.mean, method="DL")
  }, error = function(e) {
    list(b = c(NA, NA, NA, NA))
  })
  betah <- ma$b
  p <- length(vars)
  
  
  ################################################################################
  ### reconstructing IPD
  ################################################################################
  dat01ripd <- data.frame()
  arm0 <- unique(data.ma.mean[[vars[2]]])
  
  
  for(st in 1:strata){
    ad02mean <- subset(data.ma.mean, subset=strata==st)
    ad02var <- subset(data.ma.var, subset=strata==st)
    n_strata <- ad02mean$n
    
    ripd03 <- data.frame()
    x_vars <- vars[-c(1:2)]
    for(arm in arm0) {
      ripd01 <- ad02mean[rep(arm+1, each=n_strata[arm+1]),]
      for(var in x_vars) {
        ripd01[[var]] <- rnorm(n_strata[arm+1], mean=ad02mean[[var]][arm+1], sd=sqrt(ad02var[[var]][arm+1]))
      }
      
      
      ripd02 <- model.frame(formula.ma, data=ripd01)
      
      y_hat <- model.matrix(formula.ma, data=ripd02) %*% betah
      sd_y <- sqrt( max(ad02var[[vars[1]]][arm+1] - var(y_hat), 0) )
      ripd02[[1]] <- y_hat + rnorm(n_strata[arm+1], sd=sd_y)
      
      ripd02$strata <- st
      ripd02$vi <- sd_y^2
      
      ripd03 <- rbind(ripd03, ripd02)
    }
    
    dat01ripd <- rbind(dat01ripd, ripd03)
    
  }
  


  data.ipd$vi[data.ipd[vars[2]]==0] <- var(data.ipd[[vars[1]]][data.ipd[vars[2]]==0])
  data.ipd$vi[data.ipd[vars[2]]==1] <- var(data.ipd[[vars[1]]][data.ipd[vars[2]]==1])
  data.ipd$strata <- 0
  data.ipd1 <- data.ipd
  misscol <- setdiff(names(data.ipd), names(dat01ripd))
  for (col in misscol){
    data.ipd1[[col]] <- NA
  }
  data.ipd1 <- data.ipd1[names(dat01ripd)]
  if(ps_meta == 1) {
    dat01all <- rbind(data.ipd1, dat01ripd)
  }else if(ps_meta == 2) {
    dat01all <- rbind(data.ipd1, subset(dat01ripd, dat01ripd[[vars[2]]]==0))
  }else {
    stop("ps_meta need to set 1 or 2. 
         1: treatment and control groups in meta-analysis are used to calculate PS, 
         2: control group in meta-analysis is used to calculate PS.")
  }
  
  
  dat01all$set <- as.numeric(dat01all$strata==0)
  
  n1 <- sum(data.ipd[[vars[2]]])
  n0 <- sum(1-data.ipd[[vars[2]]])
  
  ################################################################################
  ### density ratio
  ################################################################################
  vars_new <- vars[-c(1:2)]
  vars_new2 <- paste("I(", vars_new, "^2)", sep="")
  
  ps_dat <- rbind(data.frame(id=1, subset(dat01all, set==1)), data.frame(id=0, dat01all))
  n_target <- sum(ps_dat$id==1)
  n_source <- sum(ps_dat$id==0)
  
  ps_formula <- reformulate(c("1", vars_new, vars_new2), response = "id")
  ps_model <- glm(ps_formula,
                  family = binomial(link = "logit"),
                  data = ps_dat)
  psval <- predict(ps_model, type="response")[ps_dat$id==0]
  dat01all$weight <- psval / (1 - psval) * (n_source/n_target)

  
  n1h <- sum(dat01all$weight[dat01all[[vars[2]]]==1])
  n0h <- sum(dat01all$weight[dat01all[[vars[2]]]==0])


  
  ################################################################################
  ### weighted regression
  ################################################################################
  
  retval <- iwlm(formula, data=dat01all, weights=dat01all$weight)

  return(retval)
}


