library(dplyr)
library(tidyr)

set.seed(202412)

source("propose_functions.R")

################################################################################
### Meta-analysis by Sanguankeo et al (2015)
################################################################################

dat01 <- read.csv("02_case study/Sanguankeo et al_2015.csv")

dat01$yi <- dat01$Mean.Difference.Random
dat01$vi <- ( (dat01$X95..CI.Upper - dat01$X95..CI.Lower) / (2 * qnorm(0.975)) )^2
dat01$x1 <- dat01$Baseline.eGFR.Statin.Mean
dat01$x1.vi <- dat01$Baseline.eGFR.Control.SD^2
dat01$x0 <- dat01$Baseline.eGFR.Control.Mean
dat01$x0.vi <- dat01$Baseline.eGFR.Statin.SD^2
dat01$n1 <- dat01$Treatment.Participants..n.
dat01$n0 <- dat01$Control.Participants..n.
dat01$fp <- dat01$Follow.up..months.


dat02 <- dat01 %>% 
  subset(Outcome=="rate of eGFR change per year") %>%
  select(Study, yi, vi, x1, x1.vi, x0, x0.vi, n1, n0, fp)
dat02

ma0 <- rma(yi=yi, vi=vi, data=dat02)
forest(ma0, slab=dat02$Study)


### create outcome and covariate of the statin and Control group

dat03 <- dat01 %>% 
  subset(Outcome=="total change in eGFR") %>%
  select(Study, yi, vi, x1, x1.vi, x0, x0.vi, n1, n0, fp) %>%
  mutate(
    z1 = yi - fp/12,
    z0 = - fp/12,
    z1.vi = n1*vi / (n1+n0),
    z0.vi = n0*vi / (n1+n0)
  )
dat03$z1 - dat03$z0
dat03$yi
dat03
ma1 <- rma(yi=yi, vi=vi, data=dat03)
forest(ma1, slab=dat03$Study)



dat04 <- rbind(
  data.frame(study=dat03$Study, yi=dat03$z1, vi=dat03$z1.vi, z=1, x=dat03$x1, x.vi=dat03$x1.vi, n=dat03$n1),
  data.frame(study=dat03$Study, yi=dat03$z0, vi=dat03$z0.vi, z=0, x=dat03$x0, x.vi=dat03$x0.vi, n=dat03$n0)
)
dat04$strata <- as.numeric(factor(dat04$study))


ma <- rma(yi=yi, vi=vi, mods= ~ z + x, data=dat04)
summary(ma)


################################################################################
### Meta-analysis without Sawara 2008
################################################################################

target <- "Sawara 2008"
dat05 <- dat04 %>% subset(study!=target)
ma2 <- rma(yi=yi, vi=vi, mods = ~ z + x, data=dat05, method="DL")
tau2.ci <- confint(ma2, type="QP")$random[1,2:3]
tau2.ci
res.ma2 <- data.frame(b=c(ma2$b, ma2$tau2), 
           se=c(ma2$se, ma2$se.tau2), 
           cil=c(ma2$ci.lb, tau2.ci[1]), 
           ciu=c(ma2$ci.ub,tau2.ci[2])
           )
res.ma2 <- round(res.ma2, 3)
rownames(res.ma2) <- c("intercept", "treatment", "baseline eGFR", "$tau^2$")
res.ma2

library(xtable)
res.ma2_latex <- xtable(res.ma2)
print(res.ma2_latex)



################################################################################
### Virtual trial mimicking Sawara 2008
################################################################################

### sample size in 1:1 allocation
dat10 <- subset(dat04, study==target)
delta.target <- dat10$yi[1]-dat10$yi[2]
sd.target <- sqrt(sum(dat10$n*dat10$vi))
n_max <- 200

### Sawara 2008 + meta-analysis
n_target <- sum(dat10$n)
dat11 <- data.frame(
  study = target,
  z = rep(1, dat10$n[1]),
  yi0 = rnorm(dat10$n[1]),
  x0 = rnorm(dat10$n[1])
) %>%
  mutate(
    yi = dat10$yi[1] + (yi0 - mean(yi0))/sd(yi0) * sqrt(dat10$n[1]*dat10$vi[1]),
    x = dat10$x[1] + (x0 - mean(x0))/sd(x0) * sqrt(dat10$x.vi[1])
  )
dat12 <- data.frame(
  study = target,
  z = rep(0, dat10$n[2]),
  yi0 = rnorm(dat10$n[2]),
  x0 = rnorm(dat10$n[2])
) %>%
  mutate(
    yi = dat10$yi[2] + (yi0 - mean(yi0))/sd(yi0) * sqrt(dat10$n[2]*dat10$vi[2]),
    x = dat10$x[2] + (x0 - mean(x0))/sd(x0) * sqrt(dat10$x.vi[2])
  )
dat13 <- rbind(dat11, dat12)

head(dat13)

lm.target1 <- lm(yi ~ z + x, data=dat13)
summary(lm.target1)
confint(lm.target1)


dat14 <- dat04 %>% subset(study!=target) %>% select(study, yi, z, x, strata, n)
head(dat14)
dat15 <- dat04 %>% subset(study!=target) %>% select(study, vi, z, x.vi, strata, n)
names(dat15) <- names(dat14)
head(dat15)
dat14$strata <- as.numeric(factor(dat14$study))
dat15$strata <- as.numeric(factor(dat15$study))


res01 <- inmass(formula = yi ~ z + x, formula.ma = yi ~ z + x, 
               data.ipd = dat13, data.ma.mean = dat14, data.ma.var = dat15, strata=4)

summary(res01)
confint(res01)



### When allocated to 10 vs. 30 participants
n1 <- n_max/2
n0 <- n_max/4
dat20 <- subset(dat04, study==target)
dat21 <- data.frame(
  study = target,
  z = rep(1, n1),
  yi0 = rnorm(n1),
  x0 = rnorm(n1)
) %>%
  mutate(
    yi = dat10$yi[1] + (yi0 - mean(yi0))/sd(yi0) * sqrt(dat10$n[1]*dat10$vi[1]),
    x = dat10$x[1] + (x0 - mean(x0))/sd(x0) * sqrt(dat10$x.vi[1])
  )
dat22 <- data.frame(
  study = target,
  z = rep(0, n0),
  yi0 = rnorm(n0),
  x0 = rnorm(n0)
) %>%
  mutate(
    yi = dat10$yi[2] + (yi0 - mean(yi0))/sd(yi0) * sqrt(dat10$n[2]*dat10$vi[2]),
    x = dat10$x[2] + (x0 - mean(x0))/sd(x0) * sqrt(dat10$x.vi[2])
  )
dat23 <- rbind(dat21, dat22)
head(dat23)

lm.target2 <- lm(yi ~ z + x, data=dat23)
summary(lm.target2)
confint(lm.target2)


res02 <- inmass_control(formula = yi ~ z + x, formula.ma = yi ~ z + x, 
               data.ipd = dat23, data.ma.mean = dat14, data.ma.var = dat15, strata=4)

summary(res02)
confint(res02)



### When allocated to 10 vs. 30 participants
n1 <- n_max/2

dat31 <- data.frame(
  study = target,
  z = rep(1, n1),
  yi0 = rnorm(n1),
  x0 = rnorm(n1)
) %>%
  mutate(
    yi = dat10$yi[1] + (yi0 - mean(yi0))/sd(yi0) * sqrt(dat10$n[1]*dat10$vi[1]),
    x = dat10$x[1] + (x0 - mean(x0))/sd(x0) * sqrt(dat10$x.vi[1])
  )

head(dat31)

## impossible when only target trial
lm.target3 <- lm(yi ~ z + x, data=dat31)
summary(lm.target3)
confint(lm.target3)


res03 <- inmass_control(formula = yi ~ z + x, formula.ma = yi ~ z + x, 
               data.ipd = dat31, data.ma.mean = dat14, data.ma.var = dat15, strata=4)

summary(res03)
confint(res03)


################################################################################
### Summarize results
################################################################################

madat <- dat03 %>%
  select(Study, yi, vi, x1, x1.vi, x0, x0.vi, n1, n0, fp)

madat$vi <- sqrt(madat$vi)
madat$x1.vi <- sqrt(madat$x1.vi)
madat$x0.vi <- sqrt(madat$x0.vi)


format_columns <- function(df, column_pairs) {
  for (pair in column_pairs) {
    col1 <- pair[1]
    col2 <- pair[2]
    new_col_name <- paste0(col1, " (", col2, ")")
    df[[new_col_name]] <- sprintf("%.1f (%.1f)", df[[col1]], df[[col2]])
  }
  return(df)
}


column_pairs <- list(
  c("yi", "vi"),
  c("x1", "x1.vi"),
  c("x0", "x0.vi")
)


formatted_madat <- format_columns(madat, column_pairs)


formatted_madat <- formatted_madat[, c("Study", "n1", "n0", "fp", "yi (vi)", "x1 (x1.vi)", "x0 (x0.vi)")]

formatted_madat

resbase <- data.frame(
  method = c("Target", "Target with InMASS", "Target(2:1)", "Target(2:1) with InMASS control", "Single-arm", "Single-arm with InMASS control"),
  n1 = c(dat10$n[1], dat10$n[1], n_max/2, n_max/2, n_max/2, n_max/2),
  n0 = c(dat10$n[2], dat10$n[2], n_max/4, n_max/4, 0, 0),
  n12 = c(dat10$n[1], n_max/2, n_max/2, n_max/2, n_max/2, n_max/2),
  n02 = c(dat10$n[2], n_max/2, n_max/4, n_max/2, 0, n_max/2)
)
resdat <- rbind(
  c(summary(lm.target1)$coefficients[2,], confint(lm.target1)[2,]),
  c(summary(res01)$coefficients[2,], confint(res01)[2,]),
  c(summary(lm.target2)$coefficients[2,], confint(lm.target2)[2,]),
  c(summary(res02)$coefficients[2,], confint(res02)[2,]),
  c(NA),
  c(summary(res03)$coefficients[2,], confint(res03)[2,])
)
res <- cbind(resbase, resdat)
res$Estimate <- round(as.numeric(res$Estimate), 1)
res$`Std. Error` <- round(as.numeric(res$`Std. Error`), 2)
res$`t value` <- round(as.numeric(res$`t value`), 2)
res$`Pr(>|t|)`[c(1,3,5)] <- round(as.numeric(res$`Pr(>|t|)`[c(1,3,5)]), 4)
res$`Pr(>|t|)` <- as.character(res$`Pr(>|t|)`)
res$`2.5 %` <- round(as.numeric(res$`2.5 %`), 2)
res$`97.5 %` <- round(as.numeric(res$`97.5 %`), 2)

res <- res %>%
  mutate(
    yi = paste(Estimate, "(", `Std. Error`, ")", sep=""),
    ci = paste("[", `2.5 %`, ", ", `97.5 %`, "]", sep="")
  ) %>%
  select(method, n1, n0, yi, ci, `t value`, `Pr(>|t|)`)
res

write.csv(formatted_madat, file="02_case study/meta-analysis-data.csv")
write.csv(res, file="02_case study/rda_case_study.csv")




