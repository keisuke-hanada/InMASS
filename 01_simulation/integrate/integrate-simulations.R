
library(dplyr)
library(ggplot2)

foldername <- "04_simulation-v1.0/"
subfolder <- c("test_1to1", "test_3to1", "test_4to0")
filename <- "/total_evaluates.csv"
fname <- paste(foldername, subfolder, filename, sep="")
dname <- c("1:1-", "3:1-", "treatment-")


### 1:1 ###
d10 <- data.frame()
d01 <- read.csv(fname[1])
d01$cmethod[d01$method=="ipd_method"] <- ""
d01$cmethod[d01$method=="propose_method"] <- "+InMASS"
d01$cformula[d01$formula==1] <- "MID"
d01$cformula[d01$formula==2] <- "ID"

d01 <- d01 %>% 
  mutate(
    Method = paste(dname[1], cformula, cmethod, sep=""),
    dname = dname[1]
  )
d10 <- rbind(d10, d01)


head(d10)


strata0 <- unique(d10$strata)
n0 <- unique(d10$n)
strata.labs <- paste(strata0, " studies", sep="")
names(strata.labs) <- strata0
n.labs <- paste("n=", n0, sep="")
names(n.labs) <- n0
model.labs <- c("Normal", "Chi-squared")
names(model.labs) <- unique(d10$model)
formula.labs <- c("Model misspesified", "Model spesified")
names(formula.labs) <- unique(d10$formula)



g1dat <- subset(d10, subset=eval=="meanfunc")
g1dat2 <- subset(d10, subset=eval=="varfunc")
g1dat$mse <- g1dat2$mean + (g1dat$mean-2)^2
g1 <- ggplot(g1dat, aes(x=n, y=mse, col=Method, shape=Method, linetype=Method)) +
  geom_line() +
  geom_point() +
  geom_abline(slope=0, intercept=0) +
  ylab("Mean Squared Error") +
  facet_grid(model ~ strata, labeller = labeller(strata = strata.labs, model = model.labs)) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1),
         shape = guide_legend(nrow = 1),
         fill  = guide_legend(nrow = 1))
g1


g2dat <- subset(d10, subset=eval=="powerfunc" & model=="make_model1")
g2 <- ggplot(g2dat, aes(x=sd-2, y=mean, group=Method, col=Method, shape=Method, linetype=Method)) +
  geom_abline(slope=0, intercept=0.05, linetype="dashed") +
  geom_abline(slope=0, intercept=0.8, linetype="dashed") +
  geom_abline(slope=0, intercept=0.9, linetype="dashed") +
  geom_line(linewidth=1) +
  ylab("Power") +
  xlab("Difference from CATE for target trial") +
  ylim(c(0,1)) +
  ggtitle("Simulate model 1: normally distributed") +
  facet_grid(n ~ strata, labeller = labeller(n = n.labs, strata = strata.labs)) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1),
         shape = guide_legend(nrow = 1),
         fill  = guide_legend(nrow = 1))
g2

g3dat <- subset(d10, subset=eval=="powerfunc" & model=="make_model2")
g3 <- ggplot(g3dat, aes(x=sd-2, y=mean, group=Method, col=Method, shape=Method, linetype=Method)) +
  geom_abline(slope=0, intercept=0.05, linetype="dashed") +
  geom_abline(slope=0, intercept=0.8, linetype="dashed") +
  geom_abline(slope=0, intercept=0.9, linetype="dashed") +
  geom_line(linewidth=1) +
  ylab("Power") +
  xlab("Difference from CATE for target trial") +
  ylim(c(0,1)) +
  ggtitle("Simulate model 2: chi-squared") +
  facet_grid(n ~ strata, labeller = labeller(n = n.labs, strata = strata.labs)) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1),
         shape = guide_legend(nrow = 1),
         fill  = guide_legend(nrow = 1))
g3

ggsave(paste(foldername, "integrate/figure-1-1to1-mse.pdf", sep="" ), g1, width=6, height=4, dpi=300)
ggsave(paste(foldername, "integrate/figure-1-1to1-power-model1.pdf", sep="" ), g2, width=8, height=8, dpi=300)
ggsave(paste(foldername, "integrate/figure-1-1to1-power-model2.pdf", sep="" ), g3, width=8, height=8, dpi=300)



g4dat <- subset(d10, subset=eval=="biasfunc" & n==40)
g4 <- ggplot(g4dat, aes(x=Method, y=mean)) +
  geom_boxplot() +
  facet_grid(model~strata)
g4



### 3:1 and only treatment ###
d10 <- data.frame()
for(i in 1:length(fname)){
  d01 <- read.csv(fname[i])
  d01$cmethod[d01$method=="ipd_method"] <- ""
  d01$cmethod[d01$method=="propose_method"] <- "+InMASS"
  d01$cformula[d01$formula==1] <- "MID"
  d01$cformula[d01$formula==2] <- "ID"
  
  d01 <- d01 %>% 
    mutate(
      Method = paste(dname[i], cformula, cmethod, sep=""),
      dname = dname[i]
    )
  d10 <- rbind(d10, d01)
}

head(d10)

d10 <- subset(d10, !Method %in% c("1:1-ID+InMASS", "1:1-MID+InMASS"))
unique(d10$Method)

d11 <- subset(d10, cformula=="ID")
g1dat <- subset(d11, subset=eval=="meanfunc")
g1dat2 <- subset(d11, subset=eval=="varfunc")
g1dat$mse <- g1dat2$mean + (g1dat$mean-2)^2
g1 <- ggplot(g1dat, aes(x=n, y=mse, col=Method, shape=Method, linetype=Method)) +
  geom_line() +
  geom_point() +
  geom_abline(slope=0, intercept=0) +
  ylab("Mean Squared Error") +
  facet_grid(model ~ strata, labeller = labeller(strata = strata.labs, model = model.labs)) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1),
         shape = guide_legend(nrow = 1),
         fill  = guide_legend(nrow = 1))
g1


g2dat <- subset(d11, subset=eval=="powerfunc" & model=="make_model1")
g2 <- ggplot(g2dat, aes(x=sd-2, y=mean, group=Method, col=Method, shape=Method, linetype=Method)) +
  geom_abline(slope=0, intercept=0.05, linetype="dashed") +
  geom_abline(slope=0, intercept=0.8, linetype="dashed") +
  geom_abline(slope=0, intercept=0.9, linetype="dashed") +
  geom_line(linewidth=1) +
  ylab("Power") +
  xlab("Difference from CATE for target trial") +
  ylim(c(0,1)) +
  ggtitle("Simulate model 1: normally distributed") +
  facet_grid(n ~ strata, labeller = labeller(n = n.labs, strata = strata.labs)) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1),
         shape = guide_legend(nrow = 1),
         fill  = guide_legend(nrow = 1))
g2

g3dat <- subset(d11, subset=eval=="powerfunc" & model=="make_model2")
g3 <- ggplot(g3dat, aes(x=sd-2, y=mean, group=Method, col=Method, shape=Method, linetype=Method)) +
  geom_abline(slope=0, intercept=0.05, linetype="dashed") +
  geom_abline(slope=0, intercept=0.8, linetype="dashed") +
  geom_abline(slope=0, intercept=0.9, linetype="dashed") +
  geom_line(linewidth=1) +
  ylab("Power") +
  xlab("Difference from CATE for target trial") +
  ylim(c(0,1)) +
  ggtitle("Simulate model 2: chi-squared") +
  facet_grid(n ~ strata, labeller = labeller(n = n.labs, strata = strata.labs)) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1),
         shape = guide_legend(nrow = 1),
         fill  = guide_legend(nrow = 1))
g3

ggsave(paste(foldername, "integrate/figure-2-controlAD-model-identify-mse.pdf", sep="" ), g1, width=6, height=4, dpi=300)
ggsave(paste(foldername, "integrate/figure-2-controlAD-model-identify-power-model1.pdf", sep="" ), g2, width=8, height=8, dpi=300)
ggsave(paste(foldername, "integrate/figure-2-controlAD-model-identify-power-model2.pdf", sep="" ), g3, width=8, height=8, dpi=300)



d11 <- subset(d10, cformula=="MID")
g1dat <- subset(d11, subset=eval=="meanfunc")
g1dat2 <- subset(d11, subset=eval=="varfunc")
g1dat$mse <- g1dat2$mean + (g1dat$mean-2)^2
g1 <- ggplot(g1dat, aes(x=n, y=mse, col=Method, shape=Method, linetype=Method)) +
  geom_line() +
  geom_point() +
  geom_abline(slope=0, intercept=0) +
  ylab("Mean Squared Error") +
  facet_grid(model ~ strata, labeller = labeller(strata = strata.labs, model = model.labs)) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1),
         shape = guide_legend(nrow = 1),
         fill  = guide_legend(nrow = 1))
g1


g2dat <- subset(d11, subset=eval=="powerfunc" & model=="make_model1")
g2 <- ggplot(g2dat, aes(x=sd-2, y=mean, group=Method, col=Method, shape=Method, linetype=Method)) +
  geom_abline(slope=0, intercept=0.05, linetype="dashed") +
  geom_abline(slope=0, intercept=0.8, linetype="dashed") +
  geom_abline(slope=0, intercept=0.9, linetype="dashed") +
  geom_line(linewidth=1) +
  ylab("Power") +
  xlab("Difference from CATE for target trial") +
  ylim(c(0,1)) +
  ggtitle("Simulate model 1: normally distributed") +
  facet_grid(n ~ strata, labeller = labeller(n = n.labs, strata = strata.labs)) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1),
         shape = guide_legend(nrow = 1),
         fill  = guide_legend(nrow = 1))
g2

g3dat <- subset(d11, subset=eval=="powerfunc" & model=="make_model2")
g3 <- ggplot(g3dat, aes(x=sd-2, y=mean, group=Method, col=Method, shape=Method, linetype=Method)) +
  geom_abline(slope=0, intercept=0.05, linetype="dashed") +
  geom_abline(slope=0, intercept=0.8, linetype="dashed") +
  geom_abline(slope=0, intercept=0.9, linetype="dashed") +
  geom_line(linewidth=1) +
  ylab("Power") +
  xlab("Difference from CATE for target trial") +
  ylim(c(0,1)) +
  ggtitle("Simulate model 2: chi-squared") +
  facet_grid(n ~ strata, labeller = labeller(n = n.labs, strata = strata.labs)) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1),
         shape = guide_legend(nrow = 1),
         fill  = guide_legend(nrow = 1))
g3

ggsave(paste(foldername, "integrate/figure-3-controlAD-model-misidentify-mse.pdf", sep="" ), g1, width=6, height=4, dpi=300)
ggsave(paste(foldername, "integrate/figure-3-controlAD-model-misidentify-power-model1.pdf", sep="" ), g2, width=8, height=8, dpi=300)
ggsave(paste(foldername, "integrate/figure-3-controlAD-model-misidentify-power-model2.pdf", sep="" ), g3, width=8, height=8, dpi=300)







### K=10 main plots ###
d10 <- data.frame()
for(i in 1:length(dname)){
  d01 <- read.csv(fname[i])
  d01$cmethod[d01$method=="ipd_method"] <- ""
  d01$cmethod[d01$method=="propose_method"] <- "+InMASS"
  d01$cmethod[d01$method=="meta_method"] <- " MA"
  d01$cformula[d01$formula==1] <- "MID"
  d01$cformula[d01$formula==2] <- "ID"
  d01 <- d01 %>% subset(strata==10) %>%
    mutate(
      Method = paste(substr(dname[i], 1, nchar(dname[i])-1), cmethod, sep=""),
      dname = dname[i]
    )
  d10 <- rbind(d10, d01)
}
d10$formula <- factor(d10$formula, levels = c(2, 1))

head(d10)


g1dat <- subset(d10, subset=eval=="meanfunc")
g1dat2 <- subset(d10, subset=eval=="varfunc")
g1dat$mse <- g1dat2$mean + (g1dat$mean-2)^2
g1 <- ggplot(g1dat, aes(x=n, y=mse, col=Method, shape=Method, linetype=Method)) +
  geom_line() +
  geom_point() +
  geom_abline(slope=0, intercept=0) +
  ylab("Mean Squared Error") +
  facet_grid(formula ~ model, labeller = labeller(formula = formula.labs, model = model.labs)) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1),
         shape = guide_legend(nrow = 1),
         fill  = guide_legend(nrow = 1))
g1


g2dat <- subset(d10, subset=eval=="powerfunc" & model=="make_model1")
g2 <- ggplot(g2dat, aes(x=sd-2, y=mean, group=Method, col=Method, shape=Method, linetype=Method)) +
  geom_abline(slope=0, intercept=0.05, linetype="dashed") +
  geom_abline(slope=0, intercept=0.8, linetype="dashed") +
  geom_abline(slope=0, intercept=0.9, linetype="dashed") +
  geom_line(linewidth=1) +
  ylab("Power") +
  xlab("Difference from CATE for target trial") +
  ylim(c(0,1)) +
  ggtitle("Simulate model 1: normally distributed") +
  facet_grid(formula ~ n, labeller = labeller(formula = formula.labs, n = n.labs)) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1),
         shape = guide_legend(nrow = 1),
         fill  = guide_legend(nrow = 1))
g2

g3dat <- subset(d10, subset=eval=="powerfunc" & model=="make_model2")
g3 <- ggplot(g3dat, aes(x=sd-2, y=mean, group=Method, col=Method, shape=Method, linetype=Method)) +
  geom_abline(slope=0, intercept=0.05, linetype="dashed") +
  geom_abline(slope=0, intercept=0.8, linetype="dashed") +
  geom_abline(slope=0, intercept=0.9, linetype="dashed") +
  geom_line(linewidth=1) +
  ylab("Power") +
  xlab("Difference from CATE for target trial") +
  ylim(c(0,1)) +
  ggtitle("Simulate model 2: chi-squared") +
  facet_grid(formula ~ n, labeller = labeller(formula = formula.labs, n = n.labs)) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1),
         shape = guide_legend(nrow = 1),
         fill  = guide_legend(nrow = 1))
g3

ggsave(paste(foldername, "integrate/main-figure-1-mse.pdf", sep="" ), g1, width=6, height=4, dpi=300)
ggsave(paste(foldername, "integrate/main-figure-2-power-model1.pdf", sep="" ), g2, width=8, height=5, dpi=300)
ggsave(paste(foldername, "integrate/main-figure-3-power-model2.pdf", sep="" ), g3, width=8, height=5, dpi=300)



g4dat <- subset(d10, subset=eval=="biasfunc" & n==40)
g4 <- ggplot(g4dat, aes(x=Method, y=mean)) +
  geom_hline(yintercept=0, linetype=2) +
  geom_boxplot() +
  ylab(expression("Bias of "~delta[T])) +
  facet_grid(formula ~ model, labeller = labeller(formula = formula.labs, model = model.labs)) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle=25, hjust=1, size=12),
        strip.text.x = element_text(size=16),
        strip.text.y = element_text(size=16),
        axis.title = element_text(size=16)) +
  guides(color = guide_legend(nrow = 1),
         shape = guide_legend(nrow = 1),
         fill  = guide_legend(nrow = 1))
g4
ggsave(paste(foldername, "integrate/main-figure-0-bias.pdf", sep="" ), g4, width=12, height=6, dpi=300)




### all bias plots ###
stratas <- c(5,10,30)
ns <- c(20,40,100)

d10 <- data.frame()
for(i in 1:length(dname)){
  d01 <- read.csv(fname[i])
  d01$cmethod[d01$method=="ipd_method"] <- ""
  d01$cmethod[d01$method=="propose_method"] <- "+InMASS"
  d01$cmethod[d01$method=="meta_method"] <- " MA"
  d01$cformula[d01$formula==1] <- "MID"
  d01$cformula[d01$formula==2] <- "ID"
  d01 <- d01 %>%
    mutate(
      Method = paste(substr(dname[i], 1, nchar(dname[i])-1), cmethod, sep=""),
      dname = dname[i]
    )
  d10 <- rbind(d10, d01)
}
d10$formula <- factor(d10$formula, levels = c(2, 1))


pdf(paste(foldername, "integrate/main-figure-0-bias-all.pdf", sep="" ), width=16, height=6)

for(st in stratas) {
  for (nval in ns) {
    
    g4dat <- subset(d10, subset=eval=="biasfunc" & strata==st & n==nval)
    g4 <- ggplot(g4dat, aes(x=Method, y=mean)) +
      geom_hline(yintercept=0, linetype=2) +
      geom_boxplot() +
      coord_cartesian(ylim=c(-2,2)) +
      labs(
        title = paste(st, " studies with n=", nval, sep=""),
        y = expression("Bias of "~delta[T])
      ) +
      facet_grid(formula ~ model, labeller = labeller(formula = formula.labs, model = model.labs)) +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle=25, hjust=1, size=12),
            strip.text.x = element_text(size=16),
            strip.text.y = element_text(size=16),
            axis.title = element_text(size=16)) +
      guides(color = guide_legend(nrow = 1),
             shape = guide_legend(nrow = 1),
             fill  = guide_legend(nrow = 1))
    print(g4)
    

  }
}
dev.off()

