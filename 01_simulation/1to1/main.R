# install.packages("jsonlite")
# install.packages("digest")

library(jsonlite)
library(digest)


simfolder <- "04_simulation-v1.0/test_1to1/"



if (!dir.exists(simfolder)) {
  dir.create(simfolder)
  dir.create(paste(simfolder, "/simdata", sep=""))
  dir.create(paste(simfolder, "/result", sep=""))
  dir.create(paste(simfolder, "/evaluation", sep=""))
  dir.create(paste(simfolder, "/figure", sep=""))
  message("フォルダを新規作成しました！")
} else {
  message("フォルダは既に存在しています。")
}




source(paste(simfolder, "/simulate_models.R", sep=""))

model_func <- c("make_model1", "make_model2")
n0 <- c(20, 40, 100)
strata0 <- c(5, 10, 30)
nsim <- 10000
seed <- 1234

update <- TRUE



if(file.exists(paste(simfolder, "/meta_simdata.csv", sep="")) & !update) {
  metadat01 <- read.csv(file=paste(simfolder, "/meta_simdata.csv", sep=""))[,-1]
}else {
  metadat01 <- data.frame()
}

for(n in n0){
  for(strata in strata0){
    message(paste("n:", n, ", strata:", strata))
    for(model in model_func){
      params <- list(n=n, sigma=1, strata=strata, nsim=nsim, model=model)
      params_json <- toJSON(params, auto_unbox = TRUE)
      unique_key <- digest(params_json, algo = "md5")
      if (sum(metadat01$key == unique_key) == 0 || update) {
        
        message(paste("Create data: ", params_json))

        set.seed(seed)
        func <- get(model)
        d01 <- func(n=n, sigma=1, strata=strata, nsim=nsim)
        
        metadat02 <- data.frame(key = unique_key, n=n, sigma=1, strata=strata, nsim=nsim, model=model)
        
        simfile_name <- paste0(simfolder, "/simdata/", unique_key, ".rds", sep="")
        saveRDS(d01, file=simfile_name)
        
        
        metadat01 <- rbind(metadat01, metadat02)
        
      }else {
        
        message(paste("Exist data: ", params_json))
        
      }
    }
    
  }
}


write.csv(metadat01, file=paste(simfolder, "/meta_simdata.csv", sep=""))

rm(d01)
gc()




source(paste(simfolder, "/analysis_methods.R", sep=""))
methods <- c("ipd_method", "propose_method", "meta_method")
formulas <- 1:2


metadat03 <- read.csv(file=paste(simfolder, "/meta_simdata.csv", sep=""))
head(metadat03)

for(l in 1:length(metadat03$X)){
  print(metadat03[l,])
  unique_key <- metadat03$key[l]
  simfile_name <- paste0(simfolder, "/simdata/", unique_key, ".rds", sep="")
  dat01 <- readRDS(simfile_name)
  dat01$params$formula <- as.formula(dat01$params$formula)
  
  
  for(i.formula in formulas){
    if(i.formula == 1){
      method_formula <- yik ~ 1 + x1k
    }else{
      method_formula <- yik ~ 1 + x1k + x2k + x1k:x2k
    }
    
    for(method in methods){
      message(method)
      func <- get(method)
      
      if (method == "meta_method") {
        
        subfname <- paste0(simfolder, "/result/", method, "_formula", i.formula, "_", unique_key, sep="")
        res01 <- func(formula=method_formula, dat=dat01, subfname=subfname)
        
      } else {

        res01 <- func(formula=method_formula, dat=dat01)
        
        resfile_name <- paste0(simfolder, "/result/", method, "_formula", i.formula, "_", unique_key, ".rds", sep="")
        saveRDS(res01, file=resfile_name)
        
      }
      
      rm(res01)
      gc()
    }

  }
  
}






source(paste(simfolder, "/evaluation_methods.R", sep=""))
evals_func <- c("meanfunc", "varfunc", "cpfunc", "powerfunc", "biasfunc")
methods <- c("ipd_method", "propose_method", "meta_method")



metadat04 <- read.csv(file=paste(simfolder, "/meta_simdata.csv", sep=""))
head(metadat04)

dat01eval <- data.frame()

for(l in 1:length(metadat04$key)){
  params <- list(n=metadat04$n[l], sigma=metadat04$sigma[l], 
                 strata=metadat04$strata[l], nsim=metadat04$nsim[l], 
                 model=metadat04$model[l])
  unique_key <- metadat04$key[l]
  x2k_file_name <- paste0(simfolder, "/result/propose_method_formula1_", unique_key, ".rds", sep="")
  dat20 <- readRDS(x2k_file_name)
  dat20$params <- NULL
  x2k_m <- as.numeric(sapply(dat20, function(x) mean(x$data$x2k[x$data$strata==0 & x$data$x1k==1])))
  rm(dat20)
  gc()

  for(i.formula in formulas){
    if(i.formula == 1){
      method_formula <- yik ~ 1 + x1k
    }else{
      method_formula <- yik ~ 1 + x1k + x2k + x1k:x2k
    }
    
    for(method in methods){
      
      if (method == "meta_method") {
        
        subfname <- paste0(simfolder, "/result/", method, "_formula", i.formula, "_", unique_key, sep="")
        resfile_name <- paste0(subfname, "_ite", 1:10000, ".rds", sep="")
        dat11 <- lapply(resfile_name, readRDS)
        dat11$params <- dat11[[1]]$params

      } else {
        
        resfile_name <- paste0(simfolder, "result/", method, "_formula", i.formula, "_", unique_key, ".rds", sep="")
        dat11 <- readRDS(resfile_name)
        
      }
      
      dat11$x2k_m <- x2k_m
      
      
      for(evalfun in evals_func){
        mes01 <- paste(l, "/", length(metadat04$key), ", method:", method, ", formula:", i.formula, ", eval:", evalfun)
        message(mes01)
        
        if ( (method != "meta_method") | (method=="meta_method" & evalfun=="biasfunc") ) {
          
          func <- get(evalfun)
          evaluate <- func(dat11, method)
          eval01 <- list(params=params, eval=evalfun, value=evaluate)
          
          
          evalfile_name <- paste0(simfolder, "/evaluation/", method, "_", evalfun, "_", unique_key, ".rds", sep="")
          saveRDS(eval01, file=evalfile_name)
          
          paramd <- as.data.frame(params)[1,]
          row.names(paramd) <- NULL
          
          
          if(evalfun %in% c("meanfunc", "varfunc", "cpfunc")){
            dat02eval <- data.frame(paramd, method=method, formula=i.formula, eval=evalfun, mean=mean(evaluate), sd=sd(evaluate), 
                                    params=params, key=unique_key)
          }else if (evalfun == "powerfunc") {
            dat02eval <- data.frame(paramd, method=method, formula=i.formula, eval=evalfun, mean=evaluate$power, sd=evaluate$x, 
                                    params=params, key=unique_key)
          }else if (evalfun == "biasfunc") {
            dat02eval <- data.frame(paramd, method=method, formula=i.formula, eval=evalfun, mean=evaluate, sd=NA, 
                                    params=params, key=unique_key)
          }
          
          dat01eval <- rbind(dat01eval, dat02eval)
          
        }
        
      }
    }
    
  }

}

head(dat01eval)

evalall_file_name <- paste0(simfolder, "/total_evaluates.csv", sep="")
write.csv(dat01eval, file=evalall_file_name)





