  #### Simulation study
  tic <- Sys.time()
  setwd("/home/friedrich79/Documents/Forschung/Goettingen/COVID-19/Simu/Simu_risk_difference/Data")
  library(data.table)
  library(Publish)
  library(MatchIt)
  library(sandwich)
  library(lmtest)
  library(parallel)
    
  #### coefficients
  N <- 40
  b_trt <- 3.128  # marginal and conditional OR coincide
  iter <- 2000
  
 # beta_trt <- c(-1.4, 1.03, 0.04, 0.16, 0.07, 0.11)
  beta_trt <- c(-2.3, 0.31, 0.03, 1.099, -0.1054, 0.1031)
  
  #beta_out <- c(-1.06, 1.88, 0.01, -24.02, -26.14, 0.69, b_trt)
  beta_out <- c(-1.06, 0.619, 0.0077, 0.9461, -1.3499, 0.0896, b_trt)
  
  simu <- function(i, ...){
    
    ## Generate covariates
    x1 <- rbinom(N, size = 1, prob = 0.5)   # sex
    x2 <- round(rnorm(N, 45, 15))  # age
    x3 <- rbinom(N, size = 2, 0.5)  # clinical status
    x3.l <- ifelse(x3 == 1, 1, 0)
    x3.u <- ifelse(x3 == 2, 1, 0)
    x4 <- round(runif(N, 0, 10))   # time since onset
    
    covs.trt <- cbind(1, x1, x2, x3.l, x3.u, x4)
    
    ### Generate treatment status
    logit_trt <- covs.trt %*% beta_trt
    pi_trt <- exp(logit_trt)/(exp(logit_trt)+1)
    
    dt <- data.table(x1 = x1, x2 = x2, x3 = x3, x4 = x4)
    dt[, trt := rbinom(N, size = 1, prob = pi_trt)]
    
    ### Generate outcome conditional on treatment
    covs2 <- cbind(covs.trt, dt$trt)
    logit_out <- covs2 %*% beta_out
    pi_out <- exp(logit_out)/(exp(logit_out)+1)
    dt[, outcome := rbinom(N, size = 1, pi_out)]
    
    
    dt[, trt := factor(trt)]
    dt[, x1 := factor(x1)]
    dt[, x3 := factor(x3)]
    
   # source("ps-methods.R")
    
    #out <- try(ps_methods(dt))
    out <- NA
    return(list(out=out, data = dt))
  }
  
  set.seed(310320)
  sim <- lapply(1:iter, simu)
  
  data <- lapply(sim, function(x) x$data)
  datname <- paste0("new_OR10_Data_N", N, ".RDS")
  
  percent <- sapply(data, function(x) prop.table(table(x$trt))["1"])
  mean(percent)
  
  saveRDS(data, file = datname)
  
  # sim <- lapply(sim, function(x) x$out)
  # index <- lapply(sim, length)
  # sim[which(index!=5)] <- NULL
  # 
  # filename <- paste0("Simu_N", N, ".RDS")
  # 
  # saveRDS(sim, file = paste0("Results/OR1-Scen1/", filename))
  # 
  Sys.time() - tic
  
