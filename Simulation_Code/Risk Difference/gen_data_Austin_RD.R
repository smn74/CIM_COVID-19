#### Simulation study
tic <- Sys.time()
#setwd("/home/friedrich79/Documents/Forschung/Goettingen/COVID-19/Simu")
library(data.table)
library(Publish)
library(parallel)
  
###
# Simulation scenario as in Austin, full PS model, 9 binary covariates

#### coefficients
N <- 1000
b_trt <- 1.032  # true RD 0.16 
iter <- 2000

# beta_{0, trt}, beta_1 - beta_6 (Notation Austin)
beta_trt <- c(-3.5, log(5), log(2), log(5), log(2), log(5), log(2))

# coefficients alpha_{0, outcome}, beta_trt, alpha_1 - alpha_6 (Notation Austin)
beta_out <- c(-5, b_trt, rep(log(5), 3), rep(log(2), 3))

simu <- function(i, ...){
  
  ## Generate covariates
  dt <- matrix(NA, nrow = N, ncol = 9)
  covs <- sapply(1:9, function(x) dt[, x] <- rbinom(N, size = 1, prob = 0.5))
  colnames(covs) <- paste0("x", 1:9)
    
  covs.trt <- cbind(1, covs[, c(1, 2, 4, 5, 7, 8)])
  
  ### Generate treatment status
  logit_trt <- covs.trt %*% beta_trt
  pi_trt <- exp(logit_trt)/(exp(logit_trt)+1)
  
  dt <- data.table(covs)
  dt[, trt := rbinom(N, size = 1, prob = pi_trt)]
  
  ### Generate outcome conditional on treatment
  covs2 <- cbind(1, dt$trt, covs[, 1:6])
  logit_out <- covs2 %*% beta_out
  pi_out <- exp(logit_out)/(exp(logit_out)+1)
  dt[, outcome := rbinom(N, size = 1, pi_out)]
  
  cols <- colnames(dt)[-NCOL(dt)]
  
  dt[, (cols):=lapply(.SD, as.factor), .SDcols = cols]
  
  #source("ps-methods_austin.R")
  
  #out <- try(ps_methods(dt))
  out <- NA
  return(list(out=out, data = dt))
}

set.seed(310320)
sim <- lapply(1:iter, simu)

data <- lapply(sim, function(x) x$data)
datname <- paste0("Austin_OR2_Data_N", N, ".RDS")

percent <- sapply(data, function(x) prop.table(table(x$trt))["1"])
mean(percent)

saveRDS(data, file = paste0("/home/friedrich79/Documents/Forschung/Goettingen/COVID-19/Simu/Simu_risk_difference/Data/", datname))

# sim <- lapply(sim, function(x) x$out)
# index <- lapply(sim, length)
# sim[which(index!=5)] <- NULL
# 
# filename <- paste0("Austin_Simu_N", N, ".RDS")
#saveRDS(sim, file = paste0("Results/OR1-Scen3/", filename))

#Sys.time() - tic
