#### Simulation study
library(data.table)
library(Publish)
library(MatchIt)
library(sandwich)
library(lmtest)
library(parallel)
  
###
# Simulation scenario as in Austin, full PS model, 9 binary covariates

#### coefficients
N <- 1000 # similar for 40, 100
b_trt <- 0  # marginal and conditional OR coincide
# 0.9707 for true OR of 2, 3.2625 for true OR of 10
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
  
  return(list(out=NA, data = dt))
}

set.seed(310320)
sim <- lapply(1:iter, simu)

data <- lapply(sim, function(x) x$data)
datname <- paste0("Austin_Data_N", N, ".RDS")
saveRDS(data, file = paste0("Results/OR1-Scen3/", datname))
