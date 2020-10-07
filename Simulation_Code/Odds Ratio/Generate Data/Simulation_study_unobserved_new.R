sim <- function(N, ...){  

#### Simulation study
  library(parallel)
  library(data.table)
  
  #### coefficients
  b_trt <- 0  # marginal and conditional OR coincide
  #b_trt <- 1.1111  # for true OR of 2
  #b_trt <- 3.4793  # for true OR of 10
  iter <- 2000
  
  beta_trt <- c(-2.3, 0.31, 0.03, 1.099, -0.1054, 0.1031, log(5))
  beta_out <- c(-1.06, 0.619, 0.0077, 0.9461, -1.3499, 0.0896, log(5), b_trt)
  
  simu <- function(i, ...){
    
    ## Generate covariates
    x1 <- rbinom(N, size = 1, prob = 0.5)   # sex
    x2 <- round(rnorm(N, 45, 15))  # age
    x3 <- rbinom(N, size = 2, 0.5)  # clinical status
    x3.l <- ifelse(x3 == 1, 1, 0)
    x3.u <- ifelse(x3 == 2, 1, 0)
    x4 <- round(runif(N, 0, 10))   # time since onset
    x5 <- rnorm(N)   # unobserved confounder, e.g. viral load at day0
    
    covs.trt <- cbind(1, x1, x2, x3.l, x3.u, x4, x5)
    
    ### Generate treatment status
    logit_trt <- covs.trt %*% beta_trt
    pi_trt <- exp(logit_trt)/(exp(logit_trt)+1)
    
    dt <- data.table(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5)
    dt[, trt := rbinom(N, size = 1, prob = pi_trt)]
    
    ### Generate outcome conditional on treatment
    covs2 <- cbind(covs.trt, dt$trt)
    logit_out <- covs2 %*% beta_out
    pi_out <- exp(logit_out)/(exp(logit_out)+1)
    dt[, outcome := rbinom(N, size = 1, pi_out)]
    
    dt[, trt := factor(trt)]
    dt[, x1 := factor(x1)]
    dt[, x3 := factor(x3)]
    
    return(list(out=NA, data = dt))
  }
  
  set.seed(310320)
  sim <- lapply(1:iter, simu)
  
  data <- lapply(sim, function(x) x$data)
  datname <- paste0("new_unobserved_Data_N", N, ".RDS")
  saveRDS(data, file = datname)
} 

Ns <- list(40, 100, 1000)
lapply(Ns, sim)