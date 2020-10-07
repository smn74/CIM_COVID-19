####
# function to calculate causal odds ratio including bootstrapped confidence intervals
#
# iter, seed: for bootstrapping

bootstrap <- function(data, iter = 1000, seed = 123, ...){
  
  # function to predict counterfactual outcomes using g-formula
  risks <- function(data, ...){
    # fit Q-model for Y~A + L
    fit <- glm(outcome ~  trt + z1+ z2 + z3 + z4 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9,
               data, family = binomial)
    
    predictions <- list()
    pot.outcome <- c("1", "0")
    for(i in 1:2){
      nd <- copy(data)
      nd[, outcome := NA]
      nd[, trt := pot.outcome[i]]
      predictions[[i]] <- predict(fit, newdata = nd, type = "response")
    
    }  
    names(predictions) <- pot.outcome
    effects <- sapply(FUN = mean, predictions)
    
    return(effects)
  }
  
  risks.original <- risks(data)
  
  # bootstrapping
  boot <- mclapply(X=1:iter, FUN=function(i, ...){
    bindex <- sample(1:NROW(data), size = NROW(data), replace = TRUE)
    bsample <- data[bindex]
    b.risks <- risks(bsample)
    return(b.risks)
  })
  
 output <- list(risks=risks.original, boot = boot)
  return(output)
}



causal.or <- function(risks, boot,  ...){
  
  Y_0 <- risks["0"]
  Y_1 <- risks["1"]
  
  causal.or <- (Y_1/(1-Y_1))/(Y_0/(1-Y_0))
  
  Y_0_b <- sapply(boot, function(x){x["0"]})  
  Y_1_b <- sapply(boot, function(x){x["1"]})
  
  ## count number of failures and remove NA
  fail <- length(which(is.na(Y_0_b)))
  if(fail != 0){
    Y_0_b <- as.numeric(Y_0_b[!is.na(Y_0_b)])
    Y_1_b <- as.numeric(Y_1_b[!is.na(Y_1_b)])
  }
  
  causal.or.b <- (Y_1_b/(1-Y_1_b))/(Y_0_b/(1-Y_0_b))
  
  q2.5 <- quantile(causal.or.b, prob = 0.025)
  q97.5 <- quantile(causal.or.b, prob = 0.975)
  
  COR <- c(q2.5, causal.or, q97.5, fail)
  names(COR) <- c("Lower", "Estimate", "Upper", "failed")
  return(COR)
}
