## calculate confidence intervals using bootstrap
###
# function to calculate causal odds ratio including bootstrapped confidence intervals
bootstrap <- function(data, form, iter = 10000, seed = 123, ...){
  
  # function to predict counterfactual outcomes for each hospital using g-formula
  risks <- function(data, form, ...){
    # fit Q-model for Y~A + L
    fit <- glm(form,
               data, family = binomial)
    
    predictions <- list()
    pot.outcome <- c(1, 0)
    for(i in 1:2){
      nd <- copy(data)
      nd[, outcome:=NA]
      nd[, trt := pot.outcome[i]]
      predictions[[i]] <- predict(fit, newdata = nd, type = "response")
      
    }  
    names(predictions) <- pot.outcome
    effects <- sapply(FUN = mean, predictions)
    
    return(effects)
  }
  
  risks.original <- risks(data, form)
  
  # bootstrapping
  boot <- mclapply(X=1:iter, FUN=function(i, ...){
    bindex <- sample(1:NROW(data), size = NROW(data), replace = TRUE)
    bsample <- data[bindex]
    b.risks <- risks(bsample, form)
    return(b.risks)
  })
  
  output <- list(risks=risks.original, boot = boot)
  return(output)
}


# function for calculating risk difference instead of OR
causal.rd <- function(risks, boot, ...){
  
  Y_0 <- risks["0"]
  Y_1 <- risks["1"]
  
  causal.rd <- Y_1 -Y_0
  
  Y_0_b <- sapply(boot, function(x){x["0"]})  
  Y_1_b <- sapply(boot, function(x){x["1"]})
  
  causal.rd.b <- Y_1_b-Y_0_b
  
  q2.5 <- quantile(causal.rd.b, prob = 0.025)
  q97.5 <- quantile(causal.rd.b, prob = 0.975)
  
  CRD <- c(q2.5, causal.rd, q97.5)
  names(CRD) <- c("Lower", "Estimate", "Upper")
  return(CRD)
}
