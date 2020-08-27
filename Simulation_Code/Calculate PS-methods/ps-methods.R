#### function for calculating the different PS methods

ps_methods <- function(dt, ...){
  fit.unadj <- glm(outcome ~ trt, data = dt, family = binomial)
  tab <- regressionTable(fit.unadj)
  tab.unadj <- tab$trt[2, ]
 
  # ps model includes all four variables 
  propmod <- glm(trt ~ x1 + x2 + x3 + x4,
                 family="binomial", data = dt)
  
  # propensity score
  dt[, ps := predict(propmod, type = "response")]
  
  ## 1.) covariate adjustment
  fit.c <- glm(outcome ~ trt + ps, data = dt, family = binomial)
  tab <- regressionTable(fit.c)
  tab.c <- tab$trt[2, ]
  
  ## 2.) Matching
  ps.mod <- matchit(trt ~ x1+ x2+x3+x4, data = dt,
                    method = "nearest", replace = FALSE)
  m.data <- match.data(ps.mod, distance = "pscore")
  # fit regression model on matched data (unconditional)
  fit.match <- glm(outcome ~ trt, data = m.data, family = binomial)
  tab <- regressionTable(fit.match)
  tab.match <- tab$trt[2, ]
  

  # calculate IPTW 
  dt[trt == 1, w:= 1/ps]
  dt[trt == 0, w:= 1/(1-ps)]
  
  # simple model using weighted data
  fit.iptw <- glm(outcome ~ trt, data = dt, weights = w, family = binomial)
  tt <- coeftest(fit.iptw, vcov = sandwich)
  
  tab.iptw <- data.frame(Variable = "", Units = "1", Missing = "", OddsRatio = exp(tt[2, 1]), 
                         Lower = exp(confint(tt))[2, 1], Upper = exp(confint(tt))[2, 2] , Pvalue = tt[2, 4])
  
  ### simple doubly robust model
  dt$z <- NA
  dt$z[dt$trt==0] <- -dt$w[dt$trt==0]
  dt$z[dt$trt==1] <- dt$w[dt$trt==1]
  
  ## CI using bootstrap
  source("causal_OR_with_Bootstrap.R")
  boots <- bootstrap(dt)
  CIboot <- causal.or(boots$risks, boots$boot)
  
  tab.iptw.double <- data.frame(Variable = "", Units = "1", Missing = "", OddsRatio = CIboot["Estimate"], 
                                Lower = CIboot["Lower"], Upper = CIboot["Upper"] , Pvalue = NA)
  
  
  table <- list(tab.unadj, tab.c, tab.match, tab.iptw, tab.iptw.double)
  names(table) <- c("unadjusted", "cov", "matching", 
                    "IPTW", "g-computation")
  #table <- as.data.table(table)
  
  return(table)
}