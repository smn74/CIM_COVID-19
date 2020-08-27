#### function for calculating the different PS methods

ps_methods.a <- function(dt, ...){
  fit.adj <- glm(outcome ~ trt + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9, data = dt, family = binomial)
  tab <- regressionTable(fit.adj)
  tab.adj <- tab$trt[2, ]
 
  # ps model includes all nine variables 
  propmod <- glm(trt ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9,
                 family="binomial", data = dt)
  
  # propensity score
  dt[, ps := predict(propmod, type = "response")]
  
  ## second doubly robust model
  qt <- quantile(log(dt$ps/(1-dt$ps)), probs = seq(0, 1, 1/5))
  dt[, logitps := log(ps/(1-ps))]
  # create dummy variables
  dt[, z1:= as.numeric(logitps < qt[2])]
  dt[, z2:= as.numeric(qt[2] < logitps & logitps < qt[3])]
  dt[, z3:= as.numeric(qt[3] < logitps & logitps < qt[4])]
  dt[, z4:= as.numeric(qt[4] < logitps & logitps < qt[5])]
  
  ## CI using bootstrap
  source("causal_OR_with_Bootstrap_austin_2ndDR.R")
  boots <- bootstrap(dt)
  CIboot <- causal.or(boots$risks, boots$boot)
  
  tab.double2 <- data.frame(Variable = "", Units = "1", Missing = "", OddsRatio = CIboot["Estimate"], 
                                Lower = CIboot["Lower"], Upper = CIboot["Upper"] , Pvalue = NA)
  
  
  
  ## g-estimation w/o DR
  # source("causal_OR_with_Bootstrap_simple_austin.R")
  # boots2 <- bootstrap(dt)
  # CIboot2 <- causal.or(boots2$risks, boots2$boot)
  # 
  # tab.gcomp <- data.frame(Variable = "", Units = "1", Missing = "", OddsRatio = CIboot["Estimate"], 
  #                         Lower = CIboot["Lower"], Upper = CIboot["Upper"] , Pvalue = NA)
  # 
  table <- list(tab.adj,  tab.double2)
  names(table) <- c("adjusted", "DR2")
  
  return(table)
}