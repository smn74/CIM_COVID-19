## applying the risk difference calculations

ps_methods <- function(dt, ...){
  dt[, trt:= as.numeric(as.character(trt))]
  
  ## 8.) AIPW
  form.ps <- "trt ~ x1 + x2 + x3 + x4"
  form.out <- "outcome ~ x1 + x2 + x3 + x4"
  source("aipw/psw-code2.R")
  tmp <- ps.aug(dt, form.ps, weight = "ATE",form.out, family="binomial" )
  ci.aipw <- c(tmp$est.risk - 1.96*tmp$std.risk, tmp$est.risk, tmp$est.risk + 1.96*tmp$std.risk)
  names(ci.aipw) <- c("Lower", "Estimate", "Upper")

  return(ci.aipw)
  
}
