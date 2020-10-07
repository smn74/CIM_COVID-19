## applying the risk difference calculations

rd.austin <- function(dt, ...){
  dt[, trt:= as.numeric(as.character(trt))]
  cols <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9")
  setDT(dt)[, (cols):= lapply(.SD, function(x) as.numeric(as.character(x))), .SDcols=cols]
  
  ## 8.) AIPW
  source("aipw/psw-code.R")
  form.ps <- "trt ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9"
  form.out <- "outcome ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9"
  tmp <- ps.aug(dt, form.ps, weight = "ATE", form.out, family="binomial" )
  ci.aipw <- c(tmp$est.risk - 1.96*tmp$std.risk, tmp$est.risk, tmp$est.risk + 1.96*tmp$std.risk)
  names(ci.aipw) <- c("Lower", "Estimate", "Upper")
  
  # table <- list(ci.unadjusted, ci.adj, ci.PScov, ci.match, ci.iptw, CIboot.rd, CIboot.rd2, ci.aipw)
  # names(table) <- c("unadjusted", "adjusted", "cov", "matching", 
  #                   "IPTW", "g-computation", "DR2", "AIPW")
  return(ci.aipw)
  
}