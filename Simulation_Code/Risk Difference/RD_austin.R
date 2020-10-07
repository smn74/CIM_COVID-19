## applying the risk difference calculations

rd.austin <- function(dt, ...){
  dt[, trt:= as.numeric(as.character(trt))]
  # ps model includes all nine variables 
  propmod <- glm(trt ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9,
                 family="binomial", data = dt)
  # propensity score
  dt[, ps := predict(propmod, type = "response")]
  
  # 1.) crude risk difference
  rd <- lm(outcome ~ trt, data = dt)
  est <- rd$coefficients["trt"]
  stderror <- coeftest(rd, vcov = vcovHC(rd, type = "HC3"))[2,2]
  ci.unadjusted <- c(est - 1.96*stderror, est, est + 1.96*stderror)
  names(ci.unadjusted) <- c("Lower", "Estimate", "Upper")
  
  
  # 2.) covariate adjustment
  lm.adj <- lm(outcome ~ trt + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9, data = dt)
  est.a <- lm.adj$coefficients["trt"]
  stderror.a <- coeftest(lm.adj, vcov = vcovHC(lm.adj, type = "HC3"))[2,2]
  ci.adj <- c(est.a - 1.96*stderror.a, est.a, est.a + 1.96*stderror.a)
  names(ci.adj) <- c("Lower", "Estimate", "Upper")
  
  
  ## 3.) PS covariate
  lm.PScov <- lm(outcome ~ trt + ps, data = dt)
  est.ps <- lm.PScov$coefficients["trt"]
  stderror.ps <- coeftest(lm.PScov, vcov = vcovHC(lm.PScov, type = "HC3"))[2,2]
  ci.PScov <- c(est.ps - 1.96*stderror.ps, est.ps, est.ps + 1.96*stderror.ps)
  names(ci.PScov) <- c("Lower", "Estimate", "Upper")
  
  ## 4.) Matching
  ps.mod <- tryCatch(matchit(trt ~ x1+ x2+x3+x4 + x5 +x6 +x7+x8+x9, data = dt,
                    method = "nearest", replace = FALSE, caliper = 0.2, distance = "linear.logit"),
                    error = function(e) NA,
                    warning = function(e) NA)
  if(!is.na(ps.mod)){
  m.data <- match.data(ps.mod)
  m.data <- as.data.table(m.data)
  # b: only treated experience the event
  b <- m.data[trt ==1 & outcome == 1, .N]
  # c: only untreated experience the event
  c <- m.data[trt ==0 & outcome == 1, .N]
  # n: number of matched pairs
  n <- (NROW(m.data)/2)
  
  diff.match <- (b-c)/n
  var <- ((b+c)-(c-b)^2/n)/n
  ci.match <- c(diff.match - 1.96*var, diff.match, diff.match + 1.96*var)
  } else {
    ci.match <- rep(NA, 3)
  }
  names(ci.match) <- c("Lower", "Estimate", "Upper")
  
  ## 5.) IPTW 
  dt[trt == 1, w:= 1/ps]
  dt[trt == 0, w:= 1/(1-ps)]
  
  msm <- svyglm(outcome ~ trt, design = svydesign(~1, weights = ~w, data = dt))
  est.msm <- msm$coefficients["trt"]
  stderror.msm <- coeftest(msm, vcov = vcovHC(msm, type = "HC3"))[2,2]
  ci.iptw <- c(est.msm - 1.96*stderror.msm, est.msm, est.msm + 1.96*stderror.msm)
  names(ci.iptw) <- c("Lower", "Estimate", "Upper")
  
  ## 6.) simple doubly robust g-computation
  dt[trt == 0, z:= -w]
  dt[trt == 1, z:= w]
  
  ## CI using bootstrap
  source("bootstrap_RD.R")
  boots <- bootstrap(dt, form = outcome ~  trt + z + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9)
  CIboot.rd <- causal.rd(boots$risks, boots$boot)
  
  ## 7.) 2nd DR
  qt <- quantile(log(dt$ps/(1-dt$ps)), probs = seq(0, 1, 1/5))
  dt[, logitps := log(ps/(1-ps))]
  
  # create dummy variables
  dt[, z1:= as.numeric(logitps < qt[2])]
  dt[, z2:= as.numeric(qt[2] < logitps & logitps < qt[3])]
  dt[, z3:= as.numeric(qt[3] < logitps & logitps < qt[4])]
  dt[, z4:= as.numeric(qt[4] < logitps & logitps < qt[5])]
  
  boots2 <- bootstrap(dt, form = outcome ~  trt + z1 + z2 + z3 + z4 +  x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9)
  CIboot.rd2 <- causal.rd(boots2$risks, boots2$boot)
  
  ## 8.) AIPW
  # source("psw-code.R")
  # form.ps <- "trt ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9"
  # form.out <- "outcome ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9"
  # tmp <- ps.aug( data = dt, form.ps = form.ps, weight = "ATE",form.outcome = form.out, family="binomial" )
  # ci.aipw <- c(tmp$est.risk - 1.96*tmp$std.risk, tmp$est.risk, tmp$est.risk + 1.96*tmp$std.risk)
  ci.aipw <- rep(NA, 3)
  names(ci.aipw) <- c("Lower", "Estimate", "Upper")
  
  table <- list(ci.unadjusted, ci.adj, ci.PScov, ci.match, ci.iptw, CIboot.rd, CIboot.rd2, ci.aipw)
  names(table) <- c("unadjusted", "adjusted", "cov", "matching", 
                    "IPTW", "g-computation", "DR2", "AIPW")
  return(table)
  
}