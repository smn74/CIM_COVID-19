#### function for calculating the different PS methods

ps_methods <- function(dt, ...){
  ## 2.) Matching
  ps.mod <- matchit(trt ~ x1+ x2+x3+x4, data = dt,
                    method = "nearest", replace = FALSE, caliper = 0.2)
  m.data <- match.data(ps.mod, distance = "pscore")
  # fit regression model on matched data (unconditional)
  fit.match <- glm(outcome ~ trt, data = m.data, family = binomial)
  tab <- regressionTable(fit.match)
  tab.match <- tab$trt[2, ]
  
  # adjusting for matched pairs
  m.data$id <- rownames(m.data)
  m.data <- as.data.table(m.data)
  rn <- ps.mod$match.matrix[!is.na(ps.mod$match.matrix),]
  m.data[id %in% names(rn), matched.id := rn]
  m.data[is.na(matched.id), class:= id]
  m.data[!is.na(matched.id), class:= matched.id]

  fit.match2 <- tryCatch(clogit(outcome ~ trt + strata(class), data = m.data),
                         error = function(e) NA, warning = function(e) NA)
  if(!is.na(fit.match2)){
    tab <- regressionTable(fit.match2)
    tab.match.clogit <- tab$trt[2, ]
  } else tab.match.clogit <- NA

  # ### GEE
  library(geepack)
  gee.matched <- geeglm(formula   = outcome ~ trt,
                        family    = binomial(link = "logit"),
                        id        = class,
                        data      = m.data,
                        corstr    = "exchangeable",
                        scale.fix = FALSE,
                        std.err = "jack"
  )
  #tab <- regressionTable(gee.matched)
  confint.geeglm <- function(object, parm, level = 0.95, ...) {
    cc <- coef(summary(object))
    mult <- qnorm((1+level)/2)
    citab <- with(as.data.frame(cc),
                  cbind(lwr=Estimate-mult*Std.err,
                        upr=Estimate+mult*Std.err))
    rownames(citab) <- rownames(cc)
    citab[parm,]
  }

  tab.match.gee <- data.frame(Units = "1", OddsRatio = exp(gee.matched$coefficients[2]),
                         Lower = exp(confint.geeglm(gee.matched, 2)[1]),
                         Upper = exp(confint.geeglm(gee.matched, 2)[2]))


  table <- list(tab.match, tab.match.clogit, tab.match.gee)
  names(table) <- c("matching unadjusted", "matching conditional", "matching gee")
  #table <- as.data.table(table)
  
  return(table)
}