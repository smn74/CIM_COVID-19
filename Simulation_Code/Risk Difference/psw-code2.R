# Code taken from PSW
ps.aug <- function( dt, form.ps, weight = "ATE",form.outcome, family="binomial"){

  dt[, x1:=as.numeric(as.character(x1))]
  # x3 must be dummy coded
  dt[, x3.l := ifelse(x3 == 1, 1, 0)]
  dt[, x3.u := ifelse(x3 == 2, 1, 0)]
  
  out.ps <- PSW:::ps.model(dat = dt, form.ps = as.formula(form.ps))
  trt.var <- as.character(terms(out.ps$fm)[[2]])
  Xname <- names(coef(out.ps$fm))[-1]
  ps.hat <- out.ps$ps.hat
  beta.hat <- as.numeric(coef(out.ps$fm))
  omega <- sapply(ps.hat, PSW:::calc.omega, weight = weight, delta = 0.002, 
                  K = K)
  Q <- dt[, trt] * ps.hat + (1 - dt[, trt]) * (1 - ps.hat)
  W <- omega/Q
  res <- list(weight = weight, ps.model = out.ps$fm, ps.hat = ps.hat, 
              W = W)
  out.outcome <- PSW:::outcome.model(dat = dt, form = as.formula(form.out), 
                               trt.var = trt.var, family = "binomial")
  out.var <- as.character(terms(out.outcome$fm1)[[2]])
  
  psw.aug.core <- function(dat, beta.hat, 
                           omega, Q, out.ps, out.outcome, 
                           trt.var, out.var, weight, 
                           family, K, delta){
    n <- nrow(dat)
    W <- omega/Q
    tmp1 <- W * dat[, trt]
    mu1.hat <- sum(tmp1 * (dat[, outcome] - out.outcome$Y.hat.m1))/sum(tmp1)
    mu2.hat <- sum(omega * out.outcome$Y.hat.m1)/sum(omega)
    tmp0 <- W * (1 - dat[, trt])
    mu3.hat <- sum(tmp0 * (dat[, outcome] - out.outcome$Y.hat.m0))/sum(tmp0)
    mu4.hat <- sum(omega * out.outcome$Y.hat.m0)/sum(omega)
    
    p1.hat <- mu1.hat + mu2.hat
    p0.hat <- mu3.hat + mu4.hat
    est.risk <- p1.hat - p0.hat
    
    alpha1.hat <- as.numeric(coef(out.outcome$fm1))
    alpha0.hat <- as.numeric(coef(out.outcome$fm0))
    beta.hat <- as.numeric(coef(out.ps$fm))
    n.alpha1 <- length(alpha1.hat)
    n.alpha0 <- length(alpha0.hat)
    n.beta <- length(beta.hat)
    n <- nrow(dat)
    Amat <- Bmat <- 0
    for (i in 1:n) {
      Xi <- unlist(c(1, dat[i, .(x1, x2 ,x3.l, x3.u, x4)]))
      Vi <- unlist(c(1, dat[i, .(x1, x2 ,x3.l, x3.u, x4)]))
      Zi <- dat[i, trt]
      Yi <- dat[i, outcome]
      Wi <- W[i]
      Qi <- Q[i]
      omegai <- omega[i]
      ei <- PSW:::calc.ps.Xbeta(Xmat = Xi, beta = beta.hat)
      ei.deriv1 <- PSW:::calc.ps.deriv1(Xmat = Xi, beta = beta.hat)
      ei.deriv2 <- PSW:::calc.ps.deriv2(Xi = Xi, beta = beta.hat)
      Wi.deriv.beta <- PSW:::calc.W.derive.beta(Zi = Zi, Xi = Xi, 
                                                omega.ei = omegai, beta.hat = beta.hat, ei = ei, 
                                                Qi = Qi, weight = weight, delta = delta, K = K)
      omegai.deriv.beta <- PSW:::calc.omega.derive.beta(Xi = Xi, 
                                                        beta.hat = beta.hat, ei = ei, weight = weight, delta = delta, 
                                                        K = K)
      tmp1 <- sum(Vi * alpha1.hat)
      tmp0 <- sum(Vi * alpha0.hat)
      m1.hat <- exp(tmp1)/(1 + exp(tmp1))
      m0.hat <- exp(tmp0)/(1 + exp(tmp0))
      m1.deriv.alpha1 <- m1.hat * (1 - m1.hat) * Vi
      m0.deriv.alpha0 <- m0.hat * (1 - m0.hat) * Vi
      s1.deriv.alpha1 <- -m1.hat * (1 - m1.hat) * outer(Vi, 
                                                        Vi)
      s0.deriv.alpha0 <- -m0.hat * (1 - m0.hat) * outer(Vi, 
                                                        Vi)
      
      this.phi.row1 <- Wi * Zi * (Yi - m1.hat - mu1.hat)
      this.phi.row2 <- omegai * (m1.hat - mu2.hat)
      this.phi.row3 <- Wi * (1 - Zi) * (Yi - m0.hat - mu3.hat)
      this.phi.row4 <- omegai * (m0.hat - mu4.hat)
      this.phi.row5 <- Zi * (Yi - m1.hat) * Vi
      this.phi.row6 <- (1 - Zi) * (Yi - m0.hat) * Vi
      this.phi.row7 <- (Zi - ei)/ei/(1 - ei) * ei.deriv1
      this.phi <- c(this.phi.row1, this.phi.row2, this.phi.row3, 
                    this.phi.row4, this.phi.row5, this.phi.row6, this.phi.row7)
      Bmat <- Bmat + outer(this.phi, this.phi)
      quad1 <- diag(c(-Wi * Zi, -omegai, -Wi * (1 - Zi), -omegai))
      quad2 <- matrix(0, nrow = n.alpha1 + n.alpha0 + n.beta, 
                      ncol = 4)
      tmp <- c(-Wi * Zi * m1.deriv.alpha1, rep(0, n.alpha0), 
               Zi * (Yi - m1.hat - mu1.hat) * Wi.deriv.beta, omegai * 
                 m1.deriv.alpha1, rep(0, n.alpha0), (m1.hat - 
                                                       mu2.hat) * omegai.deriv.beta, rep(0, n.alpha1), 
               -Wi * (1 - Zi) * m0.deriv.alpha0, (1 - Zi) * (Yi - 
                                                               m0.hat - mu3.hat) * Wi.deriv.beta, rep(0, n.alpha1), 
               omegai * m0.deriv.alpha0, (m0.hat - mu4.hat) * omegai.deriv.beta)
      quad3 <- matrix(tmp, byrow = TRUE, nrow = 4, ncol = n.alpha1 + 
                        n.alpha0 + n.beta)
      quad4.blk1 <- cbind(Zi * s1.deriv.alpha1, matrix(0, nrow = n.alpha1, 
                                                       ncol = n.alpha0 + n.beta))
      quad4.blk2 <- cbind(matrix(0, nrow = n.alpha0, ncol = n.alpha1), 
                          (1 - Zi) * s0.deriv.alpha0, matrix(0, nrow = n.alpha0, 
                                                             ncol = n.beta))
      quad4.blk3 <- cbind(matrix(0, nrow = n.beta, ncol = n.alpha1 + 
                                   n.alpha0), -ei * (1 - ei) * outer(Xi, Xi))
      quad4 <- rbind(quad4.blk1, quad4.blk2, quad4.blk3)
      this.phi.deriv <- rbind(cbind(quad1, quad3), cbind(quad2, 
                                                         quad4))
      Amat <- Amat + this.phi.deriv
    }
    Amat <- Amat/n
    Bmat <- Bmat/n
    Amat.inv <- solve(Amat)
    var.mat <- (Amat.inv %*% Bmat %*% t(Amat.inv))/n
    var.mat <- var.mat[c(1:4), c(1:4)]
    tmp.risk <- c(1, 1, -1, -1)
    var.risk <- PSW:::rowVec(tmp.risk) %*% var.mat %*% PSW:::colVec(tmp.risk)
    std.risk <- sqrt(as.numeric(var.risk))
    ans <- list(est.risk = est.risk, std.risk = std.risk)
    return(ans)
  }
  
  
  res.aug <- psw.aug.core(dat = dt, beta.hat = beta.hat, 
                          omega = omega, Q = Q, out.ps = out.ps, out.outcome = out.outcome, 
                          trt.var = trt.var, out.var = out.var, weight = weight, 
                          family = "binomial", K = K, delta = 0.002)
  
 return(res.aug)
}