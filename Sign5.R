# This is the 5th version of Sign Identification functions file

# ---- 1. General VAR and other preliminary functions ----
# ==== 1.1. Reduced VAR estimation ====
# #### 1.1.2 Function for data preparation ####
specify.var <- function(endog, exog = NULL, num.lags = 1, det.part = 'const') {
  check.ts <- function(Data) {
    temp <- complete.cases(Data)
    which.not.cc <- which(temp)
    cond <- sum(temp)/
      (which.not.cc[length(which.not.cc)] - 
         which.not.cc[1] + 1) != 1
    return(cond)
  }
  if (check.ts(endog)) {
    stop('There should be no NAs in the middle of the endog sequence')
  }
  if (!is.null(exog)) {
    if (check.ts(exog)) {
      stop('There should be no NAs in the middle of the exog sequence')
    }
  }
  if (is.null(colnames(endog))) {colnames(endog) <- paste0('v', 1:ncol(endog))}
  
  Y <- as.matrix(endog)[1:min(nrow(endog), nrow(exog)), ]
  Z <- matrix(NA, nrow = nrow(Y), ncol = 0)
  
  if ('const' %in% det.part)      { Z <- cbind(Z, rep(1, nrow(Y))) } 
  if ('trend' %in% det.part)      { Z <- cbind(Z, 1:nrow(Y)) } 
  if ('quad.trend' %in% det.part) { Z <- cbind(Z, (1:nrow(Y))^2) } 
  Z.names <- c('const', 'trend', 'quad.trend')
  Z.names <- Z.names[Z.names %in% det.part]
  
  if (!is.null(exog)) {
    Z <- cbind(Z, exog); Z.names <- c(Z.names, colnames(exog))
    exog.names <- colnames(exog)
  } else { exog.names <- NULL }
  
  if (num.lags > 0) {
    for (i in 1:num.lags) {
      lags <- rbind(matrix(NA, nrow = i, ncol = ncol(Y)), Y[1:(nrow(Y) - i), ])
      Z <- cbind(Z, lags)
      Z.names <- c(Z.names, paste0('L', i, '.', colnames(Y)))
    }
  }
  colnames(Z) <- Z.names
  
  return(list(Y = Y[complete.cases(Y) & complete.cases(Z), ],
              Z = Z[complete.cases(Y) & complete.cases(Z), ],
              num.lags = num.lags, exog.names = exog.names,
              det.part = det.part, exog.vars = exog, 
              endog.vars = endog))
}

# #### 1.1.2 OLS function ####
# - Function, made to gain consistency with other fitting functions
fit.ols.var <- function(spec.var, method = 'qr') {
  Y <- spec.var$Y; Z <- spec.var$Z
  if (method == 'qr') {
    dec <- qr(Z)
    R <- qr.R(dec)
    Phi  <- solve(R) %*% t(qr.Q(dec)) %*% Y
    M.ZZ <- solve(t(R) %*% R)
  } else if (method == 'direct') {
    Phi  <- solve(t(Z) %*% Z) %*% t(Z) %*% Y
    M.ZZ <- solve(t(Z) %*% Z)
  }
  
  resid <- Y - Z %*% Phi
  rss   <- colSums(resid^2)
  df    <- nrow(Z) - ncol(Z)
  sigma.sq <- rss/df
  Var.Phi <-lapply(1:ncol(Y), function(i) {sigma.sq[i]*M.ZZ})
  names(Var.Phi) <- colnames(Y)
  Sigma <- t(resid) %*% resid/df
  
  return(list(Phi = Phi, Var.Phi = Var.Phi, Sigma = Sigma, 
              det.part = spec.var$det.part, sigma.sq = sigma.sq,
              rss = rss, df = df, resid = resid,
              spec.var = spec.var))
}

fit.ols.restricted.var <- function(spec.var, which.restricted) {
  num.endog <- ncol(spec.var$Y)
  Phi <- vector('list', num.endog)
  names(Phi) <- colnames(spec.var$Y)
  Z.temp <- Phi
  resid <- matrix(NA, ncol = num.endog, nrow = nrow(spec.var$Y))
  df <- rep(NA,  spec.var$num.lags)
  
  for (i in 1:num.endog) {
    Z.temp[[i]] <- spec.var$Z[, is.na(which.restricted[, i])]
    Phi[[i]]    <- solve(t(Z.temp[[i]]) %*% Z.temp[[i]]) %*% 
      t(Z.temp[[i]]) %*% spec.var$Y[, i]
    resid[, i]  <- spec.var$Y[, i] - Z.temp[[i]] %*% Phi[[i]]
    df[i] <- nrow(spec.var$Y) - length(Phi[[i]])
  }
  
  rss  <- colSums(resid^2)
  sigma.sq <- rss/df
  return(list(Phi = Phi, Z.temp = Z.temp, resid = resid,
              df = df, sigma.sq = sigma.sq, spec.var = spec.var))
}

# #### 1.1.2 Rearrange fit function ####
# reorder.fit <- function(fit, new.order) {
#   num.det <- length(fit$det.part)
#   num.lags <- (nrow(fit$Phi) - num.det )/ncol(fit$Phi)
#   num.endog <- ncol(fit$Phi)
#   Phi.order <- c(1:num.det,  sapply(1:num.lags, function(l) {
#     num.det + l * num.endog + new.order
#     }))
#   
#   fit$Phi <- fit$Phi[Phi.order, new.order]
#   fit$Var.Phi <- fit$Var.Phi[new.order]
#   for (i in seq_along())
# }

# ==== 1.2. IRF, FEVD, Historical Decompositions ====
# #### 1.2.1 IRF calculation ####
calc.irf <- function(structural.params, exog.names = NULL, det.part = 'const',
                     which.shock = NULL, is.continuous = FALSE, 
                     horizon = 1e2, calc.horizon = horizon,
                     calc.long.run = TRUE, previous.irf = NULL) {
  
  # NB! Be sure that cols and rows of A have the same order
  A <- structural.params$A
  inv.A <- solve(A)
  if (!is.null(exog.names)) {
    C <- structural.params$B[exog.names, , drop = FALSE]
  }
  which.endog <- !(rownames(structural.params$B) %in% c(det.part, exog.names))
  B <- structural.params$B[which.endog, , drop = FALSE]
  
  num.lags <- nrow(B)/nrow(A)
  
  if (is.null(which.shock)) {
    which.shock <- rownames(A)
  }
  if (is.null(previous.irf)) {
    irf <- matrix(NA, nrow = horizon, ncol = nrow(A)^2)
    colnames(irf) <- c(sapply(rownames(A), function(x) {
      paste0('shock.', x, '_resp.', colnames(A))}))
    start.rep <- 1
  } else {
    irf <- previous.irf$irf
  }
  
  for (x in which.shock) {
    k <- which(colnames(A) == x) - 1
    # Which irf vector should be used
    irf.ind <- k * ncol(A) + 1:ncol(A)
    
    if (x %in% colnames(A)) {
      e0 <- matrix(as.numeric(colnames(A) == x), nrow = 1) %*% inv.A
    } else if (x %in% exog.names) {
      e0 <- matrix(as.numeric(exog.names == x), nrow = 1) %*% C %*% inv.A
    }
    
    # What the first step of IRF should be calculated?
    if (!is.null(previous.irf)) {
      which.irf.made <- which(cumsum(complete.cases(irf[, irf.ind]))/
                                cumsum(rep(1, horizon)) == 1)
      if (length(which.irf.made)) {
        start.rep <- 1
      } else {
        start.rep <- length(which.irf.made) + 1
      }
    }
    
    for (i.rep in start.rep:calc.horizon) {
      if (i.rep == 1) {
        irf[i.rep, irf.ind] <- e0
      } else {
        if (i.rep <= num.lags) {
          y.lag <- matrix(0, nrow = 1, ncol = num.lags*ncol(A))
          y.lag[, 1:(ncol(A)*(i.rep - 1))] <- 
            c(t(irf[(i.rep - 1):1, irf.ind]))
        } else {
          y.lag <- matrix(c(t(irf[i.rep - num.lags - 1  + 
                                    num.lags:1, irf.ind])), nrow = 1)
        }
        
        if (is.continuous) {
          irf[i.rep, irf.ind] <- y.lag %*% B %*% inv.A + e0
        } else {
          irf[i.rep, irf.ind] <- y.lag %*% B %*% inv.A 
        }
      }
    }
  }
  
  # Calculate LR IRFs
  if (calc.long.run) {
    sum.mat <- function(B, num.lags) {
      num.var <- nrow(B)/num.lags
      if (num.lags == 1) {
        return(B) 
      } else{
        B.sum <- B[1:num.var, ]
        for (s in 1:(num.lags - 1)) {
          B.sum <- B.sum + B[s*num.var + 1:num.var, ]
        }
        return(B.sum)
      }
    }
    
    if (x %in% exog.names) {
      lr.irf <- C %*% solve(structural.params$A - sum.mat(B, num.lags))
      rownames(lr.irf) <- paste0('shock.', rownames(C))
    } else{
      lr.irf <- solve(structural.params$A - sum.mat(B, num.lags))
      rownames(lr.irf) <- paste0('shock.', rownames(A))
    }
    colnames(lr.irf) <- paste0('resp.', colnames(A))
    return(list(irf = irf, lr.irf = lr.irf))
  } else{ return(irf) }
}

# #### 1.2.2 Forecast error variance decomposition calculation ####
rearrange.list <- function(x, type = 'main') {
  if (type == 'main') {
    n.draws <- nrow(x[[1]])
    horizon <- ncol(x[[1]])
    temp <- matrix(NA, ncol = length(x), nrow = n.draws)
    res <- vector('list', horizon)
    
    for (h in 1:horizon) {
      res[[h]] <- temp
      for ( i in 1:n.draws) {
        for (j in 1:length(x)) {
          res[[h]][i, j] <- x[[j]][i , h]
        }
      }
      colnames(res[[h]]) <- names(x)
    }
  } else if (type == 'inverse') {
    horizon <- length(x)
    n.draws <- nrow(x[[1]])
    num.shock.resp  <- ncol(x[[1]])
    
    res <- lapply(1:num.shock.resp, function(i) {
      res <- matrix(NA, nrow = n.draws, ncol = horizon)
    })
    names(res) <- colnames(x[[1]])
    
    for (i in 1:num.shock.resp) {
      for (h in 1:horizon) {
        for (j in 1:n.draws) {
          res[[i]][j, h] <- x[[h]][j, i]
        }
      }
    }
  }
  
  return(res)
}

calc.fevd <- function(irf.draws, forecast.horizon, horizon, d = NULL,
                      as.share = TRUE) {
  num.shocks <- sqrt(length(irf.draws))
  n.draws    <- nrow(irf.draws[[1]])
  sample.horizon    <- ncol(irf.draws[[1]])
  
  if (is.null(d)) {d <- rep(1, num.shocks)}
  
  irf <- rearrange.list(irf.draws)
  irf.mult <- lapply(1:(horizon + forecast.horizon), function(i) {
    res <- matrix(NA, ncol = ncol(irf[[1]]), nrow = nrow(irf[[1]]))
    colnames(res) <- colnames(irf[[1]])
    return(res)
  })
  fevd <- lapply(1:horizon, function(i) {
    res <- matrix(0, ncol = ncol(irf[[1]]), nrow = nrow(irf[[1]]))
    colnames(res) <- colnames(irf[[1]])
    return(res)
  })
  
  for (i in 1:n.draws) {
    for (h in 1:(horizon + forecast.horizon)) {
      for (s in 1:num.shocks) {
        seq <- (s - 1)*num.shocks + 1:num.shocks
        irf.mult[[h]][i, seq] <- d[s] * diag(t(irf[[h]][i , seq, drop = FALSE]) %*% 
                                    irf[[h]][i , seq, drop = FALSE])
      }
    }
  }
  
  for (h in 1:horizon) {
    for (f.h in 1:forecast.horizon) {
      if (h + f.h <= sample.horizon) {
        fevd[[h]] <- fevd[[h]] + irf.mult[[h + f.h]]
      }
    }
  }
  
  if (as.share) {
    for (h in 1:horizon) {
      for (i in 1:n.draws) {
        fevd[[h]][i, ]  <- fevd[[h]][i, ] / sum(fevd[[h]][i, ])
      }
    }
  }
  return(rearrange.list(fevd, type = 'inverse'))
}

# #### 1.2.3 Historical decomposition calculation ####
calc.hist.dec <- function(irf, resid, d = NULL, horizon = nrow(fit$spec.var$Y)) {
  
}

# ---- 2. Estimation methods ----
# ==== 2.1. Prior elicitation ====
fit.univar.ar <- function(spec.var, num.lags = NULL, det.part = NULL) {
  if (is.null(num.lags)) num.lags <- spec.var$num.lags
  if (is.null(det.part)) {
    det.part <- spec.var$det.part
  } else if (det.part == 'none') {
    det.part <- NULL
  }
  
  which.restricted <- matrix(0, nrow = ncol(spec.var$Z),
                             ncol = ncol(spec.var$Y))
  rownames(which.restricted) <- colnames(spec.var$Z)
  colnames(which.restricted) <- colnames(spec.var$Y)
  if (length(det.part) > 0) {
    which.restricted[det.part, ] <- NA
  }
  for (l in 1:num.lags) {
    diag(which.restricted[length(spec.var$det.part) + (l - 1) * ncol(spec.var$Y) +
                            1:ncol(spec.var$Y), ]) <- NA
  }
  
  return(fit.ols.restricted.var(spec.var, which.restricted))
}

define.prior.classical.sims.zha <- function(num.equations, num.lags = 1,
                                           overall.tightness = 0.5,
                                           relative.tightness = 0.2,
                                           deterministic.tightness = 1e3,
                                           lag.decay.rate = 2, 
                                           assume.random.walk = TRUE,
                                           s = NULL, det.part = 'const',
                                           exog.names = NULL,
                                           endog.names = NULL,
                                           spec.var = NULL) {
  # NB! Usually SVAR is normalized to s_r = 1, but it is not always the case;
  # If this assumption is violated, the s vector should be defined
  if(is.null(endog.names)) {endog.names <- paste0('v', 1:num.equations)}
  if (is.null(s)) {
    # s <- rep(1, length(endog.names))
    s <- sqrt(fit.univar.ar(spec.var, 1, det.part)$sigma.sq)
    }
  start <- length(det.part) + length(exog.names) 
  
  B.prior <- matrix(0, ncol = length(endog.names), nrow = 
                      length(det.part) + length(exog.names) +
                      num.lags*length(endog.names))
  rownames(B.prior) <- c(det.part, exog.names,
                         sapply(1:num.lags, function(l) {
                           paste0('L', l, '.', endog.names)
                         }))
  if (length(assume.random.walk) == 1) {
    if (assume.random.walk) {diag(B.prior[(start + 1):nrow(B.prior), ]) <- 1}
  } else {
    diag(B.prior[(start + 1):nrow(B.prior), ]) <- as.numeric(assume.random.walk)
  }
  
  colnames(B.prior) <- endog.names
  V.B.prior <- matrix(0, ncol = nrow(B.prior), nrow = nrow(B.prior))
  rownames(V.B.prior) <- rownames(B.prior)
  colnames(V.B.prior) <- rownames(B.prior)
  
  if (length(det.part) == 1) {
    V.B.prior[1, 1] <- overall.tightness * deterministic.tightness
  } else {
    diag(V.B.prior[1:length(det.part), ]) <- 
      overall.tightness * deterministic.tightness
  }
  
  for (l in 1:num.lags) {
    dim <- start + (l - 1)*length(endog.names) + 1:length(endog.names)
    diag(V.B.prior[dim, dim]) <- 
      overall.tightness * relative.tightness / (s * l^lag.decay.rate)
  }
  return(list(B.prior = B.prior, V.B.prior = V.B.prior))
}

# ==== 2.2. Estimation functions ====
specify.method <- function(inference.method) {
  if (inference.method == 'bootstrap') {
    # Arefiev's function for block bootstraping
    block.sample <- function(m, size = m, block.length = round(m^(1/3), 0)) {
      ret <- sample.int(m, size = size, replace = TRUE)
      for (i in 2:m) {
        if ((i - 1) %% block.length == 0) next
        ret[i] <- ret[i - 1] + 1
        if (ret[i] > m) ret[i] <- ret[i] - m
      }
      
      ret
    }
    
    # Residuals -> Y.sampled:
    create.endog.sample <- function(fit) {
      num.endog <- ncol(fit$Phi)
      num.lags <- (nrow(fit$Phi) - sum(rownames(fit$Phi) %in% 
                                         c(fit$spec.var$det.part, 
                                           fit$spec.var$exog.names)))/num.endog
      
      resid <- fit$resid[block.sample(nrow(fit$resid), 
                                      nrow(fit$resid) + num.lags), ]
      colnames(resid) <- colnames(fit$resid)
      
      Phi <- fit$Phi
      # Create exogenous part
      Z <- matrix(NA, nrow = nrow(resid), ncol = 0)
      
      if ('const' %in% fit$spec.var$det.part) { Z <- cbind(Z, rep(1, nrow(resid)))}
      if ('trend' %in% fit$spec.var$det.part) { Z <- cbind(Z, 1:nrow(resid))}
      if ('quad.trend' %in% fit$spec.var$det.part) { Z <- cbind(Z, 1:nrow(resid))}
      
      if (!is.null(fit$spec.var$exog.vars)) {
        Z <- cbind(Z, exog.vars)
      }
      
      sample <- matrix(NA, ncol = ncol(resid), nrow = nrow(resid))
      for (i in 1:nrow(resid)) {
        if (i == 1) {
          if (ncol(Z) == 0) {
            sample[1, ] <- resid[1, ]
            next
          } else {
            Z.temp <- cbind(Z[1, , drop = FALSE], 
                            matrix(0, ncol = num.lags*num.endog, nrow = 1 ))
          }
        } else if (i <= num.lags) {
          Z.temp <- matrix(0, nrow = 1, ncol = num.lags*num.endog)
          Z.temp[, 1:(num.endog*(i - 1))] <- 
            c(t(sample[(i - 1):1, ]))
          Z.temp <- cbind(Z[i, ], Z.temp)
        } else {
          Z.temp <- cbind(Z[i, ], matrix(c(t(sample[i - num.lags - 1  + 
                                                      num.lags:1, ])), nrow = 1))
        }
        sample[i, ] <- Z.temp %*% Phi + resid[i, ]
      }
      colnames(sample) <- colnames(fit$resid)
      return(sample)
    }
    
    # Function to create a pair (Phi, Sigma) from posterior with noniformative prior
    # at some step
    fit.step <- function(fit) {
      sample.endog <-  create.endog.sample(fit)
      spec.var <- specify.var(endog = sample.endog, exog = fit$spec.var$exog.vars, 
                              num.lags = fit$spec.var$num.lags, 
                              det.part = fit$spec.var$det.part) 
      new.fit <- try(fit.ols.var(spec.var = spec.var), silent = TRUE)
      if ("try-error" %in% class(new.fit)) { next}
      
      return(list(A = solve(chol(new.fit$Sigma)),
                  B = new.fit$Phi %*% solve(chol(new.fit$Sigma))))
    }
    
  } else if (inference.method == 'sims.zha') {
    fit.step <- function(fit, prior) {
      require(MASS)
      calc.posterior <- function(fit, prior) {
        V.B.posterior <-  solve(solve(prior$V.B.prior) +
                                  t(fit$spec.var$Z) %*% fit$spec.var$Z)
        B.posterior <- V.B.posterior %*% (
          solve(prior$V.B.prior) %*% prior$B.prior +
            t(fit$spec.var$Z) %*% fit$spec.var$Y)
        
        S.inv <- solve(t(fit$spec.var$Y) %*% fit$spec.var$Y +
                         t(prior$B.prior) %*% prior$V.B.prior %*% prior$B.prior -
                         t(B.posterior) %*% V.B.posterior %*% B.posterior )
        return(list(V.B.posterior =  V.B.posterior,
                    B.posterior =  B.posterior,
                    S.inv =  S.inv))
      }
      num.endog <- ncol(fit$spec.var$Y)
      
      # Generate A matrix
      posterior <-  calc.posterior(fit, prior)
      df.rwishart <- nrow(fit$spec.var$Y)  + num.endog + 1
      Bartlett.dec <- matrix(0, ncol = num.endog, nrow = num.endog)
      diag(Bartlett.dec) <- sqrt(sapply(1:ncol(Bartlett.dec), function(i) {
        rchisq(n = 1, df = df.rwishart - i + 1)
      }))
      Bartlett.dec[upper.tri(Bartlett.dec)] <- rnorm(num.endog*(num.endog - 1)/2)
      A <- chol(posterior$S.inv) %*% Bartlett.dec
      
      # Generate B matrix
      B <- sapply(1:num.endog, function(i) {
        mvrnorm(n = 1, mu = posterior$B.posterior %*% A[, i], 
                Sigma = posterior$V.B.posterior)
      })
      colnames(B) <- rownames(A)
      colnames(A) <- rownames(A)
      return(list(A = A, B = B))
    }
  } else if (inference.method == 'baumeister.hamilton'){
  } else if (inference.method == 'waggoner.zha'){
  }
  return(fit.step)
}


# ---- 3. Identification methods ----
# ==== 3.1 ARRW identification scheme ====
# #### 3.1.1 define preliminary functions ####
specify.functions <- function(inference.method = NULL, restrictions.type) {
  if (restrictions.type == 'only.sign') {
    main <- function(starting.mat, fit) {
      X <- matrix(rnorm(length(starting.mat)), 
                  ncol = ncol(starting.mat), nrow = nrow(starting.mat))
      new.draw <- starting.mat %*% qr.Q(qr(X))
      rownames(new.draw) <- rownames(starting.mat)
      colnames(new.draw) <- colnames(starting.mat)
      return(list(A = new.draw, 
                  B = fit$Phi %*% new.draw))
    }
    return(main)
  } else if (restrictions.type == 'zero.and.sign') {
    fit.step <- specify.method(inference.method)
    
    # Define main funciton to draw from a distribution
    main <- function(fit, 
                     zero.restrictions.on.irf = FALSE,  
                     zero.restrictions.on.A = FALSE,
                     zero.restrictions.on.B = FALSE, ...) {
      require(pracma)
      num.lags <- fit$spec.var$num.lags
      # NB! If restrictions are imposed on irf, irf should be estimated before the function call
      new.params <- fit.step(fit = fit, ...)
      
      if (zero.restrictions.on.irf) {
        # Construct starting values of IRF
        irf <- calc.irf(structural.params = new.params, 
                        exog.names = fit$spec.var$exog.names, 
                        det.part = fit$spec.var$det.part,
                        horizon = attr(zero.restrictions$irf, 'which.horizon'))
      }
      
      Q.jm1 <- matrix(0, nrow = 0, ncol = nrow(fit$Sigma))
      
      for (j in 1:nrow(fit$Sigma)) {
        
        restricted.function <- matrix(NA, nrow = 0, ncol = ncol(fit$Sigma))
        if (zero.restrictions.on.irf) {
          # restr.j  <- zero.restrictions$irf[, (j - 1)*ncol(A) + 1:ncol(A), drop = FALSE]
          if (!is.null(zero.restrictions$irf[[j]])) {
            restricted.function <- t(sapply(
              zero.restrictions$irf[[j]]$which.resp,
              function(x) {
                irf$irf[zero.restrictions$irf[[j]]$which.horizon, 
                        grepl(x, attr(zero.restrictions$irf, 'irf.colnames'))]
              }))
          }
          
          # LR restrictions
          cond.zero <- !is.na(zero.restrictions$lr.irf[j, ])
          if (any(cond.zero)) {
            restricted.function <- rbind(restricted.function,
                                         t(irf$lr.irf[, cond.zero, drop = FALSE]))
          }
        } 
        
        if (zero.restrictions.on.A) {
          # NB! j corresponds to equations in this case, due to choosing
          # not matrix Q*IRF, but A*t(Q) itself
          cond.zero <- !is.na(zero.restrictions$A[, j])
          if (any(cond.zero)) {
            restricted.function <- rbind(restricted.function,
                                         new.params$A[cond.zero, , drop = FALSE]) 
          }
        }
        
        if (zero.restrictions.on.B) {
          for (l in 1:fit$spec.var$num.lags) {
            # NB! j corresponds to equations in this case, due to choosing
            # not matrix Q*IRF, but B*t(Q) itself
            which.Bmat <- length(fit$spec.var$det.part) +
              (l - 1) * ncol(new.params$B) + 1:ncol(new.params$B)
            cond.zero <- !is.na(zero.restrictions$B[which.Bmat, j])
            if (any(cond.zero)) {
              restricted.function <- rbind(restricted.function,
                               new.params$B[which.Bmat, ][cond.zero, ]) 
            }
          }
        }
        
        R.j <- rbind(restricted.function, Q.jm1)
        if (nrow(R.j) == 0) {
          N <- diag(rep(1, ncol(fit$Sigma)))
        } else {
          N <- nullspace(R.j)
          if (is.null(N)) {
            stop('Too many restrictions are posed on the ', colnames(A)[j],
                 ' shock')
          }
        }
        
        # x.j <- matrix(rnorm(nrow(A)), ncol = 1)
        # q.j <- t(N) %*% x.j
        # q.j <- N %*% q.j/sqrt(sum(q.j^2))
        y.j <- matrix(rnorm(ncol(N)), ncol = 1)
        q.j <- N %*% y.j/sqrt(sum(y.j^2))
        
        Q.jm1 <- rbind(Q.jm1, c(q.j))
      }
      
      A <- new.params$A %*% t(Q.jm1); colnames(A) <- rownames(A)
      B <- new.params$B %*% t(Q.jm1); colnames(B) <- rownames(A)
      return(list(A = A, B= B))
    }
    
    return(main)
  }
}

create.irf.restrictions.template <- function (structural.params, horizon = 50) {
  irf <- matrix(NA, nrow = horizon, ncol = nrow(structural.params$A)^2)
  colnames(irf) <- c(sapply(rownames(structural.params$A), function(x) {
    paste0('shock.', x, '_resp.', colnames(structural.params$A))}))
  
  lr.irf <- matrix(NA, ncol = ncol(structural.params$A), 
                   nrow = nrow(structural.params$A))
  colnames(lr.irf) <- paste0('resp.', colnames(structural.params$A))
  rownames(lr.irf) <- paste0('shock.', rownames(structural.params$A))
  
  Amat <-  matrix(NA, ncol = ncol(structural.params$A), 
                  nrow = nrow(structural.params$A))
  
  colnames(Amat) <- colnames(structural.params$A)
  rownames(Amat) <- rownames(structural.params$A)
  
  Bmat <-  matrix(NA, ncol = ncol(structural.params$B), 
                  nrow = nrow(structural.params$B))
  colnames(Bmat) <- colnames(structural.params$B)
  rownames(Bmat) <- rownames(structural.params$B)
  
  return(list(irf = irf, lr.irf = lr.irf, Amat = Amat, Bmat = Bmat))
}

restrictions.irf.matrices.to.vectors <- function(restrictions, is.sign = TRUE){
  # Prepare functions
  find.shock <- function(restrictions) {
    which.shock <- colnames(restrictions$irf)[sapply(1:ncol(restrictions$irf), 
                                                     function(j) {
                                                       any(!is.na(restrictions$irf[, j]))
                                                     })]
    which.shock <- strsplit(gsub('', x = which.shock, pattern = 'shock.'), split = '_')
    return(unique(sapply(seq_along(which.shock), function(i) {
      which.shock[[i]][1]
    })))
  }
  
  find.resp <- function(restrictions) {
    which.resp <- colnames(restrictions)[sapply(1:ncol(restrictions), 
                                                function(j) {
                                                  any(!is.na(restrictions[, j]))
                                                })]
    which.resp <- strsplit(x = which.resp, split = '_')
    return(unique(sapply(seq_along(which.resp), function(i) {
      which.resp[[i]][2]
    })))
  }
  
  # 1. restructurize irf restrictions
  # 1.a Sign restrictions
  if (is.sign) {
    which.horizon <- which(sapply(1:nrow(restrictions$irf), function(i) {
      any(!is.na(restrictions$irf[i, ]))
    }))
    which.horizon <- which.horizon[length(which.horizon)]
    
    vec.restrictions <- c(restrictions$irf[1:which.horizon, ])
    attr(vec.restrictions, 'which.resp')  <- which(!is.na(vec.restrictions))
    attr(vec.restrictions, 'which.shock') <- find.shock(restrictions)
    attr(vec.restrictions, 'which.horizon') <- which.horizon
    
    restrictions$irf <- vec.restrictions
    return(restrictions)
  } else {
    n.shocks <- sqrt(ncol(restrictions$irf))
    if (ncol(restrictions$irf) %% n.shocks != 0) {
      stop('Incorrect IRF matrix specification')
    }
    shocks <- strsplit(x = colnames(restrictions$irf), split =  '_')
    shocks <- gsub(x = unique(sapply(seq_along(shocks), function(i) {
      shocks[[i]][[1]]
    })), pattern = 'shock.', replacement = '')
    
    res <- vector('list', n.shocks)
    names(res) <- shocks
    for (j in 1:n.shocks) {
      restr.j  <- restrictions$irf[, (j - 1)*n.shocks + 1:n.shocks, drop = FALSE]
      if (all(is.na(restr.j))) {
        next
      } else {
        which.horizon <- which(sapply(1:nrow(restr.j), function(i) {
          any(!is.na(restr.j[i, ]))
        }))
        which.resp <- find.resp(restr.j)
        res[[j]] <- list(which.horizon = which.horizon,
                         which.resp  = which.resp)
      }
    }
    attr(res, "which.horizon") <- max(sapply(seq_along(res), function(i) {
      if (is.null(res[[i]])) { 0 } else{
        res[[i]]$which.horizon[length(res[[i]]$which.horizon)]
      }
    }))
    
    attr(res, 'irf.colnames') <- colnames(restrictions$irf)
    restrictions$irf <- res
    return(restrictions)
  }
}

# #### 2.1.2 Main function ####
identify.arrw <- function(fit, sign.restrictions, starting.mat = NULL,
                          zero.restrictions = NULL,
                          inference.method = NULL,
                          stop.threshold = 5e8, n.draw = 3e2, 
                          irf.horizon = 50, is.continuous = FALSE,
                          spec.notification = NULL, construct.lr.irf = TRUE,
                          ...) {
  pracma::tic()
  
  check.sign <- function(x, s) {
    if (s == '+') {
      x >= 0
    } else if (s == '-') {
      x <= 0
    } else {
      stop('Sign restrictions are not correctly specified')
    }
  }
  
  if (is.null(starting.mat)) {starting.mat <- solve(chol(fit$Sigma))}
  if (is.null(zero.restrictions)) {
    restrictions.type <- 'only.sign'
  } else {
    restrictions.type <- 'zero.and.sign'
  }
  
  sign.restrictions.on.irf <- !all(is.na(sign.restrictions$irf))
  sign.restrictions.on.A <- !all(is.na(sign.restrictions$A))
  sign.restrictions.on.B <- !all(is.na(sign.restrictions$B))
  
  zero.restrictions.on.irf <- !all(is.na(zero.restrictions$irf))
  zero.restrictions.on.A <- !all(is.na(zero.restrictions$A))
  zero.restrictions.on.B <- !all(is.na(zero.restrictions$B))
  
  create.draw <- specify.functions(inference.method = inference.method,
                                   restrictions.type = restrictions.type)
  
  
  # Construct store objects
  # parameters
  A.draws <- matrix(NA, nrow = n.draw, ncol = length(starting.mat))
  colnames(A.draws) <- c(sapply(1:ncol(starting.mat), function(i) {
    paste0('eq.', colnames(starting.mat)[i], '_var.', 
           rownames(starting.mat))}))
  
  B.draws <- matrix(NA, nrow = n.draw, ncol = length(fit$Phi))
  colnames(B.draws) <- c(sapply(1:ncol(starting.mat), function(i) {
    paste0('eq.', colnames(fit$Phi)[i], '_var.', 
           rownames(fit$Phi))}))
  
  # IRFs
  irf.draws <- lapply(1:length(starting.mat), function(i) {
    return(matrix(NA, ncol = irf.horizon, nrow = n.draw))
  })
  names(irf.draws) <- c(sapply(1:ncol(starting.mat), function(i) {
    paste0('shock.', colnames(starting.mat)[i], '_resp.', 
           rownames(starting.mat))}))
  
  # LR IRFs
  if (construct.lr.irf ) {
    lr.irf.draws <- matrix(NA, nrow = n.draw, ncol = length(starting.mat))
    colnames(lr.irf.draws) <- c(sapply(1:ncol(starting.mat), function(i) {
      paste0('shock.', colnames(starting.mat)[i], '_resp.', 
             rownames(starting.mat))}))
  } else {
    lr.irf.draws < NULL
  }
  
  # Start MCMC iterations
  k <- 1; i <- 0
  while (k <= n.draw) {
    i <- i + 1 
    if (i - n.draw >= stop.threshold) { 
      stop('MCMC process terminated. Stop threshold is achieved')
    }
    if (i %% 1e3 == 0)  {
      cat(paste0(spec.notification, 'MCMC Iteration ', i, '; Stored ', k - 1, '/', n.draw, ' draws \n'))
    }
    
    # Construct a new draw
    if (restrictions.type == 'only.sign') {
      new.draw <- create.draw(starting.mat, fit)
    } else if (restrictions.type == 'zero.and.sign') {
      if (inference.method %in% c('bootstrap', 'sims.zha')) {
        new.draw <- create.draw(fit = fit,
                                zero.restrictions.on.irf =  zero.restrictions.on.irf,
                                zero.restrictions.on.A  =  zero.restrictions.on.A,
                                zero.restrictions.on.B  =  zero.restrictions.on.B, 
                                ...)
      } 
    }
    
    # Normalization
    # a. Nonnegative diagonal
    norm.cond <- diag(new.draw$A) < 0
    if (sum(norm.cond) != 0) {
      new.draw$A[, norm.cond] <- - new.draw$A[, norm.cond]
      new.draw$B[, norm.cond] <- - new.draw$B[, norm.cond]
    }
    # b. No zero elements at diagonal
    if (any(abs(diag(new.draw$A)) <= 1e-10)) { next }
    
    # Vectorize parameters draw
    
    vec.draw <- c(new.draw$A)
    
    # Prepare condition
    cond.hold <- TRUE
    if (sign.restrictions.on.irf) {
      previous.irf <- calc.irf(structural.params = new.draw, 
                                               which.shock = attr(sign.restrictions$irf, 'which.shock'),
                                               horizon = irf.horizon,
                                               calc.horizon = attr(sign.restrictions$irf, 'which.horizon'), 
                                               calc.long.run = construct.lr.irf,
                                               previous.irf = NULL,
                                               is.continuous = is.continuous)
      
      # Check sign restrictions on finite IRFs
      # It's an R trick: all(NULL) = TRUE, so there is no
      # need to specify condition on non-NA SR restrictions
      vec.irf <- previous.irf$irf[attr(sign.restrictions$irf, 'which.horizon'), ]
      cond.hold <- all(sapply(attr(sign.restrictions$irf, 'which.resp'), 
                              function(i) {
                                check.sign(x = vec.irf[i],
                                           s = sign.restrictions$irf[i])
                              }))
      
      # Check sign restrictions on LR impact
      # It's an R trick: all(NULSystem.out.println("Hello, World");L) = TRUE, so there is no
      # need to specify condition on non-NA LR restrictions
      cond.hold <- cond.hold & (all(sapply(which(!is.na(sign.restrictions$lr.irf)), 
                       function(i) {
                         check.sign(x = previous.irf$lr.irf[i],
                                    s = sign.restrictions$lr.irf[i])
                                    })))
    } 
    
    # Check restrictions on A
    if (sign.restrictions.on.A) {
      vec.restrictions <- c(sign.restrictions$Amat)
      cond.hold <- all(sapply(which(!is.na(vec.restrictions)), 
                              function(i) {
                                check.sign(x = vec.draw[i],
                                           s = vec.restrictions[i])
                              }))
      previous.irf <- NULL
    }
    
    # Check restrictions on B
    if (sign.restrictions.on.B) {
      vec.restrictions <- c(sign.restrictions$Bmat)
      vec.Bdraw <- c(new.draw$B)
      cond.hold <- all(sapply(which(!is.na(vec.restrictions)), 
                              function(i) {
                                check.sign(x = vec.Bdraw[i],
                                           s = vec.restrictions[i])
                              }))
      previous.irf <- NULL
    }
    
    if (cond.hold) {
      # Store parameters values
      A.draws[k, ] <- vec.draw
      B.draws[k, ] <- c(new.draw$B)
      
      # Calculate IRFs
      irf.new.draw <- calc.irf(structural.params = new.draw, 
                               horizon = irf.horizon,
                               calc.long.run = construct.lr.irf,
                               previous.irf = previous.irf,
                               is.continuous = is.continuous)
      
      # Store LR irf values if neccessary
      if (construct.lr.irf ) {
        lr.irf.draws[k, ] <- c(t(irf.new.draw$lr.irf))
        irf.new.draw <- irf.new.draw$irf # this is made to achieve object class consistency
      }
      
      # Store irf values
      for (x in colnames(irf.new.draw)) {
        irf.draws[[x]][k, ] <- irf.new.draw[, x]
      }
      
      k <- k + 1
    }
  }
  
  attr(irf.draws, 'n.draw')  <- n.draw
  attr(irf.draws, 'horizon') <- irf.horizon
  
  if (i == k) {
    warning('There were no draws, that did not satisfy identification conditions. Restrictions might gain no sufficient information.')
  }
  
  time.elapsed <- pracma::toc(echo = FALSE)
  cat('MCMC Process finished. Elapsed time is', floor(time.elapsed/60), 'min', time.elapsed - 60 * floor(time.elapsed/60), 'sec. \n')
  
  return(list(par.draws = list(A = A.draws, B = B.draws), 
              irf.draws = irf.draws, lr.irf.draws = lr.irf.draws,
              time.elapsed = time.elapsed))
}

# ==== 3.2 PFA method ====
# 
# # ==== 3.3 Baumeister Hamilton identification technique =====
# define.prior.baumeister.hamilton(num.equations, num.lags = 1,
#                                  overall.tightness = 0.5,
#                                  relative.tightness = 0.2,
#                                  deterministic.tightness = 1e3,
#                                  lag.decay.rate = 2, 
#                                  assume.random.walk = TRUE,
#                                  s = NULL, det.part = 'const',
#                                  exog.names = NULL,
#                                  endog.names = NULL,
#                                  spec.var = NULL, Vi = 0.1,
#                                  kappa = rep(2, num.equations),
#                                  mean.A, Sigma.A, nu.A = rep(3, num.equations)) {
#   prior.params <- define.prior.classical.sims.zha(num.equations, num.lags,
#                                                  overall.tightness,
#                                                  relative.tightness,
#                                                  deterministic.tightness,
#                                                  lag.decay.rate, 
#                                                  assume.random.walk,
#                                                  s, det.part, exog.names,
#                                                  endog.names, spec.var)
#   prior.params$A <- mean.A
#   prior.params$Sigma.A <- Sigma.A
#   
# }
identify.baumeister.hamilton <- function() {
  
}

# ==== 3.4 Geweke with Waggoner-Zha estimation ====

# ---- 4. Post-estimation functions ----
# ==== 4.1 IRF plots ====
# #### 4.1.1 Standard plot ####
plot.irf <- function(irf.draws, num.shocks = NULL, which.shock = NULL,
                     which.responce = NULL, mar=c(1,2,2.5,0),
                     draw.each.model = TRUE, horizon = NULL,
                     ...) {
  if (is.null(num.shocks)) {num.shocks <- sqrt(length(irf.draws))}
  num.resp <- length(irf.draws)/num.shocks
  old.par <- par()
  par(mfrow = c(num.shocks, num.resp),
      mar=mar, ...)
  
  for (s in seq_along(irf.draws)) {
    if (is.null(horizon)) {horizon <- ncol(irf.draws[[s]])}
    IR <- irf.draws[[s]][, 1:horizon, drop = FALSE]
    IR.med <- sapply(1:ncol(IR), function(i) { median(IR[, i])})
    IR.mode <- sapply(1:ncol(IR), function(i) { modeest::mlv(IR[, i], method = "mfv")$M})
    IR.975 <- sapply(1:ncol(IR), function(i) { quantile(IR[, i], probs = 0.975)})
    IR.025 <- sapply(1:ncol(IR), function(i) { quantile(IR[, i], probs = 0.025)})
    
    
    plot(1, type="n", xlab="", ylab="", xlim=c(0, ncol(IR) + 1), 
         ylim=c(min(IR.975, IR.025), max(IR.975, IR.025)), 
         main = names(irf.draws)[s])
    grid()
    polygon(c(1:ncol(IR), rev(1:ncol(IR))),c(IR.025 ,rev(IR.975)), 
            col=rgb(0, 0, 1, 0.2), border = NA)
    if (draw.each.model) {
      for (i in 1:nrow(IR)) {
        lines(IR[i, ], col = 'gray', lty=3)
      }
    }
    lines(IR.975, col = 'blue', lty=2)
    lines(IR.025, col = 'blue', lty=2)
    lines(IR.med, col = 'blue')
    lines(IR.mode, col = 'blue', lty = 3)
    abline(a = 0, b = 0, col = 'red')
  }
  par <- old.par
}

# #### 4.1.2 Multimodal confidence interval estimation ####
calc.multimodal.intervals <- function(x, n = 500, p.max = 0.95,
                                      step = 1e-3,  ...) {
  is.unique.int<- function(cond) {
    which.cond <- which(cond)
    res.cond <- sum(cond)/
      (which.cond[length(which.cond)] - 
         which.cond[1] + 1) == 1
    return(res.cond)
  }
  
  find.intervals <- function (cond) {
    res <- matrix(NA, 0, 2)
    for (i in which(cond)) {
      if (i == 1) {
        k <- 1
      } else if (!cond[i - 1]){
        k <- i
      }
      
      if (is.na(cond[i + 1])) {
        res <- rbind(res, c(k, i))
      } else if (!cond[i + 1]) {
        res <- rbind(res, c(k, i))
      }
    }
    return(res)
  }
  
  dens <- density(x = x, n = n, ...)
  max.freq <- max(dens$y)
  freq.loop <- seq(from = 0, to = max.freq, by = step*max.freq)
  
  p.new <- 1; cond.new <- NULL; i <- 0
  while (p.new >= p.max) {
    i <- i + 1
    p <- p.new; cond <- cond.new 
    
    cond.new <- dens$y >= freq.loop[i]
    p.new <- sum(dens$y[cond.new])/sum(dens$y)
  }
  
  check <- abs(c(p, p.new) - p.max)
  if ((check != min(check))[1]) { 
    cond <- cond.new; p <- p.new; freq  <- freq.loop[i - 1]
  } else {freq <- freq.loop[i]}
  
  x.dropped <- dens$x[cond]
  y.dropped <- dens$y[cond]
  
  if ( is.unique.int(cond)) {
    intervals <-c(which(cond)[1], which(cond)[length(which(cond))])
    x.intervals <- matrix(dens$x[intervals], ncol = 2)
  } else {
    intervals <- find.intervals(cond)
    x.intervals <- cbind(dens$x[intervals[, 1]], 
                         dens$x[intervals[, 2]])
  }
  return(list(x.dropped = x.dropped, y.dropped = y.dropped, p = p, freq = freq,
              cond = cond, 
              intervals = intervals, x.intervals = x.intervals))
}

calc.irf.intervals <- function(irf, n = 1e3, p.max = 0.95,
                               step = 1e-3, error = 1e-4, 
                               horizon = 20, ...) {
  require(LaplacesDemon)
  res <- lapply(1:length(irf), function(i) { vector('list', length = horizon)})
  names(res) <- names(irf)
  intervals <- res;
  max.n.intervals <- 1
  
  for(x in names(irf)) {
    for (i in 1:horizon) {
      if (is.multimodal(irf[[x]][, i])) {
        res[[x]][[i]] <- calc.multimodal.intervals(x = irf[[x]][, i],
                                                   n, p.max, step, ...)
        intervals[[x]][[i]] <- res[[x]][[i]]$x.intervals
        max.n.intervals <- max(max.n.intervals, nrow(intervals[[x]][[i]]))
      } else {
        res[[x]][[i]] <- NULL
        intervals[[x]][[i]] <- p.interval(obj = irf[[x]][, i], prob = p.max, ...)
      }
    }
  }
  
  for(x in names(irf)) {
    for (i in 1:horizon) {
      intervals[[x]][[i]] <- rbind(intervals[[x]][[i]], matrix(0, ncol = 2,
                                                               nrow = max.n.intervals - nrow(intervals[[x]][[i]])))
    }
  }
  return(list(res = res, intervals = intervals, max.n.intervals = max.n.intervals))
}

plot.irf.hpd <- function(irf.draws, num.shocks = NULL, which.shock = NULL,
                         which.responce = NULL, mar=c(1, 2, 2.5, 0),
                         draw.each.model = TRUE, horizon = NULL,
                         p.seq = c(0.99, 0.95, 0.9, 0.75, 0.5),
                         ...) {
  if (is.null(num.shocks)) {num.shocks <- sqrt(length(irf.draws))}
  num.resp <- length(irf.draws)/num.shocks
  old.par <- par()
  par(mfrow = c(num.shocks, num.resp),
      mar=mar, ...)
  
  irf.hpd <- lapply(p.seq, function(p) {calc.irf.intervals(irf = irf.draws, 
                                                           n = 1e3, p.max = p,
                                                           step = 1e-3, error = 1e-4, 
                                                           horizon = horizon) })
  
  for (s in seq_along(irf.draws)) {
    if (is.null(horizon)) {horizon <- ncol(irf.draws[[s]])}
    IR <- irf.draws[[s]][, 1:horizon, drop = FALSE]
    
    plot(1, type="n", xlab="", ylab="", xlim=c(0, ncol(IR) + 1), 
         ylim=c(min(unlist(irf.hpd[[which(p.seq == max(p.seq))]]$intervals[[s]])), 
                max(unlist(irf.hpd[[which(p.seq == max(p.seq))]]$intervals[[s]]))), 
         main = names(irf.draws)[s]); grid()
    
    for (i.p in seq_along(p.seq)) {
      for ( i in 1:irf.hpd[[i.p]]$max.n.intervals) {
        IR.low  <- sapply(1:horizon, function(h) {irf.hpd[[i.p]]$intervals[[s]][[h]][i, 1]})
        IR.high <- sapply(1:horizon, function(h) {irf.hpd[[i.p]]$intervals[[s]][[h]][i, 2]})
        
        polygon(c(1:ncol(IR), rev(1:ncol(IR))),c(IR.low ,rev(IR.high)), 
                col=rgb(0, 0, 1, 0.15), border = NA)
      }
    }
    abline(a = 0, b = 0, col = 'red')
  }
  par <- old.par
}

# ==== 4.2. FEVD plots ====
# ==== 4.3. Historical Decompositions plots ====

# ----- 5. Convergence Diagnostics ----