# ---- 1.  Preliminary setup and class definition ----
setOldClass("proc_time")
setClassUnion("list.or.null", c("list", "NULL", "data.frame"))
setClass('BSVAR.out', representation(draws = "list",
                                     marginal.likelihood = "list.or.null",
                                     irf   = "list.or.null",
                                     timer  = "proc_time"))

setMethod("summary", "BSVAR.out", function(object, y = NULL) {
  # Generic function for convergence check and summary statistics
  require(MCMCpack)
  # Convergence tests
  draws <- mcmc(cbind(data.frame(object@draws$b),
                 data.frame(object@draws$g)))
  heidel.welch  <- heidel.diag(draws)
  raftery.lewis <- raftery.diag(draws)$resmat
  colnames(raftery.lewis) <- paste0("raftery_", colnames(raftery.lewis))
  geweke        <- geweke.diag(draws)
  if (!is.null(y)) {
    gelman.rubin <- gelman.plot(mcmc.list(object@draws, y@draws))
  } else {
    sample1 <- 1:floor(nrow(draws)/2)
    sample2 <- (floor(nrow(draws)/2) + 1):nrow(draws)
    gelman.rubin <- gelman.diag(mcmc.list(
      mcmc(draws[sample1, ]), 
      mcmc(draws[sample2, ])))
  }
  
  convergence <- list(heidel.welch  = heidel.welch,
              raftery.lewis = raftery.lewis,
              geweke        = geweke,
              gelman.rubin  = gelman.rubin)
  
  sum.table <- cbind(summary(draws)$statistics, 
                  summary(draws)$quantiles[ ,c(1, 5)])
  colnames(sum.table)[c(5, 6)] <- paste0("q", colnames(sum.table)[c(5, 6)])
  return(list(summary = sum.table,
              convergence = convergence))
  })

# ---- 2. General functions for all SVAR samplers ----
# ==== 2.2 Restriction matrix definition ====
# 
# In notations of K. Juselius and S. Johansen R-type restriction matrix is a
# matrix that establishes a system of equations to hold: R*params = 0;
# H-type matrix is a matrix that connects vector of all parameters with the
# vector of only unrestricted parameters: all.params = H*unrestr.params.
# If there are only exclusion constraints, it is convenient to use 
# exact restriction matrix 

# #### 2.2.1 From exact restrictions to H-type #####
exact.to.H <- function(x) {
  # x is assumed to be logical
  res <- lapply(1:ncol(x), function(i) {
    cond <- x[, i]
    if (sum(cond) == 0){
      H <- diag(rep(1, nrow(x)))
      rownames(H) <- rownames(x)
      colnames(H) <- rownames(x)
    } else {
      H <- matrix(0, nrow = nrow(x), ncol = sum(!cond))
      H[!cond, ] <- diag(rep(1, sum(!cond)))
      if (!is.null(colnames(x))) {
        rownames(H) <- rownames(x)
        colnames(H) <- rownames(x)[!cond]
      }
    }
    return(H)
  })
  names(res) <- colnames(x)
  return(res)
}

# #### 2.2.3 Function to generate artificial models ####
decode.restrictions <- function(integer, n, m) {
  constraints <- matrix(1, ncol = n, nrow = m)
  for (i in 1:m) {
    for (j in 1:n) {
      if (i == j) next
      if (integer %% 2) {
        constraints[i, j] <- 0
      }
      integer <- integer %/% 2
    }
  }
  constraints
}

# ==== 2.3 Match priors on reduced-form parameters with priors on structural ones =====
structural.priors.on.unrestricted.params <-function(prior, restrictions) {
  # It is assumed that priors between equations are independent.
  # This version of function does not support non-zero mean priors for the reduced-form coefficients.
  temp <- lapply(seq_along(restrictions$B), function(i) {
    which.B.i <- (i - 1)*nrow(restrictions$B[[i]]) + 1:nrow(restrictions$B[[i]])
    which.A.i <- (i - 1)*nrow(restrictions$A[[i]]) + 1:nrow(restrictions$A[[i]])
    
    S.i <- prior$vec.A$S[which.A.i, which.A.i, drop = FALSE]
    H.i <- prior$vec.B$H[which.B.i, which.B.i, drop = FALSE]
    P.i <- prior$vec.B$P[which.B.i, which.A.i, drop = FALSE]
    
    if (ncol(restrictions$B[[i]]) != 0) { # if there is no situation when the ith column of B matrix is restricted.
      H   <- solve(t(restrictions$B[[i]]) %*% solve(H.i) %*% restrictions$B[[i]])
      P   <- H %*% t(restrictions$B[[i]]) %*% solve(H.i) %*%
        P.i %*% restrictions$A[[i]]
      S   <- solve(t(restrictions$A[[i]]) %*% solve(S.i) %*% restrictions$A[[i]] + 
                     t(restrictions$A[[i]]) %*% t(P.i) %*% solve(H.i) %*%  
                     P.i%*%restrictions$A[[i]] -
                     t(P) %*% solve(H) %*% P)
    } else { # situation when the whole ith column if B matrix is restricted
      S  <- solve(t(restrictions$A[[i]]) %*% solve(S.i) %*% restrictions$A[[i]])
      H  <- NULL
      P  <- NULL
    }
    return(list(H = H, P = P, S = S))
  })
  res <- lapply(1:3, function(i) {
    lapply(seq_along(temp), function(j) {
      temp[[j]][[i]]
    })
  })
  names(res) <- c('H', 'P', 'S')
  return(res)
}

# ==== 2.4 Define starting values ====
# Two types of starting values are considered: 
# 1. direct mle estimates,
# 2. ols estimates with random rotation matrix
define.starting.values <- function(Y, Z, restrictions, type = "both") {
  start.values.mle <- NULL; mle.exact <- NULL
  start.values.ols <- NULL; ols       <- NULL 
  if (type %in% c("both", "mle")) {
    # #### 4.3.1 Starting values using MLE estimates ####
    mle.exact <- mle(Y, Z, constraints.on.coefficients = restrictions$exact,
                     start.type = "2SLS")
    
    A.start <- lapply(1:ncol(Y), function(i) {
      cond <- !restrictions$exact$A[, i]
      mle.exact$A[1:ncol(Y), ][cond, i, drop = FALSE]
    })
    
    B.start <- lapply(1:ncol(Y), function(i) {
      cond <- !restrictions$exact$B[, i]
      - mle.exact$A[(ncol(Y) + 1):nrow(mle.exact$A), , drop = FALSE][cond, i, drop = FALSE]
    })
    
    start.values.mle <- list(b = A.start, g = B.start)
  }
  if (type %in% c("both", "ols")) {
    ols <- ols.VAR(Y, Z)
    chol.Omega <- t(chol(ols$covmat))
    
    # Define random rotation matrix
    rand.matrix <- matrix(rnorm(ncol(chol.Omega)*nrow(chol.Omega)), 
                          ncol = ncol(chol.Omega), nrow = nrow(chol.Omega))
    rotation.mat <- qr.Q(qr(t(rand.matrix*(!restrictions$exact$A))))
    
    # Define parameters
    A <- (chol.Omega %*% t(rotation.mat))
    B <- ols$A %*% rotation.mat
    
    # Define starting parameters 
    A.start <- lapply(1:ncol(Y), function(i) {
      cond <- !restrictions$exact$A[, i]
      matrix(A[cond, i], ncol = 1)
    })
    
    B.start <- lapply(1:ncol(Y), function(i) {
      cond <- !restrictions$exact$B[, i]
      matrix(B[cond, i], ncol = 1)
    })
    
    start.values.ols <- list(b = A.start, g = B.start)
  }
  
  P <- 0*rbind(restrictions$exact$A, restrictions$exact$B)
  diag(P) <- 1
  A.prior <- lapply(1:ncol(Y), function(i) {
    cond <- !restrictions$exact$A[, i]
    P[1:ncol(Y), ][cond, i, drop = FALSE]
  })
  
  B.prior <- lapply(1:ncol(Y), function(i) {
    cond <- !restrictions$exact$B[, i]
    - P[(ncol(Y) + 1):nrow(mle.exact$A), , drop = FALSE][cond, i, drop = FALSE]
  })
  
  return(list(mle = start.values.mle,
              ols = start.values.ols,
              prior = list(b = A.prior, g = B.prior),
              estim = list(mle = mle.exact,
                           ols = ols)))
}

# ==== 2.5 Construct a matrix with b draws and H-type restrictions =====
free.params.to.matrix <- function(params, restrictions) {
  res <- matrix(NA, nrow = nrow(restrictions[[1]]), ncol = length(restrictions))
  for (i in seq_along(restrictions)) {
    if (!is.null(params[[i]])) {
      res[, i] <- restrictions[[i]] %*% matrix(params[[i]], ncol = 1)
    } else {
      res[, i] <- 0
    }
  }
  return(res)
}


# ==== 2.7 Direct marginal density estimation ====
# Based on Ding, Karlsson (2013), Eklund, Karlsson (2007) and Koop et.al (2008)
# Integral of type p(y|M) = \int p(y| \theta, M)p(\theta| M) dtheta could be approximated by
# the mean \hat{p}(y| \theta, M) = R^{-1} \sum_{r = 1}{R} p(y| \theta^{(r)}, M),
# where (r) is simulation draws from p(\theta)

# #### 2.7.1 Define likelihood function ####
construct.log.lik.fun <- function(Y, Z, restrictions, type = "Waggoner-Zha") {
  n <- ncol(Y)
  nT <- nrow(Y)
  X <- as.matrix(scale(cbind(Y, Z), scale = FALSE))
  if (type == "Waggoner-Zha") {
    function(b, g) {
      A.restr <- free.params.to.matrix(b, restrictions$A)
      B.restr <- free.params.to.matrix(g, restrictions$B)
      
      er.term <- X %*% rbind(A.restr, - B.restr)
      res <- - 0.5 * nT * n * log(2 * pi) + nT * log (abs(det(A.restr))) -
        0.5 * sum(er.term ^ 2)
      return(res)
    }
  }
}

# #### 2.7.1 Construct prior density function ####
construct.prior.density <- function(prior, type = "Waggoner-Zha", nT) {
  if (type == "Waggoner-Zha") {
    function(b, g) {
      # B.restr <- free.params.to.matrix(g, restrictions$B) 
      sum(sapply(seq_along(prior$S), function(i) {
        log(dmvnorm(x = c(b[[i]]), mean = rep(0, length(b[[i]])), 
                    sigma = prior$S[[i]])) +
        ifelse(is.null(prior$P[[i]]), 0, 
               log(dmvnorm(x = c(g[[i]]), mean = c(prior$P[[i]]%*%b[[i]]), 
                           sigma = prior$H[[i]])))
        }))
    }
  } else if (type == "diffuse") {
    function(b, g) {return(1)}
  }
}

# #### 2.7.2 Construct joint density function ####
construct.joint.posterior <- function(post.params, restrictions,
                                      n, nT,
                                      type = "Waggoner-Zha") {
  if (type == "Waggoner-Zha") {
    function(b, g) {
      A.restr <- free.params.to.matrix(b, restrictions$A) 
      return(sum(sapply(1:n, function(i) {
          log(dmvnorm(x = c(b[[i]]), mean = rep(0, length(b[[i]])),
                      sigma = nT*post.params[[i]]$S)) +
            ifelse(is.null(post.params[[i]]$P), 0, 
                   log(dmvnorm(x = c(g[[i]]), mean = c(post.params[[i]]$P%*%b[[i]]),
                               sigma = post.params[[i]]$H)))
        })) + nT*log(abs(det(A.restr))) - 0.5*nT*n*log(2*pi))
    }
  }
}

construct.joint.posterior.new <- function(post.params, restrictions,
                                      n, nT,
                                      type = "Waggoner-Zha") {
  if (type == "Waggoner-Zha") {
    function(b, g) {
      A.restr <- free.params.to.matrix(b, restrictions$A)
      if (length(b) == 2) {
        return(list(b = lapply(1:2, function(i){
          nT * log(abs(det(A.restr))) +
            log(dmvnorm(x = c(b[[i]]),  mean = rep(0, length(b[[i]])),
                        sigma = sqrt(nT)*post.params[[i]]$S))
        }), 
        g = lapply(1:2, function(i) {
          ifelse(is.null(post.params[[i]]$P), 0, 
                 log(dmvnorm(x = c(g[[i]]), mean = c(post.params[[i]]$P%*%b[[i]]),
                             sigma = post.params[[i]]$H)))
        })))
      }
    }
  }
}

# #### 2.7.2 Calculate an estimation of marginal density ####
direct.marginal.density <- function(lik.fun, prior, 
                                    sample.save = 1e4, burn.in = 5e2,
                                    spec.notification = NULL) {
  # Sample from prior distribution
  num.b <- sum(sapply(seq_along(prior$S), function(i) ncol(prior$S[[i]])))
  num.g <- sum(sapply(seq_along(prior$H), function(i) {
    if (is.null(prior$H[[i]])) {
      return(0)
    } else {
      return(ncol(prior$H[[i]]))
    }
    }))
  b.draw  <- matrix(NA, ncol = num.b,  nrow = sample.save)
  g.draw  <- matrix(NA, ncol = num.g, nrow = sample.save)
  mlf.draw <- matrix(NA, ncol = 1,  nrow = sample.save)
  for (i.rep in 1:sample.save + burn.in) {
    if (i.rep %% 1e3 == 0)  {
      cat(spec.notification, 'Marg Lik Iteration', i.rep, '/', 
          sample.save + burn.in, '\n')
    }
    
    b <- lapply(seq_along(prior$S), function(i) {
      matrix(runif(n = ncol(prior$S[[i]]), min = -1, max = 1), ncol = 1)
    })

    g <- lapply(seq_along(prior$P), function(i) {
      if (is.null(prior$P[[i]])) {
        return(NULL)
      } else {
        return(matrix(runif(n = ncol(prior$H[[i]]), min = -1, max = 1), ncol = 1))
      }
    })
    
    mlf.draw[i.rep - burn.in, ] <- lik.fun(b = b, g = g)
    b.draw[i.rep - burn.in, ] <- unlist(b)
    g.draw[i.rep - burn.in, ] <- unlist(g)
  }
  
  mlf.draw <- mcmc(mlf.draw)
  colnames(mlf.draw) <- "log marg lik"
  
  return(list(mlf = mlf.draw, b = b.draw, g = g.draw))
}
# ==== 2.8 IRF function ====
irf <- function(object) {
  
}

# ---- 3. SVAR Gibbs Samplers ----
# ==== 3.1 Waggoner-Zha algorithm ====
Waggoner.Zha.gibbs.sampler <- function(Y, Z, restrictions,
                               prior, irf.settings = NULL,
                               Chib.values = NULL,
                               mle.params = NULL,
                               start.values = NULL, # set MLE estimate by default
                               burn.in = 5e2, sample.save = 1e4,
                               likelihood.normalization = TRUE,
                               spec.notification = NULL) {
  start.timer <- proc.time()
  require(MCMCpack)
  require(mvtnorm)
  require(MASS)
  # 0. Set neccessary constants
  sample.size <- sample.save + burn.in
  nT <- nrow(Y)
  n  <- ncol(Y)
  k  <- ncol(Z)
  
  # 1. Prepare neccessary functions
  # 1.1 Generate weights for free structural parameters of the contemporaneous matrix.
  generate.weights <- function(A.restr, S, chol.S, restr.mat) {
    # This function calculates the w-vectors for decomposition of structural parameters
    # Auxilliary function to generate 0 in vector if condition is not hold
    select.row <- function(x, cond) { 
      x[!cond, ] <- 0
      return(x)
    }
    
    # Generate the 1st vector of the R^(q_i) space, which produces null-space(A.restr)
    w  <- Null(A.restr)
    w1 <- t(restr.mat %*% chol.S) %*% w / sqrt(sum((t(w) %*% restr.mat %*% chol.S)^2))
    
    # Calculate other vectors separately
    w.all <- lapply(1:nrow(w1), function(i) {
      if (i == 1) {
        return(c(w1))
      } else {
        w.i   <- c(select.row(x = w1, cond = 1:nrow(w1) <= i)) # 0-elements if j >  i
        w.im1 <- c(select.row(x = w1, cond = 1:nrow(w1) < i))  # 0-elements if j >= i
        res   <- w.im1*w1[i, ]/sqrt(sum(w.i^2)*sum(w.im1^2))
        res[i] <- -sum(w.im1^2)
        return(res)
      }
    })
    return(w.all)
  }
  
  # 1.2 Function to generate bi sample of free parameters contained in A0 matrix.
  sample.bi <- function(A.restr, S, chol.S, restr.mat, nT) {
    w.all <- generate.weights(A.restr, S, chol.S, restr.mat)
    
    beta <- rep(NA, ncol(S)) # Since number of free parameters = ncol(inv.S)
    # Generate the first beta from gamma distribution
    beta[1] <- ifelse(runif(n = 1) <= 0.5, -1, 1)*
      sqrt(rgamma(n = 1, shape = 0.5*(nT + 1), scale = 2/nT))
    # Generate other betas conditional on the first beta from Normal distribution
    if (length(beta) > 1) {
      beta[2:length(beta)] <- rnorm(n = length(beta) - 1, mean = 0, sd = 1/sqrt(nT))
    }
    
    sum.w <- matrix(unlist(w.all), ncol = ncol(S), byrow = FALSE)%*%
      matrix(beta, ncol = 1)
    bi <- chol.S%*%sum.w
    return(bi)
  }
  
  # 1.3 One-step Gibbs sampling for bi
  draw.b.posterior <- function(b.prev, mle.A, restrictions, post.params, nT, n) {
    b.draws <- vector('list', n)
    
    A.restr <- free.params.to.matrix(b.prev, restrictions$H.type$A) # set starting A.restr
    for (i in 1:n) {
      b.draws[[i]] <- sample.bi(A.restr[, -i], post.params[[i]]$S, 
                                post.params[[i]]$chol.S, 
                                restrictions$H.type$A[[i]], nT)
      b.prev[[i]]  <- b.draws[[i]]
      A.restr <- free.params.to.matrix(b.prev, restrictions$H.type$A) # update A.restr subject to a new draw of bi
    }
    return(b.draws)
  }
  
  # 1.4 Construct conditional function for bi parameter
  dconditional.bi <- function(b.i, A.restr, post.params, nT,
                              restr.mat, b.other.part = NULL) {
    # Initialization
    S <- post.params$S
    chol.S <- post.params$chol.S
    
    w.all <- generate.weights(A.restr, S, chol.S, 
                              restr.mat)
    beta <- c(solve(chol.S %*% matrix(unlist(w.all), 
                    ncol = ncol(S), byrow = FALSE)) %*% b.i)
    
    return(log(abs(2*beta[1])) + log( dgamma(x = beta[1]^2, shape = 0.5*(nT + 1), scale = 2/nT) ) + 
             ifelse(length(beta) == 1, 0, sum(log(sapply(2:length(beta),
             function(i) dnorm(beta[i], mean = 0, sd = 1/sqrt(nT)))))) +
             ifelse(is.null(b.other.part), 0, b.other.part))
  }
  
  # 1.4 Calculate parameters of posterior sample 
  calc.posterior.params <- function(Y, Z, prior, restrictions) {
    nT <- nrow(Y)
    # 1.5.1 Preliminary calculations to enhance speed
    ZZ <- t(Z)%*%Z
    ZY <- t(Z)%*%Y
    YY <- t(Y)%*%Y
    
    # 1.4.2 Compute parameters for posterior functions
    post.params <- lapply(1:ncol(Y), function(i) {
      U <- restrictions$H.type$A[[i]]
      V <- restrictions$H.type$B[[i]]
      
      if (ncol(V) != 0) { # the situation when ith column of B is not fully restricted
        H <- solve(t(V) %*% ZZ %*% V + solve(prior$H[[i]]))
        P <- H %*% (t(V) %*% ZY %*% U + solve(prior$H[[i]]) %*% prior$P[[i]])
        S <- solve((t(U) %*% YY %*% U + solve(prior$S[[i]]) + 
                      t(prior$P[[i]]) %*% solve(prior$H[[i]]) %*% prior$P[[i]] -
                      t(P) %*% solve(H) %*% P)/nT)
        chol.S <- t(chol(S)) # a matrix L of the inv.S = L'L, i.e. in W-Zh notation this is T = L'
      } else { # the situation when the whole ith column of B is restricted to 0
        S <- solve((t(U) %*% YY %*% U + solve(prior$S[[i]]))/nT)
        chol.S <- t(chol(S))
        H <- NULL
        P <- NULL
      }
      return(list(H = H, P = P, S = S, chol.S = chol.S))
    })
    return(post.params)
  }
  
  # 2. Gibbs Sampling itself
  post.params <- calc.posterior.params(Y, Z, prior, restrictions)
  
  # Data.frame for b and g draws
  b.draws   <- matrix(NA, nrow = sample.save, ncol= sum(!restrictions$exact$A))
  g.draws   <- matrix(NA, nrow = sample.save, ncol = sum(!restrictions$exact$B))
  
  # Names for draws
  colnames(b.draws) <- gibbs.draw.namer(mat.restrictions = restrictions$exact$A, 
                                        name = "b")
  colnames(g.draws) <- gibbs.draw.namer(mat.restrictions = restrictions$exact$B, 
                                        name = "g")
  
  
  # Data.frame for likelihoods
  lik.fun <- construct.log.lik.fun(Y = Y, Z = Z, restrictions = restrictions$H.type, 
                                   type = "Waggoner-Zha")
  joint.posterior <-  construct.joint.posterior(post.params, restrictions$H.type, 
                                                n = n, nT = nT,
                                                type = "Waggoner-Zha")
  dprior <- construct.prior.density(prior = prior, nT = nT)

  prior.density <- matrix(1, nrow = sample.save, ncol = 1)
  lf <- matrix(1, nrow = sample.save, ncol = 1)
  dconditional.b1.sample <- matrix(1, nrow = sample.save, ncol = 1)
  
  b.prev <- start.values$b
  g.prev <- start.values$g
  
  # Extract MLE estimates. It is assumed that if there is no mle.params object defined in
  # function inputs, start values should be used for LP normalization instead.
  if (is.null(mle.params)) {
    mle.A <- free.params.to.matrix(b.prev, restrictions$H.type$A) 
  } else {
    mle.A <- free.params.to.matrix(mle.params$b, restrictions$H.type$A)
  }
  
  for (i.rep in 1:sample.size) {
    if (i.rep %% 1e3 == 0)  {
      cat(spec.notification, 'Gibbs Sampler Iteration', i.rep, '/', sample.size, '\n')
    }
    
    b.prev <- draw.b.posterior(b.prev, mle.A, restrictions, post.params, nT, n)
    if (i.rep > burn.in) {
      g.prev <- lapply(seq_along(restrictions$H.type$A), function(j) {
        if (is.null(prior$P[[j]])) {
          return(NULL)
        } else {
          return(matrix(rmvnorm(n = 1, mean = post.params[[j]]$P %*% matrix(b.prev[[j]], 
                 ncol = 1), sigma = post.params[[j]]$H), ncol = 1))
        }
      })
      
      # Likelihood preserving normalization
      A.restr <- free.params.to.matrix(b.prev, restrictions$H.type$A) 
      if (likelihood.normalization) {
        for (i in 1:n) {
          normalize <- ifelse(sum(solve(A.restr)[i, ]*mle.A[, i]) < 0, -1, 1)
          b.prev[[i]] <- normalize*b.prev[[i]]
          g.prev[[i]] <- normalize*g.prev[[i]]
        }
      }
      
      # 3. Joint posterior density function
      # 3.1 Chib method   
      # 3.2 Direct method
      lf[i.rep - burn.in, ] <- lik.fun(b = b.prev, g = g.prev)
      
      # 3.3 Prior density 
      prior.density[i.rep - burn.in, ] <- dprior(b = b.prev, g = g.prev)
      
      # Conditional b1 density
      b.dcon <- b.prev; b.dcon[[1]] <- Chib.values$b[[1]]
      A.restr.temp <- free.params.to.matrix(b.dcon, restrictions$H.type$A)

      b.other.part <- NULL
      dconditional.b1.sample[i.rep - burn.in, ] <- dconditional.bi(b.i = b.dcon[[1]], 
                              A.restr = A.restr.temp[, - 1, drop = FALSE],
                              post.params = post.params[[1]], 
                              restr.mat = restrictions$H.type$A[[1]],
                              nT = nT, b.other.part = b.other.part)
      # If there is ... nononono

      b.draws[i.rep - burn.in, ] <- unlist(b.prev)
      g.draws[i.rep - burn.in, ] <- unlist(g.prev)
    }
  }
  
  # 4. Marginalized likelihood
  # 4.1 Chib's method
  if (n == 2) {
    # Chib's 3 blocks method
    # 4.1.1 Calculate dconditional for the first b:
    b.star <- Chib.values$b
    g.star <- Chib.values$g
    
    dconditional.b1 <- log(mean(exp(dconditional.b1.sample)))
    
#     4.1.2 Calculate dconditional for the second b:
    A.restr <- free.params.to.matrix(params = b.star, 
                                     restrictions = restrictions$H.type$A)
#     dconditional.b2 <- dconditional.bi(b = b.star, i = 2, post.params = post.params,
#                                       nT = nT, restr.mat = restrictions$H.type$A, 
#                                       q.i = length(b.star[[1]]))
#     b.other.part <- - 0.5 * t(b.star[[1]]) %*% solve(post.params[[1]]$S) %*% b.star[[1]]
    b.other.part <- NULL
    dconditional.b2 <- dconditional.bi(b.i = b.star[[2]],
                    A.restr = A.restr[, - 2, drop = FALSE],
#                     A.restr = A.restr,
                    post.params = post.params[[2]], nT = nT,
                    restr.mat = restrictions$H.type$A[[2]],
                    b.other.part = b.other.part)

    # 4.1.3 Calculate dconditional for g:
    dconditional.g <- sum(sapply(1:n, function(i) {
      ifelse(is.null(post.params[[i]]$P), 0, 
             log(dmvnorm(x = c(g.star[[i]]), mean = c(post.params[[i]]$P %*% b.star[[i]]),
                         sigma = post.params[[i]]$H)))
    }))
    
    # Marginal likelihood evaluation
    lf.my <- lik.fun(b = b.star, g = g.star)
    prior.dens <- dprior(b = b.star, g = g.star)
    dposterior <- dconditional.b2 + dconditional.b1 + dconditional.g
    
    chib.exact <- lf.my + prior.dens - dposterior
  } else {
    chib.exact <- NULL # This version does not calculate more than 3-block GS densities.
  }
  
  # 4.2 BIC method
  var.bic <- lik.fun(b = start.values$b, g = start.values$g) +
    (sum(!restrictions$exact$A) + sum(!restrictions$exact$B))*log(nT)

  # 5. Save all draws as MCMC objects
  b.draws.mcmc <- mcmc(b.draws); colnames(b.draws.mcmc) <- colnames(b.draws)
  g.draws.mcmc <- mcmc(g.draws); colnames(g.draws.mcmc) <- colnames(g.draws)

  return(new(Class = "BSVAR.out", draws = list(b = b.draws.mcmc,
                                               g = g.draws.mcmc),
                                  marginal.likelihood = list(dlikelihood = lf,
                                                             dprior = prior.dens,
                                                             dconditional.b1 = dconditional.b1,
                                                             dconditional.b2 = dconditional.b2,
                                                             dconditional.g  = dconditional.g,
                                                             chib.exact = chib.exact,
                                                             dposterior = dposterior,
                                                             BIC = var.bic,
                                                             lf.my = lf.my),
                                  irf = NULL,
                                  timer = proc.time() - start.timer))
}


# ---- 4. Bayesian Model Averaging functions ----
# ==== 4.1 Draw one restriction ====
decode.restrictions <- function(integer, n, m) {
  constraints <- matrix(1, ncol = n, nrow = m)
  for (i in 1:m) {
    for (j in 1:n) {
      if (i == j) next
      if (integer %% 2) {
        constraints[i, j] <- 0
      }
      integer <- integer %/% 2
    }
  }
  constraints
}

# ==== 4.2 Calculate Bayes Factor ====
construct.bayes.factors <- function(log.marg.lik, method = "Chib", 
                                    weights = NULL, prob = .975) {
  if (method == "direct") {
    const <- mean(log.marg.lik[[2]])
    mean  <- mean(exp(log.marg.lik[[1]] - const))/mean(exp(log.marg.lik[[2]] - const))
    
    max <- quantile(exp(log.marg.lik[[1]] - const), prob = min(prob, 1 - prob))/
      quantile(exp(log.marg.lik[[2]] - const), prob = max(prob, 1 - prob))
    min <- quantile(exp(log.marg.lik[[1]] - const), prob = max(prob, 1- prob))/
      quantile(exp(log.marg.lik[[2]] - const), prob = min(prob, 1 - prob))
    
    return(list(max = max, min = min, mean = mean))
  }
  if (method == "Chib") {
    return(exp(log.marg.lik[[1]] - log.marg.lik[[2]]))
  }
  if (method == "importance") {
    if (is.null(weights)) stop("Weights must be specified")
    const <- mean((log.marg.lik[[2]]))
  }
  if (method == "harmonic_mean") {
  }
}

# ==== 4.3 Simulation loops ====
direct.draw.model <- function(ind.model, A, B, prior, 
                              eps, Z, calc.direct.marg = FALSE,
                              name.digits = 2, notification = NULL) {
  
  # 0. Initialization
  n.obs <- nrow(Z)
  n <- ncol(A)
  num.exog <- ncol(Z)
  m <- 2^(n*(n - 1 + num.exog))
  
  # 1. Generate restrictions 
  true.model <-  decode.restrictions(ind.model, n, n + num.exog)
  A.restr <- true.model[1:n, ]
  B.restr <- true.model[n + 1:num.exog, , drop = FALSE]
  
  # 2 Generate artificial data set
  Y <- (Z %*% (B.restr*B) + eps) %*% solve(A.restr*A)

  # 3. Simulate under prior model.
  estimated.prior.models <- lapply(1:m - 1, function(k) {
    # ##### 6.4.0. Define restrictions  for prior model #####
    model <- decode.restrictions(k, n, n + num.exog)
    A.restr <- model[1:n, ] == 0
    B.restr <- model[n + 1:num.exog, , drop = FALSE] == 0
    
    # #### 3.1 Set restrictions ####
    exact.restrictions <- list(A = A.restr,
                               B = B.restr)
    # 3.1.1 Define H-type restrictions ####
    Htype.restrictions <- lapply(seq_along(exact.restrictions), 
                                 function(i) exact.to.H(exact.restrictions[[i]]))
    names(Htype.restrictions) <- c('A', 'B')
    # 3.1.2 Combine restrictions on A and on B ####
    restrictions <- list(exact  = exact.restrictions,
                         H.type = Htype.restrictions)
    # #### 3.2 Set structural prior ####
    str.prior <- structural.priors.on.unrestricted.params(prior, restrictions$H.type)
    
    # #### 3.3 Define start values
    start.values <- define.starting.values(Y, Z, restrictions, type = "both")
    
    # #### 3.4 Scale in order to compare with Nickolay Gennadievich's approach ####
    # (optional, if deterministic term is considered, turn off this code)
    Y <- as.matrix(scale(Y, scale = FALSE)) 
    Z <- as.matrix(scale(Z, scale = FALSE)) 
    # #### 3.5 Calculate BSVAR ####
    test.passed <- FALSE; ind.test <- 1
    while (!test.passed) {
      BSVAR <- Waggoner.Zha.gibbs.sampler(Y = Y, Z = Z, 
                                          restrictions = restrictions, 
                                          prior = str.prior, start.values = start.values$mle,
                                          Chib.values = start.values$mle,
                                          spec.notification = paste0(notification, "Testing Model ", 
                                          sprintf(paste0("%0", name.digits, "d"), k + 1)))
      if (k %in% c(0,  12)) {
        H.test <- heidel.diag(BSVAR@marginal.likelihood$dlikelihood)
      } else {
        H.test <- heidel.diag(mcmc(cbind(data.frame(BSVAR@draws$b), 
                               data.frame(BSVAR@draws$g))))
      }
      
      test.passed <- all(as.logical(c(H.test[, c(1, 4)])))
      if (ind.test >= 3) {
        warning(paste0("Convergence is not achieved for model M", 
                       ind.model + 1, "testing model M", k ))
        test.passed <- TRUE
      }
      ind.test <- ind.test + 1
    }
    if (calc.direct.marg) {
      lik.fun <- construct.log.lik.fun(Y = Y, Z = Z,
                                       restrictions = restrictions$H.type,
                                       type = "Waggoner-Zha")
      marg <- direct.marginal.density(lik.fun, str.prior, spec.notification = paste("Testing Model", k + 1))
    } else marg <- NULL
    
    return(list(BSVAR = BSVAR, direct.marg = marg, mle = start.values$estim$mle))
    })
  names(estimated.prior.models) <- c(paste0("M", sprintf(paste0("%0", name.digits, "d"), 
                                                         1:m)))
  return(estimated.prior.models)
}

# ==== 4.4 Summary of model output ====
# #### 4.4.1 Coefficients diagnostics plots ####
# CAUTION: Very slow 
plot.coefficients <- function(output, ind.model, true.coeff, file) {
  require(ggplot2)
  require(gridExtra)
  
  restrictions <- output[[ind.model]]$mle$constraints.on.coefficients
  mle.estim <- output[[ind.model]]$mle$A
  if (nrow(mle.estim) > ncol(mle.estim)) {
    mle.estim[(ncol(mle.estim) + 1):nrow(mle.estim), ] <-
      (- 1)*mle.estim[(ncol(mle.estim) + 1):nrow(mle.estim), ]
  }
  
  plot.env <- new.env()
  plot.env$S <- cbind(data.frame(output[[ind.model]]$BSVAR@draws$b),
             data.frame(output[[ind.model]]$BSVAR@draws$g))
  plot.env$mle.estim <- mle.estim
  plot.env$true.coeff <- true.coeff
  
  dummy.frame <- data.frame(X = -2:2, Y = seq(from = 0, to = 1, length = 5))
  plotlist <- vector("list", nrow(mle.estim)*ncol(mle.estim))
  names(plotlist) <- paste0("c", 1:(nrow(mle.estim)*ncol(mle.estim)))
  
  pdf(file, width = 14, height = 7)
  grid.newpage()
  pushViewport(
    viewport(x = unit(0.5, "npc"), y = unit(0.5, "npc"),
             layout = grid.layout(ncol(restrictions), nrow(restrictions)),
             width = unit(1, "npc"), # aspect ratio preserved
             height = unit(0.5, "npc")))
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  k <- 1
  for (i in 1:nrow(mle.estim)) {
    for (j in 1:ncol(mle.estim)) {
      if (restrictions[i, j]) {
        command <- paste0('plot.temp <- ggplot(S, 
          aes_string(x = paste0(ifelse(i >  ncol(true.coeff), "g", "b"),', i - 
          ifelse(i > ncol(mle.estim), ncol(mle.estim), 0), ',',  j, '))) + 
          geom_density(fill = "blue", col = "blue", alpha = 0.2) + 
          geom_vline(xintercept = mean(S[, paste0(ifelse(i >  ncol(true.coeff), 
          "g", "b"),', i - ifelse(i > ncol(mle.estim), ncol(mle.estim), 0), ',', 
          j, ')]), col = "blue", linetype = "dashed") +
          geom_vline(xintercept = mle.estim[', i, ',',  j, '], col = "red", linetype = "dashed") + 
          geom_vline(xintercept = true.coeff[', i, ',',  j, '], col = "red") +
          xlab(paste0(ifelse(i >  ncol(true.coeff), "g", "b"),', i - ifelse(i > ncol(mle.estim),
          ncol(mle.estim), 0), ',',  j, ')) +
          theme_bw()')
        eval(parse(text = command), envir = plot.env)
        plot.temp <- plot.env$plot.temp
        print(plot.temp, vp = vplayout(j, i))
        k <- k + 1
      } else {
        command <- paste0('plot.temp <- ggplot(dummy.frame, aes(x = X, y = Y)) + 
          geom_blank() +
          geom_vline(xintercept = 0, col = "red") +
          ylab("density") +
          xlab(paste0(ifelse(i > ncol(true.coeff), "g", "b"),', i, ',',  j, '))  + 
          theme_bw()')
        eval(parse(text = command))
        print(plot.temp, vp = vplayout(j, i))
      }
    }
  }
  dev.off()
}


# ==== Prior model probabilities ====
define.prior.model.probabilities <- function(vertex.priors, n, k, 
                                             is.prob.equal.zero = TRUE) {
  m <- 2^(n*(n - 1 + k))
  sapply(1:m, function(ind.model) {
    if (is.prob.equal.zero) {
      model <- decode.restrictions(ind.model - 1, n, n + k)
    } else {
      model <- 1 - decode.restrictions(ind.model - 1, n, n + k)
    }
    
    prod(sapply(seq_along(model), function(i) {
      vertex.priors[i]^(model[i])*(1 - vertex.priors[i])^(1 - model[i])
      }))
    })
}

# ==== Posterior model probabilities ====
calc.posterior.prob <- function(i, mean.likelihoods, prior.probabilities = NULL) {
  if (is.null(prior.probabilities)) {
    prior.probabilities <- rep(0.5, length(mean.likelihoods))
  }
  if (mean.likelihoods[i] < 1e-50) {
    return(0)
  } else {
    return(mean.likelihoods[i] * prior.probabilities[i]/
             sum(sapply(1:length(mean.likelihoods), function(j) {
               mean.likelihoods[j]*prior.probabilities[j]})))
  }
}

# ==== Function to summary and plot results ====
summary.and.plot.prior.models <- function(output, filename, wd = getwd(),
                                          graph.format = ".pdf",
                                          plot.densities = FALSE, 
                                          plot.diagnostics = FALSE,
                                          marg.type = "Chib",
                                          A = NULL, ind.true.model = NULL) {
  require(ggplot2)
  setwd(wd)
  m <- length(output)
  # Likelihood function estimation from Gibbs Sampler
  mean.log.likelihoods <- sapply(seq_along(output), function(i) {
    if (marg.type == "direct") {
      mean(output[[i]]$BSVAR@marginal.likelihood$dlikelihood)
    } else if (marg.type == "Chib") {
      output[[i]]$BSVAR@marginal.likelihood$chib.exact
    } else if (marg.type == "BIC") {
      output[[i]]$BSVAR@marginal.likelihood$BIC
    }
  })
  cat("Log likelihoods\n")
  print(mean.log.likelihoods)
  
  # Maximum likelihood estimation
  max.lf <- sapply(seq_along(output), function(i) {
    - output[[i]]$mle$optimization$objective
  })
  
  # Comparative plot 
  plot.lf <- data.frame(lf = c(mean.log.likelihoods, max.lf))
  plot.lf$type <- c(rep("From Gibbs", m), rep("Max Lik", m))
  plot.lf$model <- names(output)
  
  lf.gs.plot <- ggplot(plot.lf, aes(x = model, y = lf, fill = type)) +
    ylab("Likelihoods") + 
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw()
  ggsave(lf.gs.plot, file = paste0("likplot_", filename, graph.format))
  
  # Calculate posterior probabilities using draws from the Gibbbs Sampler
  const <- max(mean.log.likelihoods)
  mean.likelihoods <- sapply(seq_along(output), function(i) {
    if (marg.type == "direct") {
      mean(exp(output[[i]]$BSVAR@marginal.likelihood$dlikelihood - const))
    } else if (marg.type == "Chib") {
      exp(output[[i]]$BSVAR@marginal.likelihood$chib.exact - const)
    } else if (marg.type == "BIC") {
      mean(exp(output[[i]]$BSVAR@marginal.likelihood$BIC - const))
    }
  })
  
  cat("Likelihoods\n")
  print(mean.likelihoods)
  
  # Bad parametrization. Check inference/
  cond <- complete.cases(mean.likelihoods)&(mean.likelihoods < 1e20)
  if (any(!cond)) {
    stop(paste0("Marginal likelihoods of the following models are not correctly calculated: \n",
                   seq_along(mean.likelihoods )[cond], "\n results of these models have been truncated"))
#     mean.likelihoods[!cond] <- 0
  }
  
  prob.gs <- sapply(seq_along(mean.likelihoods), function(i) {
    calc.posterior.prob(i, mean.likelihoods)
  }); names(prob.gs) <- names(output)
  
  # View on screen to make the user feel sad
  cat(paste0(filename, " posterior probabilities: \n"))
  print(prob.gs)
  
  # Plot posterior probabilities
  prob.data <- data.frame(prob.gs); colnames(prob.data) <- "P"
  prob.data$model <- paste0("M", c(paste0(0, 1:9), 10:16))
  prob.data$pc <- ifelse(prob.gs> 0.001, as.character(signif(prob.gs, 2)),
                         "")
  prob.data$moralization <- c(rep("C1", 3), "C2", "C1", "C3", "C1", "C5",
                                   "C1", "C1", "C4", "C6", rep("C7", 3), "C8")
  
  
  prob.plot <- ggplot(prob.data) +
    ylab("Posterior probability")+
    ggtitle(paste0("True model: M", sprintf("%02d", ind.true.model))) + 
    geom_bar(aes(x = model, y = P, label = pc, fill = moralization), stat = "identity") + 
    geom_text(aes(x = model, y = P, label = pc), size = 3, vjust = 0) +
    theme(axis.text.x=element_blank(), axis.ticks=element_blank(),
          axis.title.x=element_blank(), legend.title=element_blank(),
          axis.title.y=element_blank(), text = element_text(size = 10)) +
    scale_fill_hue(l=40, guide = FALSE) +
    theme_bw()
  ggsave(prob.plot, file = paste0("posterior_prob_plot_", filename, graph.format))
  
  # Plot densities (optional because of high CPU memory consumption)
  if (plot.diagnostics) {
    for (i in seq_along(output)) {
      require(ggmcmc)
      K <- ggs(mcmc(cbind(data.frame(output[[i]]$BSVAR@draws$b,
               data.frame(output[[i]]$BSVAR@draws$g)))))
      
      ggmcmc(K, file = paste0("diag_", i, "_", filename, ".pdf"), 
             plot = c("ggs_traceplot()","ggs_compare_partial()", "ggs_running()",
             "ggs_autocorrelation()", "ggs_crosscorrelation()", 
             "ggs_geweke()"))
      
      lik <- output[[i]]$BSVAR@marginal.likelihood$dlikelihood
      lik <- data.frame(lik)
      
#       # Open new device
      lik.plot <- ggplot(data.frame(lik), aes(x = lik)) + 
        geom_density(fill = "blue", col = "blue", alpha = 0.2) +
        geom_vline(xintercept = max.lf[[i]], col = "red") +
        xlab("likelihood") +
        theme_bw()
      ggsave(lik.plot, file = paste0("lik", i, filename, ".pdf"))
      
      # Save coefficients of the true model  
    }
  }
  # Plot coefficients densities
  if (plot.densities) {
    plot.coefficients(output = output, ind.model = ind.true.model,
                      true.coeff = A, 
                      file = paste0(filename, "_densplot.pdf"))
  }
  # Save in file 
  save(file = paste0("variables_", filename, ".RData"), 
             list = c("output", "prob.gs", "max.lf"))
  return(list(prob.gs = prob.gs, max.lf = max.lf))
}

# ==== 4.2 Tables creator ====
prob.results.to.latex.tab <- function(c.vec, num.models, wd) {
  require(xtable)
  old.wd <- getwd()
  res <- vector("list", length(c.vec))
  k <- 1
  for (c in c.vec) {
    temp <- matrix(NA, num.models, num.models)
    for (ind.true.model in 1:num.models) {
      setwd(old.wd)
      load(paste0("variables_c", 
                  ifelse(c >= 1, sprintf("%02d", c),
                         paste0("0dot", c)), "_RRM", 
                  sprintf("%02d", ind.true.model),
                  ".RData"))
      temp[ind.true.model, ] <- prob.gs
    }
    
    colnames(temp) <- paste0("M", sprintf("%02d", 1:num.models))
    rownames(temp) <- paste0("M", sprintf("%02d", 1:num.models))
    
    setwd(wd)
    temp.table <- xtable(temp)
    res[[k]] <- temp.table
    print(temp.table,  file=paste0("Prob2_45_", ifelse(c > 1, sprintf("%02d", c),
                                          paste0("0dot", c*10)), ".tex"),
          size = "footnotesize")
    k <- k + 1
  }
  setwd(old.wd)
  return(res)
}
                                 

# ---- 5. MC3 procedure ----
# ==== 5.1 Function to define model order, neighbours and parents ====
define.order <- function(n, k, type = "naive") {
  # Return: adjacency matrix of model choice graph
  m <- 2^(n*(n - 1 + k))
  all.models <- lapply(1:m, function(i) {
    decode.restrictions(i - 1, n, n + k)
    })
  adjacency <- matrix(0, m, m)
  for (i in 1:m) {
    for (j in (1:m)[-i]) {
      if (type == "naive") {
        cond <- sum(abs(all.models[[i]] - all.models[[j]])) == 1
      } else if (type == "moralized.by.concentration") {
        # Indian codewriter mode on
        params.i <- all.models[[i]] %*% t(all.models[[i]]) != 0
        params.j <- all.models[[j]] %*% t(all.models[[j]]) != 0
        cond <- sum(abs(params.i[upper.tri(params.i)] -
                          params.j[upper.tri(params.j)])) == 1
      } 
      if (cond) {
        adjacency[i, j] <- 1
        adjacency[j, i] <- 1
      }
    }
  }
  return(adjacency != 0)
}

# ==== 5.2 MC3 procedure function ====
mc3.bsvar <- function(Y, Z, prior, 
                      prior.probabilities = NULL, start.model = NULL,
                      order = NULL, burn.in = 1e2, sample.save = 1e3,
                      calc.direct.marg = FALSE,
                      name.digits = 2, spec.notification = NULL,
                      wd = getwd(), marg.type = "Chib") {
  # 0. Initialization
  n.obs <- nrow(Z)
  n <- ncol(A)
  k <- ncol(Z)
  m <- 2^(n*(n - 1 + k))
  posterior <- matrix(NA, ncol = 1, nrow = sample.save)
  
  sample.size <- burn.in + sample.save
  if (is.null(prior.probabilities) ) {
    prior.probabilities <- rep(0.5, m)
  }
  if (is.null(start.model)) {
    ind.model <- 1
  } else ind.model <- start.model
  
  if (is.null(order)) {
    order <- define.order(n, k, type = "naive") 
  }
  
  for (i.rep in 1:sample.size) {
    if (i.rep %% 1e3 == 0)  {
      cat(spec.notification, 'MC3 Iteration', i.rep, '/', sample.size, '\n')
    }
    # ##### 3. Propose model ####
    # ##### Randomly select model from neighbourhood ####
    if (i.rep == 1) {
      ind.model.proposal <- ind.model
    } else {
      ind.model.proposal <- (1:m)[order[ind.model, ]][floor(runif(1)*
                                  sum(order[ind.model, ])) + 1]
    }
    
    # ##### 5. Find posterior probability ####
    all.files <- list.files(wd)
    # ##### 6. If there is no calculated probability, calculate it #####
    if (!(paste0(spec.notification, "mc3_posterior_M", sprintf("%02d", ind.model.proposal), ".RData") %in% all.files)) {
      # ##### 4. Define restrictions  for prior model #####
      model <- decode.restrictions(ind.model.proposal - 1, n, n + k)
      
      A.restr <- model[1:n, ] == 0
      B.restr <- model[(n + 1):nrow(model), , drop = FALSE] == 0
      
      # #### 4.1 Set restrictions ####
      exact.restrictions <- list(A = A.restr,
                                 B = B.restr)
      # 4.1.1 Define H-type restrictions ####
      Htype.restrictions <- lapply(seq_along(exact.restrictions), 
                                   function(i) exact.to.H(exact.restrictions[[i]]))
      names(Htype.restrictions) <- c('A', 'B')
      # 4.1.2 Combine restrictions on A and on B ####
      restrictions <- list(exact  = exact.restrictions,
                           H.type = Htype.restrictions)
      # #### 4.2 Set structural prior ####
      str.prior <- structural.priors.on.unrestricted.params(prior, restrictions$H.type)
      
      # #### 4.3 Define start values
      start.values <- define.starting.values(Y, Z, restrictions, type = "both")
      # #### 3.4 Scale in order to compare with Nickolay Gennadievich's approach ####
      # (optional, if deterministic term is considered, turn off this code)
      Y <- as.matrix(scale(Y, scale = FALSE)) 
      Z <- as.matrix(scale(Z, scale = FALSE)) 
      # #### 3.5 Calculate BSVAR ####
      BSVAR <- Waggoner.Zha.gibbs.sampler(Y = Y, Z = Z, 
                                          restrictions = restrictions, 
                                          prior = str.prior, start.values = start.values$mle,
                                          spec.notification = paste0(spec.notification, "Proposal Model ", 
                                          sprintf(paste0("%02d"), ind.model.proposal + 1)))
      output <- list(BSVAR = BSVAR, mle = start.values$estim$mle)
      save(file = paste0(spec.notification, "mc3_posterior_M", sprintf("%02d", ind.model.proposal), ".RData"),
                 list = "output")
    } else load(paste0(spec.notification, "mc3_posterior_M", sprintf("%02d", ind.model.proposal), ".RData"))
    # #### 3.6 Extract model posterior probability ####
    if (marg.type == "Chib"){
      proposed.log.posterior.probability <- output$BSVAR@marginal.likelihood$chib.exact
    } else if (marg.type == "direct") {
      proposed.log.posterior.probability <- mean(output$BSVAR@marginal.likelihood$dlikelihood)
    }
   
    # Define acceptance ratio of the MCMCMH
    if (i.rep == 1) {
      state.log.posterior.probability <- proposed.log.posterior.probability
    } else  {
      if (is.null(prior.probabilities)) {
        acceptance.ratio <- min(1, exp(proposed.log.posterior.probability -
                                           state.log.posterior.probability))
      } else {
        prior.proposal.subject.to.state <- prior.probabilities[ind.model.proposal]/
          prior.probabilities[order[ind.model, ]]
        prior.state.subject.to.proposal <- prior.probabilities[ind.model]/
          prior.probabilities[order[ind.model.proposal]]
        acceptance.ratio <- min(1, prior.state.subject.to.proposal * 
                   exp(proposed.log.posterior.probability - state.log.posterior.probability) *
                   prior.probabilities[ind.model.proposal]/(prior.probabilities[ind.model]*
                   prior.proposal.subject.to.state ))
      }
      # Accept proposal with defined probability
      if (runif(1) <= acceptance.ratio) ind.model <- ind.model.proposal
      
      # Save results 
      if (i.rep > burn.in) {
        posterior[i.rep - burn.in, ] <- ind.model
      }
    }
    # #### In order to clean space on a disk ####
    rm(list = "output")
  }
  return(posterior)
}

mc3.draw.model <- function(ind.model, A, B, prior, 
                              eps, Z, notification = "MC3") {
  start.timer <- proc.time()
  # 0. Initialization
  n.obs <- nrow(Z)
  n <- ncol(A)
  k <- ncol(Z)
  m <- 2^(n*(n - 1 + k))
  
  # 1. Generate restrictions 
  true.model <-  decode.restrictions(ind.model, n, n + k)
  A.restr <- true.model[1:2, ]
  B.restr <- true.model[3, , drop = FALSE]
  
  # 2 Generate artificial data set
  Y <- (Z %*% (B.restr*B) + eps) %*% solve(A.restr*A)
  
  # 3 Start MC3 loop
  mc.draws <- mc3.bsvar(Y, Z, prior, spec.notification =  
          paste0("M_true", sprintf("%02d", ind.model + 1), "_"))
  timer <- proc.time() - start.timer
  return(list(mc.draws = mc.draws, timer = timer))
}
