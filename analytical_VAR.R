library(MCMCpack)
setwd("/home/shere-khan/Dropbox/Research/bvar_codes/myBVAR")

Y.raw <- read.table("Yraw.dat", header = FALSE)
colnames(Y.raw) <- c("oh", "ah", "uh")

# ---- 0. Generate function ----
# matrix.to.vector
matrix.to.vector <- function(A, ch = ".to.") {
  res <- matrix(c(A), ncol = 1)
  rownames(res) <- c(sapply(colnames(A), function(c) {
    paste0(c, ch, rownames(A))
  }))
  return(res)
}
# generate lags
lag.gen <- function(Y, lag.order = 1) {
  Y <- as.matrix(Y)
  res <- matrix(NA, ncol = lag.order*ncol(Y), nrow = nrow(Y))
  
  for (i in 1:lag.order) {
    res[(lag.order + 1):nrow(Y), (i - 1)*ncol(Y) + 1:ncol(Y)] = 
      Y[(1 + lag.order - i):(nrow(Y) - i), 1:ncol(Y)]
  }
  colnames(res) <- c(sapply(1:lag.order, 
                     function(i) paste0("L", i, '.', colnames(Y))))
  return(as.data.frame(res))
}

var.ols <- function(Y, X) {
  nT <- nrow(Y)
  A  <- solve(t(X)%*%X)%*%(t(X)%*%Y)
  
  SSE       <- t(Y - X%*%A)%*%(Y - X%*%A)
  sigma.sq     <- SSE/(nT - k)
  return(list(A        = A, 
              vec.A    = matrix.to.vector(A),
              SSE      = SSE,
              sigma.sq = sigma.sq,
              k        = nrow(A),
              n        = ncol(A),
              nT       = nT))
}

# Generate prior function
construct.Minnesota.prior <- function(Y, A.prior, ols, 
             use.lag.directly = FALSE,
             hyper.par = list(a.bar1 = 0.5,
                              a.bar2 = 0.5,
                              a.bar3 = 1,
                              a.bar4 = 100)) {
      lag.order <- nrow(A.prior)/ncol(A.prior)
      k <- ols$k
      # 1. Construct single AR process and  generate sigma.sq 
      sigma.sq <- sapply(1:ncol(Y), function(i) {
        lag.Y <- as.matrix(na.omit(lag.gen(Y[, i, drop = FALSE], lag.order)))
        Y     <- as.matrix(Y[(lag.order+ 1):nrow(Y), i])
        alpha <- solve(t(lag.Y)%*%lag.Y)%*%(t(lag.Y)%*%Y)
        return(sum((Y - lag.Y%*%alpha)^2)/(nrow(Y) - lag.order + 1))
      })
      names(sigma.sq) <- colnames(Y)
      
      # 2. Define standart deviations for parameters priors.
      V <- 0*A.prior
      own.lag <- !(lower.tri(A.prior)|upper.tri(A.prior))
      
      if ("const" %in% rownames(A.prior)) {
        V["const", ] <- hyper.par$a.bar4*sigma.sq # variance for coefficient for constant 
      }
      if ("trend" %in% rownames(A.prior)) {
        V["trend", ] <- hyper.par$a.bar4*sigma.sq # variance for coefficient for trend
      }
      
      # 3. Define variance for own lags of variables
      V[own.lag] <- hyper.par$a.bar1/(ifelse(use.lag.directly, 
                          1:3, ceiling(1:3/ncol(Y)))^(2*hyper.par$a.bar3))
      
      # 4. Define variance for lags of other variables
      for (i in 1:1:ncol(V)) {
        for (j in 1:ncol(V)) {
          for (p in 1:lag.order) {
            if ((p == 1)&(i == j)) {
            } else {
              V[(p - 1)*ncol(Y) + i, j] <- hyper.par$a.bar2*
                sigma.sq[j]/(sigma.sq[i]*(p)^(2*hyper.par$a.bar3))
            }
          }
        }
      }
      
      # 5. Create vectorization of A prior
      vec.A.prior = matrix.to.vector(A.prior)
      
      # 6. Diagonalize variance-covariance matrix
      V.prior <- diag(c(V))
      colnames(V.prior) <- rownames(vec.A.prior)
      rownames(V.prior) <- rownames(vec.A.prior)
      
      return(list(A.prior     = A.prior,
                  vec.A.prior = vec.A.prior,
                  V.prior     = V.prior,
                  sigma.sq    = ols$sigma.sq))
}

# Construct posterior
construct.posterior <- function(X, Y, prior, prior.type, ols) {
  if (prior.type == "Noninformative") {
    # Posterior of coeff|Data ~ Multi-T(kronecker(SSE, solve(t(X)%*%(X))), A.ols, nT - k)
    V.post     <- solve(t(X)%*%X)
    vec.A.post <- ols$vec.A
    A.post     <- ols$A
    
    # Posterior of Sigma|Data ~ inv-Wishart(SSE, T - k)
    S.post <- SSE
    v.post <- nT - k
    
    # Mean and variance of posterior 
    coeff.mean <- A.post
    coeff.var  <- kronecker(S.post, V.post)/(v.post - ols$n - 1)
  } else if (prior.type == "Minnesota") {
    # For coeff
    V.post <- solve(solve(prior$V) + kronecker(solve(prior$sigma.sq)))
    A.post <- V.post%*%(solve(V .prior)%*%A.prior + 
                          kronecker(solve(prior$sigma.sq), t(X)%*%X)%*%ols.A)
    
    A.mean <- A.post
  } else if (prior.type == "Normal-Wishart") {
    # For coeff
    V.post <- solve(solve(prior$V) + t(X)%*%X)
    A.post <- V.post %*% (solve(prior$V)%*%prior$A + t(X)%*%X%*%ols$A)
    vec.A.post <- matrix.to.vector(A.post)
    
    # For Sigma
    S.post <- ols$SSE + prior$S + t(ols$A)%*%t(X)%*%X%*%ols$A + 
      t(prior$A)%*%solve(V.prior)%*%prior$A - t(A.post)%*%
      (solve(prior$V) + t(X)%*%X)%*%A.post
    v.post <- ols$nT + prior$v 
  }
  return(list(coeff.mean = coeff.mean, coeff.var = coeff.var))
}

# ---- 1. Specification of VAR model ----
deterministic.terms      <- "const"
lag.order                <- 2
forecasting              <- TRUE
forecast.method          <- "direct"
h                        <- 4
m                        <- nrow(Y.raw)
n                        <- ncol(Y.raw)

prior.type <- "Minnesota"

# ---- 2. Data handling ----
# ==== 2.1 Generate Y in training sample and Y in prediction sample ==== 
if (forecasting) {
  if (h <= 0) stop("Incorrect forecasting length")
  
  if (forecast.method == "direct") {
    Y.ts <- Y.raw[(h + 1):m, , drop = FALSE]
    Y.ps <- Y.raw[2:(m - h), , drop = FALSE]
    m    <-  m - h - 1
  } else if (forecast.method == "iterated") {
    Y.ts <- Y.raw
    Y.ps <- Y.raw
  }
} else {
  Y.ts <- Y.raw
  Y.ps <- Y.raw
}

# ==== 2.2 Generate lagged Y
X.ts <- lag.gen(Y.ps, lag.order)
X.ts <- na.omit(X.ts)

# ==== 2.3 Generate deterministic terms ====
if ("const" %in% deterministic.terms)  X.ts[["const"]] <- rep(1, nrow(X.ts))
if ("trend" %in% deterministic.terms)  X.ts[["trend"]] <- 1:nrow(X.ts)

nT3 <- nrow(X.ts)
k   <- ncol(X.ts)

Z.ts <- kronecker(diag(rep(1, n)), as.matrix(X.ts))
Y.ts <- Y.ts[(lag.order + 1):m, ]
nT   <- m - lag.order

# ==== 2.4 Forecasting setup ====
if (forecasting) {
  if (forecast.method == "direct") {
    Y <- as.matrix(Y.ts[1:(nrow(Y.ts) - 1), ])
    X <- as.matrix(X.ts[1:(nrow(X.ts) - 1), ])
    Z <- kronecker(diag(rep(1, n)), X)
    nT <- nT - 1
  } else {
    Y <- as.matrix(Y.ts[1:(nrow(Y.ts) - h), ])
    X <- as.matrix(X.ts[1:(nrow(X.ts) - h), ])
    Z <- kronecker(diag(rep(1, n)), X)
    nT <- nT - h
  }
} else {
  Y <- as.matrix(Y.ts)
  X <- as.matrix(X.ts)
  Z <- as.matrix(Z.ts)
}

# ---- 3. Priors -----
# ==== 3.1 OLS Estimation ====
ols <- var.ols(Y, X)

# ==== 3.2 Prior hyperparameters for BVAR ====
# ##### 3.3 Prior on coefficients #####
if (prior.type == "Minnesota") {
  A.prior <- 0*ols$A
#   diag(A.prior) <- 1
  prior.coeff <- construct.Minnesota.prior(Y = Y.raw, A.prior = A.prior, 
                                     ols = ols, use.lag.directly = FALSE,
                                     hyper.par = list(a.bar1 = 0.5,
                                                      a.bar2 = 0.5,
                                                      a.bar3 = 1,
                                                      a.bar4 = 100))
} else if (prior.type == "Normal-Wishart") {
  # Non-informative priors
  prior.coeff <- vector("list", length = 4)
  names(prior.coeff)  <- c("A.prior", "vec.A.prior",
                    "V.prior",     "sigma.sq")
  prior.coeff$A.prior <- 0*ols$A
  
  prior.coeff$V.prior <- diag(rep(10, nrow(ols$vec.A)))
  colnames(prior.coeff$V.prior) <-  rownames(ols$vec.A)
  rownames(prior.coeff$V.prior) <-  rownames(ols$vec.A)
} else if (prior.type == "Noninformative") {
  prior.coeff <- NULL
}

# #### 3.4 Wishart priors ####
prior.concentration        <- vector("list", 2)
names(prior.concentration) <- c("v", "S") 
prior.concentration$v <- n

prior.concentration$S <- diag(rep(10, n))
colnames(prior.concentration$S) <- colnames(Y)
rownames(prior.concentration$S) <- colnames(Y)

# ---- 4. Posteriors ----
library(vars)
data(Canada)
mclist <- mcmc.list(mcmc(Canada[, 1]), mcmc(Canada[, 2]))
D <- gelman.diag(mclist)
gelman.plot(mclist)


prior <- function(x) {
  a <- 1
  b <- 4
  return(a + b)
}
