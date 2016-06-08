# ---- 1. Data handling ----
# matrix.to.vector
matrix.to.vector <- function(A, sep = c("c_", ".r_")) {
  res <- matrix(c(A), ncol = 1)
  rownames(res) <- c(sapply(colnames(A), function(c) {
    paste0(sep[1], c, sep[2], rownames(A))
  }))
  return(res)
}

# vector.to.matrix 
vector.to.matrix <- function(x, sep = c("c_", ".r_"),
                             nrow = NULL, ncol = NULL) {
  split <- strsplit(x = gsub(x = rownames(x),
                            pattern = sep[1], 
                            replacement = ''),
                   split = sep[2])
  names.col <- unique(sapply(seq_along(split), function(i) split[[i]][1]))
  names.row <- unique(sapply(seq_along(split), function(i) split[[i]][2]))
  if ((is.null(ncol))&(is.null(nrow))) {
    ncol <- length(names.col)
    nrow <- length(names.row)
  } else if (is.na(ncol)) {
    ncol <- 1
  } else if (is.na(nrow)) {
    nrow <- 1
  }
  res <- matrix(x, nrow = nrow, ncol = ncol)
  colnames(res) <- names.col
  rownames(res) <- names.row
  
  return(res)
}

# generate lags
lag.gen <- function(Y, lag.order = 1) {
  Y <- as.matrix(Y)
  res <- matrix(NA, ncol = lag.order*ncol(Y), nrow = nrow(Y))
  
  for (i in 1:lag.order) {
    res[(i + 1):nrow(Y), (i - 1)*ncol(Y) + 1:ncol(Y)] = 
      Y[1:(nrow(Y) - i), 1:ncol(Y)]
  }
  colnames(res) <- c(sapply(1:lag.order, 
                            function(i) paste0("L", i, '.', colnames(Y))))
  return(as.data.frame(res))
}

# Simplified OLS function
ols.VAR <- function(Y, X, method = "direct") {
  nT <- nrow(Y)
  if (method == "qr") {
    QR <- qr(X)
    A <- solve(qr.R(QR))%*%t(qr.Q(QR))%*%Y
  } else if (method == "direct") {
    A  <- solve(t(X)%*%X)%*%(t(X)%*%Y)
  }
  
  SSE       <- t(Y - X%*%A)%*%(Y - X%*%A)
  sigma.sq     <- SSE/(nT - k)
  return(list(A        = A, 
              SSE      = SSE,
              sigma.sq = sigma.sq,
              k        = nrow(A),
              n        = ncol(A),
              nT       = nT))
}

# ---- 2. Prior specification ----
construct.Minnesota.prior <- function(Y, A.prior, ols, 
                                      use.lag.directly = FALSE,
                                      hyper.par = list(pi1 = 0.5,
                                                       pi2 = 0.5,
                                                       pi3 = 1,
                                                       pi4 = 200)) {
  lag.order <- (nrow(A.prior) - sum(c("const", "trend") %in% rownames(A.prior)))/
    ncol(A.prior)
  k <- ols$k
  # 1. Construct single AR process and  generate sigma.sq 
  sigma.sq <- sapply(1:ncol(Y), function(i) {
    Yi       <- as.matrix(Y[(lag.order + 1):nrow(Y), i])
    lag.Yi   <- as.matrix(na.omit(lag.gen(Y[, i, drop = FALSE], lag.order)))
    alpha    <- solve(t(lag.Yi)%*%lag.Yi)%*%(t(lag.Yi)%*%Yi)
    return(sum((Yi - lag.Yi%*%alpha)^2)/(nrow(Yi) - lag.order + 1))
  })
  names(sigma.sq) <- colnames(Y)
  
  # 2. Define standart deviations for parameters priors.
  V <- 0*A.prior
  own.first.lag <- !(lower.tri(A.prior)|upper.tri(A.prior))
  
  if ("const" %in% rownames(A.prior)) {
    V["const", ] <- hyper.par$pi1 * hyper.par$pi4 * sigma.sq # variance for coefficient for constant 
  }
  if ("trend" %in% rownames(A.prior)) {
    V["trend", ] <- hyper.par$pi1 * hyper.par$pi4 * sigma.sq # variance for coefficient for trend
  }
  
  # 3. Define variance for own lags of variables
  V[own.first.lag] <- hyper.par$pi1/((lag.order)^(2*hyper.par$pi3))
  
  # 4. Define variance for lags of other variables
  for (i in 1:1:ncol(V)) {
    for (j in 1:ncol(V)) {
      for (p in 1:lag.order) {
        if ((p == 1)&(i == j)) {
        } else {
          V[(p - 1)*ncol(Y) + i, j] <- hyper.par$pi2*hyper.par$pi1*
            sigma.sq[j]/(sigma.sq[i]*(lag.order)^(2*hyper.par$pi3))
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
  
  return(list(A     = A.prior,
              vec.A = vec.A.prior,
              V     = V.prior,
              sigma.sq    = ols$sigma.sq))
}

# ==== 2.2 Restriction matrix definition ====
# 
# In notations of K. Juselius and S. Johansen R-type restriction matrix is a
# matrix that establishes a system of equations to hold: R*params = 0;
# H-type matrix is a matrix that connects vector of all parameters with the
# vector of only unrestricted parameters: all.params = H*unrestr.params.
# If there are only exclusion constraints, it is convenient to use 
# exact restriction matrix 

# #### 2.2.1 From R-type to H-type matrix #####
R.to.H <- function(R){
  # ACHTUNG!!! This function doesn't work with combined restrictions
  # This function turns R matrix into its H counterpart
  # 1. Preliminary preparation
  names.restr <- colnames(R)
  
  r.restr <- length(R[, 1])
  k <- length(R[1, ])
  par.free <- k - r.restr
  
  H <- matrix(0, k, par.free)
  rownames(H) <- names.restr
  colnames(H) <- rep("smth", par.free)
  
  N.restr <- c()
  unrestr <- c()
  E.restr <- list()
  EE.restr <- list()
  
  n.n <- 1
  n.e <- 1
  n.ee <- 1
  n.u <- 1
  
  # 2. Classificate restrictions
  # E.restr <- euality restriction of type: b11 = b12
  # N.restr <- restriction of equality to zero: b11 = 0
  # unrestr <- no restriction at all
  for (i in 1:r.restr){
    if (sum(R[i, ]) == 1){
      N.restr[n.n] <- which(R[i, ] == 1)
      n.n <- n.n + 1
    } else if (sum(R[i, ]) == 2){
      E.restr[[n.e]] <- which(R[i, ] == 1)
      n.e <- n.e + 1
    } else if ((sum(R[i, ]) == 0)&(!all(R[i,] == 0))){
      temp <- rep(1, 2)
      temp[1] <- which(R[i, ] == 1)
      temp[2] <- - which(R[i, ] == -1)
      
      E.restr[[n.e]] <- temp
      n.e <- n.e + 1
    }
  }
  
  # 2.2 Check for unrestricted parameters
  for (i in 1:k){
    if (all(R[, i] == 0)){
      unrestr[n.u] <- i
      n.u <- n.u + 1
    }
  }
  
  # 3. Reproduce the same restrictions in R-type matrix
  # 3.1 Create N-type restrictions (zero-row at i-th element)
  if (length(N.restr) > 0){
    for (j in N.restr){
      H[j, ] <- rep(0, length(H[j, ]))
    }
  } 
  
  # 3.2 Create E-type restriction
  k.res <- 1
  if (length(E.restr) > 0){
    for (j in 1:length(E.restr)) {
      if (E.restr[[j]][1]*E.restr[[j]][2] > 0) { # c(1, 1) type
        H[E.restr[[j]], ] <- matrix(0, 2, par.free)
        H[E.restr[[j]], k.res] <- c(1, -1)
        
        if (length(names.restr) > 0){
          colnames(H)[k.res] <- paste(names.restr[E.restr[[j]][1]], 
                                      ".minus.", names.restr[E.restr[[j]][2]], sep = "")
        }
        
        k.res <- k.res + 1
      } else if (E.restr[[j]][1]*E.restr[[j]][2] < 0) { # c(1, -1) type
        H[abs(E.restr[[j]]), ] <- matrix(0, 2, par.free)
        H[abs(E.restr[[j]]), k.res] <- rep(1, 2)
        
        if (length(names.restr) > 0){
          colnames(H)[k.res] <- paste(names.restr[abs(E.restr[[j]])[1]], 
                                      ".plus.", names.restr[abs(E.restr[[j]])[2]], sep = "")
        }
        
        k.res <- k.res + 1
      }
    }
  }
  
  # 3.3 Create unrestricted restrictions
  if (length(unrestr) > 0) {
    for (j in unrestr){
      H[j, ] <- rep(0, par.free)
      if (length(k.res) > 0){
        H[j, k.res] <- 1
      }
      
      if (length(names.restr) > 0){
        colnames(H)[k.res] <- names.restr[j]
      }
    }
  } 
  
  return(H)
}

# #### 2.2.2 From R-type  to H-type #####
H.to.R <- function(H){
  # This function returns a t(R) matrix of restriction from H matrix
  par.free <- length(H[1, ])
  k <- length(H[, 1])
  N.restr <- k - par.free
  
  Minus.restr <- c()
  Plus.restr <- c()
  
  R <- matrix(0, N.restr, k)
  colnames(R) <- rownames(H)
  k.r <- 1
  
  # 1. Restrictions of equality
  for (i in seq_along(H[1, ])){ # handle by cols in one row
    if (all(H[, i] == 0)){
      warning("I found the Null-column in the H matrix and deleted it. Drop it from your H matrix immediately!")
    } else {
      Plus.restr <- which(H[, i] == 1)
      Minus.restr <- which(H[, i] == -1) 
      
      num.equal <- sum(abs(H[ ,i]))
      
      if (num.equal == 1){
        # This is the case of absolutely free parameter
        # Don't do anything
      } else {
        if (length(Plus.restr) > 1){
          for (j in 2:length(Plus.restr)){
            R[k.r, c(Plus.restr[j - 1], Plus.restr[j])] <- c(1, -1)
            k.r <- k.r + 1
          }
        }
        
        if (length(Minus.restr) > 1){
          for (j in 2:length(Minus.restr)){
            R[k.r, c(Minus.restr[j - 1], Minus.restr[j])] <- c(1, -1)
            k.r <- k.r + 1
          }
        }
        
        if ((length(Plus.restr) > 0)&(length(Minus.restr)) > 0){
          R[k.r, c(Plus.restr[1], Minus.restr[1])] <- c(1, 1)
          k.r <- k.r + 1
        }
      }
    }
  }
  
  # 2. Restrictions of equality to zero.
  for (i in seq_along(H[, 1])){ # handle by rows in one column
    if (all(H[i, ] == 0)){
      R[k.r, ] <- rep(0, k)
      R[k.r, i] <- 1
      k.r <- k.r + 1
    }
  }
  return(R)
}

# #### 2.2.3 From exact restrictions to H-type #####
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

# ==== 2.3 Match priors on reduced-form parameters with priors on structural ones =====
structural.priors.on.unrestricted.params <-function(prior, restrictions) {
  # It is assumed that priors between equations are independent.
  # This version of function does not support non-zero mean priors for the reduced-form coefficients.
  temp <- lapply(seq_along(restrictions$B), function(i) {
    which.B.i <- (i - 1)*nrow(restrictions$B[[i]]) + 1:nrow(restrictions$B[[i]])
    which.A.i <- (i - 1)*nrow(restrictions$A[[i]]) + 1:nrow(restrictions$A[[i]])
    
    S.i <- prior$vec.A$S[which.A.i, which.A.i]
    H.i <- prior$vec.B$H[which.B.i, which.B.i]
    P.i <- prior$vec.B$P[which.B.i, which.A.i]
    
    H   <- solve(t(restrictions$B[[i]]) %*% solve(H.i) %*% restrictions$B[[i]])
    P   <- H %*% t(restrictions$B[[i]]) %*% solve(H.i) %*%
      P.i %*% restrictions$A[[i]]
    S   <- solve(t(restrictions$A[[i]]) %*% solve(S.i) %*% restrictions$A[[i]] + 
                   t(restrictions$A[[i]]) %*% t(P.i) %*% solve(H.i) %*%  
                   P.i%*%restrictions$A[[i]] -
                   t(P) %*% solve(H) %*% P)
    list(H = H, P = P, S = S)
    })
  res <- lapply(1:3, function(i) {
    lapply(seq_along(temp), function(j) {
      temp[[j]][[i]]
      })
    })
  names(res) <- c('H', 'P', 'S')
  return(res)
}

# ==== 2.4 Construct a matrix with b draws and H-type restrictions =====
free.params.to.matrix <- function(params, restrictions) {
  res <- matrix(NA, nrow = nrow(restrictions[[1]]), ncol = length(restrictions))
  for (i in seq_along(restrictions)) {
    res[, i] <- restrictions[[i]] %*% matrix(params[[i]], ncol = 1)
  }
  return(res)
}

# ---- 3. Posterior functions ----
analytical.posterior <- function(X, Y, prior, ols) {
  if (prior$type == "Noninformative") {
    # Posterior of coeff|Data ~ Multi-T(kronecker(SSE, solve(t(X)%*%(X))), A.ols, nT - k)
    V.post     <- solve(t(X)%*%X)
    vec.A.post <- ols$vec.A
    A.post     <- ols$A
    
    # Posterior of Sigma|Data ~ inv-Wishart(SSE, T - k)
    S.post <- SSE
    v.post <- nT - k
    
    # Mean and variance of posterior 
    A.mean <- A.post
    A.var  <- kronecker(S.post, V.post)/(v.post - ols$n - 1)
  } else if (prior$type == "Minnesota") {
    # For coeff
    V.post <- solve(solve(prior$coeff$V) + kronecker(solve(prior$coeff$sigma.sq)))
    A.post <- V.post%*%(solve(prior$V)%*%prior$coeff$A + 
                          kronecker(solve(prior$coeff$sigma.sq), t(X)%*%X)%*%ols$A)
    
    A.mean <- A.post
    A.var  <- NULL
    S.post <- NULL
    v.post <- NULL
  } else if (prior$type == "Normal-Wishart") {
    # For coeff
    V.post <- solve(solve(prior$coeff$V) + t(X)%*%X)
    A.post <- V.post %*% (solve(prior$coeff$V)%*%prior$coeff$A + t(X)%*%X%*%ols$A)
    vec.A.post <- matrix.to.vector(A.post)
    
    # For Sigma
    S.post <- ols$SSE + prior$covmat$S + t(ols$A)%*%t(X)%*%X%*%ols$A + 
      t(prior$A)%*%solve(prior$coeff$V)%*%prior$coeff$A - t(A.post)%*%
      (solve(prior$coeff$V) + t(X)%*%X)%*%A.post
    v.post <- ols$nT + prior$covmat$v 
    A.var <- (1/(v.post - ols$n - 1))*kronecker(S.post, V.post)
  } else if (prior$type == "Minnesota-Wishart") {
    V.post <- solve(solve(prior$coeff$V) + kronecker(solve(prior$coeff$sigma.sq)))
    A.post <- V.post%*%(solve(prior$V)%*%prior$coeff$A + 
                          kronecker(solve(prior$coeff$sigma.sq), t(X)%*%X)%*%ols$A)
    
    A.mean <- A.post
    
    # For Sigma
    S.post <- ols$SSE + prior$covmat$S + t(ols$A)%*%t(X)%*%X%*%ols$A + 
      t(prior$A)%*%solve(prior$coeff$V)%*%prior$coeff$A - t(A.post)%*%
      (solve(prior$coeff$V) + t(X)%*%X)%*%A.post
    v.post <- ols$nT + prior$covmat$v 
    A.var <- (1/(v.post - ols$n - 1))*kronecker(S.post, V.post)
  }
  return(list(coeff  = list(mean = coeff.mean, var = coeff.var),
              covmat = list(S.post = S.post, v.post = v.post)))
}

# ---- 4. Reduced form ----
VAR.gibbs.sampler <- function(Y, X, burn.in, sample.save, 
                              prior, repfor, lag.order,
                              irf.settings) {
  require(MCMCpack)
  # 1. Set prliminary values and create empty objects
  # 1.1 Main values
  n  <- ncol(Y)
  k  <- ncol(Z)
  nT <- nrow(Y)
  h <- irf.settings$horizon
  Z <- kronecker(diag(rep(1, n)), X)
  
  sample.size <- sample.save + burn.in
  bigj <- cbind(diag(rep(1, n)), matrix(0, nrow = n, ncol = (lag.order - 1)*n))
  
  irf.list <- lapply(1:irf.settings$horizon, function(i) {
    res <- matrix(NA, nrow = n*n, ncol = sample.save)
    rownames(res) <- c(sapply(colnames(Y), function(c) {
      paste0("resp.", c, ":imp.", colnames(Y))
    }))
    return(res)
  })
  
  # 1.2 Extract elements from prior object
  vec.A.prior <- matrix.to.vector(prior$coeff$A)
  vec.Sigma   <- matrix.to.vector(prior$covmat$S)
  # 1.3 Empty objects for Gibbs Sampling
  vec.A.post.list <- matrix(NA, ncol = nrow(vec.A.prior), nrow = sample.save)
  rownames(vec.A.post.list) <- colnames(vec.A.prior)
  vec.Sigma.post.list  <- matrix(NA, ncol = nrow(vec.Sigma), nrow = sample.save)
  rownames(vec.Sigma.post.list) <- colnames(vec.Sigma)
  
  Sigma <- prior$coeff$sigma.sq
  for (i.rep in 1:sample.size) {
    if (i.rep %% 1e3 == 0)  {
      cat('Iteration', i.rep, '/', sample.size, '\n')
    }
    
    # Posterior of A.post | Sigma.post, Data
    Variance <- kronecker(solve(Sigma), diag(rep(1, nT)))
    #                   sparseMatrix(i = 1:nT, j = 1:nT, x = rep(1, nT)))
    V.post     <- solve(solve(prior$coeff$V) + t(Z)%*%Variance%*%Z)
    vec.A.post <- V.post%*%(solve(prior$coeff$V)%*%prior$coeff$vec.A +
                              t(Z)%*%Variance%*%matrix.to.vector(Y))
    A.post <- vector.to.matrix(vec.A.post)
    
    # Posterior of Sigma.post | A.post, Data ~ iW(inv(S), v.post)
    v.post <- nT + prior$covmat$v
    S.post <- prior$covmat$S + t(Y - X%*%A.post)%*%(Y - X%*%A.post)
    Sigma  <- solve(rwish(v = v.post, S = solve(S.post)))
    
    # Out of the burn-in sample 
    if (i.rep > burn.in) {
      # IRF simulation
      biga <- matrix(0, n*lag.order, n*lag.order)
      for (j in 1:(lag.order - 1)) {
        biga[(j*n + 1:(j*n)), (j*n + 1:(j*n))] <- diag(rep(1, n))
      }
      atemp <- A.post[2:nrow(A.post), ]
      atemp <- matrix.to.vector(atemp)
      splace <- 0
      for (l in 1:lag.order) {
        for (s in 1:n) {
          biga[s, ((l - 1)*n + 1):(l*n)] <- atemp[splace + 1:n, 1]
          splace <- splace + n
        }
      }
      # St dev matrix for SVAR
      STCO <- chol(Sigma)
      bigai <-biga
      
      imp.resp <- vector("list", irf.settings$horizon)
      for (k in seq_along(imp.resp)) {
        imp.resp[[k]] <- bigj%*%bigai%*%t(bigj)%*%STCO
        colnames(imp.resp[[k]]) <- paste0("resp.", colnames(Y))
        rownames(imp.resp[[k]]) <- paste0("imp.", colnames(Y))
        bigai <- bigai%*%biga
      }
      
      # Save all draws
      vec.A.post.list[i.rep - burn.in, ]     <- vec.A.post
      vec.Sigma.post.list[i.rep - burn.in, ] <- matrix.to.vector(Sigma)
      for (k in 1:irf.settings$horizon) {
        irf.list[[k]][, i.rep - burn.in] <- c(matrix.to.vector(imp.resp[[k]],
                                                         sep = c("", ":")))
      }
    }
  }
  
  return(list(vec.A.post      = vec.A.post.list,
              vec.Sigma.post  = vec.Sigma.list,
              irf.post        = irf.list))
}

# ---- 5. Structural form ----
# ==== 5.1 Priors and restrictions handling  ====

# ==== 5.2 Waggoner-Zha algorithm ====
SVAR.gibbs.sampler <- function(Y, Z, restrictions,
                              prior, irf.settings,
                              start.values = NULL,
                              burn.in = 5e2, sample.save = 1e4,
                              likelihood.normalization = FALSE) {
  require(mvtnorm)
  require(MASS)
  nT <- nrow(Y)
  # 1. Prepare neccessary functions
  # 1.1 Generate weights for free structural parameters of the contemporaneous matrix.
  generate.weights <- function(A.restr, S, chol.S, restr.mat) {
    # This function calculates the w-vectors for decomposition of structural parameters
    # Auxilliary function to generate 0 in vector if condition is not hold
   select.row <- function(x, cond){ 
     x[!cond, ] <- 0
     return(x)
   }
   
   # Generate the 1st vector of the R^(q_i) space
   qr.A.restr <- qr(A.restr)
   Q <- qr.Q(qr.A.restr)
   w <- Q[, ncol(Q), drop = FALSE] # this w solves the problem: w'%*%t.A.restr = 0, see ...
   w1 <- t(restr.mat%*%chol.S)%*%w / sqrt(c(t(w)%*%restr.mat%*%S%*%t(restr.mat)%*%w) )
   
   # Calculate other vectors separately
   w.all <- lapply(1:nrow(w1), function(i) {
     if (i == 1) {
       return(c(w1))
     } else {
       w.i   <- c(select.row(x = w1, cond = 1:nrow(w1) <= i)) # 0-elements j >  i
       w.im1 <- c(select.row(x = w1, cond = 1:nrow(w1) < i))  # 0-elements j >= j
       res   <- w.im1*w1[i, ]/sqrt(sum(w.i^2)*sum(w.im1^2))
       res[i] <- -sum(w.im1^2)
       return(res)
     }
   })
   return(w.all)
  }
  
  # 1.2 Function to generate dk sample of free parameters contained in A0 matrix.
  sample.bk <- function(A.restr, S, restr.mat, nT) {
    chol.S <- t(chol(S)) # a matrix L of the inv.S = L'L, i.e. in W-Zh notation this is T = L'
    w.all <- generate.weights(A.restr, S, chol.S, restr.mat)
    
    beta <- rep(NA, ncol(S)) # Since number of free parameters = ncol(inv.S)
    # Generate the first beta from gamma distribution
    beta[1] <- ifelse(runif(n = 1) <= 0.5, -1, 1)*
      rgamma(n = 1, shape = 0.5*(nT + 1), scale = 0.5*nT)
    # Generate other betas conditional on the first beta from Normal distribution
    if (length(beta) > 1) {
      beta[2:length(beta)] <- rnorm(n = length(beta) - 1, mean = 0, sd = 1/sqrt(nT))
    }
    
    sum.w <- matrix(unlist(w.all), ncol = ncol(S), byrow = FALSE)%*%
      matrix(beta, ncol = 1)
    bk <- chol.S%*%sum.w
    return(bk)
  }
  
  # 1.5 Calculate parameters of posterior sample 
  calc.posterior.params <- function(Y, Z, prior, restrictions) {
    nT <- nrow(Y)
    # 1.5.1 Preliminary calculations to enhance speed
    ZZ <- t(Z)%*%Z
    ZY <- t(Z)%*%Y
    YY <- t(Y)%*%Y
    
    # 1.5.2 Compute parameters for posterior functions
    post.params <- lapply(1:ncol(Y), function(i) {
      U <- restrictions$H.type$A[[i]]
      V <- restrictions$H.type$B[[i]]
      
      H <- solve(t(V)%*%ZZ%*%V + solve(prior$H[[i]]))
      P <- H%*%(t(V)%*%ZY%*%U + solve(prior$H[[i]])%*%prior$P[[i]])
      S <- solve(t(U)%*%YY%*%U + solve(prior$S[[i]]) + 
                       t(prior$P[[i]])%*%solve(prior$H[[i]])%*%prior$P[[i]] -
                       t(P)%*%solve(H)%*%P)/nT
      return(list(H = H, P = P, S = S))
    })
    return(post.params)
  }
  # 2. One-step function for gibbs sampler realization
  posterior.sampling <- function(b.prev, g.prev, post.params, prior, 
                                 restrictions, burn = FALSE) {
    nT <- nrow(Y)
    # 1.4.1 Sample b.k values
    A.all.restr <- free.params.to.matrix(b.prev, restrictions$H.type$A)
    b <- lapply(seq_along(restrictions$H.type$A), function(k) {
      sample.bk(A.restr   = A.all.restr[, -k, drop = FALSE], 
                S         = post.params[[k]]$S, 
                restr.mat = restrictions$H.type$A[[k]], 
                nT) 
      }) 
    if (!burn) {
      g <- lapply(seq_along(restrictions$H.type$A), function(k) {
        rmvnorm(n = 1, mean = post.params[[k]]$P%*%matrix(b.prev[[k]], ncol = 1), 
                sigma = post.params[[k]]$H)
        })
    } else g <- g.prev
    return(list(b = b, g = g))
  }
  
  # 3. Start doing something already
  # 3.1 Initialization
  sample.size <- sample.save + burn.in
  
  b.prev <- start.values$b
  g.prev <- start.values$g
  # 3.2 Calculate posterior parameters
  post.params <- calc.posterior.params(Y, Z, prior, restrictions)
  
  # 3.3. Start sample
  for (i.rep in 1:sample.size) {
    if (i.rep %% 1e3 == 0)  {
      cat('Iteration', i.rep, '/', sample.size, '\n')
    }
    
    sample <- posterior.sampling(b.prev, g.prev, post.params, prior, 
                       restrictions, burn = i.rep <= burn.in)
    b.prev <- sample$b
    g.prev <- sample$g

  }
  
}