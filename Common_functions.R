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

# Namer for matrix elements 
gibbs.draw.namer <- function(mat.restrictions, name = "X") {
  unlist(lapply(1:ncol(mat.restrictions), function(i) {
    if (sum(!mat.restrictions[, i]) != 0) {
      return(paste0(name, 
            (1:nrow(mat.restrictions))[!mat.restrictions[, i]], i))
    } else {
      return(NULL)
    }
    }))
}


# ---- 2. Frequentist's estimation methods ----
# ==== 2.1 Simplified OLS function ====
ols.VAR <- function(Y, X, method = "direct") {
  k  <- nrow(A)
  nT <- nrow(Y)
  if (method == "qr") {
    QR <- qr(X)
    A <- solve(qr.R(QR))%*%t(qr.Q(QR))%*%Y
  } else if (method == "direct") {
    A  <- solve(t(X)%*%X)%*%(t(X)%*%Y)
  }
  
  SSE          <- t(Y - X%*%A)%*%(Y - X%*%A)
  sigma.sq     <- SSE/(nT - k)
  return(list(A        = A, 
              covmat   = t(Y)%*%Y/nT,
              SSE      = SSE,
              sigma.sq = sigma.sq,
              k        = k,
              n        = ncol(A),
              nT       = nT))
}

# ==== 2.2 Simplified 2SLS function ====
# Since there is a high sensitivity of MLE results to the start parameter points
# in case when the true parameter much greater than 1, 
# we propose to use 2SLS in order to find start values.

# #### 2.2.1 Function to generate random parameters ####
generate.random.parameters <- function(restrictions) {
  num.endog <- ncol(restrictions)
  num.exog <- nrow(restrictions) - num.endog
  S <- matrix(rnorm(num.exog^2), ncol = num.exog)
  A <- restrictions * 0
  for (i in 1:(num.endog + num.exog)) {
    for (j in 1:num.endog) {
      if (i == j) {
        A[i, j] <- runif(1)
      } else if (restrictions[i, j] > 0) {
        A[i, j] <- rnorm(1)
      }
    }
  }
  list('A' = A, 'S' = as.matrix(S), 
       'restrictions' = restrictions,
       'num.exogenous' = num.exog,
       'num.endogenous' = num.endog)
}

# #### 2.2.2 Function to generate random data from random parameters ####
generate.random.data <- function(model, num.observations) {
  Z <-  matrix(rnorm(num.observations * model$num.exogenous), ncol = model$num.exogenous) %*% solve(model$S)
  er <- matrix(rnorm(num.observations * model$num.endogenous), ncol = model$num.endogenous)
  A0 <- model$A[1:model$num.endogenous, , drop = FALSE]
  A1 <- model$A[(model$num.endogenous + 1):nrow(model$A), , drop = FALSE]
  Y <- (Z %*% A1 + er) %*% solve(A0)
  list('Y' = Y, 'Z' = Z, 'X' = cbind(Y, Z), 'residuals' = er)
}

# #### 2.2.3 2SLS  function ####
tsls <- function (data, restrictions) {
  parents <- function(data, node.index, restrictions) {
    as.matrix(data$X[ ,setdiff(which(restrictions[ , node.index] > 0), node.index), drop = FALSE])
  }
  instruments <- function(data, residuals = matrix(NA, ncol = ncol(data$Y), nrow = nrow(data$Y))) {
    ret <- cbind(data$Z, residuals)
    ret[ , colMeans(is.na(ret)) == 0, drop = FALSE]    
  }
  
  estimated.residuals <- matrix(NA, ncol = ncol(data$Y), nrow = nrow(data$Y))
  
  estimate.equation <- function(data, equation.index, restrictions, estimated.residuals, tolerance = 1e-5) {
    par <- parents(data, node.index = equation.index, restrictions = restrictions)
    instr <- instruments(data, residuals = estimated.residuals)
    par <- scale(par, scale = FALSE)
    instr <- scale(instr, scale = FALSE)
    node <- data$Y[ ,equation.index]
    Pi <- as.matrix(solve(as.matrix(t(instr) %*% instr)) %*% t(instr) %*% par)
    pi <- as.matrix(solve(as.matrix(t(instr) %*% instr)) %*% t(instr) %*% node)
    test.matrix <- as.matrix(t(Pi) %*% Pi)
    if (ncol(par) == 0) {
      is.identified <- TRUE
    } else {
      is.identified <- (sum(abs(eigen(test.matrix)$value) > tolerance) >= ncol(par))
    }
    if (!is.identified)
      return(list('is.identified' = FALSE))
    if (ncol(par) > 0) {
      coef <- - solve(t(Pi) %*% Pi) %*% t(Pi) %*% pi
    } else {
      coef <- integer(0)
    }
    
    est.residuals <- node - par %*% coef
    estimated.equation <- restrictions[ ,equation.index]
    counter <- 0
    for (i in 1:length(estimated.equation)) {
      if (i == equation.index) {
        estimated.equation[i] <- 1
        next
      } 
      if (estimated.equation[i] == 0) next
      counter <- counter + 1
      estimated.equation[i] <- coef[counter]
    }
    list('estimated.equation' = estimated.equation  / sd(data$X%*%estimated.equation),
         'is.identified' = TRUE)
  }
  
  A <- restrictions * NA
  estimated.residuals <- matrix(NA, ncol = ncol(data$Y), nrow = nrow(data$Y))
  is.identified <- rep(FALSE, ncol(restrictions))
  new.node.has.been.identified <- TRUE
  while (new.node.has.been.identified) {
    new.node.has.been.identified <- FALSE
    for (i in 1:length(is.identified)){
      if (is.identified[i]) next
      equation <- estimate.equation(data, i, restrictions, estimated.residuals)
      if (equation$is.identified) {
        new.node.has.been.identified <- TRUE
        A[ ,i] <- equation$estimated.equation
        estimated.residuals[ ,i] <- data$X %*% A[ ,i]
        is.identified[i] <- TRUE
      }
    }
    if (all(is.identified)) break
  }
  list('is.identified' = is.identified, 'A' = A)
}

# Nickolay Gennadyevich's MLE function
mle <- function(Y, Z, 
                constraints.on.coefficients,
                is.exact.constraints = TRUE,
                irf.length = 20,
                tolerance = 1e-7,
                maxiter = 1000,
                theta.start = NULL,
                start.type = "chol"){
  
  require(nloptr)
  constraints.on.coefficients <- rbind(!constraints.on.coefficients$A,                                     
                                         !constraints.on.coefficients$B)
  n <- ncol(constraints.on.coefficients)
  if ((ncol(Y)) != n || 
        (ncol(Y) + ncol(Z) != nrow(constraints.on.coefficients)))
    stop('Incompatible parameters')
  constraints <- list('A' = constraints.on.coefficients,
                      'num.A' = sum(constraints.on.coefficients))
  
  matrices.to.vector <- function(A, constraints) {
    vec <- rep(0, constraints$num.A)
    counter <- 0
    for (i in 1:nrow(constraints$A)) {
      for (j in 1:ncol(constraints$A)) {
        if (constraints$A[i, j] > 0) {
          counter <- counter + 1
          vec[counter] <- A[i, j]
        }
      }
    }
    vec
  }
  
  vector.to.matrices <- function(vec, constraints) {
    A <- constraints$A * 0
    counter <- 0
    for (i in 1:nrow(constraints$A)) {
      for (j in 1:ncol(constraints$A)) {
        if (constraints$A[i, j] > 0) {
          counter <- counter + 1
          A[i, j] <- vec[counter]
        }
      }
    }
    return(A)
  }
  
  construct.nll <- function(X, constraints) {
    X <- as.matrix(scale(X, scale = FALSE))
    n <- ncol(constraints$A)
    m <- nrow(X)
    N <- nrow(constraints$A)
    const <- 0.5 * n * log(2 * pi) 
    nll <- function(Params.vec) {
      # negative log likelihood
      if (length(Params.vec) != constraints$num.A) 
        stop('The supplied vector of parameters is not compatible with the imposed constraints')
      A <- vector.to.matrices(Params.vec, constraints)
      
      er.term <- X %*% A  # Orghogonal structural shocks
      
      A0 <- A[1:n, 1:n]
      objective <-  (const - log (abs(det(A0)))) * m + 0.5 * sum(er.term ^ 2)
      
      M <- matrix(0, ncol = n, nrow = N)
      for (i in 1:N) {
        for (j in 1:n) {
          M[i, j] <- sum(X[, i] * er.term[ ,j])
        }
      }
      grad.A <-  M - (m * rbind(t(solve(A0)), matrix(0, ncol = n, nrow = nrow(constraints$A)-n))) 
      nonorthogonal.er.term <- X %*% A      
      
      grad <- matrices.to.vector(grad.A, constraints)
      
      
      return(list('objective' = objective, 
                  'gradient' = grad))
    }
    return(nll)
  }
  
  
  
  X <- cbind(Y, Z)
  X <- X[complete.cases(X), ]
  #M <- t(chol(empirical.concentration(X)))[ ,1:n] # starting point for minimization of the nll
  nll <- construct.nll(X, constraints)
  if (is.null(theta.start)){
    if (start.type == "chol") {
      A.start <- t(chol(solve(t(X)%*%X)))[ ,1:n]
    } else if (start.type == "2SLS") {
      model <- generate.random.parameters(constraints.on.coefficients)
      rdata <- generate.random.data(model, nrow(Y))
      tsls.result <- tsls(rdata, restrictions = constraints.on.coefficients)
      if (all(tsls.result$is.identified)) {
        A.start <- tsls.result$A
      } else {
        warning("The model is not identified. Cholesky transformation is used instead of 2SLS estimation.")
        A.start <- t(chol(solve(t(X)%*%X)))[ ,1:n]
      }
    }
    theta.start <- matrices.to.vector(A.start, constraints)
  } else {
    theta.start <-  matrices.to.vector(theta.start$A, constraints)    
  }
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel"  = tolerance )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel"  = tolerance,
                "maxeval"   = maxiter,
                "local_opts" = local_opts )
  out <- nloptr( x0 = theta.start,
                 eval_f = nll,
                 opts = opts)
  
  
  solution <- vector.to.matrices(out$solution, constraints)
  A <- solution
  residuals <- as.matrix(scale(X, scale = FALSE)) %*% A 
  
#   irf <- function(A, S, imp.index, var.list, irf.length) {
#     A <- A %*% S
#     len <- irf.length
#     n <- ncol(A)
#     A0.inv <- solve(A[1:n, 1:n])
#     num.lags <- nrow(A) %/% n - 1
#     resp <- matrix(0, ncol = n, nrow = len)
#     for (i in 1:len) {
#       E <- rep(0, n) 
#       if (i == 1) {
#         E[imp.index] <- 1 * A[imp.index, imp.index]
#       }
#       for (j in 1:num.lags) {
#         if (i - j > 0)
#           E <- E - resp[i - j, ] %*% A[(n * j + 1):(n * (j + 1)), ]
#       }
#       resp[i, ] <- E %*% A0.inv
#     }
#     resp <- as.data.frame(resp)
#     #     for (i in 1:n) {
#     #       names(resp)[i] <- paste('imp: ', var.list[imp.index], ' resp: ', var.list[i], sep = '')
#     #     }
#     
#     resp
#   }
#   resp <- data.frame("Period" = 1:irf.length)
#   for (i in 1:n) {
#     resp <- cbind(resp, irf(A, S, imp.index = i, var.list, irf.length))
#   }
#   coefs <- A
#   for (i in 1:n) {
#     coefs[ ,i] <- - coefs[ ,i] / A[i, i]
#   }
  
  list('Y' = Y, 'Z' = Z,"is.identified" = NA, "A" = A, 
       "constraints.on.coefficients" = constraints.on.coefficients,
       "optimization" = out) 
}

# ---- 3. R Cookbook function for multiple plots ----
# CAUTION! Very slow!!!
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
