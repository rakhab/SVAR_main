a <- 0
sigma <- 1
sigma.X <- 1

alpha = 0.5

num.observations <- 300

X <- rnorm(num.observations, mean = 0, sd = sigma.X)
Y <- a * X + rnorm(num.observations, mean = 0, sd = sigma)
plot(X, Y, pch = 20, col = 'blue')
abline(lm(Y ~ X), col = 'red')

X <- X - mean(X)
Y <- Y - mean(Y)

log.likelihood <- function(a, X, Y) {
  residuals <- Y - a * X
  -0.5 * length(residuals) * log(2 * pi * sigma^2) - sum(residuals^2) / (2 * sigma^2)
}

construct.integrand <- function(X, Y) {
  function(a) {
    sapply(a, function (a) exp(log.likelihood(a, X, Y) - log.likelihood(0, X, Y)))
  }
}

integrand <- construct.integrand(X, Y)

pdf.data.M2 <- integrate(integrand, lower = -10, upper = 10)$value

posterior <- alpha / (alpha + (1 - alpha) * pdf.data.M2)
posterior


# ---- 1. Normal model case ----
calculate.posterior.prob.normal <- function(a.seq, sigma, sigma.X, 
                                     num.observations, num.iter,
                                     alpha) {
  log.likelihood <- function(a, X, Y) {
    residuals <- Y - a * X
    -0.5 * length(residuals) * log(2 * pi * sigma^2) - sum(residuals^2) / (2 * sigma^2)
  }
  
  construct.integrand <- function(X, Y) {
    function(a) {
      sapply(a, function (a) exp(log.likelihood(a, X, Y) - log.likelihood(0, X, Y)))
    }
  }
  
  construct.posterior.prob <- function(sigma, sigma.X, num.observations,
                                       alpha) { 
    function(a){
      X <- rnorm(num.observations, mean = 0, sd = sigma.X)
      Y <- a * X + rnorm(num.observations, mean = 0, sd = sigma)
      
      X <- X - mean(X)
      Y <- Y - mean(Y)
      
      integrand <- construct.integrand(X, Y)
      pdf.data.M2 <- integrate(integrand, lower = -10, upper = 10)$value
      
      posterior <- alpha / (alpha + (1 - alpha) * pdf.data.M2)
      return(posterior)
    }
  }
  
  construct.mean.posterior.prob <- function(sigma, sigma.X, num.iter,
                                            num.observations, alpha) {
    posterior.prob <- construct.posterior.prob(sigma, sigma.X, 
                                               num.observations, alpha)
    function(a) {
      sapply(1:num.iter, function(i) posterior.prob(a))
    }
  }
  seq.posterior.prob <- construct.mean.posterior.prob(sigma, sigma.X, num.iter,
                                                       num.observations, alpha)
  
  res <- lapply(a.seq, function(a) {
    cat(paste0("a = ", a, "\n"))
    seq.posterior.prob(a)
    })
  
  names(res) <- paste0("param", a.seq)
  class(res) <- "iterated.posterior.prob"
  return(res)
}

PP <- calculate.posterior.prob.normal(a.seq = seq(from = -1, to = 1, length = 301), 
                         sigma = 1, 
                         sigma.X = 1, 
                         num.observations = 1e3, 
                         num.iter = 1e3,
                         alpha = 0.5)

summary.iterated.posterior.prob <- function(x, pc = c(0.01, 0.05, 0.95, 0.99)) {
  res <- vector("list", 4 + length(pc))
  names(res) <- c("mean", "median", "max", "min", paste0("p", pc))
  f <- c(mean, median, max, min)
  
  for (k in seq_along(f)) {
    res[[k]] <- sapply(seq_along(x), function(i) f[[k]](x[[i]]))
    names(res[[k]]) <- names(x)
  }
  
  for (p in seq_along(pc)) {
    res[[4 + p]] <- sapply(seq_along(x), function(i) quantile(x[[i]], pc[p]))
    names(res[[4 + p]]) <- names(x)
  }
  
  res[["param"]] <- as.numeric(gsub(pattern = "param", replacement = "", x = names(x)))
  return(as.data.frame(res))
}
summary(PP)


my.plot.iterated.posterior.prob  <- function(x, ...) {
  require(ggplot2)
  x.sum <- summary(x)
  ggplot(x.sum, aes(x = param)) +
    coord_cartesian(ylim=c(0,1)) +
    geom_ribbon(aes(ymin = p0.01, ymax = p0.99),
                fill = "blue",alpha = 0.1) + 
    ylab("Posterior Probabilities") +
    geom_line(aes(y = p0.05), col = "blue", linetype = "dashed") +
    geom_line(aes(y = p0.95), col = "blue", linetype = "dashed") +
    geom_line(aes(y = median), col = "red") +
    geom_line(aes(y = mean), col = "red", linetype = "dotted") +
    theme_bw()
}

my.plot.iterated.posterior.prob(PP)

# ---- 2. Partial Correlations Case ----
rho <- seq(from = -1, to = 1, length = 1e3)
nu <- 300
plot( (1 - rho^2)^(nu/2 - 1), type = "l")

calculate.posterior.prob.wish <- function(a.seq, sigma, sigma.X, 
                                            num.observations, num.iter,
                                            alpha) {
  log.likelihood <- function(a, X, Y) {
    residuals <- Y - a * X
    -0.5 * length(residuals) * log(2 * pi * sigma^2) - sum(residuals^2) / (2 * sigma^2)
  }
  
  omega <- seq(from = -2, to = 2, length = 100)
  int <- omega^(-1)*(omega + 1/omega - 2*r*rho)^(k - n)
  plot(x = omega, y = int, type = "l")
  
  
  construct.integrand <- function(X, Y) {
    vec <- rbind(Y, X)
    C   <- vec%*%t(vec)
    rho <- C[1, 2]/sqrt(C[1, 1]*C[2, 2])
    nu <- length(X)
    
    function(r) {
      integrand0 <- function(omega) {
        sapply(omega, function(omega) {
         omega^(-1)*(omega + 1/omega - 2*r*rho)^(1 - nu)
        })
      }
      
      sapply(r, function (r) {
        calc.integrand0 <- integrate(integrand0, lower = 0, upper = 10)$value
      }
    }
  }
  
  construct.posterior.prob <- function(sigma, sigma.X, num.observations,
                                       alpha) { 
    function(a){
      X <- rnorm(num.observations, mean = 0, sd = sigma.X)
      Y <- a * X + rnorm(num.observations, mean = 0, sd = sigma)
      
      X <- X - mean(X)
      Y <- Y - mean(Y)
      
      integrand <- construct.integrand(X, Y)
      pdf.data.M2 <- integrate(integrand, lower = -10, upper = 10)$value
      
      posterior <- alpha / (alpha + (1 - alpha) * pdf.data.M2)
      return(posterior)
    }
  }
  
  construct.mean.posterior.prob <- function(sigma, sigma.X, num.iter,
                                            num.observations, alpha) {
    posterior.prob <- construct.posterior.prob(sigma, sigma.X, 
                                               num.observations, alpha)
    function(a) {
      sapply(1:num.iter, function(i) posterior.prob(a))
    }
  }
  seq.posterior.prob <- construct.mean.posterior.prob(sigma, sigma.X, num.iter,
                                                      num.observations, alpha)
  
  res <- lapply(a.seq, function(a) {
    cat(paste0("a = ", a, "\n"))
    seq.posterior.prob(a)
  })
  
  names(res) <- paste0("a", a.seq)
  class(res) <- "iterated.posterior.prob"
  return(res)
}