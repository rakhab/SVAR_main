rm(list = ls())

library(MCMCpack)
library(mvtnorm)
library(MASS)
library(ggplot2)
library(gridExtra)
library(igraph)
# ---- 1. Simulation study ----
setwd('/home/shere-khan/Dropbox/Research/bvar_codes/myBVAR')
source('Common_functions.R')
source('SVAR_gibbs.R')
# Prior elicitation: 
# Due to Lindesley-Bartlett paradox we don't use diffuse priors
c <- 20
scale <- 1

# ==== 1.1 Generate artificial data sets
n.obs <- 300
#Z <- rmvnorm(n = n.obs, mean = rep(0, 2), sigma = diag(rep(1, 2)))
# colnames(Z) <- c("L1.oh", "L2.uh")
Z <- matrix(rnorm(n = n.obs, mean = 0, sd = 1), ncol = 1)
colnames(Z) <- "Z"

A <- scale * matrix(c(1, 0.5, 0.5, 1), ncol = 2, nrow = 2, byrow = TRUE) 
diag(A) <- 1
colnames(A) <- c("oh", "uh")
rownames(A) <- c("oh", "uh")

B <-  scale * matrix(c(0.5, 0.3), ncol = 2, nrow = 1, byrow = TRUE)
colnames(B) <- c("oh", "uh")
# rownames(B) <- c("L1.oh", "L2.uh")
rownames(B) <- "Z"
eps <- matrix(rnorm(n = 2*n.obs, mean = 0, sd = 1), ncol = 2)

# ==== 1.2 Check if there are unexpected zero values in concentration matrix ====
check.params <- lapply(0:15, function(ind.model){
  R <- decode.restrictions(ind.model, 2, 3)
  P <- rbind(A, -B)*R
  return(P %*% t(P))
  })

scale <- 1

B <-  scale * matrix(c(0.5, 0.3), ncol = 2, nrow = 1, byrow = TRUE)
res <- vector("list", 4)
k <- 1
for (a21 in c(0.05, 0.1, 0.2, 0.4)) {
  A <- scale * matrix(c(1, 0.5, a21, 1), ncol = 2, nrow = 2, byrow = TRUE) 
  diag(A) <- 1
  colnames(A) <- c("oh", "uh")
  rownames(A) <- c("oh", "uh")
  
  # ==== 1.2 Check if there are unexpected zero values in concentration matrix ====
  R <- decode.restrictions(5, 2, 3) # ind.model = 5 <=> M06 in paper's notation
  P <- rbind(A, -B)*R
  res[[k]] <- P %*% t(P) 
  k <- k + 1
}


# Check reduced form parameters
check.params2 <- lapply(0:15, function(ind.model){
  R <- decode.restrictions(ind.model, 2, 3)
  Phi <- (B*R[3, 1:2, drop = FALSE]) %*% solve(A * R[1:2, 1:2])
  Sigma <- (A * R[1:2, 1:2]) %*% t(A * R[1:2, 1:2])
  return(list(Phi = Phi, Sigma = Sigma))
  })


# ==== 1.3 General prior setting on unrestricted SVAR parameters ====
prior <- list(vec.A = list(S = NULL), vec.B = list(P = NULL, H = NULL))

prior$vec.B$P <- 0*matrix(0, nrow = ncol(Z)*ncol(eps), ncol = ncol(eps)*ncol(eps))
prior$vec.B$H <- diag(rep(c, nrow(prior$vec.B$P)))
prior$vec.A$S <- diag(rep(c, ncol(eps)^2))

# ==== 1.4 Start loop models ====
set.seed(1234) # for replication of results
# #### 1.4.1 Small parameter values ####
setwd('/home/shere-khan/Dropbox/Research/final/test2')
timer <-vector("list", 16)
for (ind.true.model in 0:15) {
  start.timer <- proc.time()
  out <- direct.draw.model(ind.model = ind.true.model, 
                    A = A, B = B, prior = prior, 
                    eps = eps, Z = Z, calc.direct.marg = FALSE,
                    notification = paste0("True model: ", sprintf("%02d", 
                    ind.true.model + 1), "; "))
  # Chib method
  summary.and.plot.prior.models(output = out, filename =
                                paste0("c", ifelse(c >= 1, sprintf("%02d", c),
                                                   paste0("0dot", c*10)),
                                "_M", sprintf("%02d", ind.true.model + 1)),
                                wd = getwd(),
                                graph.format = ".pdf",
                                plot.densities = TRUE,
                                plot.diagnostics = FALSE,
                                marg.type = "Chib",
                                A = rbind(A, B),
                                ind.true.model = ind.true.model + 1) 
#   # Direct method
  summary.and.plot.prior.models(output = out, filename =
                                  paste0("c", ifelse(c >= 1, sprintf("%02d", c),
                                                     paste0("0dot", c*10)),
                                "_direct_M", sprintf("%02d", ind.true.model + 1)),
                                wd = getwd(),
                                graph.format = ".pdf",
                                plot.densities = FALSE,
                                plot.diagnostics = FALSE,
                                marg.type = "direct",
                                A = rbind(A, B),
                                ind.true.model = ind.true.model + 1) 
  # To lower CPU consumption at each iteration
  rm(list = "out")
  timer[[ind.true.model + 1]] <- proc.time() - start.timer
}

R <- prob.results.to.latex.tab(c.vec = 20, 
                               num.models = 16, 
                               wd = "/home/shere-khan/Dropbox/Research/texts") 

# #### 1.4.3 Try MC3 procedure ####
for (ind.true.model in 0:15) {
  mc3.draws <- mc3.draw.model(ind.model = ind.true.model, A, B, prior, 
                             eps, Z, notification = "MC3") 
  save(file = paste0("MC3_M", sprintf("%2d", ind.true.model), ".RData"),
       list = "mc3.draws")
  rm(list = "mc3.draws")
}

timer <-vector("list", 16)
for (ind.true.model in 0:15) {
  start.timer <- proc.time()
  print(ind.true.model)
  load(paste0("variables_c20_M", sprintf("%02d", ind.true.model + 1), ".RData"))
  summary.and.plot.prior.models(output = output, filename =
                                  paste0("c20_M", sprintf("%02d", ind.true.model + 1)),
                                wd = getwd(),
                                graph.format = ".pdf",
                                plot.densities = FALSE,
                                plot.diagnostics = FALSE,
                                marg.type = "Chib",
                                A = rbind(A, B),
                                ind.true.model = ind.true.model + 1) 
  timer[[ind.true.model + 1]] <- proc.time() - start.timer
}

K <- matrix(NA, ncol = 16, nrow = 8)
rownames(K) <- c("Log Chib estimate", "Log likelihood", "Max.lik", "Log Posterior", "Log prior",
                 "Conditional b1", "Conditional b2", "Conditional g")
for (ind.true.model in 0:15) {
  start.timer <- proc.time()
  K[, ind.true.model + 1] <- c(output[[ind.true.model + 1]]$BSVAR@marginal.likelihood$chib.exact,
                               output[[ind.true.model + 1]]$BSVAR@marginal.likelihood$lf.my,
                               output[[ind.true.model + 1]]$mle$optimization$objective,
                               output[[ind.true.model + 1]]$BSVAR@marginal.likelihood$dposterior,
                               output[[ind.true.model + 1]]$BSVAR@marginal.likelihood$dprior,
                               output[[ind.true.model + 1]]$BSVAR@marginal.likelihood$dconditional.b1,
                               output[[ind.true.model + 1]]$BSVAR@marginal.likelihood$dconditional.b2,
                               output[[ind.true.model + 1]]$BSVAR@marginal.likelihood$dconditional.g)
}

