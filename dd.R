A <- matrix(rnorm(9, 0, 1), 3)
A.r <- A
A.r[1, ]    <- 0
A.r[2, 2:3] <- 0
A.r[3, 3]   <- 0

A1 <- cbind(rep(0, 3), A.r[, c(2, 3)])
A1 <-  A.r[, c(2, 3)]
QR1 <- qr(A1)

Q1 <- qr.Q(QR1)
R1 <- qr.R(QR1)

Q1%*%R1 - A1

w <- Q1[, ncol(Q1), drop = FALSE]

t(Q1)%*%Q1
t(w)%*%A1

azaza <- sapply(1:1e3, function(i) {
  A <- matrix(rf(9, 200, 50), 3)
  A.r <- A
  A.r[1, ]    <- 0
  A.r[2, 2:3] <- 0
  A.r[3, 3]   <- 0
  
  A1 <-  A.r[, c(2, 3)]
  QR1 <- qr(A1)
  
  Q1 <- qr.Q(QR1)
  R1 <- qr.R(QR1)
  
  return(t(w)%*%A1)
})

sum(azaza)

