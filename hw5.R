library(expm)
# 1) 
econX <- read.table('econ B.txt')
x.vect <- as.matrix(econX)
n <- length(x.vect)
std.time <- seq(from=0, to=1, length.out=n)
Z.mat <- as.matrix(cbind(rep(x=1, times=n), std.time))
colnames(Z.mat) <- c("z0", "z1")

reg.summary <- function(x.vect, Z.mat) {
  n <- dim(Z.mat)[1]
  r <- dim(Z.mat)[2]
  ZpZ <- t(Z.mat) %*% Z.mat
  ZpX <- t(Z.mat) %*% x.vect
  ZpZ.inv <- solve(ZpZ)
  b.hat <- ZpZ.inv %*% ZpX
  mu <- Z.mat %*% b.hat
  res <- as.numeric(x.vect - mu)
  SSE <- sum(res^2)
  AIC <- n*log(2*pi*SSE/n) + n + 2*(r+1)
  result <- list(b.hat=b.hat, mu=mu, res=res, SSE=SSE, AIC=AIC)
  return(result)
}

ls.stat <- reg.summary(x.vect, Z.mat)
AIC <- ls.stat$AIC
print(AIC)           # 405.9867

# 2) 
model.parms <- list()
model.parms$x.vect <- x.vect
model.parms$Z.mat <- Z.mat
comp.parms <- list()
comp.parms$max.lag <- 1
comp.parms$n.samp <- 10000
result <- run.dsamp.reg.wn(model.parms, comp.parms)
result$stat$DIC            # 402.1446

# 3) 
max.lag <-25
ls.stat <- reg.summary(x.vect, Z.mat)
acf.seq <- get.acf.ar1(res=ls.stat$res, max.lag)
ls.stat <- reg.summary.acorr(x.vect, Z.mat, acf.seq, p.acf=2)
AIC <- ls.stat$AIC
AIC                              # 338.716

# 4) 
model.parms <- list()
model.parms$x <- x.vect
model.parms$x0 <- 10
comp.parms <- list()
comp.parms$n.samp <- 10000
comp.parms$m <- 1
comp.parms$max.kappa <- 1000
comp.parms$x0.rad <- 37.8
comp.parms$max.lag <- 1
result <- run.mcmc.diff.wn(model.parms, comp.parms)
result$stat$WAIC                  # 325.678

# 5) 
sez <- read.table("seasonal B.txt")
x.vect <- as.matrix(sez)
omega <- 5/45
n <- length(x.vect)
time.seq <- seq(from=0, to=1, length.out=n)
Z.mat <- matrix(data=0, nrow=n, ncol=5)
Z.mat[,1] <- rep(x=1, times=n)
Z.mat[,2] <- time.seq
Z.mat[,3] <- time.seq^2
Z.mat[,4] <- cos(2*pi*omega*time.seq)
Z.mat[,5] <- sin(2*pi*omega*time.seq)
ls.stat <- reg.summary(x.vect, Z.mat[,1:4])
AIC <- ls.stat$AIC
AIC
# 350.3021 for 1A
# 307.3012 for 2A 
# 306.6867 for 3A ----- lowest
# 308.278 for 4A

# 6) 
max.lag <- 25
ls.stat <- reg.summary(x.vect, Z.mat[,1:4])
acf.seq <- get.acf.ar1(res=ls.stat$res, max.lag)
ls.stat <- reg.summary.acorr(x.vect, Z.mat[,1:4], acf.seq, p.acf=2)
AIC <- ls.stat$AIC
AIC
# 326.535 for 1B 
# 308.0537 for 2B ------- lowest
# 308.1618 for 3B 
# 309.8597 for 4B

# 7)
delta <- 0.4
lambda <- 0.6
sigw <- 6.5
x0 <- 19
#ewma.sim <- function(n, sigw, lambda, delta=0, x0=0) {
  w0 <- rnorm(n=1, mean=0, sd=sigw)
  wn.samp <- rnorm(n=n, mean=0, sd=sigw)
  ewma.samp <- numeric(length=n)
  ewma.samp[1] <- delta + x0 + wn.samp[1] - lambda*w0
  for (t in 2:n) {
    ewma.samp[t] <- delta + ewma.samp[t-1] + wn.samp[t] - lambda*wn.samp[t-1]
  }
  ewma.ts <- ts(ewma.samp)
  return(ewma.ts)
}
ewma.ts <- c(23.8, 37.9, 23.0, 27.3, 31.8, 37.5, 33.0)
n <- length(ewma.ts)
#ewma.ts <- ewma.sim(n=n, sigw=sigw, lambda=lambda, delta=delta, x0=x0)
predx.measr <- numeric(length=n)
predx.measr[1] <- delta + x0
for (t in 2:n) {
  predx.measr[t] <- (1-lambda)*ewma.ts[t-1] + lambda*predx.measr[t-1] + delta
}
onestep <- (1-lambda)*ewma.ts[n] + lambda*predx.measr[n] + delta

max.m <- 2
predx <- numeric(length=max.m)
predx[1] <- onestep
for (m in 2:max.m) {
  predx[m] = (m-1)*delta + predx[1]
}
predx                   # 33.86833

# 8) 
pred.err <- numeric(length=max.m)
for (m in 1:max.m) {
  pred.err[m] <- (1 + (m-1)*(1-lambda)^2)*sigw^2
}
pred.err                # 49.01

# 9) 
x.vect <- c(35.3, 39.0, 36.3, 36.6, 30.4, 24.1, 26.4, 37.4, 42.3, 36.1, 30.8)
p <- 1
d <- 0
q <- 1
arma.tsa <- arima(x.vect, order=c(p,d,q))
arma.tsa$coef   # 0.1894630    

# 10)
sqrt(diag(arma.tsa$var.coef))   # 0.3252693  

# 11) 
sqrt(arma.tsa$sigma2)          # 3.59674

# 12)
arma.tsa$aic                  # 68.75239
