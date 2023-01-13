# 1)
x.ts <- read.table('astronomy A.txt', header=FALSE)
x.ts <- x.ts$V1
n = length(x.ts)
m = 2
n.lag <- (2*m)+1
wgts <- rep(x=1, times=n.lag)/n.lag
x.hat.ts <- ma.smooth(x.ts, wgts)
x.hat.ts[32]          # 54.86

# 2) 
wgts <- Gauss.wgts(m)
x.hat.ts <- ma.smooth(x.ts, wgts)
x.hat.ts[12]         # 50.82538

# 3) 
x.vect <- as.matrix(x.ts)
Z.mat <- cbspline.regmat(n, nspline=8, ratio=1.3)
result <- reg.summary(x.vect, Z.mat)
result$mu[13]      # 44.35398

# 4) 
time.seq <- 1:n
omega <- 5/46      # 46 is n
r <- dim(Z.mat)[2]
Z.mat <- cbind(Z.mat, matrix(data=0, nrow=n, ncol=2))
Z.mat[,r+1] <- cos(2*pi*omega*time.seq)
Z.mat[,r+2] <- sin(2*pi*omega*time.seq)
r <- r + 2
result <- reg.summary(x.vect, Z.mat)
result$mu[13]           # 53.01345

# 5)
wgts <- Gauss.wgts(m)
x.hat.ts <- reg.smooth(x.ts, wgts)
x.hat.ts[23]            # 55.65723

# 6) 
alpha = c(1.15, 0.55)
Er2.form <- alpha[1] / (1 - alpha[2])        
Er2.form         # 2.555556

# 7) 
Er4.form <- 3*alpha[1]^2*(1+alpha[2]) / ((1 - alpha[2])*(1-3*alpha[2]^2))
kappa.form <- Er4.form / (Er2.form)^2             
kappa.form           # 22.62162

# 8) 
n = 10^5
zcut = 2.17
arch.data <- arch.sim(n=n, alpha=alpha)
cutoff <- zcut*sqrt(Er2.form)
table(abs(arch.data$ts) > cutoff) / n           # 0.03401  
