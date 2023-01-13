#1) 
# xt=3.0cos(2πωt − 0.5)
# A= 3.0, phi = -0.5
# U1 = Acos(phi)
u1 = 3 * cos(-0.5)
u1    # 2.632748

#2)
phi <- c(-1, -0.3)
theta <- c(-0.4)
sigw <- 5
omega <- 0.26

imag.const <- sqrt(as.complex(-1))
p <- length(phi)
q <- length(theta)
n.val <- length(omega)
val.seq <- numeric(length=n.val)
for (i.val in 1:n.val) {
  numer <- Mod(1 + sum(theta*exp(-2*(1:q)*pi*imag.const*omega[[i.val]])))
  denom <- Mod(1 - sum(phi*exp(-2*(1:p)*pi*imag.const*omega[[i.val]])))
  val.seq[i.val] <- (sigw*numer/denom)^2
}
val.seq           # 22.72351

#3) 
phi <- c(0.75, -0.9)
sigw <- 7.5
omega.seq <- seq(0,0.5,by=0.0001)
spdens <- ar.spdens(omega.seq, phi, sigw)          # 6a line 202
max(spdens)
plot(omega.seq, spdens, type="l", main="", xlab="frequency", ylab="")

index3 <- which(spdens == max(spdens))
ans <- 0.0001 * index3
ans      # 0.1853

#4) 
x.ts <- read.table('environ C.txt', header=FALSE)
x.ts <- x.ts$V1
n <- length(x.ts)
x.dft <- fft(x.ts) / sqrt(n)
x.dft
#k/n = 8/59 #so looking for (k+1)th entry which is 9th entry in the dft vector
x.dft[9]
add <- abs(Re(x.dft[9])) + abs(Im(x.dft[9]))
add  # 9.640275

# 5) 
x.pgram <- Mod(x.dft)^2
x.pgram[9]          # 47.03155

#6) 
m <-4
L <- 2*m + 1
wgts <- rep(x=1/L, times=L)
x.spgram <- smooth.pgram(x.pgram, wgts)     #smooth.pgram in 6b line 59
x.spgram[9]           # 90.02003

# 7) 
# 2 iterations
# line 179 6b
L <- 9
base.wgts <- c(0.35, 0.45, 0.2)
n.iter <- 2
wgts <- Daniell.wgts(base.wgts, n.iter)
x.spgram <- smooth.pgram(x.pgram, wgts)
x.spgram[9]     # 122.5231

#8)
df <- 2/sum(wgts^2)
df          # 10.47466

#9) 
m <- 4
L <- 2*m + 1
omega <- 8/59
ir.lag <- -m:m
pt.funct <- power.transfer(omega, wgts, ir.lag)        # line 382 6b
pt.funct       # 0.2043917


#10)
y.ts <- filter.ts(x.ts, wgts, ir.lag)             # line 395 6b
y.ts[13]        # 4.227036

