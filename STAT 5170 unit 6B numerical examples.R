# =================================================================================
# STAT 5170: Applied Time Series
# Numerical examples for part B of learning unit 6
# =================================================================================
#
# =================================================================================
# Smoothing the periodogram
# =================================================================================

# Recall the example from part A of the learning unit in which we explored the the fish population count data next alongside a time series that measures the Southern Oscillation Index (SOI), both of which are found in the "astsa" package.

data(soi, package="astsa")
data(rec, package="astsa")

# Plots of the data are generated as follows.

par(mfrow = c(2, 1), mai=c(0.50,1.00,0.1,0.1)) # mai=c(bottom, left, top, right)
ts.plot(soi, ylab="SOI", main="")
ts.plot(rec, ylab="count", main="")
par(mfrow = c(1, 1))

# Corresponding periodograms are calculated after centering with respect the a time series' sample mean and scaling by a factor that induces the results of a parallel analysis calculated elsewhere.

n <- length(rec)
soi.dft <- fft(soi) / sqrt(n)
soi.pgram <- Mod(soi.dft)^2
rec.dft <- fft(rec) / sqrt(n)
rec.pgram <- Mod(rec.dft)^2

# To remove the influcence of the sample mean, the first entry of each periodogram is set to zero.

soi.pgram[1] <- 0
rec.pgram[1] <- 0

# This is equivalent to centering each time series by subtracting its sample mean. (Can you see why?)

# Plots of the periodograms are generated as follows. Frequencies are translated so that one unit is a yearly frequency of omega=1/12, two units is a half-year frequency of omega=2/12, etc.

n.yr <- 453/12
freq.seq <- 1:ceiling(6*n.yr)
freq.idx <- freq.seq+1
par(mfrow = c(2, 1), mai=c(0.50,1.00,0.1,0.1)) # mai=c(bottom, left, top, right)
plot(freq.seq/n.yr, soi.pgram[freq.idx], type="l", lty=1, lwd=1, ylab="Periodogram", xlab="Frequency")
abline(v=1/4, lty="dotted")
plot(freq.seq/n.yr, rec.pgram[freq.idx], type="l", lty=1, lwd=1, ylab="Periodogram", xlab="Frequency")
abline(v=1/4, lty="dotted")
par(mfrow = c(1, 1))

# ---------------------------------------------------------------------------------
# Unweighted smoothing
# ---------------------------------------------------------------------------------

# The following user-defined function can be used for general smoothing, unweighted or weighted. The argument 'wgts' is a sequence of weights that sum to one. Its length, L, is assumed to follow from the the formula
#
# L = 2*m+1
#
# which implies that L must be odd.

smooth.pgram <- function(raw.pgram, wgts) {
	n.freq <- length(raw.pgram)
	L <- length(wgts)
	m <- 0.5*(L-1)
	new.pgram <- numeric(length=n.freq)
	for (i.freq in 1:m) {
		idx.seq <- (m-i.freq+2):L
		wgts0 <- wgts[idx.seq] / sum(wgts[idx.seq])
		new.pgram[i.freq] <- sum(wgts0*raw.pgram[1:(i.freq+m)])
	}
	for (i.freq in (1+m):(n.freq-m)) {
		new.pgram[i.freq] <- sum(wgts*raw.pgram[(i.freq-m):(i.freq+m)])
	}
	for (i.freq in (n.freq-m+1):n.freq) {
		idx.seq <- 1:(m+(n.freq-i.freq)+1)
		wgts0 <- wgts[idx.seq] / sum(wgts[idx.seq])
		new.pgram[i.freq] <- sum(wgts0*raw.pgram[(i.freq-m):n.freq])
	}
	return(new.pgram)
}

# For unweighted smoothing with m=4, hence L=2*m+1=9, define the weights as follows

m <- 4
L <- 2*m + 1
wgts <- rep(x=1/L, times=L)
wgts

# Smoothed periodograms and corresponding plots are generated as follows:

soi.spgram <- smooth.pgram(soi.pgram, wgts)
rec.spgram <- smooth.pgram(rec.pgram, wgts)

n.yr <- 453/12
freq.seq <- 1:ceiling(6*n.yr)
freq.idx <- freq.seq+1
par(mfrow = c(2, 1), mai=c(0.50,1.00,0.1,0.1)) # mai=c(bottom, left, top, right)
plot(freq.seq/n.yr, soi.spgram[freq.idx], type="l", lty=1, lwd=1, ylab="Unweighted smoothing", xlab="Frequency")
abline(v=1/4, lty=2)
plot(freq.seq/n.yr, rec.spgram[freq.idx], type="l", lty=1, lwd=1, ylab="Unweighted smoothing", xlab="Frequency")
abline(v=1/4, lty=2)
par(mfrow = c(1, 1))

# ---------------------------------------------------------------------------------
# Confidence intervals
# ---------------------------------------------------------------------------------

# The confidence interval calculation carried out in part A of the learning unit may be repeated on the smoothed periodogram.

# With a smoothed periodogram, chi-square quantiles are calculated using a special formula for degrees of freedom:

df <- 2/sum(wgts^2)

# When the weights are uniform, as they are in this example, this value is identical to the value 2L, the length of the weights.

c(2L, df)

# Confidence intervals for the periodogram associated with the one-year (omega=1/12) cycles and four-year (omega=1/48) cycles are calculated as follows. 

# The value associated with the one year cycle found at element 39 of the periodogram vector.

i.val <- 39
 
# Confidence intervals are calculated as follows

alpha <- 0.05
cutoff <- qchisq(p=c(alpha/2, 1-alpha/2), df=df)

# SOI time series
val <- soi.spgram[i.val]
lo.bnd <- df*val/cutoff[2]
hi.bnd <- df*val/cutoff[1]
print(paste("Point estimate: ", sprintf("%6.4f", val), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", lo.bnd), ", ", sprintf("%6.4f", hi.bnd), "]", sep=""))
print("Log transform:")
print(paste("Point estimate: ", sprintf("%6.4f", log(val)), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", log(lo.bnd)), ", ", sprintf("%6.4f", log(hi.bnd)), "]", sep=""))

# Fish population count time series
val <- rec.spgram[i.val]
lo.bnd <- df*val/cutoff[2]
hi.bnd <- df*val/cutoff[1]
print(paste("Point estimate: ", sprintf("%6.4f", val), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", lo.bnd), ", ", sprintf("%6.4f", hi.bnd), "]", sep=""))
print("Log transform:")
print(paste("Point estimate: ", sprintf("%6.4f", log(val)), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", log(lo.bnd)), ", ", sprintf("%6.4f", log(hi.bnd)), "]", sep=""))

# The value associated with the four year cycle found at element 10 of the periodogram vector.

i.val <- 10

# The associated confidence intervals are calculated as follows

# SOI time series
val <- soi.spgram[i.val]
lo.bnd <- df*val/cutoff[2]
hi.bnd <- df*val/cutoff[1]
print(paste("Point estimate: ", sprintf("%6.4f", val), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", lo.bnd), ", ", sprintf("%6.4f", hi.bnd), "]", sep=""))
print("Log transform:")
print(paste("Point estimate: ", sprintf("%6.4f", log(val)), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", log(lo.bnd)), ", ", sprintf("%6.4f", log(hi.bnd)), "]", sep=""))

# Fish population count time series
val <- rec.spgram[i.val]
lo.bnd <- df*val/cutoff[2]
hi.bnd <- df*val/cutoff[1]
print(paste("Point estimate: ", sprintf("%6.4f", val), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", lo.bnd), ", ", sprintf("%6.4f", hi.bnd), "]", sep=""))
print("Log transform:")
print(paste("Point estimate: ", sprintf("%6.4f", log(val)), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", log(lo.bnd)), ", ", sprintf("%6.4f", log(hi.bnd)), "]", sep=""))

# ---------------------------------------------------------------------------------
# Daniell kernels
# ---------------------------------------------------------------------------------

# A recursive algorithm that generates the weights of a Daniell kernel is coded into a user-defined function as follows. The algorithm begins with an initial set of weights and takes it through a give number of iterations through a refinement step.

Daniell.wgts <- function(base.wgts, n.iter) {
	for (iter in 0:n.iter) {
		if (iter == 0) {
			L2 <- length(base.wgts)
			wgts <- base.wgts
		} else {
			L2 <- 2*L1 - 1
			wgts <- rep(x=0, times=L2)
			for (i.wgt in 1:L1) {
				wgts[i.wgt:(i.wgt+L1-1)] <- wgts[i.wgt:(i.wgt+L1-1)] + last.wgts[i.wgt]*last.wgts
			}
		}
		last.wgts <- wgts
		L1 <- L2
	}
	return(wgts)
}

# The weight calculation is demonstrated as follows. For initial exploration, start with an initial set of uniform weights of length L = 2*m + 1 = 3, which corresponds to the setting m=1. It is then passed through a single iteration of the refinement algorithm to produce a set of weights of length L = 2*m + 1 = 5, which corresponds to the setting m=2. Starting with an initial set of uniform weights produces what is called an "ordinary" Daniell kernel.

base.wgts.A <- c(1, 1, 1)/3
n.iter <- 1
wgts.A <- Daniell.wgts(base.wgts.A, n.iter)
wgts.A

# The following repeats the calculation using a different set of initial weights. This initial setting produces what is called an "modified" Daniell kernel.

base.wgts.B <- c(1, 2, 1)/4
wgts.B <- Daniell.wgts(base.wgts.B, n.iter)
wgts.B

# Notice that in all cases, the weights sum to one:

sum(base.wgts.A)
sum(wgts.A)
sum(base.wgts.B)
sum(wgts.B)

# Corresponding plots are generated as follows.

par(mfrow = c(2, 2), mai=c(0.50,1.00,0.1,0.1)) # mai=c(bottom, left, top, right)
disp.wgts <- base.wgts.A
L <- length(disp.wgts)
m <- 0.5*(L - 1)
plot(c(-m, m), c(0, 0), type="l", lty=2, lwd=1, xlim=c(-m-0.25, m+0.25), ylim=c(0, 1.1*max(disp.wgts)), col="white", ylab="weight", xlab="", main="base weights (0 iter)")
for (i.wgt in 1:L) {
	lines((i.wgt-m-1)*c(1, 1), c(0, disp.wgts[i.wgt]), lty=1, lwd=3, col="blue")
}
abline(h=0, lty=2)
disp.wgts <- base.wgts.B
L <- length(disp.wgts)
m <- 0.5*(L - 1)
plot(c(-m, m), c(0, 0), type="l", lty=2, lwd=1, xlim=c(-m-0.25, m+0.25), ylim=c(0, 1.1*max(disp.wgts)), col="white", ylab="weight", xlab="", main="base weights (0 iter)")
for (i.wgt in 1:L) {
	lines((i.wgt-m-1)*c(1, 1), c(0, disp.wgts[i.wgt]), lty=1, lwd=3, col="blue")
}
abline(h=0, lty=2)
disp.wgts <- wgts.A
L <- length(disp.wgts)
m <- 0.5*(L - 1)
plot(c(-m, m), c(0, 0), type="l", lty=2, lwd=1, xlim=c(-m-0.25, m+0.25), ylim=c(0, 1.1*max(disp.wgts)), col="white", ylab="weight", xlab="", main=paste("final weights (", n.iter, " iter)", sep=""))
for (i.wgt in 1:L) {
	lines((i.wgt-m-1)*c(1, 1), c(0, disp.wgts[i.wgt]), lty=1, lwd=3, col="blue")
}
abline(h=0, lty=2)
disp.wgts <- wgts.B
L <- length(disp.wgts)
m <- 0.5*(L - 1)
plot(c(-m, m), c(0, 0), type="l", lty=2, lwd=1, xlim=c(-m-0.25, m+0.25), ylim=c(0, 1.1*max(disp.wgts)), col="white", ylab="weight", xlab="", main=paste("final weights (", n.iter, " iter)", sep=""))
for (i.wgt in 1:L) {
	lines((i.wgt-m-1)*c(1, 1), c(0, disp.wgts[i.wgt]), lty=1, lwd=3, col="blue")
}
abline(h=0, lty=2)
par(mfrow = c(1, 1))

# Suppose instead the initial set of weights is

base.wgts <- c(1, 2, 2, 2, 2, 2, 1)/12

# The length of this set of weights is L = 2*m + 1 = 7, which corresponds to m=3. Passing it through one iteration of the refinement step produces a set of weights having length L = 2*m + 1 = 13, which corresponds to m=6.

n.iter <- 1
wgts <- Daniell.wgts(base.wgts, n.iter)
wgts

# The pattern here is that each iteration of the refinement step takes a set of weights having length L1 to a new set of weights having length L2 = 2*L1 - 1. In terms of the parameter m, the corresponding relationship is that a set of weights with that parameter at m=m1 expands to a set of weights with the parameter at m=m2=2*m1.

# A plot of the inital and final weights is generated by the following code.

par(mfrow = c(2, 1), mai=c(0.50,1.00,0.1,0.1)) # mai=c(bottom, left, top, right)
disp.wgts <- base.wgts
L <- length(disp.wgts)
m <- 0.5*(L - 1)
plot(c(-m, m), c(0, 0), type="l", lty=2, lwd=1, xlim=c(-m-0.25, m+0.25), ylim=c(0, 1.1*max(disp.wgts)), col="white", ylab="weight", xlab="", main="base weights (0 iter)")
for (i.wgt in 1:L) {
	lines((i.wgt-m-1)*c(1, 1), c(0, disp.wgts[i.wgt]), lty=1, lwd=3, col="blue")
}
abline(h=0, lty=2)
disp.wgts <- wgts
L <- length(disp.wgts)
m <- 0.5*(L - 1)
plot(c(-m, m), c(0, 0), type="l", lty=2, lwd=1, xlim=c(-m-0.25, m+0.25), ylim=c(0, 1.1*max(disp.wgts)), col="white", ylab="weight", xlab="", main=paste("final weights (", n.iter, " iter)", sep=""))
for (i.wgt in 1:L) {
	lines((i.wgt-m-1)*c(1, 1), c(0, disp.wgts[i.wgt]), lty=1, lwd=3, col="blue")
}
abline(h=0, lty=2)
par(mfrow = c(1, 1))

# ---------------------------------------------------------------------------------
# Smoothing with Daniell kernels
# ---------------------------------------------------------------------------------

# The following code produces smoothed periodograms and corresponding plots of the fish population count and SOI time series using the modified Daniell kernel specified above.

soi.spgram <- smooth.pgram(soi.pgram, wgts)
rec.spgram <- smooth.pgram(rec.pgram, wgts)

n.yr <- 453/12
freq.seq <- 1:ceiling(6*n.yr)
freq.idx <- freq.seq+1
par(mfrow = c(2, 1), mai=c(0.50,1.00,0.1,0.1)) # mai=c(bottom, left, top, right)
plot(freq.seq/n.yr, soi.spgram[freq.idx], type="l", lty=1, lwd=1, ylab="Weighted smoothing", xlab="Frequency")
abline(v=1/4, lty=2)
plot(freq.seq/n.yr, rec.spgram[freq.idx], type="l", lty=1, lwd=1, ylab="Weighted smoothing", xlab="Frequency")
abline(v=1/4, lty=2)
par(mfrow = c(1, 1))

# The peaks in these smoothed periodograms are much sharper than those produced using uniform smoothing.

# ---------------------------------------------------------------------------------
# Confidence intervals
# ---------------------------------------------------------------------------------

# The following code repeats the previous confidence intervals calculations on the smoothed periodograms with weighted smoothing. The special formula for degrees of freedom is  

df <- 2/sum(wgts^2)
df

# Don't be thrown of by this being a non-integer value for degrees of freedom.

# Confidence intervals for the periodogram associated with the one-year (omega=1/12) cycles and four-year (omega=1/48) cycles are calculated as follows. 

# The value associated with the one year cycle found at element 39 of the periodogram vector.

i.val <- 39
 
# Confidence intervals are calculated as follows

alpha <- 0.05
cutoff <- qchisq(p=c(alpha/2, 1-alpha/2), df=df)

# SOI time series
val <- soi.spgram[i.val]
lo.bnd <- df*val/cutoff[2]
hi.bnd <- df*val/cutoff[1]
print(paste("Point estimate: ", sprintf("%6.4f", val), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", lo.bnd), ", ", sprintf("%6.4f", hi.bnd), "]", sep=""))
print("Log transform:")
print(paste("Point estimate: ", sprintf("%6.4f", log(val)), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", log(lo.bnd)), ", ", sprintf("%6.4f", log(hi.bnd)), "]", sep=""))

# Fish population count time series
val <- rec.spgram[i.val]
lo.bnd <- df*val/cutoff[2]
hi.bnd <- df*val/cutoff[1]
print(paste("Point estimate: ", sprintf("%6.4f", val), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", lo.bnd), ", ", sprintf("%6.4f", hi.bnd), "]", sep=""))
print("Log transform:")
print(paste("Point estimate: ", sprintf("%6.4f", log(val)), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", log(lo.bnd)), ", ", sprintf("%6.4f", log(hi.bnd)), "]", sep=""))

# The value associated with the four year cycle found at element 10 of the periodogram vector.

i.val <- 10

# The associated confidence intervals are calculated as follows

# SOI time series
val <- soi.spgram[i.val]
lo.bnd <- df*val/cutoff[2]
hi.bnd <- df*val/cutoff[1]
print(paste("Point estimate: ", sprintf("%6.4f", val), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", lo.bnd), ", ", sprintf("%6.4f", hi.bnd), "]", sep=""))
print("Log transform:")
print(paste("Point estimate: ", sprintf("%6.4f", log(val)), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", log(lo.bnd)), ", ", sprintf("%6.4f", log(hi.bnd)), "]", sep=""))

# Fish population count time series
val <- rec.spgram[i.val]
lo.bnd <- df*val/cutoff[2]
hi.bnd <- df*val/cutoff[1]
print(paste("Point estimate: ", sprintf("%6.4f", val), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", lo.bnd), ", ", sprintf("%6.4f", hi.bnd), "]", sep=""))
print("Log transform:")
print(paste("Point estimate: ", sprintf("%6.4f", log(val)), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", log(lo.bnd)), ", ", sprintf("%6.4f", log(hi.bnd)), "]", sep=""))

# =================================================================================
# Linear filters
# =================================================================================

# The following user-defined function calculates the power transfer function from a given impulse response function, at a given frequency. The lag values associated with the impulse response function must also be provided. 

power.transfer <- function(omega, ir.funct, ir.lag) {
	imag.const <- sqrt(as.complex(-1))
	n.lag <- length(ir.lag)
	fr.funct <- 0
	for (i.lag in 1:n.lag) {
		fr.funct <- fr.funct + ir.funct[i.lag]*exp(-2*pi*imag.const*omega*ir.lag[i.lag])
	}
	pt.val <- Mod(fr.funct)^2
	return(pt.val)
}

# The function below applies the linear-filter specifications to a given time series. Data values at times before or after the times of the measured time series are set to the sample mean.

filter.ts <- function(x.ts, ir.funct, ir.lag) {
	n <- length(x.ts)
	mu.hat <- mean(x.ts)
	val.seq <- numeric(length=length(ir.lag))
	y.ts <- x.ts
	for (t in 1:n) {
		idx.seq <- t-ir.lag
		sub.idx.in <- which((idx.seq >= 1) & (idx.seq <= n))
		val.seq[sub.idx.in] <- x.ts[idx.seq[sub.idx.in]]
		sub.idx.out <- which((idx.seq < 1) | (idx.seq > n))
		val.seq[sub.idx.out] <- mu.hat
		y.ts[t] <- sum(ir.funct*val.seq)
	}
	return(y.ts)
}

# ---------------------------------------------------------------------------------
# High-pass filter
# ---------------------------------------------------------------------------------

# A high-pass filter that is equivalent to differencing is defined by the following specifications

ir.lag <- c(0, 1)
ir.funct <- c(1, -1)

# The code below produces a plot of the impulse response function.

min.lag <- min(ir.lag); max.lag <- max(ir.lag); rng.lag <- max.lag-min.lag
min.funct <- min(ir.funct); max.funct <- max(ir.funct); rng.funct <- max.funct-min.funct
plot(c(min.lag, max.lag), c(0, 0), type="l", lty=2, lwd=1, xlim=c(min.lag-0.1*rng.lag, max.lag+0.1*rng.lag), ylim=c(min.funct-0.1*rng.funct, max.funct+0.1*rng.funct), col="white", ylab="weight", xlab="", main="impulse response function")
for (i.lag in 1:length(ir.lag)) {
	lines(ir.lag[i.lag]*c(1, 1), c(0, ir.funct[i.lag]), lty=1, lwd=3, col="blue")
}
abline(h=0, lty=2)

# A plot of the corresponding power transfer function is generated as follows

omega.seq <- seq(from=0, to=0.5,by=0.01)
pt.funct <- power.transfer(omega.seq, ir.funct, ir.lag)
plot(omega.seq, pt.funct, type="l", lty=1, lwd=1, xlab="frequency", ylab="", main="high-pass filter")

# The following code generates a plot of the SOI time series and a plot of the same time series passed through the filter.

y.ts <- filter.ts(soi, ir.funct, ir.lag)
par(mfrow = c(2, 1), mai=c(0.50,1.00,0.1,0.1)) # mai=c(bottom, left, top, right)
ts.plot(soi, type="l", lty=1, lwd=1, ylab="raw value", xlab="time")
ts.plot(y.ts, type="l", lty=1, lwd=1, ylab="filtered value", xlab="time")
par(mfrow = c(1, 1))

# The next set of code produces plots of periodograms for the SOI time series and the same time series passed through the filter.

n.yr <- 453/12
freq.seq <- 1:ceiling(6*n.yr)
freq.idx <- freq.seq+1
pt.funct <- power.transfer(freq.seq/453, ir.funct, ir.lag)
par(mfrow = c(2, 1), mai=c(0.50,1.00,0.1,0.1)) # mai=c(bottom, left, top, right)
plot(freq.seq/n.yr, soi.pgram[freq.idx], type="l", lty=1, lwd=1, ylab="raw value", xlab="Frequency")
abline(v=1/4, lty="dotted")
plot(freq.seq/n.yr, pt.funct*soi.pgram[freq.idx], type="l", lty=1, lwd=1, ylab="filtered value", xlab="Frequency")
abline(v=1/4, lty="dotted")
par(mfrow = c(1, 1))

# ---------------------------------------------------------------------------------
# Low-pass filter
# ---------------------------------------------------------------------------------

# A low-pass filter is defined from a modified Daniell kernel specified as follows

m <- 6
L <- 2*m + 1
ir.lag <- -m:m
ir.funct <- rep(x=0.5/m, times=L)
ir.funct[c(1,L)] <- 0.25/m

# The code below produces a plot of the impulse response function.

min.lag <- min(ir.lag); max.lag <- max(ir.lag); rng.lag <- max.lag-min.lag
min.funct <- min(ir.funct); max.funct <- max(ir.funct); rng.funct <- max.funct-min.funct
plot(c(min.lag, max.lag), c(0, 0), type="l", lty=2, lwd=1, xlim=c(min.lag-0.1*rng.lag, max.lag+0.1*rng.lag), ylim=c(min.funct-0.1*rng.funct, max.funct+0.1*rng.funct), col="white", ylab="weight", xlab="", main="impulse response function")
for (i.lag in 1:length(ir.lag)) {
	lines(ir.lag[i.lag]*c(1, 1), c(0, ir.funct[i.lag]), lty=1, lwd=3, col="blue")
}
abline(h=0, lty=2)

# A plot of the corresponding power transfer function is generated as follows
omega.seq <- seq(from=0, to=0.5,by=0.01)
pt.funct <- power.transfer(omega.seq, ir.funct, ir.lag)
plot(omega.seq, pt.funct, type="l", lty=1, lwd=1, xlab="frequency", ylab="", main="low-pass filter")

# The following code generates a plot of the SOI time series and a plot of the same time series passed through the filter.

y.ts <- filter.ts(soi, ir.funct, ir.lag)
par(mfrow = c(2, 1), mai=c(0.50,1.00,0.1,0.1)) # mai=c(bottom, left, top, right)
ts.plot(soi, type="l", lty=1, lwd=1, ylab="raw value", xlab="time")
ts.plot(y.ts, type="l", lty=1, lwd=1, ylab="filtered value", xlab="time")
par(mfrow = c(1, 1))

# The next set of code produces plots of periodograms for the SOI time series and the same time series passed through the filter.

n.yr <- 453/12
freq.seq <- 1:ceiling(6*n.yr)
freq.idx <- freq.seq+1
pt.funct <- power.transfer(freq.seq/453, ir.funct, ir.lag)
par(mfrow = c(2, 1), mai=c(0.50,1.00,0.1,0.1)) # mai=c(bottom, left, top, right)
plot(freq.seq/n.yr, soi.pgram[freq.idx], type="l", lty=1, lwd=1, ylab="raw value", xlab="Frequency")
abline(v=1/4, lty="dotted")
plot(freq.seq/n.yr, pt.funct*soi.pgram[freq.idx], type="l", lty=1, lwd=1, ylab="filtered value", xlab="Frequency")
abline(v=1/4, lty="dotted")
par(mfrow = c(1, 1))
