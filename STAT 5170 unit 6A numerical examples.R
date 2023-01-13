# =================================================================================
# STAT 5170: Applied Time Series
# Numerical examples for part A of learning unit 6
# =================================================================================
#
# =================================================================================
# Aggregating simple trigonometric time series
# =================================================================================

# This example illustrated how simple trigonometric time series may be aggregated to form a time series that exhibits complex patterns.

# ---------------------------------------------------------------------------------
# Fixed time series
# ---------------------------------------------------------------------------------

# Start by specifying a time sequence 

n <- 100
time.seq <- 1:n

# Now specify three simple trigonometric time series.

# The first time series has frequency 6/100, and Fourier coefficients 2 and 3.

k <- 6
omega <- k/n
U1 <- 2
U2 <- 3
x1 <- U1*cos(2*pi*omega*time.seq) + U2*sin(2*pi*omega*time.seq)

# The second time series has frequency 10/100, and Fourier coefficients 4 and 5.

k <- 10
omega <- k/n
U1 <- 4
U2 <- 5
x2 <- U1*cos(2*pi*omega*time.seq) + U2*sin(2*pi*omega*time.seq)

# The third time series has frequency 40/100, and Fourier coefficients 6 and 7.

k <- 40
omega <- k/n
U1 <- 6
U2 <- 7
x3 <- U1*cos(2*pi*omega*time.seq) + U2*sin(2*pi*omega*time.seq)

# The code below generates separate plots of the three simple time series.

par(mfrow = c(3, 1), mai=c(0.3,0.6,0.1,0.2)) # mai=c(bottom, left, top, right)
plot(time.seq, x1, xlab="", ylab="x1", type="l")
plot(time.seq, x2, xlab="", ylab="x2", type="l")
plot(time.seq, x3, xlab="", ylab="x3", type="l")
par(mfrow = c(1, 1))

# The time series are aggregated by adding them together.

x <- x1 + x2 + x3

# The following code generates a plot of the aggregated time series

par(mfrow = c(3, 1), mai=c(0.3,0.6,0.1,0.2)) # mai=c(bottom, left, top, right)
plot(time.seq, x, xlab="", ylab="aggregated", type="l")
par(mfrow = c(1, 1))

# ---------------------------------------------------------------------------------
# Random time series
# ---------------------------------------------------------------------------------

# The example above is now repeated using random specifications of the frequencies and Fourier coefficients

# The first time series selects k randomly between 0 and 24 in calculating the frequencey omega = k/n.

k <- sample.int(n=25, size=1)-1
omega <- k/n
U1 <- runif(n=1, min=0, max=10)
U2 <- runif(n=1, min=0, max=10)
x1 <- U1*cos(2*pi*omega*time.seq) + U2*sin(2*pi*omega*time.seq)

print(paste("x1: k=", k, ", U1=", sprintf("%4.2f", U1), ", U2=", sprintf("%4.2f", U2), sep=""))
print(paste("x1: frequency=", sprintf("%4.2f", omega), ", amplitude=", sprintf("%4.2f", sqrt(U1^2 + U2^2)), sep=""))

# The second time series selects k randomly between 5 and 34 in calculating the frequencey omega = k/n.

k <- 4 + sample.int(n=30, size=1)
omega <- k/n
U1 <- runif(n=1, min=0, max=10)
U2 <- runif(n=1, min=0, max=10)
x2 <- U1*cos(2*pi*omega*time.seq) + U2*sin(2*pi*omega*time.seq)

print(paste("x2: k=", k, ", U1=", sprintf("%4.2f", U1), ", U2=", sprintf("%4.2f", U2), sep=""))
print(paste("x2: frequency=", sprintf("%4.2f", omega), ", amplitude=", sprintf("%4.2f", sqrt(U1^2 + U2^2)), sep=""))

# The third time series selects k randomly between 15 and 49 in calculating the frequencey omega = k/n.

k <- 14 + sample.int(n=35, size=1)
omega <- k/n
U1 <- runif(n=1, min=0, max=10)
U2 <- runif(n=1, min=0, max=10)
x3 <- U1*cos(2*pi*omega*time.seq) + U2*sin(2*pi*omega*time.seq)

print(paste("x3: k=", k, ", U1=", sprintf("%4.2f", U1), ", U2=", sprintf("%4.2f", U2), sep=""))
print(paste("x3: frequency=", sprintf("%4.2f", omega), ", amplitude=", sprintf("%4.2f", sqrt(U1^2 + U2^2)), sep=""))

# Separate plots of the three simple time series and the aggregated time series are generated as follows.

x <- x1 + x2 + x3
par(mfrow = c(4, 1), mai=c(0.3,0.6,0.1,0.2)) # mai=c(bottom, left, top, right)
plot(time.seq, x1, xlab="", ylab="x1", type="l")
plot(time.seq, x2, xlab="", ylab="x2", type="l")
plot(time.seq, x3, xlab="", ylab="x3", type="l")
plot(time.seq, x, xlab="", ylab="aggregated", type="l")
par(mfrow = c(1, 1))

# =================================================================================
# Spectral density of an MA(1) time series
# =================================================================================

# This example explores the spectral density of moving-average time series. We will make use of the user-defined function we have seen before that simulates such time series.

ma.sim <- function(n, theta, sigw, w0=NA) {
	q <- length(theta)
	wn.samp <- rnorm(n=n, mean=0, sd=sigw)
	if (is.na(w0)) {
		pastw <- rep(x=0, times=q)
	} else {
		pastw <- c(w0, rep(x=0, times=q-length(w0)))
	}
	ma.samp <- numeric(length=n)
	for (t in 1:n) {
		ma.samp[t] <- sum(theta*pastw) + wn.samp[t]
		pastw <- c(wn.samp[t], pastw[1:q-1])
	}
	ma.ts <- ts(ma.samp)
	return(ma.ts)
}

# ---------------------------------------------------------------------------------
# Fixed time series model
# ---------------------------------------------------------------------------------

# The following code generates a plot of a simulated MA(1) time series with moving-average parameter theta1=0.5 alongside a plot of its spectral density.

theta1 <- 0.5
n <- 200
time.seq <- 1:n
sigw <- 1
omega.seq <- seq(from=0, to=0.5,by=0.01)

par(mfrow = c(2, 1), mai=c(0.9,0.1,0.1,0.1)) # mai=c(bottom, left, top, right)
ma.ts <- ma.sim(n, theta1, sigw)
plot(ma.ts, type="l", main="", xlab="time", ylab="")
spdens <- sigw^2*(1 + theta1^2 + 2*theta1*cos(2*pi*omega.seq))
plot(omega.seq, spdens, type="l", main="", xlab="frequency", ylab="")
par(mfrow = c(1, 1))

# ---------------------------------------------------------------------------------
# Random time series model
# ---------------------------------------------------------------------------------

# The example above is now repeated using random specifications of the moving-average parameter

theta1 <- runif(n=1, min=-1, max=1)
par(mfrow = c(2, 1), mai=c(0.9,0.1,0.25,0.1)) # mai=c(bottom, left, top, right)
ma.ts <- ma.sim(n, theta1, sigw)
plot(ma.ts, type="l", xlab="time", ylab="", main=paste("theta1=", sprintf("%4.2f", theta1), sep=""))
spdens <- sigw^2*(1 + theta1^2 + 2*theta1*cos(2*pi*omega.seq))
plot(omega.seq, spdens, type="l", main="", xlab="frequency", ylab="")
par(mfrow = c(1, 1))

# =================================================================================
# Complex numbers in R
# =================================================================================

# R supports functionality for working with complex numbers. The imaginary constant, i, is specified as follows

imag.const <- sqrt(as.complex(-1))

# With this, the spectral density of an MA time series model may be calculated using the built-in function Mod() to find the modulus of a complex number, which is squared in the formula for the spectral density.

# Observe that the numerical value resulting from the MA(1) spectral density formula from the previous example is duplicated by the subsequent formula that reflects the general formula provided in the notes, and that makes use of complex arithematic.

omega <- runif(n=1, min=0, max=0.5)

sigw^2*(1 + theta1^2 + 2*theta1*cos(2*pi*omega))
sigw^2*Mod(1 + theta1*exp(-2*pi*imag.const*omega))^2

# The following code uses this functionality to define functions for calculating the spectral density of a general MA or AR time series.

sigw^2*Mod(1 + theta1*exp(-2*pi*imag.const*omega))^2

ma.spdens <- function(omega, theta, sigw) {
	imag.const <- sqrt(as.complex(-1))
	q <- length(theta)
	n.val <- length(omega)
	val.seq <- numeric(length=n.val)
	for (i.val in 1:n.val) {
		val.seq[i.val] <- (sigw*Mod(1 + sum(theta*exp(-2*(1:q)*pi*imag.const*omega[[i.val]]))))^2
	}
	return(val.seq)
}

ar.spdens <- function(omega, phi, sigw) {
	imag.const <- sqrt(as.complex(-1))
	p <- length(phi)
	n.val <- length(omega)
	val.seq <- numeric(length=n.val)
	for (i.val in 1:n.val) {
		val.seq[i.val] <- (sigw/Mod(1 - sum(phi*exp(-2*(1:p)*pi*imag.const*omega[[i.val]]))))^2
	}
	return(val.seq)
}

# The following code repeates the previous example involving an MA(1) time series using the associated user-defined function. A demonstration for autoregressive time series appears subsequently.

theta1 <- runif(n=1, min=-1, max=1)
n <- 200
time.seq <- 1:n
sigw <- 1
omega.seq <- seq(0,0.5,by=0.01)

par(mfrow = c(2, 1), mai=c(0.9,0.1,0.25,0.1)) # mai=c(bottom, left, top, right)
ma.ts <- ma.sim(n, theta1, sigw)
plot(ma.ts, type="l", xlab="time", ylab="", main=paste("theta1=", sprintf("%4.2f", theta1), sep=""))
spdens <- ma.spdens(omega.seq, theta1, sigw)
plot(omega.seq, spdens, type="l", main="", xlab="frequency", ylab="")
par(mfrow = c(1, 1))

# =================================================================================
# Spectral density of an AR(2) time series
# =================================================================================

# This example explores the spectral density of autoregressive time series. We will make use of the user-defined function we have seen before that simulates such time series.

ar.sim <- function(n, phi, sigw, x0=NA) {
	p <- length(phi)
	wn.samp <- rnorm(n=n, mean=0, sd=sigw)
	if (is.na(x0)) {
		pastx <- rep(x=0, times=p)
	} else {
		pastx <- c(x0, rep(x=0, times=p-length(x0)))
	}
	ar.samp <- numeric(length=n)
	for (t in 1:n) {
		ar.samp[t] <- sum(phi*pastx) + wn.samp[t]
		pastx <- c(ar.samp[t], pastx[1:p-1])
	}
	ar.ts <- ts(ar.samp)
	return(ar.ts)
}

# ---------------------------------------------------------------------------------
# Fixed time series model
# ---------------------------------------------------------------------------------

# The following code generates a plot of a simulated AR(2) time series, with autoregressive parameters specified below, alongside a plot of its spectral density.

phi1 <- 1.0
phi2 <- -0.9
phi <- c(phi1, phi2)
n <- 75
time.seq <- 1:n
sigw <- 1
omega.seq <- seq(0,0.5,by=0.01)

par(mfrow = c(2, 1), mai=c(0.9,0.1,0.1,0.1)) # mai=c(bottom, left, top, right)
ar.ts <- ar.sim(n, phi, sigw)
plot(ar.ts, type="l", main="", xlab="time", ylab="")
spdens <- ar.spdens(omega.seq, phi, sigw)
plot(omega.seq, spdens, type="l", main="", xlab="frequency", ylab="")
par(mfrow = c(1, 1))

# =================================================================================
# Aliasing
# =================================================================================

# The following code recreates the graphic in the notes that illustrates the aliasing concept.

# Observe that the simple trogonometric time series depicted by a solid line has frequency omega=1/4; that depicted by a dashed line has frequency omega=3/4.

time.seq <- seq(from=0, to=8, by=0.05)
plot(time.seq, cos(2*pi*time.seq*1/4), type="l", lty=1, lwd=1, xlab="time", ylab="y", main="", axes=F)
lines(time.seq, cos(2*pi*time.seq*3/4), lty=2, lwd=1)
points(c(0:8), cos(2*pi*c(0:8)/4), pch=19, cex=1)
abline(h=0)
axis(1, at=c(1,2,3,4,5,6,7))
axis(1)
axis(2)
box()

# =================================================================================
# Star brightness
# =================================================================================

# The star brightness time series is included in the "astsa" package.

data(star, package="astsa")
ts.plot(star, xlab="", ylab="brightness")
n <- length(star)

# R includes a built-in function called "fft" that calculates a rescaling of the discrete Fourier transform of a time series using the fast Fourier transform algorithm. To calcualte the DFT, divide the output of the fft() function by the square-root of the time-series length.

star.dft <- fft(star) / sqrt(n)

# This produces a vector of length n,each of whose entry is the DFT value at one of the fundamental Fourier frequencies, starting with frequency 0. That is, the output corresponds to the fundamental Fourier frequencies k/n for k=0, ..., n-1.

# The the inverse DFT is calculated by applying the fft() function to the DFT vector, with the inverse=TRUE option declared, and dividing the output by the square-root of the time-series length.

star.idft <- fft(star.dft, inverse=TRUE) / sqrt(n)

# As it is applied in the code above, the original time series is recovered.

t <- sample.int(n=n, size=1)
t
star[t]
star.idft[t]
Re(star.idft[t])

# Note: the Re() function returns the the real part of a complext number. The function Im() returns the imaginary part.

# The periodogram is calculated as the squared modulus of the DFT.

star.pgram <- Mod(star.dft)^2

# Besides the frequency k/n for k=0, which corresponds to the sample mean, the two periodicities with highest strength appear in entries 22 and 26 of the vector output.

max.idx <- 50
ord.seq <- order(-star.pgram[1:max.idx])
hi1.idx <- ord.seq[2]
hi2.idx <- ord.seq[3]
c(hi1.idx, hi2.idx)

# Since the first entry of the vector output is for k=0, these values correspond to the values k=21 and k=25. The corresponding fundamental Fourier frequencies are k/n=21/600 and k/n=25/600.

# The corresponding periods are n/k=600/21=28.6 and n/k=600/25=24 days.

# A plot of the time series and periodogram (omitting k=0) are generated by the code below.

par(mfrow = c(2, 1), mai=c(0.50,1.00,0.1,0.1)) # mai=c(bottom, left, top, right)
ts.plot(star, xlab="", ylab="brightness")
freq.seq <- (1:n-1)/n
plot(freq.seq[2:max.idx], star.pgram[2:max.idx], type="h", lty=1, lwd=3, ylab="Periodogram", xlab="Frequency")
text(x=0.050, y=7000, labels="24 day cycle")
text(x=0.027, y=9000, labels="29 day cycle")
par(mfrow = c(1, 1))

# =================================================================================
# Southern Oscillation Index (SOI) time series
# =================================================================================

# The next example examines the fish population count data next alongside a time series that measures the Southern Oscillation Index (SOI). Both of these are found in the "astsa" package.

data(soi, package="astsa")
data(rec, package="astsa")

# The data sets are measured at the same times, and have the same length.

length(soi)
length(rec)

# The following code generates a plot the two time series.

par(mfrow = c(2, 1), mai=c(0.50,1.00,0.1,0.1)) # mai=c(bottom, left, top, right)
ts.plot(soi, ylab="SOI", main="")
ts.plot(rec, ylab="count", main="")
par(mfrow = c(1, 1))

# ---------------------------------------------------------------------------------
# Periodograms
# ---------------------------------------------------------------------------------

# The two time series' periodograms are calculated by the following code. Each time series is centered by subtracting its sample mean, and also scaled by a factor that allows its results to match an existing set of calculations made elsewhere.

n <- length(rec)
soi.dft <- fft(soi) / sqrt(n)
soi.pgram <- Mod(soi.dft)^2
rec.dft <- fft(rec) / sqrt(n)
rec.pgram <- Mod(rec.dft)^2

# The code below generates a plot of the periodograms, with frequencies translated so that one unit is a yearly frequency of omega=1/12, two units is a half-year frequency of omega=2/12, etc.

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
# Confidence intervals
# ---------------------------------------------------------------------------------

# Confidence intervals for the periodogram associated with the one-year (omega=1/12) cycles, four-year (omega=1/48) cycles, and those of any other frequency are are calculated as follows. 

# The value associated with the one year cycle found at element 39 of the periodogram vector.

i.val <- 39

# This value corresponds to k=38, for which the corresponding frequency is omega=k/n=0.0839. The period is n/k=11.92, or about 12 months.

# The maximum value of each time series's periodogram is at this frequency.

soi.pgram[i.val]
max(soi.pgram)
 
rec.pgram[i.val]
max(rec.pgram)
 
# Confidence intervals are calculated as follows

alpha <- 0.05
cutoff <- qchisq(p=c(alpha/2, 1-alpha/2), df=2)

# SOI time series
val <- soi.pgram[i.val]
lo.bnd <- 2*val/cutoff[2]
hi.bnd <- 2*val/cutoff[1]
print(paste("Point estimate: ", sprintf("%6.4f", val), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", lo.bnd), ", ", sprintf("%6.4f", hi.bnd), "]", sep=""))
print("Log transform:")
print(paste("Point estimate: ", sprintf("%6.4f", log(val)), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", log(lo.bnd)), ", ", sprintf("%6.4f", log(hi.bnd)), "]", sep=""))

# Fish population count time series
val <- rec.pgram[i.val]
lo.bnd <- 2*val/cutoff[2]
hi.bnd <- 2*val/cutoff[1]
print(paste("Point estimate: ", sprintf("%6.4f", val), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", lo.bnd), ", ", sprintf("%6.4f", hi.bnd), "]", sep=""))
print("Log transform:")
print(paste("Point estimate: ", sprintf("%6.4f", log(val)), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", log(lo.bnd)), ", ", sprintf("%6.4f", log(hi.bnd)), "]", sep=""))

# The value associated with the four year cycle found at element 9 of the periodogram vector.

i.val <- 10

# This value corresponds to k=9, for which the corresponding frequency is omega=k/n=0.0199. The period is n/k=50.33, or a little over 48 months.

# The associated confidence intervals are calculated as follows

# SOI time series
val <- soi.pgram[i.val]
lo.bnd <- 2*val/cutoff[2]
hi.bnd <- 2*val/cutoff[1]
print(paste("Point estimate: ", sprintf("%6.4f", val), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", lo.bnd), ", ", sprintf("%6.4f", hi.bnd), "]", sep=""))
print("Log transform:")
print(paste("Point estimate: ", sprintf("%6.4f", log(val)), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", log(lo.bnd)), ", ", sprintf("%6.4f", log(hi.bnd)), "]", sep=""))

# Fish population count time series
val <- rec.pgram[i.val]
lo.bnd <- 2*val/cutoff[2]
hi.bnd <- 2*val/cutoff[1]
print(paste("Point estimate: ", sprintf("%6.4f", val), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", lo.bnd), ", ", sprintf("%6.4f", hi.bnd), "]", sep=""))
print("Log transform:")
print(paste("Point estimate: ", sprintf("%6.4f", log(val)), sep=""))
print(paste(sprintf("%2.0f", (1-alpha)*100), "% confidence interval: [", sprintf("%6.4f", log(lo.bnd)), ", ", sprintf("%6.4f", log(hi.bnd)), "]", sep=""))


