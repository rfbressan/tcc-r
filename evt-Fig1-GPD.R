# Simula uma GPD e sua densidade para alguns valores de qsi e beta

library(fExtremes)
setwd("C:\\Users\\Rafael\\Documents\\UDESC\\TCC\\TCC Latex svn\\R codes")

x <- seq(0.1, 6, 0.10)

gminus2 <- c(1, dgpd(x, -0.5))
g5 <- c(1, dgpd(x, 0.5))
g0 <- c(1, dgpd(x, 0))
Gminus2 <- c(0, pgpd(x, -0.5))
G5 <- c(0, pgpd(x, 0.5))
G0 <- c(0, pgpd(x, 0))

x <- c(0, x)

op <- par(mfrow=c(1,2))
plot(x, Gminus2, type="l", lty = 3, lwd = 1.5, ylab = "G(x)", main="c.d.f")
lines(x, G5, lwd = 1.5)
lines(x, G0, lty = 2, lwd = 1.5)
plot(x, gminus2, type="l", lty = 3, lwd = 1.5, ylab = "g(x)", main="p.d.f")
lines(x, g5, lwd = 1.5)
lines(x, g0, lty = 2, lwd = 1.5)

par(op)
