###############################################################################
## Exemplo do capítulo sobre Fatos Estilizados
## Gráfico ACF e ACF dos retornos absolutos
## Dados da Swap de juros 360d, cambio spot e ibovespa
###############################################################################
library(xts)
library(PerformanceAnalytics)
library(fBasics)

setwd("C:\\Users\\Rafael\\Documents\\UDESC\\TCC\\TCC Latex svn\\TCC-R-codes")

irate <- as.xts(read.zoo(read.csv2("JuroFut.csv"), format = "%d/%m/%Y", FUN = as.Date))
xrate <- as.xts(read.zoo(read.csv2("Cambio.csv"), format = "%d/%m/%Y", FUN = as.Date))
ibov <- as.xts(read.zoo(read.csv("Ibovespa.csv"), format = "%Y-%m-%d", FUN = as.Date))

prices <- cbind(1000/(1+irate[,1]/100), 1000/(1+irate[,2]/100), 1000/(1+irate[,3]/100))
iretornos <- na.omit(Return.calculate(prices, method = "log"))
xretornos <- na.omit(Return.calculate(xrate, method = "log"))
ibovretornos <- na.omit(Return.calculate(ibov[,"Adj.Close"], method = "log"))
names(iretornos) <- c("Retornos 30d", "Retornos 180d", "Retornos 360d")
names(xretornos) <- "Retornos Câmbio"
names(ibovretornos) <- "Retornos Ibovespa"

op <- par(mfrow=c(3,2))
#Juros
acf(iretornos[,3], lag.max=15, main="", lwd=2, xlab="")
grid(nx=NA, ny=NULL, col = "black")
legend("topleft", legend = "a)", bty = "n", xjust = 0, yjust = 0)
acf(abs(iretornos[,3]), lag.max=15, main="", lwd=2, xlab="")
grid(nx=NA, ny=NULL, col = "black")
#Cambio
acf(xretornos, lag.max=15, main="", lwd=2, xlab="")
grid(nx=NA, ny=NULL, col = "black")
legend("topleft", legend = "b)", bty = "n", xjust = 0, yjust = 0)
acf(abs(xretornos), lag.max=15, main="", lwd=2, xlab="")
grid(nx=NA, ny=NULL, col = "black")
#Ibovespa
acf(ibovretornos, lag.max = 15, main="", lwd = 2)
grid(nx=NA, ny=NULL, col = "black")
legend("topleft", legend = "c)", bty = "n", xjust = 0, yjust = 0)
acf(abs(ibovretornos), lag.max = 15, main="", lwd = 2)
grid(nx=NA, ny=NULL, col = "black")
par(op)