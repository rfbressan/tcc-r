###############################################################################
## Exemplo do capítulo sobre Fatos Estilizados
## Gráfico QQ dos retornos extremos
## Dados da VALE5
###############################################################################
library(xts)
library(PerformanceAnalytics)
library(fBasics)
library(fExtremes)

setwd("C:\\Users\\Rfbre\\Documents\\UDESC\\TCC\\TCC Latex svn\\TCC-R-codes")

precos <- as.xts(read.zoo(read.csv("VALE5.csv"), format = "%Y-%m-%d", FUN = as.Date))
retornos <- na.omit(Return.calculate(precos[,"Adj.Close"], method = "log"))
# As 100 maiores perdas
extremedata <- retornos[order(coredata(retornos))[1:100]]
spaces <- as.numeric(diff(time(extremedata)))

op <- par(mfrow=c(1,2))
# Grafico dos retornos
plot(extremedata, type = "h", main = "", major.format="%m/%Y", xlab = "", ylab = "Perdas")
# Grafico QQ com relacao a uma distribuicao exponencial
qqplot(qexp(ppoints(length(spaces))), spaces, xlab = "Quantis exponenciais", ylab = "Dados ordenados")
qqline(spaces, distribution = function(p) qexp(p), col = 2)

par(op)