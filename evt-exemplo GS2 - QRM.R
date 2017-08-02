###############################################################################
## Exemplo do capítulo EVT sobre estimação através de processos pontuais das medidas de risco
## o método processo pontual auto-excitavel.
## Utilização do pacote QRM
## Dados da GS entre 05/05/1999 a 04/10/2016
###############################################################################

library(QRM)
library(PerformanceAnalytics)

setwd("C:\\Users\\Rfbre\\Documents\\UDESC\\TCC\\TCC Latex svn\\TCC-R-codes")

prices <- as.xts(read.zoo(read.csv("evt-exemplo GS.csv"), format = "%Y-%m-%d", FUN = as.Date))
returns <- na.omit(Return.calculate(prices, method = "log"))
losses <- -1*returns
losses <- losses[,"Adj.Close"]
names(losses) <- "losses"
tslosses <- as.timeSeries.xts(losses)

ppdata <- extremalPP(tslosses, threshold = quantile(losses, 0.9))
ppmodel <- fit.seMPP(ppdata, model = "Hawkes", mark.influence = TRUE)

####################################################################################################
## Teste
data("sp500")
l <- -returns(sp500)
lw <- window(l, start = "1995-12-31", end = end(l))
testpp <- extremalPP(lw, threshold = quantile(lw, 0.9))
testmodel <- fit.seMPP(testpp, model = "Hawkes")
traceback()
