###############################################################################
## Exemplo do capítulo sobre EVT de medidas condicionais de risco utilizando
## o método processo pontual auto-excitavel.
## Dados da GS entre 05/05/1999 a 04/10/2016
## Utilizacao do pacote SAPP
###############################################################################
library(xts)
library(PerformanceAnalytics)
library(xtable)
library(SAPP)
library(PtProcess)
library(fExtremes)
library(gPdtest)
library(plotly)
library(boot)

setwd("C:\\Users\\Rafael\\Documents\\UDESC\\TCC\\TCC Latex svn\\TCC-R-codes")

prices <- as.xts(read.zoo(read.csv("evt-exemplo GS.csv"), format = "%Y-%m-%d", FUN = as.Date))
returns <- 100*na.omit(Return.calculate(prices, method = "log"))
losses <- -1*returns
losses <- losses[,"Adj.Close"]
names(losses) <- "losses"
u <- quantile(losses, 0.9)
indextimes <- as.numeric(index(losses)-index(losses[1,]))

## Vamos primeiro estimar os parametros da GPD
gpd <- gpdFit(as.timeSeries(losses[,"losses"]), u=u, type = "mle")
xi <- gpd@fit$par.ests[1]
beta <- gpd@fit$par.ests[2]

parameters <- c(0.05, 3.1, 1.5, 0.1, 1.1) # Valores iniciais dos parametros para lambdag
stime <- indextimes[1]
etime <- indextimes[length(indextimes)]

## Estimando o modelo ETAS
etasmodel <- etasap(indextimes, coredata(losses), threshold = u, parami = parameters, tstart = stime,
                    zte = etime, approx = 1, plot = TRUE)

############################################################################################################
## Com os parametros estimados mu, K, c, alpha e p em etasmodel$param, calcular lambdag através do pacote
## PtProcess

dados <- data.frame(indextimes, coredata(losses))
names(dados) <- c("time", "magnitude")
dados <- dados[which(dados$magnitude>u),]-u
params <- c(etasmodel$param[[1]], etasmodel$param[[2]], etasmodel$param[[4]], etasmodel$param[[3]],
            etasmodel$param[[5]])
## Calculo de lambdag atraves do pacote PtProcess com os parametros estimados por SAPP
## ######## ALGUM PROBLEMA COM ESCALA DE LAMBDA. EH PRECISO REDUZIR SEUS VALORES CALCULADOS
########### PARA AS ESTIMACOES DE VaR E ES FAZEREM MAIS SENTIDO
######################################################################################################
lambdag <- etas_gif(dados, evalpts = indextimes, params = params)

############################################################################################################
## Calculo de VaR e ES com base em lambdag

VaRESpp <- function(u, beta, xi, alpha=0.99, lambda){
  VaRpp <- u+(beta/xi)*(((1-alpha)/lambda)^(-xi)-1)
  ESpp <- ((VaRpp)/(1-xi))+((beta-xi*u)/(1-xi))
  return(cbind(VaRpp, ESpp))
}

vares <- VaRESpp(u, beta, xi, alpha=0.99, lambdag)

## Calculo de VaR e ES não condicionais, com base na GPD apenas
urisk <- gpdRiskMeasures(gpd, prob = 0.99)

## Serie xts com dados das perdas, perdas em excesso e lambdag
ts <- cbind(losses, indextimes, (losses[,"losses"]-u)-((losses[,"losses"]-u)<0)*(losses[,"losses"]-u),
            lambdag, vares)
names(ts) <- c("losses", "indextimes", "excess", "intensity", "VaR", "ES")

## Grafico comparando a intesidade estimada e as perdas em excesso
op <- par(mfrow=c(2,1))
plot(ts[,"intensity"], xlab = "", ylab = "Intensidade",
     main = "Processo Pontual de Auto-Excitação", major.format="%m/%Y")
plot(ts[,"excess"], type="h", xlab = "Data", ylab = "Perdas em excesso",
     main = "", major.format="%m/%Y")
par(op)

## Grafico do VaR e ES sobre as perdas efetivas

plot(ts[,"losses"], type="h", xlab = "Data", ylab = "Perdas",
     main = "", major.format="%m/%Y")
lines(ts[,"VaR"], col="red")
lines(ts[,"ES"], col="blue")
abline(h=urisk[c(2,3)], col=c("red", "blue"))

