###############################################################################
## Artigo a ser proposto na disciplina MEFCA. Baseado no método de McNeil2000
## Filtar os retornos diários de vários índices primeiramente por um modelo 
## ARMA-GARCH para em seguida ajustar uma distribuição GPD aos resíduos.
## Com esta modelagem é possível estimar os valores de VaR e ES
###############################################################################

library(fGarch)
library(fExtremes)
library(QRM)
library(rugarch)
library(timeSeries)
library(xts)
library(PerformanceAnalytics)
library(xtable)

start <- as.Date("2005-08-31")
end <- as.Date("2016-08-31")
backstart <- as.Date("2016-09-01")

tb <- read.csv("artigo-BVSP.csv", stringsAsFactors = FALSE)
tb <- tb[-which(tb$Adj.Close == "null"),] # remove linhas com "null"
tb[,2:7] <- lapply(tb[,2:7], as.numeric)
prices <- as.xts(read.zoo(tb, format = "%Y-%m-%d", FUN = as.Date))

returns <- 100*na.omit(Return.calculate(prices$Adj.Close, method = "log"))
losses <- -1*returns

###################################################################################
## Lembrando, o modelo das perdas é AR(1) e a volatilidade é GARCH(1,1)
## L_t=mu_t+e_t      mu_t=phi*mu_t-1
## e_t=sigma_t*z_t   sigma^2_t=alpha0+alpha1*e^2_t-1+beta*sigma^2t-1
## O modelo retorno 5 parametros:
## mu = valor do intercepto da equacao da media, nosso modelo=0 mas pode não ser
## ar1 = phi do modelo
## omega = alpha0 do modelo
## alpha1 = alpha1 do modelo
## beta1 = beta do modelo
###################################################################################

# Lfit1 usa método de fGarch
# Lfit2 usa método de rugarch
Lfit1 <- garchFit(formula = ~arma(1,0)+garch(1,1), data = losses[,"Adj.Close"],
                 include.mean = F, algorithm = "lbfgsb+nm")

ruspec <- ugarchspec(mean.model = list(armaOrder = c(1,0)))
Lfit2 <- ugarchfit(ruspec, losses["2009-01/","Adj.Close"], solver = "hybrid")
show(Lfit2)
op <- par(mfrow=c(2,3))
plot(Lfit2, which = 1)
plot(Lfit2, which = 4)
plot(Lfit2, which = 5)
plot(Lfit2, which = 9)
plot(Lfit2, which = 10)
plot(Lfit2, which = 11)
par(op)

matcoef <- Lfit2@fit$matcoef
dimnames(matcoef) <- list(c("$\\mu$", "$\\phi_1$", "$\\omega$", "$\\alpha_1$", "$\\beta_1$"),
                          c("Estimativa", "Erro Padr\\~ao", "Valor t", "Pr(>|t|)"))

print.xtable(xtable(matcoef, caption = "Par\\^ametros estimados para o modelo AR-GARCH de AAPL",
                    label = "tab:tabevtAAPL", digits = 4),
             file = "..\\tables\\tabevtAAPL.tex", sanitize.text.function = function(x) {x})

#########################################################################################
## Modelo EVT para os residuos padronizados
## Os residuos podem ser retirados atraves do metodo residuals com a opcao "standardize=T"
##

zt <- as.timeSeries(residuals(Lfit2, standardize=T))

mrlPlot(zt) # Mean Residual Life para escolher treshold u
# u por volta de 1.5% parece ser um valor adequado
# como o quantil 95% eh 1.62% vamos manter o quantil que eh o padrao
# da funcao gpdFit

# Contagem do numero de excessos apenas para verificar se eh um numero
# razoavelmente alto
sum(zt>quantile(zt, 0.95)) # 97, OK

evtfit <- gpdFit(zt)
op <- par(mfrow=c(2,2))
plot(evtfit, which='all')
par(op)

# E por fim calcula as medidas de risco para os residuos zt
risks <- gpdRiskMeasures(evtfit, prob = 0.99) # Medidas sem intervalo de conf.
tail <- gpdTailPlot(evtfit)
varci <- gpdQPlot(tail)
esci <- gpdSfallPlot(tail)
risktable <- rbind(varci, esci)
dimnames(risktable) <- list(c("$z_{.99}$", "$s_{.99}$"), c("Inf", "Estimativa", "Sup"))
print.xtable(xtable(risktable, caption = "Valores de $z_{.99}$ e $s_{.99}$ encontrados e seus
                    respectivos intervalos de confian\\c ca a 95\\%",
                    label = "tab:tabevtAAPL2", align = c("r", "r", "c", "r")),
             file = "..\\tables\\tabevtAAPL2.tex", sanitize.text.function = function(x) {x})

