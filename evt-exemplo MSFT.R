###############################################################################
## Exemplo do capítulo sobre EVT de medidas condicionais de risco utilizando
## o método processo pontual auto-excitavel.
## Dados da MSFT entre 03/01/1997 a 29/12/2006
###############################################################################
library(xts)
library(PerformanceAnalytics)
library(xtable)
library(PtProcess)
library(SAPP)
library(fExtremes) # Para utilizar as funcoes [dr]gpd
library(gPdtest)

setwd("C:\\Users\\Rafael\\Documents\\UDESC\\TCC\\TCC Latex svn\\R codes")

prices <- as.xts(read.zoo(read.csv("evt-exemplo MSFT2.csv"), format = "%Y-%m-%d", FUN = as.Date))
returns <- 100*na.omit(Return.calculate(prices, method = "log"))
losses <- -1*returns

#############################################################################################################
## Pode-se usar os pacotes PtProcess (mais geral) ou SAPP (função etasap)
## Utilizando o pacote PtProcess
## Definindo o modelo de processo pontual marcado - mpp
## Nosso modelo sera um ETAS e a distribuição das marcas uma distribuicao de Pareto com qsi > 0


# Primeiro vamos testar alguns dados e conhecer o pacote
data("Ogata")
# E um DF com magnitudes acima de um treshold e os tempos destas observações. Equivale a sequencia
# (Tj, Lj-u)
plot(Ogata$time, Ogata$magnitude + 3.5, type="h")
t <- seq(0,800,0.1)
x <- etas_gif(Ogata, t, c(0.02, 1, 0.2, 1, 1.5))
plot(t, x, type="l")
# x vira um vetor de numeros com a intensidade base de um processo ETAS

#############################################################################################################
# Quero um DF com tempos e valores de perdas acima de u=1.5%
# E outro DF com todo o historico do processo Ht
u <- 1.5
Lu <- losses[which(losses[,"Adj.Close"]>u),"Adj.Close"]
Tj <- as.numeric(index(Lu)-index(losses[1,]))
Lj <- coredata(Lu) # Eh o valor original da perda, dado que ela seja maior que u
Luj <- Lj - u # Eh a perda em excesso, que sera a entrada para as densidades lambdag e dgpd_mark
Datadf <- data.frame(Luj, time = Tj) # tem que ajustar o nome da coluna de Luj
colnames(Datadf)[1] <- "magnitude"
#Ht <- data.frame(time=as.numeric(index(losses)-index(losses[1,])), magnitude=coredata(losses[,"Adj.Close"]))# tem que ajustar o nome da coluna
#colnames(Ht)[2] <- "magnitude"
#evalpts <- Ht[,"time"] # valores onde sera avaliada a funcao intensidade
plot(Datadf$time, Datadf$magnitude, type="h")
#############################################################################################################
## Agora que ja temos o DF de dados e os pontos de avaliacao, temos de especificar o modelo MPP
## A funcao intensidade de base (lambdag) sabemos que será uma ETAS
## A densidade das marcas eh uma GPD com parametros qsi e beta
## Precisamos especificar as funcoes densidade e simulacao das marcas, dgpd_mark e rgpd_mark
## o vetor de parametros contem TODOS os parametros das duas densidades, sendo os primeiros
## pertencentes ao modelo ETAS (5), apos vem os parametros da GPD (2), Logo
## params = mu, A, alpha, c, p, xi, beta
dgpd_mark <- function(x, data, params){
# x deve conter um DF com tempos e magnitudes. Pode ser apenas uma parte de Datadf
# data deve ser um DF com o historico do processo (Ht = Datadf)
# params eh um vetor com os parametros das distribuicoes lambdag e gpd
  y <- dgpd(x[,"magnitude"], xi = params[6], mu = 0, beta = params[7], log = TRUE)
  return(as.numeric(y))
}

rgpd_mark <- function(ti, data, params){
  y <- rgpd(1, xi = params[6], mu = 0, beta = params[7])
  return(list(magnitude = y))
}

#############################################################################################################
## Como aparentemente as funcoes acima nao funcionaram
## Vamos criar outras
dgP_mark <- function(x, data, params){
  g <- (1/params[7])*(1+((params[6]*x[,"magnitude"])/params[7]))^(-1-(1/params[6]))
  return(log(g))
}

rgP_mark <- function(ti, data, params){
  y <- rgp(1, params[6], params[7])
  return(list(magnitude = y))
}
#############################################################################################################
## Mapeamento dos parametros
gmap <- expression(params[1:5])
mmap <- expression(params)

## O objeto MPP propriamente dito
Dados <- Datadf
#Dados$magnitude <- Dados$magnitude - 4.95 # Apenas quando for Phuket
parameters <- c(0.05, 1, 0.5, 0.5, 1.1, 0.5, 1) # Valores iniciais dos parametros
TT <- c(Dados$time[1], Dados$time[nrow(Dados)])
sepp <- mpp(data = Dados, gif = etas_gif, marks = list(dgP_mark, rgP_mark), params = parameters,
            gmap = gmap, mmap = mmap, TT = TT)

## Mapeamento para o otimizador
nlmmap <- function(y, p){
  y$params <- p
  return(y)
}

expmap <- function(y, p){
  y$params <- exp(p)
  return(y)
}
## Como os parametros do modelo devem ser todos positivos, podemos tirar o log deles e passar
## para o otimizador para que este trabalhe com todos os reais
initial <- parameters
## A funcao negativa logaritimica de verossimilhanca negloglik deve ser minimizada
# Utilizar dois passos. Primeiro optim e depois nlm
sepphat <- optim(initial, neglogLik, gr = NULL, object = sepp, pmap = NULL,
                 method = "L-BFGS-B", lower=c(0,0,0,0.01,1,0,0),
                 control = list(trace = 1, maxit = 500))
if(sepphat$convergence == 0){
  # Então a convergencia foi bem sucedida
  seppmodel <- expmap(sepp, sepphat$par)
}else{
  # Agora pode-se usar os parametros estimados para serem valores iniciais de uma nova otimizacao
  initial <- sepphat$par
  sepphat <- nlm(neglogLik, initial, object = sepp, pmap = expmap,
                 print.level = 2, iterlim = 500, typsize = initial)
  seppmodel <- expmap(sepp, sepphat$estimate)
}

## Dado a convergencia do modelo podemos plotar a funcao intensidade de base
lambdag <- etas_gif(Dados, evalpts = c(TT[1]:TT[2]), params = seppmodel$params)
op <- par(mfrow=c(2,1))
plot(TT[1]:TT[2], lambdag, type = "l", xlab = "Período")
plot(Dados$time, Dados$magnitude, type="h", xlab = "Período")
par(op)

## Analise da adequacao do fit
plot(residuals(seppmodel), xlab = "Evento Número", ylab = "Tempo Transformado", pty = "s")

# Vetor theta com os parametros conforme modelo na monografia
theta <- with(seppmodel, c(params[1], params[2]*params[4]^params[5], params[3], params[4], params[5]-1,
                           params[6], params[7]))
