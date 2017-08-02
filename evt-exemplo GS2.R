###############################################################################
## Exemplo do capítulo sobre EVT de medidas condicionais de risco utilizando
## o método processo pontual auto-excitavel.
## Dados da GS entre 05/05/1999 a 04/10/2016
###############################################################################
library(xts)
library(PerformanceAnalytics)
library(xtable)
library(PtProcess)
library(fExtremes)
library(gPdtest)
library(plotly)
library(boot)

setwd("C:\\Users\\Rafael\\Documents\\UDESC\\TCC\\TCC Latex svn\\TCC-R-codes")

prices <- as.xts(read.zoo(read.csv("evt-exemplo GS.csv"), format = "%Y-%m-%d", FUN = as.Date))
returns <- na.omit(Return.calculate(prices, method = "log"))
losses <- -1*returns
losses <- losses[,"Adj.Close"]
names(losses) <- "losses"
#############################################################################################################
## Pode-se usar os pacotes PtProcess (mais geral) ou SAPP (função etasap)
## Utilizando o pacote PtProcess
## Definindo o modelo de processo pontual marcado - mpp
## Nosso modelo sera um ETAS e a distribuição das marcas uma distribuicao de Pareto com qsi > 0
#############################################################################################################
# Quero um DF com tempos e valores de perdas acima de u=1.5%
u <- quantile(losses, 0.9)
Lu <- losses[which(losses[,"losses"]>u),"losses"]
Tj <- as.numeric(index(Lu)-index(losses[1,]))
Lj <- coredata(Lu) # Eh o valor original da perda, dado que ela seja maior que u
Luj <- Lj - u # Eh a perda em excesso, que sera a entrada para as densidades lambdag e dgP_mark
Datadf <- data.frame(Luj, time = Tj/365) # tem que ajustar o nome da coluna de Luj
colnames(Datadf)[1] <- "magnitude"
#Ht <- data.frame(time=as.numeric(index(losses)-index(losses[1,])), magnitude=coredata(losses[,"Adj.Close"]))# tem que ajustar o nome da coluna
#colnames(Ht)[2] <- "magnitude"
#evalpts <- Ht[,"time"] # valores onde sera avaliada a funcao intensidade
plot(Datadf$time, Datadf$magnitude, type="h")
#############################################################################################################
## Vamos primeiro estimar os parametros da GPD
gpd <- gpdFit(as.timeSeries(losses[,"losses"]), u=u, type = "mle")
xi <- gpd@fit$par.ests[1]
beta <- gpd@fit$par.ests[2]
#############################################################################################################
## Agora que ja temos o DF de dados e os pontos de avaliacao, temos de especificar o modelo MPP
## A funcao intensidade de base (lambdag) sabemos que será uma ETAS
## A densidade das marcas eh uma GPD com parametros qsi e beta
## Precisamos especificar as funcoes densidade e simulacao das marcas, dgP_mark e rgP_mark
## o vetor de parametros contem TODOS os parametros das duas densidades, sendo os primeiros
## pertencentes ao modelo ETAS (5), apos vem os parametros da GPD (2), Logo
## params = mu, A, alpha, c, p, xi, beta

#############################################################################################################
dgP_mark <- function(x, data, params){
  g <- (1/beta)*(1+((xi*x[,"magnitude"])/beta))^(-1-(1/xi))
  return(log(g))
}

rgP_mark <- function(ti, data, params){
  y <- rgp(1, beta, xi)
  return(list(magnitude = y))
}
#############################################################################################################
## Mapeamento dos parametros
gmap <- expression(params[1:5])
#mmap <- expression(params)
## Mapeamento para o otimizador
nlmmap <- function(y, p){
  y$params <- p
  return(y)
}
expmap <- function(y, p){
  y$params <- exp(p)
  return(y)
}
## O objeto MPP propriamente dito 

  Dados <- Datadf
  parameters <- c(0.05, 3.1, 1.5, 0.1, 1.1) # Valores iniciais dos parametros para lambdag
  TT <- c(Dados$time[1], Dados$time[nrow(Dados)])
  sepp <- mpp(data = Dados, gif = etas_gif, marks = list(dgP_mark, rgP_mark),
              params = parameters, gmap = gmap, mmap = NULL, TT = TT)
  
  ## Como os parametros do modelo devem ser todos positivos, podemos tirar o log deles e passar
  ## para o otimizador para que este trabalhe com todos os reais
  initial <- log(parameters)
  ## A funcao negativa logaritimica de verossimilhanca negloglik deve ser minimizada
  # Utilizar dois passos. Primeiro optim e depois nlm
  sepphat <- optim(initial, neglogLik, gr = NULL, object = sepp, pmap = expmap,
                   control = list(trace = 0, maxit = 100))
  if(sepphat$convergence == 0){
    # Então a convergencia foi bem sucedida
    seppmodel <- expmap(sepp, sepphat$par)
  }else{
    # Agora pode-se usar os parametros estimados para serem valores iniciais de uma nova otimizacao
    initial <- sepphat$par
    sepphat <- nlm(neglogLik, initial, object = sepp, pmap = expmap,
                   hessian = TRUE, print.level = 0, iterlim = 500, typsize = initial)
    seppmodel <- expmap(sepp, sepphat$estimate)
  }
  theta <- with(seppmodel, c(params[1]/365, params[2]/365, params[3],
                             params[4], params[5], xi, beta))
  names(theta) <- c("tau", "psi", "delta", "gamma", "rho", "xi", "beta")
  
  ## Estimando os erros padroes dos parametros de lambdag
  H <- sepphat$hessian
  CovM <- solve(H) # Matriz de covariancias eh a inversa da Hessiana
  # SE e^params = sqrt(se^2*e^2*params)
  setheta <- sqrt(exp(2*sepphat$estimate)*diag(CovM))
  setheta <- c(setheta[1]/365, setheta[2]/365, setheta[3:5])

## Dado a convergencia do modelo podemos plotar a funcao intensidade de base
pts <- seq(0, as.numeric(index(losses[nrow(losses),])-index(losses[1,]))/365,
           length.out = nrow(losses))
lambdag <- etas_gif(Dados, evalpts = pts, params = seppmodel$params)/365
exc <- (losses[,"losses"]-u)-((losses[,"losses"]-u)<0)*(losses[,"losses"]-u)

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
ts <- cbind(losses, pts, exc, lambdag, vares)
names(ts) <- c("losses", "indextimes", "excess", "intensity", "VaR", "ES")

op <- par(mfrow=c(2,1))
plot(ts[,"intensity"], xlab = "", ylab = "Intensidade",
     main = "Processo Pontual de Auto-Excitação", major.format="%m/%Y")
plot(ts[,"excess"], type="h", xlab = "Data", ylab = "Perdas em excesso",
     main = "", major.format="%m/%Y")
par(op)

## Grafico do VaR e ES sobre as perdas efetivas

plot(ts[-1,"losses"], type="h", xlab = "Data", ylab = "Perdas",
     main = "Medidas Condicionais de Risco",
     ylim = c(0, 0.22), major.format="%m/%Y") # A primeira perda nao tem VaR para comparar
lines(ts[-length(ts[,"VaR"]),"VaR"], col="red", lwd = 2) # O ultimo VaR so sera comparado com a proxima perda que ainda nao ocorreu
lines(ts[-length(ts[,"ES"]),"ES"], col="blue", lwd = 2)  # O ultimo ES so sera comparado com a proxima perda que ainda nao ocorreu
abline(h=urisk[c(2,3)], col=c("red", "blue"), lty = c(2,4), lwd = 2)
legend("topright", legend = c("VaR incondicional", "ES incondicional"), col = c("red", "blue"),
       lty = c(2, 4), lwd = 2, cex = 1.2)

## Analise da adequacao do fit
plot(residuals(seppmodel), xlab = "Evento Número", ylab = "Tempo Transformado", pty = "s")
abline(0, 1)


