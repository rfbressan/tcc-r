## Teste do pacote PtProcess
library(PtProcess)
data("Phuket")
Phuket$magnitude <- Phuket$magnitude - 4.95

## vou utilizar minhas distribuicoes de marcas GPD
dgP_mark <- function(x, data, params){
  g <- (1/params[7])*(1+((params[6]*x[,"magnitude"])/params[7]))^(-1-(1/params[6]))
  return(log(g))
}

rgP_mark <- function(ti, data, params){
  y <- rgp(1, params[6], params[7])
  return(list(magnitude = y))
}

TT <- c(0, 1827) #Eh todo o span de tempo
params <- c(0.05, 3.1, 1.3, 0.02, 1.1, 0.5, 1)

x <- mpp(data = Phuket, gif = etas_gif, marks = list(dgP_mark, rgP_mark), params = params, TT = TT,
         gmap = expression(params[1:5]), mmap = expression(params))

## Mapeamento para o otimizador
## Como os parametros do modelo devem ser todos positivos, podemos tirar o log deles e passar
## para o otimizador para que este trabalhe com todos os reais
expmap <- function(y, p){
  y$params <- exp(p)
  return(y)
}
initial <- log(params)
## A funcao negativa logaritimica de verossimilhanca negloglik deve ser minimizada
# Utilizar dois passos. Primeiro optim e depois nlm
z <- optim(initial, neglogLik, object = x, pmap = expmap, control = list(trace = 1, maxit = 100))
initial <- z$par
z <- nlm(neglogLik, initial, object = x, pmap = expmap, print.level = 2, iterlim = 500, typsize = initial)
x0 <- expmap(x, z$estimate)

lambdag <- etas_gif(Phuket, c(0:1827), x0$params)
plot(c(0:1827), lambdag, type = "l")
