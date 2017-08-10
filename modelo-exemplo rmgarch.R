###############################################################################
## Teste do funcionamento de copulas com os pacotes rmgarch + spd + copula
## rmgarch + spd para modelar cada um dos ativos, onde as distribuicoes margi
## nais dos residuos do garch serao uma semi-parametrica.
## rmgarch tambem modela a copula t, mas nao as arquimedianas, o pacote copula 
## pode gerar estas a partir dos dados das margens obtidos com rmgarch.
## com a copula pode-se gerar "n x d" simulacoes de retornos 1 passo a frente
## e 10 passos a frente. Então retirar os quantis para estimar VaR e ES
## Extensao multivariavel do modelo proposto por McNeil e Frey 2000
## 
###############################################################################

library(rmgarch) # Carrega rugarch e parallel
library(copula)
library(xts)
library(psych) # para o comando pairs.panel se quiser utilizar


#setwd("C:\\Users\\rfbre\\Documents\\UDESC\\TCC\\TCC Latex svn\\TCC-R-codes")

################################################################################################################
## Aquisicao de dados e calculo do xts com os retornos
################################################################################################################
## Utilizar dados de MSFT e GE para o teste
#p_msft <- as.tbl(read.csv("evt-exemplo MSFT2.csv")[, c(1, 7)]) # Apenas datas e Fechamento. Transforma em tibble
#p_ge <- as.tbl(read.csv("evt-exemplo GE.csv")[, c(1, 7)]) # Apenas datas e Fechamento. Transforma em tibble

# Apenas datas e Fechamento.
p_msft <- as.xts(read.zoo(read.csv("evt-exemplo MSFT2.csv")[, c(1, 7)], format = "%Y-%m-%d", FUN = as.Date))
p_ge <- as.xts(read.zoo(read.csv("evt-exemplo GE.csv")[, c(1, 7)], format = "%Y-%m-%d", FUN = as.Date))

# Calcula os retornos logaritmicos em termos percentuais
r_msft <- na.omit(100*diff(log(p_msft)))
r_ge <- na.omit(100*diff(log(p_ge)))
ret <- merge(r_msft, r_ge, join = "inner") # Retornos dos 2 ativos em um xts apenas

################################################################################################################
## Definicoes Gerais
################################################################################################################
d <- ncol(ret)            # Dimensao da copula = numero de ativos
n_ahead <- 10             # Numero de dias a frente a serem simulados
m_sim <- 5                # Numero de simulacoes a serem feitas para cada ativo/dia a frente. Monte Carlo
n_roll <- 12              # Numero de dias a serem rolados na janela movel
# Ordem ARMA(p,q) e GARCH(m,s)
p <- 1
q <- 1
m <- 1
s <- 1
wport <- rep(1/d, d)      # Portfolio de pesos iguais
labels <- c("MSFT", "GE") # Nome dos ativos para apresentar nos graficos

# Perfil de tempo do algoritimo
tic <- Sys.time()
################################################################################################################
## Modelo GARCH
################################################################################################################
# Especifica o modelo AR-GARCH univariavel. Mesmo modelo para todos os ativos
# ghyp porque esta estah no MDA de uma GPD, a qual desejamos para utilizar EVT no modelo.
# poderiam ser inovacoes normais e os parametros estimados do ar-garch nao seriam enviesados.
# GARCH padrao # ARMA(1, 1) para a media com constante # Assumir generalizada hiperbolica, que eh MDA de GPD
uspec <- ugarchspec(variance.model = list(model = "sGARCH",
                                          garchOrder = c(m, s)),
                    mean.model = list(armaOrder = c(p, q),
                                      include.mean = TRUE),
                    distribution.model = "ghyp")

# Modelo da copula. Aqui se pode utilizar uma copula t ou normal. Pode-se tambem especificar que as
# margens seguem distribuicoes semi-parametricas (spd) e a implementacao eh feita com GPD nas caudas
t_cspec <- cgarchspec(uspec = multispec( replicate(d, uspec)),
                      distribution.model = list(copula = "mvt", 
                                                method = "Kendall",
                                                transformation = "spd"))

# Adequa os dados ao modelo. Deixa-se 252 observacoes de fora para rodar simulacoes depois.
# eval.se para estimar o erro padrao dos parametros
# spd.control para passar parametros da distribuicao semi-parametrica
# Abrir um cluster para programacao paralela
if(Sys.info()["sysname"]=="Windows"){
  cl <- makePSOCKcluster(detectCores())
}else{
  cl <- makeForkCluster(detectCores())
}
fit_tc <- cgarchfit(t_cspec, data = ret, out.sample = n_roll, 
                    fit.control = list(eval.se=FALSE),
                    spd.control = list(lower = 0.1, upper = 0.9, type = "mle", kernel = "normal"),
                    cluster = cl)
# Encerra o cluster. Nunca se esquecer!!
stopCluster(cl)
rm(cl)

## Checagens de convergencia do ajuste
if(fit_tc@mfit$convergence != 0) stop("cGARCHfit nao convergiu")


# Rolling Simulation ------------------------------------------------------
t_orig <- nrow(ret) - n_roll
simMu <- simS <- filtMu <- filtS <- matrix(NA, ncol = d, nrow = n_roll)
simC <- filtC <- array(NA, dim = c(d, d, n_roll))
colSd <- function(x) apply(x, 2, "sd")
specx <- t_cspec

for(i in 1:d) specx@umodel$fixed.pars[[i]] <-  as.list(fit_tc@model$mpars[fit_tc@model$midx[,i]==1,i])
setfixed(specx) <- as.list(fit_tc@model$mpars[fit_tc@model$midx[,d+1]==1, d+1])

# Filtra todos os dados para resgatar as matrizes prereturns, presigma e preresidual
filtro_c <- cgarchfilter(specx, ret, 
                       filter.control = list(n.old = t_orig),
                       spd.control = list(lower = 0.1, upper = 0.9, type = "mle", kernel = "normal"))

presigmas <- matrix(tail(sigma(filtro_c), n_roll+s), ncol = d)
preresiduals <- matrix(tail(residuals(filtro_c), n_roll+s), ncol = d)
prereturns <- matrix(tail(coredata(ret), n_roll+s), ncol = d) # Retornos observados

for(i in 1:n_roll){
  # if(i==1){
  #   presigma = matrix(tail(sigma(fit_tc), s), ncol = d)
  #   # arma = c(p,q) therefore need p lags
  #   prereturns = matrix(unlist(Dat[(t_orig-p):t_orig, ]), ncol = d, nrow = p)
  #   preresiduals = matrix(tail(residuals(fit_tc),2), ncol = d, nrow = max(c(q,m)))
  #   
  #   tmp = cgarchfilter(specx, Dat[1:(t_orig+1), ], filter.control = list(n.old = t_orig))
  #   filtMu[i,] = tail(fitted(tmp), 1)
  #   filtS[i,] = tail(sigma(tmp), 1)
  #   filtC[,,i] = last(rcov(tmp))[,,1]
  # } else{
  #   presigma = matrix(tail(sigma(tmp), 2), ncol = 3)
  #   # arma = c(2,1) therefore need 2 lags
  #   prereturns = matrix(unlist(Dat[(t_orig+i-2):(t_orig+i-1), ]), ncol = 3, nrow = 2)
  #   preresiduals = matrix(tail(residuals(tmp),2), ncol = 3, nrow = 2)
  #   
  #   tmp = cgarchfilter(specx, Dat[1:(t_orig+i), ], filter.control = list(n.old = t_orig))			
  #   filtMu[i,] = tail(fitted(tmp), 1)
  #   filtS[i,] = tail(sigma(tmp), 1)
  #   filtC[,,i] = last(rcov(tmp))[,,1]
  # }
  sim3 = cgarchsim(fit_tc, n.sim = 1, m.sim = 10000, startMethod = "sample", prereturns = prereturns,
                   presigma = presigma, preresiduals = preresiduals)
  simx = t(sapply(sim3@msim$simX, FUN = function(x) x[1,]))
  simMu[i,] = colMeans(simx)
  simS[i,] = colSd(simx)
  simC[,,i] = cov(simx)
  print(i)
  
  # CHECK:
  # X[t+1] = mu[t+1] + e[t+1]
  # sim3@msim$simX[[i]][1,] - (sim3@msim$simZ[,,i]*sqrt(diag(sim3@msim$simH[[i]][,,1])))
  # is equal to filtMu[i,]
  if(i < 3 ){
    print(all.equal(sim3@msim$simX[[2]][1,] - (sim3@msim$simZ[,,2]*sqrt(diag(sim3@msim$simH[[2]][,,1]))), 
                    filtMu[i,]))
    print(all.equal(sim3@msim$simX[[10000]][1,] - (sim3@msim$simZ[,,10000]*sqrt(diag(sim3@msim$simH[[10000]][,,1]))), 
                    filtMu[i,]))
  }
}
