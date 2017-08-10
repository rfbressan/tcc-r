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
# simMu <- simS <- filtMu <- filtS <- matrix(NA, ncol = d, nrow = n_roll)
# simC <- filtC <- array(NA, dim = c(d, d, n_roll))
# colSd <- function(x) apply(x, 2, "sd")
specx <- t_cspec

for(i in 1:d) specx@umodel$fixed.pars[[i]] <-  as.list(fit_tc@model$mpars[fit_tc@model$midx[,i]==1,i])
setfixed(specx) <- as.list(fit_tc@model$mpars[fit_tc@model$midx[,d+1]==1, d+1])

# Filtra todos os dados para resgatar as matrizes prereturns, presigma e preresidual
filtro_c <- cgarchfilter(specx, ret, 
                       filter.control = list(n.old = t_orig),
                       spd.control = list(lower = 0.1, upper = 0.9, type = "mle", kernel = "normal"))
# Matrizes de condicoes iniciais
# Devem conter o numero de observacoes a serem simuladas + a quantidade necessaria para as CIs
presigmas <- matrix(tail(sigma(filtro_c), n_roll+max(c(m, s, p, q))), ncol = d)
preresiduals <- matrix(tail(residuals(filtro_c), n_roll+max(c(m, s, p, q))), ncol = d)
prereturns <- matrix(tail(coredata(ret), n_roll+max(c(m, s, p, q))), ncol = d) # Retornos observados

# pre-aloca os arrays com as simulacoes 1 passo a frente e n_ahead passos a frente
simx1 <- array(dim = c(m_sim, d, n_roll))
simx_ahead <- array(dim = c(m_sim, d, n_roll))

for(i in 1:n_roll){
  sim_roll <-  cgarchsim(fit_tc, n.sim = n_ahead, m.sim = m_sim, startMethod = "sample",
                   prereturns = prereturns[i:(max(c(m, s, p, q))+i-1), ,drop = FALSE],
                   presigma = presigmas[i:(max(c(m, s, p, q))+i-1), ,drop = FALSE], 
                   preresiduals = preresiduals[i:(max(c(m, s, p, q))+i-1), ,drop = FALSE])
  simx1[,,i] <-  t(sapply(sim_roll@msim$simX, FUN = function(x) x[1,]))
  # simx_ahead deve ser o retorno CUMULATIVO do dia 1 ate n_ahead.
  simx_ahead[,,i] <-  t(sapply(sim_roll@msim$simX, FUN = function(x) apply(x, 2, sum)))
  # CHECK:
  # X[t+1] = mu[t+1] + e[t+1]
  # sim3@msim$simX[[i]][1,] - (sim3@msim$simZ[,,i]*sqrt(diag(sim3@msim$simH[[i]][,,1])))
  # is equal to filtMu[i,]
  # if(i < 3 ){
  #   print(all.equal(sim3@msim$simX[[2]][1,] - (sim3@msim$simZ[,,2]*sqrt(diag(sim3@msim$simH[[2]][,,1]))), 
  #                   filtMu[i,]))
  #   print(all.equal(sim3@msim$simX[[10000]][1,] - (sim3@msim$simZ[,,10000]*sqrt(diag(sim3@msim$simH[[10000]][,,1]))), 
  #                   filtMu[i,]))
  
}
# Arrays simx possuem as seguintes dimensoes [m_sim, d, n_roll]
# Grafico com os caminhos das simulacoes
plot(1:n_roll, simx1[1,1,], type = "l", ylim = c(min(simx1[,1,]), max(simx1[,1,])), col = "blue")
for(i in 2:m_sim) lines(1:n_roll, simx1[i,1,], col = i)
lines(1:n_roll, coredata(ret)[(t_orig+1):nrow(ret), 1], col = "black")


# Composicao do portfolio -------------------------------------------------
# array de dimensao [m_sim, n_roll] com os retornos simulados de 1 passo a frente
port_simx1 <- apply(simx1, c(1,3), `%*%`, wport)
port_simx_ahead <- apply(simx_ahead, c(1,3), `%*%`, wport)

