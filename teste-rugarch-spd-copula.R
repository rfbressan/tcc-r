###############################################################################
## Teste do funcionamento de copulas com os pacotes rugarch + spd + copula
## rugarch + spd para modelar cada um dos ativos, onde as distribuicoes margi
## nais dos residuos do garch serao uma semi-parametrica
## com a copula pode-se gerar "n x d" simulacoes de retornos 1 passo a frente
## e 10 passos a frente. Então retirar os quantis para estimar VaR e ES
## Extensao multivariavel do modelo proposto por McNeil e Frey 2000
## 
###############################################################################

library(rugarch)
#library(parallel) # rugarch carrega parallel
library(spd)
library(copula)
library(psych) # para o comando pairs.panel se quiser utilizar
library(xts)
library(tidyverse)

setwd("C:\\Users\\rfbre\\Documents\\UDESC\\TCC\\TCC Latex svn\\TCC-R-codes")

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
m_sim <- 5              # Numero de simulacoes a serem feitas para cada ativo/dia a frente. Monte Carlo
wport <- rep(1/d, d)      # Portfolio de pesos iguais
labels <- c("MSFT", "GE") # Nome dos ativos para apresentar nos graficos

# Perfil de tempo do algoritimo
tic <- Sys.time()
################################################################################################################
## Modelo GARCH
################################################################################################################
# Especifica o modelo AR-GARCH univariavel. Mesmo modelo para todos os ativos
uspec <- ugarchspec(variance.model = list(model = "sGARCH",
                                          garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(1, 0),
                                      include.mean = TRUE),
                    distribution.model = "ghyp")

# Faz o fit para todos os modelos
# Esta estimacao parece ser a que leva mais tempo. Existe o parametro cluster na funcao multifit
# que pode ser usada para fazer computacao paralela a partir do pacote parallel.
# Deve-se criar um cluster com a funcao makeCluster e repassa-lo a multifit. Entao ao final do fit
# nao esquecer de fechar o cluster.
cl <- makeCluster(detectCores()) # No meu windows 10 core i7 clusters = 4
gfit.list <- multifit(multispec(replicate(d, uspec)),
                      ret,
                      cluster = cl)

# Finaliza o cluster criado
stopCluster(cl)
rm(cl)

# Testa a convergencia de cada ativo
ex_conv <- function(x) {return(x@fit$convergence)} # Extrai os valores de convergencia
if(sum(unlist(lapply(gfit.list@fit, ex_conv)))) stop("Algum modelo GARCH nao convergiu")

# Coleta os residuos padronizados de cada ativo e guarda em um data.frame por colunas
resid_st <- residuals(gfit.list, standardize = TRUE)

################################################################################################################
## Modelo de Copula
################################################################################################################

# A partir dos residuos padronizados, estimar uma distribuicao semi-parametrica para cada um dos ativos
#spd_fit.list <- apply(resid_st, 2, spdfit, type = "mle", kernelfit = "normal")
spd_fit.list <- as_tibble(resid_st) %>% 
  map(~spdfit(.x, type = "mle", kernelfit = "normal"))

# E entao gerar a distribuicao pseudo-uniforme
# esta eh uma distribuicao uniforme, com valores entre 0 e 1 que mantem a relacao de dependencia 
# original dos dados
## AQUI FOI APLICADA A TEORIA DO VALOR EXTREMO uma vez que esta distribuicao levou em conta os
## valores de cauda modelados com GPD
u_spd <- spd_fit.list %>% 
  map_df(~pspd(.x@data, .x, linear = TRUE))

# Verificar a correlacao entre as series
pairs.panels(u_spd, method = "kendall",
             labels = labels) # Distribuicoes uniformes e correlacao Pearson de 0.5 para MSFT e GE
                                        # Kendall de 0.35

# Com esta distribuicao pseudo-uniforme modelar a Copula
# Vamos estimar uma Copula-t
t_cop <- tCopula(dim = d)
t_cop_fit <- fitCopula(t_cop, u_spd, method = "ml")

# E agora é possível gerar dados aleatórios, a partir da copula
# para residuos padronizados que serao entao embutidos na distribuicao SPD e posteriormente no
# modelo GARCH

# set.seed(999) # Para reprodutibilidade dos dados
# u_t_sim <- array(dim = c(n_ahead, d, m_sim)) # Pre-aloca o array de simulacoes no espaco [0,1]
# for(i in 1:m_sim){
#   u_t_sim[,,i] <- rCopula(n_ahead, t_cop_fit@copula) # Array de dimensao [n x d x m] no espaco [0,1]
# }

## Outra forma de estruturar os dados com tibble e nested data frames
u_tibble <- function(nahead, copfit, labels) {
  # Retorna um tibble com d+1 colunas, n_ahead e d colunas de simulacoes de ut para cada ativo
  tib <- as_tibble(data.frame(cbind(1:nahead, 
                                    rCopula(nahead, copfit@copula))))
  colnames(tib) <- c("n_ahead", labels)
  return(tib)
}

# Aninha os data frames. Cada linha correspondente a 1 m_sim, tera um tibble resultante de uma simulacao
u_t_sim <- tibble(m_sim = 1:m_sim) %>% 
  mutate(roll = m_sim %>% 
           map(~u_tibble(n_ahead, t_cop_fit, labels)))

# Verificar a estrutura de dependencia dos u simulados
# pairs.panels(u_t_sim[,,1], method = "kendall",
#              labels = labels)       # Parece OK

# A partir dos u simulados, obter os z simulados
q_gen <- function(j, spd, u){
  qspd(u[,j], spd[[j]], linear = TRUE)
}

z_t_sim <- array(dim = c(n_ahead, d, m_sim)) # Pre-aloca o array de simulacoes no espaco [-inf, inf]
for(i in 1:m_sim){
  z_t_sim[,,i] <- sapply(1:d, q_gen, spd = spd_fit.list, u = u_t_sim[,,i]) # [n x d x m] espaço [-inf, inf] mas 
                                                                          # com media 0 e var 1 nas colunas
}

# Verificar a estrutura de dependencia dos u simulados
# pairs.panels(z_t_sim, method = "kendall",
#              labels = labels)       # Parece OK


################################################################################################################
## Simulacao
################################################################################################################

# A matriz distfit deve ter ordem [n.sim x m.sim]
# Objeo sim sera uma lista de dimensao igual ao numero de ativos, cada qual com informacoes sobre a simulacao
# realizada.
sim <- lapply(1:d, function(j)
  ugarchsim(gfit.list@fit[[j]], n.sim = n_ahead, m.sim = m_sim, startMethod = "sample",
            custom.dist = list(name = "sample",
                               distfit = z_t_sim[,j,])))

# Vamos retirar para cada ativo os valores de retorno simulados e calcular ao longo dos dias de simulacao
# n_ahead o retorno cumulativo. Como os retornos estao em termos logaritimicos, basta somar diretamente para o 
# retorno acumulado ser o retorno commposto.

# ATENCAO: AGORA OS ATIVOS ESTARAO NA TERCEIRA DIMENSAO DO ARRAY E NAO MAIS NA SEGUNDA. POIS FOI ASSIM QUE 
# A FUNCAO DE SIMULACAO RETORNOU
cumret_sim <- array(dim = c(n_ahead, m_sim, d)) # Pre-aloca o array, COM AS DIMENSOES ALTERADAS
for(i in 1:d){
  cumret_sim[,,i] <- apply(fitted(sim[[i]]), 2, cumsum)
}

# Agora calcula-se os retornos cumulativos do portfolio, para cada uma das simulacoes.
# Terminaremos com uma matriz de dimensao [n_ahead x m_sim] para os retornos do portfolio.
port_cumret <- apply(cumret_sim, c(1,2), `%*%`, wport) # Portfolio cumulative returns

# Grafico com os caminhos das simulacoes
plot(1:n_ahead, port_cumret[,1], type = "l", ylim = c(min(port_cumret), max(port_cumret)))
for(i in 2:m_sim) lines(1:n_ahead, port_cumret[,i], col = i)

elapsed <- Sys.time()-tic
cat("Elapsed time:", elapsed, "\n")

### Remove todas as variaveis da memoria!!
### Fica pesado demais depois da simulacao
rm(list = ls())

################################################################################################################
## Simulacao rolling
################################################################################################################

## Pesquisar sobre o metodo multifilter, com os argumentos n.old e out.sample.
## se este metodo retornar as matrizes sigma, retornos e residuos para cada um dos dias
## entao talvez seja possivel utilizar estas matrizes para fazer a simulacao de monte carlo
## em uma janela movel.
## multiforecast serve para fazer previsao um passo a frente e ir rolando. tem opcao n_ahead

