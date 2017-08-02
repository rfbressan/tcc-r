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
m_sim <- 500              # Numero de simulacoes a serem feitas para cada ativo/dia a frente. Monte Carlo
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
                                          garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(1, 1),
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
fit_tc <- cgarchfit(t_cspec, data = ret, out.sample = 252, 
                    fit.control = list(eval.se=FALSE),
                    spd.control = list(lower = 0.1, upper = 0.9, type = "mle", kernel = "normal"),
                    cluster = cl)
# Encerra o cluster. Nunca se esquecer!!
stopCluster(cl)
rm(cl)

## Checagens de convergencia do ajuste
if(fit_tc@mfit$convergence != 0) stop("cGARCHfit nao convergiu")

## Metodos sobre cgarchfit
fit2resid <- residuals(fit2)  # residuals() retira os residuos do modelo como um xts.
                              # o mesmo que fit2@model$residuals

fitted2 <- fitted(fit2) # fitted() retira as estimativas da media. resid + fitted = observacao
                        # o mesmo que fit2@model$mu

sigma2 <- sigma(fit2) # sigma() retira os desvios padrao calculados no modelo GARCH
                      # o mesmo que fit2@model$sigma
    
stdresid2 <- fit2resid / sigma2 # residuos padronizados do modelo. (zt). O vetor stdresid em fit2@mfit
                                # não parece bater com os valores acima.

Tin <- dim(Dat)[1]-100

# Filtragem dos dados ?? Parece ser rodar o modelo ajustado com os dados da propria amostra e alem
specx2 <- spec2
for(i in 1:d) specx2@umodel$fixed.pars[[i]] <- as.list(fit2@model$mpars[fit2@model$midx[,i]==1,i])
setfixed(specx2)<-as.list(fit2@model$mpars[fit2@model$midx[, d-1] == 1, d-1])

# including n.old filters based on the full assumptions of the fitted model
filt2a <- cgarchfilter(specx2, data = Dat[1:(Tin+100), ], filter.control  = list(n.old = Tin))
# without using n.old
filt2b <- cgarchfilter(specx2, data = Dat[1:(Tin+100), ])

filt2c<- cgarchfilter(specx2, data = Dat[1:(Tin), ])

options(width = 120)
zz <- file("test3c2.txt", open="wt")
sink(zz)
# What happens is:
# The initial variance used to start the garch iteration is based on the
# full data set provided (T+100) unless explicitly using n.old. Therefore,
# the first few values of the covariance are different for filt2b which used
# an additional 100 points to estimate the initial variance....
print(all.equal(rcov(fit2)[,,1], rcov(filt2a)[,,1]))
print(all.equal(rcov(fit2)[,,1], rcov(filt2b)[,,1]))
# ... as T grow, the impact of the initial values decays and the results is
# the same for methods using either n.old and without
print(all.equal(rcov(fit2)[,,Tin], rcov(filt2a)[,,Tin]))
print(all.equal(rcov(fit2)[,,Tin], rcov(filt2b)[,,Tin]))
sink(type="message")
sink()
close(zz)