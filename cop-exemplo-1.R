
# Exemplo 1 Capitulo Copulas ----------------------------------------------
library(xtable)

x <- rnorm(10)
y <- rnorm(10)

z <- exp(x)
t <- exp(2*y)

mean(x)
mean(y)
cov(x, y)
cor(x, y)
var(x)
var(y)
jpeg("figcopexe1.jpeg", width = 830, height = 830)
op <- par(mfrow = c(2,1))
plot(x,y, pch = 19, main = "(a)")
plot(z,t, pch = 19, main = "(b)")
par(op)
dev.off()

#i <- seq_along(x)
tab1.df <- rbind(x, y, z, t)
row.names(tab1.df) <- c("X", "Y", "Z", "T")
tab1 <- xtable(tab1.df, caption = "Dados do exemplo \\ref{exe:copexe1}.",
                digits = 2,
                label = "tab:copexe1",
                auto = TRUE)
print.xtable(tab1, 
             file = "tabcopexe1.tex",
             caption.placement = "top",
             table.placement = "ht")

tab2.df <- rbind(rank(x), rank(y))
tab3.df <- rbind(rank(z), rank(t))
all.equal(tab2.df, tab3.df)
row.names(tab2.df) <- c("R", "S")
tab2 <- xtable(tab2.df, caption = "Postos das observações de X e Y do exemplo \\ref{exe:copexe1}.",
               digits = 2,
               label = "tab:copexe2",
               auto = TRUE)
print.xtable(tab2, 
             file = "tabcopexe2.tex",
             caption.placement = "top",
             table.placement = "ht")
jpeg("figcopexe2.jpeg", width = 830, height = 830)
op <- par(mfrow = c(2,1))
plot(rank(x), rank(y), pch = 19, main = "(a)", xlab = "R", ylab = "S")
plot(rank(z), rank(t), pch = 19, main = "(b)", xlab = "R*", ylab = "S*")
par(op)
dev.off()
