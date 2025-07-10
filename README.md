rm(list = ls())
setwd("C:/Users/diretorio") #diretorio

dados <- read.table("PZT1_00000001.SySHM", header = FALSE, skip = 1) #documento
FREQ <- dados$V1
IMPED <- dados$V2

model <- lm(IMPED ~ FREQ)
p <- length(coef(model))
n <- length(IMPED)

##############################################
# boxplot por ciclos de 30
##############################################

block_size <- 30
n_blocks <- floor(length(FREQ) / block_size)

freq_medianas <- numeric(n_blocks)
imped_medianas <- numeric(n_blocks)

blocos_imped_orig <- list()
blocos_imped_clean <- list()
outliers_por_bloco <- numeric(n_blocks)

outliers_pos_blocos <- list() 

for (i in 1:n_blocks) {
  idx_ini <- (i - 1) * block_size + 1
  idx_fim <- i * block_size
  
  bloco_freq <- FREQ[idx_ini:idx_fim]
  bloco_imped <- IMPED[idx_ini:idx_fim]
  
  blocos_imped_orig[[i]] <- bloco_imped
  
  # Limites boxplot
  Q1 <- quantile(bloco_imped, 0.25)
  Q3 <- quantile(bloco_imped, 0.75)
  IQR <- Q3 - Q1
  lower_limit <- Q1 - 1.5 * IQR
  upper_limit <- Q3 + 1.5 * IQR
  idx_outliers <- which(bloco_imped < lower_limit | bloco_imped > upper_limit)
  
  # Posições reais no vetor original (índices completos)
  pos_outliers_reais <- idx_outliers + (i - 1) * block_size
  outliers_pos_blocos[[i]] <- pos_outliers_reais
  
  outliers_por_bloco[i] <- length(idx_outliers)
  
  if (length(idx_outliers) > 0) {
    bloco_imped_clean <- bloco_imped[-idx_outliers]
  } else {
    bloco_imped_clean <- bloco_imped
  }
  
  blocos_imped_clean[[i]] <- bloco_imped_clean
  
  freq_medianas[i] <- mean(bloco_freq)
  imped_medianas[i] <- median(bloco_imped_clean)
}


##############################################
# Outliers DFFITS
##############################################

inf_measures <- influence.measures(model)
infmat <- inf_measures$infmat

dffits_values <- abs(infmat[, 3])
cutoff_dffits <- 2 * sqrt(p / n)
out_dffits <- which(dffits_values > cutoff_dffits)

##############################################
# Cook's Distance
##############################################

cooksD <- abs(infmat[, 5])
cutoff_cooks <- 1
out_cooks <- which(cooksD > cutoff_cooks)

##############################################
# Leverage (hii)
##############################################

hii <- abs(infmat[, 6])
cutoff_hii <- 2 * p / n
out_hii <- which(hii > cutoff_hii)


##############################################
# gráficos comparativos
##############################################

# original
par(mfrow = c(1,1))
plot(FREQ, IMPED,
     main = "FREQ vs IMPED (Original)",
     xlab = "Frequência", ylab = "Impedância",
     pch = 20)

# boxplot blocos de 30 (medianas
par(mfrow = c(1,1))
plot(freq_medianas, imped_medianas,
     main = "FREQ vs IMPED ( Gráfico impedâncias pós outliers)",
     xlab = "Frequência (Média por Bloco)",
     ylab = "Impedância (Mediana por Bloco)",
     type = "b",
     pch = 20)

# boxplots por bloco - originais
par(mfrow = c(1,1))
boxplot(blocos_imped_orig,
        main = "Boxplots Originais por Bloco",
        xlab = "Bloco", ylab = "Impedância")

# boxplots por bloco - limpos
par(mfrow = c(1,1))
boxplot(blocos_imped_clean,
        main = "Boxplots Sem Outliers por Bloco",
        xlab = "Bloco", ylab = "Impedância")




# DFFITS
par(mfrow = c(1,1))
if (length(out_dffits) > 0) {
  plot(FREQ[-out_dffits], IMPED[-out_dffits],
       main = "FREQ vs IMPED (Sem DFFITS Outliers)",
       xlab = "Frequência", ylab = "Impedância",
       pch = 20)
} else {
  plot(FREQ, IMPED,
       main = "FREQ vs IMPED (DFFITS → Nenhum Outlier)",
       xlab = "Frequência", ylab = "Impedância",
       pch = 20)
}

# Cook's Distance
par(mfrow = c(1,1))
if (length(out_cooks) > 0) {
  plot(FREQ[-out_cooks], IMPED[-out_cooks],
       main = "FREQ vs IMPED (Sem Cook's Outliers)",
       xlab = "Frequência", ylab = "Impedância",
       pch = 20)
} else {
  plot(FREQ, IMPED,
       main = "FREQ vs IMPED (Cook → Nenhum Outlier)",
       xlab = "Frequência", ylab = "Impedância",
       pch = 20)
}

# Leverage
par(mfrow = c(1,1))
if (length(out_hii) > 0) {
  plot(FREQ[-out_hii], IMPED[-out_hii],
       main = "FREQ vs IMPED (Sem Leverage Outliers)",
       xlab = "Frequência", ylab = "Impedância",
       pch = 20)
} else {
  plot(FREQ, IMPED,
       main = "FREQ vs IMPED (Leverage → Nenhum Outlier)",
       xlab = "Frequência", ylab = "Impedância",
       pch = 20)
}

##############################################
# DFIITS com linha de corte
##############################################


par(mfrow = c(1,1))
plot(dffits_values, main = "DFFITS", ylab = "DFFITS", pch = 20)
abline(h = cutoff_dffits, col = "red", lty = 2)
abline(h = -cutoff_dffits, col = "red", lty = 2)
legend("topright", legend = paste0("|DFFITS| >", round(cutoff_dffits, 2)),
       col = "red", lty = 2, bty = "n")
points(out_dffits, dffits_values[out_dffits], col = "blue", pch = 19)


##############################################
# cook com linha de corte
##############################################

par(mfrow = c(1,1))
plot(cooksD, main = "Cook's Distance", ylab = "Cook's D", pch = 20)
abline(h = cutoff_cooks, col = "red", lty = 2)
legend("topright", legend = paste0("Cook's D > ", cutoff_cooks),
       col = "red", lty = 2, bty = "n")
points(out_cooks, cooksD[out_cooks], col = "blue", pch = 19)


##############################################
# leverage com linha de corte
##############################################

par(mfrow = c(1,1))
plot(hii, main = "Leverage (hii)", ylab = "hii", pch = 20)
abline(h = cutoff_hii, col = "red", lty = 2)
legend("bottomright", legend = paste0("hii >", round(cutoff_hii, 4)),
       col = "red", lty = 2, bty = "n")
points(out_hii, hii[out_hii], col = "blue", pch = 19)








##############################################
# resumo total
##############################################


cat("\n### RESUMO ###\n")
cat("Outliers DFFITS:", length(out_dffits), "\n")
if (length(out_dffits) > 0) {
  cat("Índices DFFITS:", out_dffits, "\n")
}

cat("Outliers Cook:", length(out_cooks), "\n")
if (length(out_cooks) > 0) {
  cat("Índices Cook:", out_cooks, "\n")
}

cat("Outliers Leverage:", length(out_hii), "\n")
if (length(out_hii) > 0) {
  cat("Índices Leverage:", out_hii, "\n")
}

cat("Total Outliers Boxplot Blocos:", sum(outliers_por_bloco), "\n")
cat("Outliers por Bloco:", outliers_por_bloco, "\n")

#print as posições dos outliers
outliers_filtrados <- outliers_pos_blocos[sapply(outliers_pos_blocos, length) > 0]
todos_outliers <- unlist(outliers_filtrados)
cat("Outliers por boxplot:", paste(todos_outliers, collapse = ", "), "\n")
