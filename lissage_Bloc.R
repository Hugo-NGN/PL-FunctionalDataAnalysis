library(tidyverse)
library(fda)
source("./utils/lissage.R")

data <- read.csv("./data/var_SV_2018-01-01_00H_nan-forced_depth.csv", sep=";")

taille_courbes <- apply(data, MARGIN = 1, function(x) which(is.na(x))[1]) - 1

data <- data[taille_courbes == 48, ]

data <- data[1:100, ]

data <- data[,-1]
data <- data[, -48]

data <- as.data.frame(lapply(data, as.numeric))



l_grille = 10^seq(-3, 5, length.out = 1000)

D <- 30

sub <- gsub("^X", "", colnames(data))
colnames(data) <- sub

z <- as.numeric(sub)

data_lisse <- spline_lissage_bloc_quantile(data, l_grille, D, z)


fd_obj <- data_lisse$fd_obj

plot(fd_obj[[1]], col = 1, lty = 1, main = "Courbes fonctionnelles (bloc 47)", ylab = "Célérité", xlab = "Profondeurs")

for (i in 2:length(fd_obj)) {
  lines(fd_obj[[i]], col = i, lty = 1)
}