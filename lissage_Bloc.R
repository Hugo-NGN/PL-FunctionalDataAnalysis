library(tidyverse)
library(fda)
source("./utils/lissage.R")
source("./utils/prepocess.R")

#charge les donnees
data <- read.csv("./data/var_SV_2018-01-01_00H_nan-forced_depth.csv", sep=";")
data <-preprocess(data)

sub <- gsub("^X", "", colnames(data))
colnames(data) <-  sub

#definition des variables utiles pour le lissage
l_grille = 10^seq(-3, 5, length.out = 1000)
D <- 30
z <- as.numeric(colnames(data))

#lissage B-spline
data_lisse <- spline_lissage_bloc_quantile(data, l_grille, D, z)

#Resultats du lissage: objets fonctionnels
fd_obj <- data_lisse$fd_obj

#Plot des courbes lissees
plot(fd_obj[[1]], col = 1, lty = 1, main = "Courbes fonctionnelles (bloc 47)", ylab = "Célérité", xlab = "Profondeurs")
for (i in 2:length(fd_obj)) {
  lines(fd_obj[[i]], col = i, lty = 1)
}

#







