library(tidyverse)
library(fda)
source("./utils/lissage.R")
source("./utils/prepocess.R")
source("./utils/distances_fonctionnelles2.R")

# ------------------------- Chargement les donnees  ----------------------------


data <- read.csv("./data/var_SV_2018-01-01_00H_nan-forced_depth.csv", sep=";")
data <-preprocess(data)

sub <- gsub("^X", "", colnames(data))
colnames(data) <-  sub

# -------------------------------- Lissage  ------------------------------------
#definition des variables utiles pour le lissage
l_grille = 10^seq(-3, 5, length.out = 1000)
D <- 30
z <- as.numeric(colnames(data))

#lissage B-spline
fd_obj <- spline_lissage_bloc_quantile(data, l_grille, D, z)

## ------------------------ Affichage courbes lisses ---------------------------
plot(fd_obj[[1]], col = 1, lty = 1, main = "Courbes fonctionnelles (bloc 47)", ylab = "Célérité", xlab = "Profondeurs")
for (i in 2:length(fd_obj)) {
  lines(fd_obj[[i]], col = i, lty = 1)
}

# ----------------------------  Calcul matrice D0  -----------------------------  
fine_grid <- seq(min(as.numeric(colnames(data))),
                 max(as.numeric(colnames(data))),
                 length.out = 10000)

D0_matrix <- calculate_D0_matrix(fd_obj, fine_grid)

# ----------------------------  Calcul matrice D1  -----------------------------

fine_grid <- seq(min(as.numeric(colnames(data))),
                 max(as.numeric(colnames(data))),
                 length.out = 10000)

D1_matrix <- calculate_D1_matrix(fd_obj, fine_grid)

# ----------------------------  Calcul matrice Dp  -----------------------------

fine_grid <- seq(min(as.numeric(colnames(data))), max(as.numeric(colnames(data))), length.out = 10000)

omega <- 0.5

Dp_matrix <- calculate_Dp_matrix(fd_obj, fine_grid, omega)

