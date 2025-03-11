library(tidyverse)
library(fda)
library(stats)
library(factoextra)
source("./utils/lissage.R")
source("./utils/prepocess.R")
source("./utils/distances_fonctionnelles.R")

# ------------------------- Chargement des donnees  ----------------------------


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

saveRDS(fd_obj, file = "./data/fdata.rds")
## ------------------------ Affichage courbes lisses ---------------------------
plot(fd_obj[[1]], col = 1, lty = 1, main = "Courbes fonctionnelles (bloc 47)",
     ylab = "Célérité", xlab = "Profondeurs")

for (i in 2:length(fd_obj)) {
  lines(fd_obj[[i]], col = i, lty = 1)
}

# ------------------------- Calcul matrices des distances ----------------------
## ---------------------------- Calcul matrice D0  ----------------------------- 
fine_grid <- seq(min(as.numeric(colnames(data))),
                 max(as.numeric(colnames(data))),
                 length.out = 10000)

D0_matrix <- calculate_D0_matrix(fd_obj, fine_grid)

## ---------------------------- Calcul matrice D1  -----------------------------

fine_grid <- seq(min(as.numeric(colnames(data))),
                 max(as.numeric(colnames(data))),
                 length.out = 10000)

D1_matrix <- calculate_D1_matrix(fd_obj, fine_grid)

## ---------------------------- Calcul matrice Dp ------------------------------

fine_grid <- seq(min(as.numeric(colnames(data))),
                 max(as.numeric(colnames(data))),
                 length.out = 10000)

omega <- 0.5

Dp_matrix <- calculate_Dp_matrix(fd_obj, fine_grid, omega)


# ------------------------------- Clustering -----------------------------------
## --------------------------------- CAH ---------------------------------------
hc <- hclust(as.dist(D0_matrix), method = "complete")

groups <- cutree(hc, k = 3)

# Visualiser le dendrogramme avec les groupes
fviz_dend(hc, k = 3,
                cex = 0.6,
                hang = -1,
                rect = TRUE,
                main = "Dendrogramme de la CAH")
print(groups)

group_colors <- c("red", "blue", 'green')

plot(fd_obj[[1]], col = group_colors[groups[1]], lty = 1,
     main = "Courbes fonctionnelles par groupe",
     ylab = "Célérité", xlab = "Profondeurs")

for (i in 2:length(fd_obj)) {
  lines(fd_obj[[i]], col = group_colors[groups[i]], lty = 1)
}

legend("topright", legend = paste("Groupe", 1:3), col = group_colors, lty = 1)

## ------------------------------ K MEANS --------------------------------------
