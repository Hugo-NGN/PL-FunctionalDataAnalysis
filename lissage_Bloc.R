library(tidyverse)
library(fda)
library(stats)
library(factoextra)
library(plotly)
source("./utils/lissage.R")
source("./utils/prepocess.R")
source("./utils/distances_fonctionnelles.R")
source("./utils/derive_fonctionnelle.R")


# ------------------------- Chargement des donnees  ----------------------------


data <- read.csv("./data/var_SV_2018-01-01_00H_nan-forced_depth.csv", sep=";")
data <-preprocess(data, extract_n_data = 300)

sub <- gsub("^X", "", colnames(data))
colnames(data) <-  sub

# -------------------------------- Lissage  ------------------------------------
#definition des variables utiles pour le lissage
l_grille = 10^seq(-3, 5, length.out = 1000)
D <- 30
z <- as.numeric(colnames(data))

#lissage B-spline
fd_obj <- spline_lissage_bloc_quantile(data, l_grille, D, z)

# sauvegarde du lissage
#saveRDS(fd_obj, file = "./data/fdata.rds")

#fd_obj <- readRDS("./data/fdata.rds")
## ------------------------ Affichage courbes lisses ---------------------------
plot(fd_obj[[1]], col = 1, lty = 1, main = "Profils lissés (bloc 47)",
     ylab = "Célérité", xlab = "Profondeurs")

for (i in 2:length(fd_obj)) {
  lines(fd_obj[[i]], col = i, lty = 1)
}


## ----------------------- Affichage courbes derivees --------------------------
fine_grid <- seq(min(as.numeric(colnames(data))),
                 max(as.numeric(colnames(data))),
                 length.out = 10000)

deriv_list = list()
for (i in seq(1, length(fd_obj))){
  deriv_list[[i]] <- eval.fd(fine_grid, fd_obj[[i]], Lfdobj=1)
}

#saveRDS(deriv_list, file = "./data/fdata_deriv.rds")

#deriv_list <- readRDS("./data/fdata_deriv.rds")

plot(fine_grid, deriv_list[[1]], type ="l", col = 1, lty = 1, main = "Dérivés des profils(bloc 47)",
     ylab = "Dérivés de la célérité/profondeur", xlab = "Profondeurs")

for (i in 2:length(fd_obj)) {
  lines(deriv_list[[i]], col = i, lty = 1)
}


# ------------------------- Calcul matrices des distances ----------------------
## ---------------------------- Calcul matrice D0  ----------------------------- 
fine_grid <- seq(min(as.numeric(colnames(data))),
                 max(as.numeric(colnames(data))),
                 length.out = 1000)

### ----------------------- calcul D0 sequentiellement ------------------------- 
#D0_matrix <- calculate_D0_matrix(fd_obj, fine_grid)
# saveRDS(D0_matrix, file="./data/D0_matrix.rds")
#D0_matrix <- readRDS("./data/D0_matrix.rds")

### -------------------- calcul D0 avec parallelisation ------------------------
D0_matrix <- calculate_D0_matrix_parallel(fd_obj, fine_grid)
#saveRDS(D0_matrix, file="./data/D0_matrix_n1000.rds")
#D0_matrix <- readRDS("./data/D0_matrix.rds")

plot_ly(z = ~D0_matrix, type = "surface") %>%
  layout(
    scene = list(
      xaxis = list(title = "Profil"),
      yaxis = list(title = "Profil"),
      zaxis = list(title = "Distance D0")
    )
  )

## ---------------------------- Calcul matrice D1  -----------------------------

fine_grid <- seq(min(as.numeric(colnames(data))),
                 max(as.numeric(colnames(data))),
                 length.out = 1000)
### ----------------------- calcul D1 sequentiellement -------------------------
# D1_matrix <- calculate_D1_matrix(fd_obj, fine_grid)
# saveRDS(D1_matrix, file="./data/D1_matrix.rds")
#D1_matrix <- readRDS("./data/D1_matrix.rds")

### -------------------- calcul D1 avec parallelisation ------------------------
D1_matrix <- calculate_D1_matrix_parallel(fd_obj, fine_grid)
#saveRDS(D1_matrix, file="./data/D1_matrix_n1000.rds")
#D1_matrix <- readRDS("./data/D1_matrix.rds")

plot_ly(z = ~D1_matrix, type = "surface") %>%
  layout(
    scene = list(
      xaxis = list(title = "Profil"),
      yaxis = list(title = "Profil"),
      zaxis = list(title = "Distance D1")
    )
  )

## ---------------------------- Calcul matrice Dp ------------------------------

fine_grid <- seq(min(as.numeric(colnames(data))),
                 max(as.numeric(colnames(data))),
                 length.out = 1000)

omega <- 0.5
### ----------------------- calcul Dp sequentiellement -------------------------
# Dp_matrix <- calculate_Dp_matrix(fd_obj, fine_grid, omega)
# saveRDS(Dp_matrix, file="./data/D0_matrix.rds")
#Dp_matrix <- readRDS("./data/Dp_matrix_omega05.rds")

### -------------------- calcul Dp avec parallelisation ------------------------
Dp_matrix <- calculate_Dp_matrix_parallel(fd_obj, fine_grid, omega)
#saveRDS(Dp_matrix, file="./data/Dp_matrix_omega05_n1000.rds")
#Dp_matrix <- readRDS("./data/Dp_matrix_omega05.rds")

plot_ly(z = ~Dp_matrix, type = "surface") %>%
  layout(
    scene = list(
      xaxis = list(title = "Profil"),
      yaxis = list(title = "Profil"),
      zaxis = list(title = "Distance Dp")
    )
  )
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


silhouette_score <- silhouette(groups, as.dist(D0_matrix))
avg_silhouette_width <- mean(silhouette_score[, 3])

## ------------------------------ K MEANS --------------------------------------
