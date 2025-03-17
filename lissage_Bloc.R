library(fda)
library(stats)
library(plotly)
library(cluster)
library(tidyverse)
library(factoextra)
source("./utils/lissage.R")
source("./utils/preprocess.R")
source("./utils/clustering.R")
source("./utils/derive_fonctionnelle.R")
source("./utils/distances_fonctionnelles.R")

# ------------------------- Chargement des donnees  ----------------------------


data <- read.csv("./data/var_SV_2018-01-01_00H_nan-forced_depth.csv", sep=";")
data <-preprocess(data, extract_n_data = 200)

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
#system.time(D0_matrix <- calculate_D0_matrix(fd_obj, fine_grid))
# saveRDS(D0_matrix, file="./data/D0_matrix.rds")
#D0_matrix <- readRDS("./data/D0_matrix.rds")

### -------------------- calcul D0 avec parallelisation ------------------------
system.time(D0_matrix <- calculate_D0_matrix_parallel(fd_obj, fine_grid))
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
system.time(D1_matrix <- calculate_D1_matrix_parallel(fd_obj, fine_grid))
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
system.time(Dp_matrix <- calculate_Dp_matrix_parallel(fd_obj, fine_grid, omega, STANDARDIZE = TRUE))
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
### ------------------------------ CAH  D0 -------------------------------------
cah_silhouette_opti_D0 <- cah_optimal_silhouette(D0_matrix, fd_obj)

hc_D0 <- cah_silhouette_opti_D0$hc

plot(hc_D0,
     main = "Dendrogramme de la CAH D0",
     sub = "",
     xlab = "",
     hang = -1,
     cex = 0.6)

rect.hclust(hc_D0, k = cah_silhouette_opti_D0$k_optimal, border = "green")

### ------------------------------ CAH  D1 -------------------------------------
cah_silhouette_opti_D1 <- cah_optimal_silhouette(D1_matrix, fd_obj)

hc_D1 <- cah_silhouette_opti_D1$hc

plot(hc_D1,
     main = "Dendrogramme de la CAH D1",
     sub = "",
     xlab = "",
     hang = -1,
     cex = 0.6)

rect.hclust(hc_D1, k = cah_silhouette_opti_D1$k_optimal, border = "green")

### ------------------------------ CAH  Dp -------------------------------------
cah_silhouette_opti_Dp <- cah_optimal_silhouette(Dp_matrix, fd_obj)

hc_Dp <- cah_silhouette_opti_Dp$hc

plot(hc_Dp,
     main = "Dendrogramme de la CAH Dp",
     sub = "",
     xlab = "",
     hang = -1,
     cex = 0.6)

rect.hclust(hc_Dp, k = cah_silhouette_opti_Dp$k_optimal, border = "green")












## ------------------------------ K MEANS --------------------------------------
# Initialiser un vecteur pour stocker les coefficients de silhouette moyens
silhouette_scores <- numeric(9)

# Boucle sur les valeurs de k de 2 à 10
for (k in 2:10) {
  # Appliquer k-means
  km_result <- kmeans(D1_matrix, centers = k, nstart = 25)
  
  # Calculer le coefficient de silhouette
  silhouette_score <- silhouette(km_result$cluster, as.dist(D0_matrix))
  avg_silhouette_width <- mean(silhouette_score[, 3])
  
  # Stocker le coefficient de silhouette moyen
  silhouette_scores[k - 1] <- avg_silhouette_width
}

# Tracer la courbe du coefficient de silhouette en fonction de k
plot(2:10, silhouette_scores, type = "b",
     xlab = "Nombre de clusters (k)",
     ylab = "Coefficient de silhouette moyen",
     main = "Coefficient de silhouette en fonction de k")

# Appliquer k-means avec k = 3
km_result <- kmeans(D0_matrix, centers = 2, nstart = 25)

# Définir les couleurs pour chaque groupe
group_colors <- c("red", "blue")

# Visualiser les groupes
plot(fd_obj[[1]], col = group_colors[km_result$cluster[1]], lty = 1,
     main = "Courbes fonctionnelles par groupe",
     ylab = "Célérité", xlab = "Profondeurs")

for (i in 2:length(fd_obj)) {
  lines(fd_obj[[i]], col = group_colors[km_result$cluster[i]], lty = 1)
}

# Ajouter une légende
legend("topright", legend = paste("Groupe", 1:2), col = group_colors, lty = 1)

# ------------------------------ FIN -------------------------------------------

