library(fda)
library(stats)
library(mclust)
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
data <-preprocess(data)

sub <- gsub("^X", "", colnames(data))
colnames(data) <-  sub


## ---------------------------- Visualisation ----------------------------------

x_values <- as.numeric(colnames(data))
# Tracer les courbes
matplot(x_values, t(data), type = "l", lty = 1, col = rainbow(1000), xaxt = "n", ylab ="Célérité", xlab="Profondeur")
title("Courbes des profils")
axis(1, at = x_values, labels = colnames(data), las = 2)


# -------------------------------- Lissage  ------------------------------------
#definition des variables utiles pour le lissage
l_grille = 10^seq(-3, 5, length.out = 1000)
D <- 30
z <- as.numeric(colnames(data))

#lissage B-spline
fd_list <- spline_lissage_bloc_quantile(data, l_grille, D, z)

# sauvegarde du lissage
#saveRDS(fd_list, file = "./data/fdata.rds")

#fd_list <- readRDS("./data/fdata.rds")
## ------------------------ Affichage courbes lisses ---------------------------
plot(fd_list[[1]], col = 1, lty = 1, main = "Profils lissés (bloc 47)",
     ylab = "Célérité", xlab = "Profondeurs")

for (i in 2:length(fd_list)) {
  lines(fd_list[[i]], col = i, lty = 1)
}


## ----------------------- Affichage courbes derivees --------------------------
fine_grid <- seq(min(as.numeric(colnames(data))),
                 max(as.numeric(colnames(data))),
                 length.out = 10000)

deriv_list = list()
for (i in seq(1, length(fd_list))){
  deriv_list[[i]] <- eval.fd(fine_grid, fd_list[[i]], Lfdobj=1)
}

#saveRDS(deriv_list, file = "./data/fdata_deriv.rds")

#deriv_list <- readRDS("./data/fdata_deriv.rds")

plot(fine_grid, deriv_list[[1]], type ="l", col = 1, lty = 1,
     main = "Dérivés des profils (bloc 47)",
     ylab = "Dérivés de la célérité/profondeur",
     xlab = "Profondeurs")

for (i in 2:length(fd_list)) {
  lines(deriv_list[[i]], col = i, lty = 1)
}


# ------------------------- Calcul matrices des distances ----------------------
## ---------------------------- Calcul matrice D0  ----------------------------- 
fine_grid <- seq(min(as.numeric(colnames(data))),
                 max(as.numeric(colnames(data))),
                 length.out = 1000)

### ----------------------- calcul D0 sequentiellement ------------------------- 
#system.time(D0_matrix <- calculate_D0_matrix(fd_list, fine_grid))
# saveRDS(D0_matrix, file="./data/D0_matrix.rds")
#D0_matrix <- readRDS("./data/D0_matrix.rds")

### -------------------- calcul D0 avec parallelisation ------------------------
system.time(D0_matrix <- calculate_D0_matrix_parallel(fd_list, fine_grid))
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
# D1_matrix <- calculate_D1_matrix(fd_list, fine_grid)
# saveRDS(D1_matrix, file="./data/D1_matrix.rds")
#D1_matrix <- readRDS("./data/D1_matrix.rds")

### -------------------- calcul D1 avec parallelisation ------------------------
system.time(D1_matrix <- calculate_D1_matrix_parallel(fd_list, fine_grid))
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
# Dp_matrix <- calculate_Dp_matrix(fd_list, fine_grid, omega)
# saveRDS(Dp_matrix, file="./data/D0_matrix.rds")
#Dp_matrix <- readRDS("./data/Dp_matrix_omega05.rds")

### -------------------- calcul Dp avec parallelisation ------------------------

#system.time(Dp_matrix <- calculate_Dp_matrix_parallel(fd_list, fine_grid, omega, STANDARDIZE = FALSE))
system.time(Dp_matrix_stdz <- calculate_Dp_matrix_parallel(fd_list, fine_grid, omega = 0.5, STANDARDIZE = TRUE))

#saveRDS(Dp_matrix_stdz, file="./data/Dp_matrix_omega05_n1000.rds")
#Dp_matrix <- readRDS("./data/Dp_matrix_omega05.rds")

plot_ly(z = ~Dp_matrix_stdz, type = "surface") %>%
  layout(
    scene = list(
      xaxis = list(title = "Profil"),
      yaxis = list(title = "Profil"),
      zaxis = list(title = "Distance Dp")
    )
  )


# ------------------------------- Clustering -----------------------------------

## ------------------------------ BASELINE -------------------------------------

# Matrice de distance euclidienne des profils
D_eucli_matrix <- as.matrix(dist(data, method = "euclidean"))

k_optimal <- kmeans_optimal_k(D_eucli_matrix)

kmeans_eucli <- kmeans(data, centers = k_optimal, nstart = 25)

# calcul du score de Davies Bouldin
baseline_db_score <- davies.bouldin(data, kmeans_eucli$cluster)

# calcul du coefficient de silhouette
baseline_silhouette_score <- mean(silhouette(kmeans_eucli$cluster, as.dist(D_eucli_matrix))[, 3])


## --------------------------------- CAH ---------------------------------------
### ------------------------------ CAH  D0 -------------------------------------


cah_silhouette_opti_D0 <- cah_optimal_silhouette(D0_matrix, fd_list, method ="complete")

hc_D0 <- cah_silhouette_opti_D0$hc

plot(hc_D0,
     main = "Dendrogramme de la CAH D0",
     sub = "",
     xlab = "",
     hang = -1,
     cex = 0.6)

rect.hclust(hc_D0,
            k = cah_silhouette_opti_D0$k_optimal,
            border = "green")

### ------------------------------ CAH  D1 -------------------------------------


cah_silhouette_opti_D1 <- cah_optimal_silhouette(D1_matrix, fd_list, method ="complete")

hc_D1 <- cah_silhouette_opti_D1$hc

plot(hc_D1,
     main = "Dendrogramme de la CAH D1",
     sub = "",
     xlab = "",
     hang = -1,
     cex = 0.6)

rect.hclust(hc_D1, k = cah_silhouette_opti_D1$k_optimal, border = "green")

### ------------------------------ CAH  Dp -------------------------------------


cah_silhouette_opti_Dp <- cah_optimal_silhouette(Dp_matrix_stdz, fd_list, method = "complete")

hc_Dp <- cah_silhouette_opti_Dp$hc

plot(hc_Dp,
     main = "Dendrogramme de la CAH Dp",
     sub = "",
     xlab = "",
     hang = -1,
     cex = 0.6)

rect.hclust(hc_Dp, k = cah_silhouette_opti_Dp$k_optimal, border = "green")


## ------------------------------ K MEANS --------------------------------------

# choix du k pour maximiser le coeff de Silhouette
k_D0 <- kmeans_optimal_k(D0_matrix)
k_D1 <- kmeans_optimal_k(D1_matrix)
k_Dp <- kmeans_optimal_k(Dp_matrix_stdz)

# clustering avec kmeans
kmeans_D0 <- kmeans_fd(D0_matrix, k_D0, fd_list)
kmeans_D1 <- kmeans_fd(D1_matrix, k_D1, fd_list)
kmeans_Dp <- kmeans_fd(Dp_matrix_stdz, k_Dp, fd_list)


# Calcul et affichage des scores pour D0
silhouette_D0 <- mean(silhouette(kmeans_D0$km_result$cluster, as.dist(D0_matrix))[, 3])
db_score_D0 <- davies.bouldin(as.matrix(D0_matrix), kmeans_D0$km_result$cluster)
cat("Score de Silhouette moyen pour D0:", silhouette_D0, "\n")
cat("Score de Davies-Bouldin pour D0:", db_score_D0, "\n")

# Calcul et affichage des scores pour D1
silhouette_D1 <- mean(silhouette(kmeans_D1$km_result$cluster, as.dist(D1_matrix))[, 3])
db_score_D1 <- davies.bouldin(as.matrix(D1_matrix), kmeans_D1$km_result$cluster)
cat("Score de Silhouette moyen pour D1:", silhouette_D1, "\n")
cat("Score de Davies-Bouldin pour D1:", db_score_D1, "\n")

# Calcul et affichage des scores pour Dp
silhouette_Dp <- mean(silhouette(kmeans_Dp$km_result$cluster, as.dist(Dp_matrix_stdz))[, 3])
db_score_Dp <- davies.bouldin(as.matrix(Dp_matrix_stdz), kmeans_Dp$km_result$cluster)
cat("Score de Silhouette moyen pour Dp:", silhouette_Dp, "\n")
cat("Score de Davies-Bouldin pour Dp:", db_score_Dp, "\n")

## -------------------------------- CAH + KMEANS ------------------------------


hybride_classif_D0 <- cah_kmeans(D0_matrix, fd_list, method="complete")
hybride_classif_D1 <- cah_kmeans(D1_matrix, fd_list, method="complete")
hybride_classif_Dp <- cah_kmeans(Dp_matrix_stdz, fd_list, method="complete")


## ------------------- Comparaison classification ARI --------------------------


ARI_D0 <- compute_ARI(kmeans_D0$km_result$cluster, cah_silhouette_opti_D0$cluster)
ARI_D1 <- compute_ARI(kmeans_D1$km_result$cluster, cah_silhouette_opti_D1$cluster)
ARI_Dp <- compute_ARI(kmeans_Dp$km_result$cluster, cah_silhouette_opti_Dp$cluster)


# ------------------------------ FIN -------------------------------------------
