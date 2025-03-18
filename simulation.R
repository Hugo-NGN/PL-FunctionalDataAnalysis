# Simulation name : fda - dataset : growth 
library(fda)
library(cluster)
library(tidyverse)

source("./utils/clustering.R")
source("./utils/distances_fonctionnelles.R")

# -------------------------- Données Simulation --------------------------------


data_growth <- fda::growth


# ---------------------------- VISUALISATION -----------------------------------
## ----------------------------- Courbes M -------------------------------------
matplot(data_growth$age,
        data_growth$hgtm,
        type="l",
        lty=1,
        xlab="Age",
        ylab="Height")

title("Courbes de croissance (M)")

## ----------------------------- Courbes F -------------------------------------
matplot(data_growth$age,
        data_growth$hgtf,
        type="l",
        lty=1,
        xlab="Age",
        ylab="Height")

title("Courbes de croissance (F)")

## ---------------------------- Comparaison M/F --------------------------------

data_combined <- cbind(data_growth$hgtm, data_growth$hgtf)

colors <- c(rep("blue", ncol(data_growth$hgtm)), rep("red", ncol(data_growth$hgtf)))

matplot(data_growth$age, data_combined, type="l", lty=1, col=colors, xlab="Age", ylab="Height")

title("Courbes de croissance par genre")

legend("topright", legend=c("M", "F"), col=c("blue", "red"), lty=1)



# ---------------------------- LISSAGE B-SPLINE --------------------------------
# Base de spline
basis <- create.bspline.basis(rangeval=range(data_growth$age), nbasis=7)

# lissage par sexe
fd_m <- smooth.basis(data_growth$age, data_growth$hgtm, basis)$fd
fd_f <- smooth.basis(data_growth$age, data_growth$hgtf, basis)$fd


# fusion en une liste d'objet fd
fd_list <- list()
for (i in seq_len(ncol(data_growth$hgtm))) {
  fd_list[[i]] <- smooth.basis(data_growth$age, data_growth$hgtm[, i], basis)$fd
}
for (i in seq_len(ncol(data_growth$hgtf))) {
  fd_list[[length(fd_list) + 1]] <- smooth.basis(data_growth$age, data_growth$hgtf[, i], basis)$fd
}

# Visualisation 
plot(fd_m, col="blue", lwd=1, main="Courbes de croissance lissées (M)")
plot(fd_f, col="red", lwd=1, main="Courbes de croissance lissées (F)")

# Comparaison des courbes lissées
#plot(fd_m, col=colors, lwd=1, main="Courbes de croissance lissées par genre")
#plot(fd_f, add=TRUE, col=colors, lwd=2)
#legend("topright", legend=c("M", "F"), col=c("blue", "red"), lty=1)

# ---------------------------- DERIVEES DES COURBES ---------------------------

# Calcul des dérivées
fd_m_prime <- eval.fd(fdobj = fd_m,evalarg = data_growth$age,  Lfdobj=1)
fd_f_prime <- eval.fd(fdobj = fd_f, evalarg = data_growth$age, Lfdobj=1)

# Visualisation
fd_m_prime %>% matplot(type="l",
                       main = "Dérivées des courbes de croissance (M)",
                       xlab="Age",
                       ylab="Dérivée")

fd_f_prime %>% matplot(type="l",
                       main = "Dérivées des courbes de croissance (F)",
                       xlab="Age",
                       ylab="Dérivée")


# ---------------------------- MATRICE DES DISTANCES ---------------------------
fine_grid <- seq(min(data_growth$age), max(data_growth$age), length.out = 1000)

D0_matrix <- calculate_D0_matrix_parallel(fd_list, fine_grid=fine_grid)
D1_matrix <- calculate_D1_matrix_parallel(fd_list, fine_grid=fine_grid)
Dp_matrix <- calculate_Dp_matrix_parallel(fd_list, fine_grid=fine_grid, omega= 0.5)



# ---------------------------- CLUSTERING -------------------------------------
## ---------------------------- BASELINE ---------------------------------------
data <- rbind(t(data_growth$hgtm), t(data_growth$hgtf))

# Matrice de distance euclidienne des profils
D_eucli_matrix <- as.matrix(dist(data, method = "euclidean"))


k_optimal <- kmeans_optimal_k(as.dist(D_eucli_matrix))

kmeans_eucli <- kmeans(D_eucli_matrix, centers = k_optimal, nstart = 25)

# calcul du score de Davies Bouldin
baseline_db_score <- davies.bouldin(data, kmeans_eucli$cluster)

# calcul du coefficient de silhouette
baseline_silhouette_score <- mean(silhouette(kmeans_eucli$cluster, as.dist(D_eucli_matrix))[, 3])

## ----------------------- K-MEANS fonctionnel ---------------------------------

### ----------------------------- D0 -------------------------------------------
k_optimal_D0 <- kmeans_optimal_k(D0_matrix)
kmeans_D0 <- kmeans(D0_matrix, centers = k_optimal_D0, nstart = 25)
db_score_D0_kmean <- davies.bouldin(data, kmeans_D0$cluster)
silhouette_score_D0_kmean <- mean(silhouette(kmeans_D0$cluster, as.dist(D0_matrix))[, 3])

### ----------------------------- D1 -------------------------------------------

k_optimal_D1 <- kmeans_optimal_k(D1_matrix)
kmeans_D1 <- kmeans(D1_matrix, centers = k_optimal_D1, nstart = 25)
db_score_D1_kmean <- davies.bouldin(data, kmeans_D1$cluster)
silhouette_score_D1_kmean <- mean(silhouette(kmeans_D1$cluster, as.dist(D1_matrix))[, 3])

### ----------------------------- Dp -------------------------------------------

k_optimal_Dp <- kmeans_optimal_k(Dp_matrix)
kmeans_Dp <- kmeans(Dp_matrix, centers = 2, nstart = 25)
db_score_Dp_kmean <- davies.bouldin(data, kmeans_Dp$cluster)
silhouette_score_Dp_kmean <- mean(silhouette(kmeans_Dp$cluster, as.dist(Dp_matrix))[, 3])

## ----------------------- CAH fonctionnelle -----------------------------------
### ----------------------------- D0 -------------------------------------------
cah_silhouette_opti_D0 <- cah_optimal_silhouette(D0_matrix, fd_list, method ="complete")
#cah_silhouette_opti_D0_wardD <- cah_optimal_silhouette(D0_matrix, fd_list, method ="ward.D")
#cah_silhouette_opti_D0_wardD2 <- cah_optimal_silhouette(D0_matrix, fd_list, method ="ward.D2")


hc_D0 <- hclust(as.dist(D0_matrix), method = "complete")

plot(hc_D0,
     main = "Dendrogramme de la CAH D0",
     sub = "",
     xlab = "",
     hang = -1,
     cex = 0.6)

rect.hclust(hc_D0,
            k = 2,
            border = "green")

### ----------------------------- D1 -------------------------------------------
cah_silhouette_opti_D1 <- cah_optimal_silhouette(D1_matrix, fd_list, method ="complete")
#cah_silhouette_opti_D1_wardD <- cah_optimal_silhouette(D1_matrix, fd_list, method ="ward.D")
#cah_silhouette_opti_D1_wardD2 <- cah_optimal_silhouette(D1_matrix, fd_list, method ="ward.D2")

hc_D1 <- hclust(as.dist(D1_matrix), method = "complete")

plot(hc_D1,
     main = "Dendrogramme de la CAH D1",
     sub = "",
     xlab = "",
     hang = -1,
     cex = 0.6)

rect.hclust(hc_D1, k = 2, border = "green")

### ----------------------------- Dp -------------------------------------------
cah_silhouette_opti_Dp <- cah_optimal_silhouette(Dp_matrix, fd_list, method = "complete", force_k=2)

hc_Dp <- hclust(as.dist(Dp_matrix), method = "complete")

plot(hc_Dp,
     main = "Dendrogramme de la CAH Dp",
     sub = "",
     xlab = "",
     hang = -1,
     cex = 0.6)

rect.hclust(hc_Dp, k = 2, border = "green")

## ---------------------- CAH + K-MEANS fonctionnel ----------------------------

hybride_classif_D0 <- cah_kmeans(D0_matrix, fd_list, cut_tree = 20)
hybride_classif_D1 <- cah_kmeans(D1_matrix, fd_list, cut_tree = 20)
hybride_classif_Dp <- cah_kmeans(Dp_matrix, fd_list,cut_tree = 20, kmeans_k = 2)

# ---------------------------- RESULTATS ---------------------------------------
# Vrais labels (sexe: H/F)
true_labels <- c(rep("M", ncol(data_growth$hgtm)), rep("F", ncol(data_growth$hgtf)))

## ---------------------------- BASELINE ---------------------------------------
ari_baseline <- adjustedRandIndex(true_labels, kmeans_eucli$cluster)
cat("ARI Baseline:", ari_baseline, "\n")

## ---------------------------- K-MEANS fonctionnel ----------------------------
ari_D0_kmean <- adjustedRandIndex(true_labels, kmeans_D0$cluster)
cat("ARI K-means D0:", ari_D0_kmean, "\n")

ari_D1_kmean <- adjustedRandIndex(true_labels, kmeans_D1$cluster)
cat("ARI K-means D1:", ari_D1_kmean, "\n")

ari_Dp_kmean <- adjustedRandIndex(true_labels, kmeans_Dp$cluster)
cat("ARI K-means Dp:", ari_Dp_kmean, "\n")

## ---------------------------- CAH fonctionnelle ------------------------------
# Extraire les clusters de la CAH
get_clusters_from_hclust <- function(hc, k) {
  return(cutree(hc, k))
}

clusters_D0_cah <- get_clusters_from_hclust(hc_D0, 2)
ari_D0_cah <- adjustedRandIndex(true_labels, clusters_D0_cah)
cat("ARI CAH D0:", ari_D0_cah, "\n")

clusters_D1_cah <- get_clusters_from_hclust(hc_D1, 2)
ari_D1_cah <- adjustedRandIndex(true_labels, clusters_D1_cah)
cat("ARI CAH D1:", ari_D1_cah, "\n")

clusters_Dp_cah <- get_clusters_from_hclust(hc_Dp, 2)
ari_Dp_cah <- adjustedRandIndex(true_labels, clusters_Dp_cah)
cat("ARI CAH Dp:", ari_Dp_cah, "\n")

## ---------------------- CAH + K-MEANS fonctionnel ----------------------------

ari_D0_hybride <- adjustedRandIndex(true_labels, hybride_classif_D0$km_result$cluster)
cat("ARI CAH + K-means D0:", ari_D0_hybride, "\n")

ari_D1_hybride <- adjustedRandIndex(true_labels, hybride_classif_D1$km_result$cluster)
cat("ARI CAH + K-means D1:", ari_D1_hybride, "\n")

ari_Dp_hybride <- adjustedRandIndex(true_labels, hybride_classif_Dp$km_result$cluster)
cat("ARI CAH + K-means Dp:", ari_Dp_hybride, "\n")









