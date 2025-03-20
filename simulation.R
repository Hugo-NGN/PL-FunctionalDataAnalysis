library(fda)
library(mclust)
library(quantmod)
source("./utils/clustering.R")
source("./utils/distances_fonctionnelles.R")

#------------------------------ Paramètres -------------------------------------

n <- 100  # Nombre de trajectoires
T <- 1    # Intervalle de temps
dt <- 0.01  # Pas de temps
time <- seq(0, T, by = dt)
mu <- 1  # Drift
sigma <- 0.2  # Volatilité

#------------------------------ Simulation -------------------------------------

## ------------------------ Mouvement Brownien Géométrique ---------------------
# Générer 100 mouvements browniens géométriques
set.seed(123)  # Pour la reproductibilité
gbm <- matrix(nrow = length(time), ncol = n)

for (i in 1:n) {
  # Générer un mouvement brownien standard
  brownian_path <- cumsum(rnorm(length(time) - 1, mean = 0, sd = sqrt(dt)))
  brownian_path <- c(0, brownian_path)
  
  # Transformer en mouvement brownien géométrique
  gbm_path <- exp((mu - 0.5 * sigma^2) * time + sigma * brownian_path) - 1
  
  # Stocker la trajectoire
  gbm[, i] <- gbm_path
}

## ---------------------------- Pont Brownien ----------------------------------
# Générer 100 ponts browniens
set.seed(123)  # Pour la reproductibilité
p_brownien <- matrix(nrow = length(time), ncol = n)

for (i in 1:n) {
  # Générer un mouvement brownien standard
  brownian_path <- cumsum(rnorm(length(time) - 1, mean = 0, sd = sqrt(dt)))
  brownian_path <- c(0, brownian_path)
  
  # Transformer en pont brownien
  brownian_bridge <- brownian_path - time * brownian_path[length(time)]
  
  # Stocker la trajectoire
  p_brownien[, i] <- brownian_bridge
}

## -------------------- concaténation des données ------------------------------
data <- cbind(gbm, p_brownien)
true_labels <- c(rep("gbm", ncol(gbm)), rep("p", ncol(p_brownien)))

## ---------------------- Visualisation des trajectoires -----------------------
plot(time, gbm[, 1], type = "l", col = rgb(0, 0, 1, alpha = 0.1),
     ylim = range(gbm, p_brownien, na.rm = TRUE),
     main = "Trajectoires de Mouvements Browniens Géométriques et Ponts Browniens",
     xlab = "Temps", ylab = "Valeur")

for (i in 2:n) {
  lines(time, gbm[, i], col = rgb(0, 0, 1, alpha = 0.1), lwd=1)
}

for (i in 1:n) {
  lines(time, p_brownien[, i], col = rgb(1, 0, 0, alpha = 0.1), lwd=1)
}

legend("topright", legend = c("Mouvement Brownien Géométrique", "Pont Brownien"),
       col = c(rgb(0, 0, 1, alpha = 0.5), rgb(1, 0, 0, alpha = 0.5)), lty = 1, lwd = 3)

# -------------------------------- Lissage -------------------------------------

basis <- create.bspline.basis(rangeval = range(time), nbasis = 10)

fd_gbm <- list()
for (i in 1:ncol(gbm)) {
  fd_gbm[[i]] <- smooth.basis(time, gbm[, i], basis)$fd
}

fd_p <- list()
for (i in 1:ncol(p_brownien)) {
  fd_p[[i]] <- smooth.basis(time, p_brownien[, i], basis)$fd
}

fd_list <- list()
for (i in 1:ncol(data)) {
  fd_list[[i]] <- smooth.basis(time, data[, i], basis)$fd
}

## ------------------------ Visualisation des courbes lissées -------------------

plot(fd_list[[1]], col = rgb(0, 0, 1, alpha = 0.1), ylab = "Valeur",
     main = "Objets Fonctionnels",
     ylim = range(-2,2)
)

for (i in 2:length(fd_list)) {
  lines(fd_list[[i]], col = rgb(0, 0, 1, alpha = 0.1))
}

for (i in 1:length(fd_p)) {
  lines(fd_p[[i]], col = rgb(1, 0, 0, alpha = 0.1))
}

legend("topright", legend = c("Mouvement Brownien Géométrique", "Pont Brownien"),
       col = c(rgb(0, 0, 1, alpha = 0.5), rgb(1, 0, 0, alpha = 0.5)), lty = 1, lwd = 2)

## --------------------- Calcul et visualisation des dérivées ------------------

# Calcul des dérivées
deriv_gbm <- lapply(fd_gbm, deriv.fd)
deriv_p <- lapply(fd_p, deriv.fd)

# Visualisation des dérivées
plot(deriv_gbm[[1]], col = rgb(0, 0, 1, alpha = 0.1), ylab = "Dérivée",
     main = "Dérivées des Trajectoires",
     ylim = range(unlist(lapply(deriv_gbm, eval.fd, time)),
                  unlist(lapply(deriv_p, eval.fd, time))))

for (i in 2:length(deriv_gbm)) {
  lines(deriv_gbm[[i]], col = rgb(0, 0, 1, alpha = 0.1))
}

for (i in 1:length(deriv_p)) {
  lines(deriv_p[[i]], col = rgb(1, 0, 0, alpha = 0.1))
}

legend("topright", legend = c("Dérivée Mouvement Brownien Géométrique", "Dérivée Pont Brownien"),
       col = c(rgb(0, 0, 1, alpha = 0.5), rgb(1, 0, 0, alpha = 0.5)), lty = 1, lwd = 2)

# ------------------------------- CLASSIFICATION -------------------------------

## ------------------------------ Baseline model -------------------------------

D_eucli_matrix <- as.matrix(dist(t(data), method = "euclidean"))

k_optimal <- kmeans_optimal_k(as.dist(D_eucli_matrix))
kmeans_eucli <- kmeans(D_eucli_matrix, centers = k_optimal, nstart = 25)

kmean_baseline_db_score <- davies.bouldin(t(data), kmeans_eucli$cluster)
kmean_baseline_silhouette_score <- mean(silhouette(kmeans_eucli$cluster, as.dist(D_eucli_matrix))[, 3])
cat("Score de Silhouette moyen pour le modèle de base:", kmean_baseline_silhouette_score, "\n")
cat("Score de Davies-Bouldin pour le modèle de base:", kmean_baseline_db_score, "\n")
cat("ARI Baseline:", adjustedRandIndex(true_labels, kmeans_eucli$cluster), "\n")

hc_baseline <- cah_optimal_silhouette(D_eucli_matrix, fd_list, method = "ward.D2", SIMULATION = TRUE)$hc
ari_hc_baseline <- adjustedRandIndex(true_labels, cutree(hc_baseline, 2))
cat("ARI CAH Baseline:", ari_hc_baseline, "\n")


hybride_baseline <- cah_kmeans(D_eucli_matrix, fd_list, cut_tree = 20, SIMULATION = TRUE)
ari_hybride_baseline <- adjustedRandIndex(true_labels, hybride_baseline$km_result$cluster)
cat("ARI CAH + K-means Baseline:", ari_hybride_baseline, "\n")

## ------------------------- Matrices de distances -----------------------------
fine_grid <- time

# Rq: les fonctions avec parrellelisation doivent avoir comme argument des objet 
# nommés "fd_list" et "fine_grid" spécifiquement

D0_matrix <- calculate_D0_matrix_parallel(fd_list, fine_grid)
D1_matrix <- calculate_D1_matrix_parallel(fd_list, fine_grid)
Dp_matrix <- calculate_Dp_matrix_parallel(fd_list, fine_grid, omega = 0.5, STANDARDIZE = TRUE)

## ------------------------- Clustering hiérarchique ---------------------------
cah_silhouette_opti_D0 <- cah_optimal_silhouette(D0_matrix, fd_list, method = "complete", SIMULATION = TRUE)
hc_D0 <-  hclust(as.dist(D0_matrix), method = "complete")

cah_silhouette_opti_D1 <- cah_optimal_silhouette(D1_matrix, fd_list, method = "complete", SIMULATION = TRUE)
hc_D1 <- hclust(as.dist(D1_matrix), method = "complete")

cah_silhouette_opti_Dp <- cah_optimal_silhouette(Dp_matrix, fd_list, method = "complete", SIMULATION = TRUE)
hc_Dp <- hclust(as.dist(Dp_matrix), method = "complete")

## ------------------------------ K-MEANS --------------------------------------
k_optimal_D0 <- kmeans_optimal_k(D0_matrix)
k_optimal_D1 <- kmeans_optimal_k(D1_matrix)
k_optimal_Dp <- kmeans_optimal_k(Dp_matrix)

kmeans_D0 <- kmeans(D0_matrix, centers = k_optimal_D0, nstart = 25)
kmeans_D1 <- kmeans(D1_matrix, centers = k_optimal_D1, nstart = 25)
kmeans_Dp <- kmeans(Dp_matrix, centers = k_optimal_Dp, nstart = 25)


# Calcul et affichage des scores pour D0
silhouette_D0 <- mean(silhouette(kmeans_D0$cluster, as.dist(D0_matrix))[, 3])
db_score_D0 <- davies.bouldin(as.matrix(D0_matrix), kmeans_D0$cluster)
cat("Score de Silhouette moyen pour D0:", silhouette_D0, "\n")
cat("Score de Davies-Bouldin pour D0:", db_score_D0, "\n")

# Calcul et affichage des scores pour D1
silhouette_D1 <- mean(silhouette(kmeans_D1$cluster, as.dist(D1_matrix))[, 3])
db_score_D1 <- davies.bouldin(as.matrix(D1_matrix), kmeans_D1$cluster)
cat("Score de Silhouette moyen pour D1:", silhouette_D1, "\n")
cat("Score de Davies-Bouldin pour D1:", db_score_D1, "\n")

# Calcul et affichage des scores pour Dp
silhouette_Dp <- mean(silhouette(kmeans_Dp$cluster, as.dist(Dp_matrix))[, 3])
db_score_Dp <- davies.bouldin(as.matrix(Dp_matrix), kmeans_Dp$cluster)
cat("Score de Silhouette moyen pour Dp:", silhouette_Dp, "\n")
cat("Score de Davies-Bouldin pour Dp:", db_score_Dp, "\n")


## ------------------------- classification hybride ----------------------------
hybride_classif_D0 <- cah_kmeans(D0_matrix, fd_list, cut_tree = 20, SIMULATION = TRUE)
hybride_classif_D1 <- cah_kmeans(D1_matrix, fd_list, cut_tree = 20, SIMULATION = TRUE)
hybride_classif_Dp <- cah_kmeans(Dp_matrix, fd_list,cut_tree = 20, kmeans_k = 2, SIMULATION = TRUE)

# ------------------------------- Résultats -------------------------------------
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

ari_D0_kmean <- adjustedRandIndex(true_labels, kmeans_D0$cluster)
cat("ARI K-means D0:", ari_D0_kmean, "\n")

ari_D1_kmean <- adjustedRandIndex(true_labels, kmeans_D1$cluster)
cat("ARI K-means D1:", ari_D1_kmean, "\n")

ari_Dp_kmean <- adjustedRandIndex(true_labels, kmeans_Dp$cluster)
cat("ARI K-means Dp:", ari_Dp_kmean, "\n")

ari_D0_hybride <- adjustedRandIndex(true_labels, hybride_classif_D0$km_result$cluster)
cat("ARI CAH + K-means D0:", ari_D0_hybride, "\n")

ari_D1_hybride <- adjustedRandIndex(true_labels, hybride_classif_D1$km_result$cluster)
cat("ARI CAH + K-means D1:", ari_D1_hybride, "\n")

ari_Dp_hybride <- adjustedRandIndex(true_labels, hybride_classif_Dp$km_result$cluster)
cat("ARI CAH + K-means Dp:", ari_Dp_hybride, "\n")






# ----------------------------- SCORE MAXIMUMS ---------------------------------

numeric_labels <- as.numeric(factor(true_labels))


silhouette_avg <- mean(silhouette(numeric_labels, D_eucli_matrix)[, 3])

# Calcul du score de Davies-Bouldin
db_score <- davies.bouldin(t(data), numeric_labels)

cat("Score de Silhouette moyen:", silhouette_avg, "\n")
cat("Score de Davies-Bouldin:", db_score, "\n")


# --------------------------------- FIN ----------------------------------------
