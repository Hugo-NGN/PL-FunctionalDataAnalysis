library(cluster)
library(factoextra)



# --------------------------------- CAH ----------------------------------------
## --------------------------- CAH avec seuil ----------------------------------
# Effectuer le clustering hiérarchique
cah_with_treshold <- function(D_matrix, treshold = 10){
  hc <- hclust(as.dist(D_matrix), method = "complete")
  
  
  # Définir un seuil de hauteur pour la coupure
  height_threshold <- 15  # Ajustez cette valeur en fonction de votre dendrogramme
  
  # Effectuer le clustering avec le seuil de hauteur
  groups_optimal <- cutree(hc, h = height_threshold)
  
  # Afficher le nombre de clusters résultant
  num_clusters <- length(unique(groups_optimal))
  print(paste("Le nombre de clusters avec le seuil de hauteur est:", num_clusters))
  
  # Visualiser le dendrogramme avec les groupes optimaux
  fviz_dend(hc, k = num_clusters,
            cex = 0.6,
            hang = -1,
            rect = TRUE,
            main = "Dendrogramme de la CAH")
  
  # Afficher les groupes
  print(groups_optimal)
  
  # Visualiser les courbes fonctionnelles par groupe
  group_colors <- rainbow(num_clusters)
  
  plot(fd_obj[[1]], col = group_colors[groups_optimal[1]], lty = 1,
       main = "Courbes fonctionnelles par groupe",
       ylab = "Célérité", xlab = "Profondeurs")
  
  for (i in 2:length(fd_obj)) {
    lines(fd_obj[[i]], col = group_colors[groups_optimal[i]], lty = 1)
  }
  
  legend("topright",
         legend = paste("Groupe", 1:num_clusters),
         col = group_colors,
         lty = 1)
  
}

## ---------------------- CAH avec silhouette opti -----------------------------
library(cluster)
library(factoextra)

cah_optimal_silhouette <- function(D_matrix, fd_obj) {
  # Calculer la silhouette moyenne pour différents nombres de clusters
  silhouette_scores <- sapply(2:5, function(k) {
    hc <- hclust(as.dist(D_matrix), method = "complete")
    groups <- cutree(hc, k = k)
    silhouette_avg <- mean(silhouette(groups, as.dist(D_matrix))[, 3])
    return(silhouette_avg)
  })
  
  # Trouver le nombre de clusters qui maximise la silhouette moyenne
  optimal_k <- which.max(silhouette_scores) + 1  
  
  # Effectuer le clustering avec le nombre optimal de clusters
  hc_optimal <- hclust(as.dist(D_matrix), method = "complete")
  groups_optimal <- cutree(hc_optimal, k = optimal_k)
  
  # Calculer les coefficients de silhouette
  sil_info <- silhouette(groups_optimal, as.dist(D_matrix))
  sil_avg <- mean(sil_info[, 3])
  
  # Calculer le score de Davies-Bouldin
  db_score <- davies.bouldin(D_matrix, groups_optimal)
  
  # Afficher les coefficients de silhouette et le score de Davies-Bouldin
  print(paste("Silhouette moyenne pour", optimal_k, "clusters:", round(sil_avg, 2)))
  print(paste("Score de Davies-Bouldin pour", optimal_k, "clusters:", round(db_score, 2)))
  
  # Visualiser le dendrogramme avec les groupes optimaux
  fviz_dend(hc_optimal, k = optimal_k,
            cex = 0.6,
            hang = -1,
            rect = TRUE,
            main = "Dendrogramme de la CAH")
  
  # Visualiser les courbes fonctionnelles par groupe
  group_colors <- rainbow(optimal_k)
  
  plot(fd_obj[[1]], col = group_colors[groups_optimal[1]], lty = 1,
       main = "Courbes fonctionnelles par groupe",
       ylab = "Célérité", xlab = "Profondeurs")
  
  for (i in 2:length(fd_obj)) {
    lines(fd_obj[[i]], col = group_colors[groups_optimal[i]], lty = 1)
  }
  
  legend("topright", legend = paste("Groupe", 1:optimal_k), col = group_colors, lty = 1)
  
  # Retourner les informations de silhouette et le score de Davies-Bouldin
  return(list(hc = hc_optimal, sil_info = sil_info, avg_silhouette = sil_avg, 
              db_score = db_score, groups = groups_optimal, k_optimal = optimal_k))
}



# ------------------------- Davies Bouldin score -------------------------------
davies.bouldin <- function(D_matrix, clusters) {
  library(clusterSim)
  db <- index.DB(D_matrix, clusters, centrotypes = "centroids")$DB
  return(db)
}

