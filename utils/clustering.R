library(fpc)
library(cluster)
library(factoextra)
# ------------------------- Davies Bouldin score -------------------------------
davies.bouldin <- function(D_matrix, clusters) {
  library(clusterSim)
  db <- index.DB(D_matrix, clusters, centrotypes = "centroids")$DB
  return(db)
}



# --------------------------------- CAH ----------------------------------------
## --------------------------- CAH avec seuil ----------------------------------
# Effectuer le clustering hiérarchique
cah_with_treshold <- function(D_matrix, height_threshold = 10){
  hc <- hclust(as.dist(D_matrix), method = "complete")
  
  # Effectue le clustering avec le seuil
  groups_optimal <- cutree(hc, h = height_threshold)
  
  num_clusters <- length(unique(groups_optimal))
  
  print(paste("Le nombre de clusters avec le seuil de hauteur est:",
              num_clusters))
  
  fviz_dend(hc, k = num_clusters,
            cex = 0.6,
            hang = -1,
            rect = TRUE,
            main = "Dendrogramme de la CAH")
  
  
  # affiche les courbes fonctionnelles par groupe
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
cah_optimal_silhouette <- function(D_matrix, fd_obj) {
  # Calculer la silhouette moyenne pour différents nombres de clusters
  silhouette_scores <- sapply(2:5, function(k) {
    hc <- hclust(as.dist(D_matrix), method = "complete")
    groups <- cutree(hc, k = k)
    silhouette_avg <- mean(silhouette(groups, as.dist(D_matrix))[, 3])
    return(silhouette_avg)
  })
  
  # Trouve le nombre de clusters qui maximise la silhouette moyenne
  optimal_k <- which.max(silhouette_scores) + 1  
  
  # Effectue le clustering avec le nombre optimal de clusters
  hc_optimal <- hclust(as.dist(D_matrix), method = "complete")
  groups_optimal <- cutree(hc_optimal, k = optimal_k)
  
  # Calcule les coefficients de silhouette
  sil_info <- silhouette(groups_optimal, as.dist(D_matrix))
  sil_avg <- mean(sil_info[, 3])
  
  # Calcule le score de Davies-Bouldin
  db_score <- davies.bouldin(D_matrix, groups_optimal)
  
  # Affiche les coefficients de silhouette et le score de Davies-Bouldin
  print(paste("Silhouette moyenne pour",
              optimal_k, "clusters:",
              round(sil_avg, 2)))
  
  print(paste("Score de Davies-Bouldin pour",
              optimal_k, "clusters:",
              round(db_score, 2)))
  
  # Visualisation du dendrogramme de CAH
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
  
  legend("topright",
         legend = paste("Groupe", 1:optimal_k),
         col = group_colors,
         lty = 1)
  
  
  return(list(hc = hc_optimal, sil_info = sil_info, avg_silhouette = sil_avg, 
              db_score = db_score, groups = groups_optimal, k_optimal = optimal_k))
}



# --------------------------------- K MEANS ------------------------------------
kmeans_optimal_k <- function(D_matrix){
  silhouette_scores <- numeric(9)
  
  # Boucle sur les valeurs de k de 2 à 10
  for (k in 2:10) {
    # Appliquer k-means
    km_result <- kmeans(D_matrix, centers = k, nstart = 25)
    
    # Calculer le coefficient de silhouette
    silhouette_score <- silhouette(km_result$cluster, as.dist(D_matrix))
    avg_silhouette_width <- mean(silhouette_score[, 3])
    
    # Stocker le coefficient de silhouette moyen
    silhouette_scores[k - 1] <- avg_silhouette_width
  }
  
  # Tracer la courbe du coefficient de silhouette en fonction de k
  plot(2:10, silhouette_scores, type = "b",
       xlab = "Nombre de clusters (k)",
       ylab = "Coefficient de silhouette moyen",
       main = "Coefficient de silhouette en fonction de k")
  
  # Retourne le nombre optimal de k
  return(which.max(silhouette_scores) + 1)
}


kmeans_fd <- function(D_matrix, k, fd_obj){
  
  km_result <- kmeans(D_matrix, centers = k, nstart = 25)
  
  # Calcul du coefficient de silhouette moyen
  sil_info <- silhouette(km_result$cluster, as.dist(D_matrix))
  sil_avg <- mean(sil_info[, 3])
  
  # Calcul du score de Davies-Bouldin
  db_score <- davies.bouldin(D_matrix, km_result$cluster)
  
  # Afficher les scores
  print(paste("Silhouette moyenne pour", k, "clusters:", round(sil_avg, 2)))
  print(paste("Score de Davies-Bouldin pour", k, "clusters:", round(db_score, 2)))
  
  # Visualiser les groupes
  group_colors <- rainbow(k)
  plot(fd_obj[[1]], col = group_colors[km_result$cluster[1]], lty = 1,
       main = "Courbes fonctionnelles par groupe",
       ylab = "Célérité", xlab = "Profondeurs")
  
  for (i in 2:length(fd_obj)) {
    lines(fd_obj[[i]], col = group_colors[km_result$cluster[i]], lty = 1)
  }
  
  # Ajouter une légende
  legend("topright", legend = paste("Groupe", 1:k), col = group_colors, lty = 1)
  
  return(list(km_result = km_result, sil_info = sil_info, avg_silhouette = sil_avg, db_score = db_score))
}


