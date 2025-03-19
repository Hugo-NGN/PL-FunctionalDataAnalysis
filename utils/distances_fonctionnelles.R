library(fda)
library(parallel)
# ------------------ VERSIONS SEQUENTIELLES (non parallelisees) ----------------
## -------------------------------- D0 -----------------------------------------
calculate_D0_matrix <- function(fd_list, fine_grid) {
  n <- length(fd_list)
  D0_matrix <- matrix(0, nrow = n, ncol = n)
  
  # Largeur des intervalles pour la somme de Riemann
  delta <- diff(range(fine_grid)) / (length(fine_grid) - 1)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      # Évaluer les courbes sur une grille fine
      fd_i_vals <- eval.fd(fine_grid, fd_list[[i]])
      fd_j_vals <- eval.fd(fine_grid, fd_list[[j]])
      
      # Calculer la différence au carré
      diff_squared <- (fd_i_vals - fd_j_vals)^2
      
      # Approcher l'intégrale par la somme de Riemann
      integral <- sum(diff_squared) * delta
      
      # Calculer la distance D0
      D0_matrix[i, j] <- sqrt(integral)
      D0_matrix[j, i] <- D0_matrix[i, j]
    }
  }
  
  return(D0_matrix)
}


## -------------------------------- D1 -----------------------------------------
calculate_D1_matrix <- function(fd_list, fine_grid){
  fd_list_deriv = list()
  for (i in seq(1, length(fd_list))){
    fd_list_deriv[[i]] <- eval.fd(fine_grid, fd_list[[i]], Lfdobj=1)
  }
  
  n <- length(fd_list_deriv)
  D1_matrix <- matrix(0, nrow = n, ncol = n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        # Calculer la différence des coordonnées
        diff <- fd_list_deriv[[i]] - fd_list_deriv[[j]]
        # Sommer les différences au carré - approximation de l'intégrale de D1 par une somme de Riemann 
        D1_matrix[i, j] <- sum(diff^2)
      }
    }
  }
  return(D1_matrix)
}

## -------------------------------- Dp -----------------------------------------

calculate_Dp_matrix <- function(fd_list, fine_grid, omega) {
  
  D0_matrix <- calculate_D0_matrix(fd_list, fine_grid)
  D1_matrix <- calculate_D1_matrix(fd_list, fine_grid)
  
  Dp_matrix <- sqrt((1 - omega) * D0_matrix^2 + omega * D1_matrix^2)
  
  return(Dp_matrix)
  
}

calculate_Dp_matrix_givenD0D1 <- function(D0_matrix, D1_matrix, omega){
  Dp_matrix <- sqrt((1-omega)*D0_matrix^2 + omega*D1_matrix^2)
  return(Dp_matrix)
} 


# ------------------------ VERSIONS AVEC PARALLELISATION -----------------------
## -------------------------------- D0-par -------------------------------------

calculate_D0_matrix_parallel <- function(fd_list, fine_grid) {
  n <- length(fd_list)
  D0_matrix <- matrix(0, nrow = n, ncol = n)
  
  # Créer un cluster avec le nombre de cœurs disponibles
  no_cores <- detectCores() - 1  # Laisser un cœur libre
  cl <- makeCluster(no_cores)
  
  # Exporter les variables nécessaires aux workers
  clusterExport(cl, varlist = c("fd_list", "fine_grid", "eval.fd"))
  
  # Fonction pour calculer une ligne de la matrice
  compute_row <- function(i) {
    n <- length(fd_list)
    row_values <- numeric(n)
    
    # Recalculer delta à l'intérieur de la fonction (évite l'erreur)
    delta <- diff(range(fine_grid)) / (length(fine_grid) - 1)
    
    for (j in (i+1):n) {
      fd_i_vals <- eval.fd(fine_grid, fd_list[[i]])
      fd_j_vals <- eval.fd(fine_grid, fd_list[[j]])
      
      diff_squared <- (fd_i_vals - fd_j_vals)^2
      integral <- sum(diff_squared) * delta
      
      row_values[j] <- sqrt(integral)
    }
    return(row_values)
  }
  
  # Appliquer la fonction en parallèle
  result_list <- parLapply(cl, 1:(n-1), compute_row)
  
  # Fermer le cluster
  stopCluster(cl)
  
  # Remplir la matrice de résultats
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      D0_matrix[i, j] <- result_list[[i]][j]
      D0_matrix[j, i] <- D0_matrix[i, j]  
    }
  }
  
  return(D0_matrix)
}

## -------------------------------- D1-par -------------------------------------
calculate_D1_matrix_parallel <- function(fd_list, fine_grid) {
  n <- length(fd_list)
  D1_matrix <- matrix(0, nrow = n, ncol = n)
  
  # Créer un cluster
  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores)
  
  # Exporter les variables nécessaires
  clusterExport(cl, varlist = c("fd_list", "fine_grid", "eval.fd"))
  
  # Fonction pour calculer une ligne de la matrice
  compute_row <- function(i) {
    fd_list_deriv_i <- eval.fd(fine_grid, fd_list[[i]], Lfdobj = 1)
    row_values <- numeric(n)
    
    for (j in 1:n) {
      if (i != j) {
        fd_list_deriv_j <- eval.fd(fine_grid, fd_list[[j]], Lfdobj = 1)
        diff <- fd_list_deriv_i - fd_list_deriv_j
        row_values[j] <- sum(diff^2)  # Approximation par somme de Riemann
      }
    }
    
    return(row_values)
  }
  
  # Appliquer en parallèle
  result_list <- parLapply(cl, 1:n, compute_row)
  
  # Fermer le cluster
  stopCluster(cl)
  
  # Remplir la matrice
  for (i in 1:n) {
    D1_matrix[i, ] <- result_list[[i]]
  }
  
  return(D1_matrix)
}

## -------------------------------- Dp-par -------------------------------------

calculate_Dp_matrix_parallel <- function(fd_list, fine_grid, omega, STANDARDIZE = TRUE) {
  
  D0_matrix <- calculate_D0_matrix_parallel(fd_list, fine_grid)
  D1_matrix <- calculate_D1_matrix_parallel(fd_list, fine_grid)
  
  if (STANDARDIZE){
  # Standardisation de D0 et D1 avec min-max
  D0_min <- min(D0_matrix)
  D0_max <- max(D0_matrix)
  D0_standardized <- (D0_matrix - D0_min) / (D0_max - D0_min)
  
  D1_min <- min(D1_matrix)
  D1_max <- max(D1_matrix)
  D1_standardized <- (D1_matrix - D1_min) / (D1_max - D1_min)
  
  Dp_matrix <- sqrt((1 - omega) * D0_standardized^2 + omega * D1_standardized^2)
  } else {
    Dp_matrix <- sqrt((1 - omega) * D0_matrix^2 + omega * D1_matrix^2)
  }
  return(Dp_matrix)
}


calculate_Dp_matrix_byD0D1 <- function(D0_matrix, D1_matrix, omega, STANDARDIZE = FALSE){
  if (!STANDARDIZE){
    return(sqrt((1-omega)*D0_matrix^2 + omega*D1_matrix^2))  
  } else {
    D0_min <- min(D0_matrix)
    D0_max <- max(D0_matrix)
    D0_standardized <- (D0_matrix - D0_min) / (D0_max - D0_min)
    
    D1_min <- min(D1_matrix)
    D1_max <- max(D1_matrix)
    D1_standardized <- (D1_matrix - D1_min) / (D1_max - D1_min)
    
    return(sqrt((1 - omega) * D0_standardized^2 + omega * D1_standardized^2))
  }
  
}










