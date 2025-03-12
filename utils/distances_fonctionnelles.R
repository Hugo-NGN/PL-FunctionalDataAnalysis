library(fda)

# -------------------------------- D0 ------------------------------------------
calculate_D0_matrix <- function(fd_obj, fine_grid) {
  n <- length(fd_obj)
  D0_matrix <- matrix(0, nrow = n, ncol = n)
  
  # Largeur des intervalles pour la somme de Riemann
  delta <- diff(range(fine_grid)) / (length(fine_grid) - 1)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      # Évaluer les courbes sur une grille fine
      fd_i_vals <- eval.fd(fine_grid, fd_obj[[i]])
      fd_j_vals <- eval.fd(fine_grid, fd_obj[[j]])
      
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


# -------------------------------- D1 ------------------------------------------
calculate_D1_matrix <- function(fd_obj, fine_grid){
  fd_obj_deriv = list()
  for (i in seq(1, length(fd_obj))){
    fd_obj_deriv[[i]] <- eval.fd(fine_grid, fd_obj[[i]], Lfdobj=1)
  }
  
  n <- length(fd_obj_deriv)
  D1_matrix <- matrix(0, nrow = n, ncol = n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        # Calculer la différence des coordonnées
        diff <- fd_obj_deriv[[i]] - fd_obj_deriv[[j]]
        # Sommer les différences au carré - approximation de l'intégrale de D1 par une somme de Riemann 
        D1_matrix[i, j] <- sum(diff^2)
      }
    }
  }
  return(D1_matrix)
}


# -------------------------------- Dp ------------------------------------------

calculate_Dp_matrix <- function(fd_obj, fine_grid, omega) {
  
  D0_matrix <- calculate_D0_matrix(fd_obj, fine_grid)
  D1_matrix <- calculate_D1_matrix(fd_obj, fine_grid)
  
  Dp_matrix <- sqrt((1 - omega) * D0_matrix^2 + omega * D1_matrix^2)
  
  return(Dp_matrix)
  
}

