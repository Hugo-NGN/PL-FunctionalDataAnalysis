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
# calculate_D1_matrix <- function(fd_obj, fine_grid) {
#   n <- length(fd_obj)
#   D1_matrix <- matrix(0, nrow = n, ncol = n)
#   
#   # Largeur des intervalles pour la somme de Riemann
#   delta <- diff(range(fine_grid)) / (length(fine_grid) - 1)
#   
#   # Vérifier si fd_obj est bien une liste d'objets fd
#   if (!all(sapply(fd_obj, is.fd))) {
#     stop("fd_obj doit être une liste d'objets fonctionnels fd")
#   }
#   
#   # Calculer les dérivées des objets fd
#   fd_deriv_list <- lapply(fd_obj, derive_fd)
#   
#   for (i in 1:(n-1)) {
#     for (j in (i+1):n) {
#       # Évaluer les dérivées des courbes sur une grille fine
#       fd_i_deriv_vals <- eval.fd(fine_grid, fd_deriv_list[[i]])
#       fd_j_deriv_vals <- eval.fd(fine_grid, fd_deriv_list[[j]])
#       
#       # Calculer la différence au carré des dérivées
#       deriv_diff_squared <- (fd_i_deriv_vals - fd_j_deriv_vals)^2
#       
#       # Approcher l'intégrale par la somme de Riemann
#       integral <- sum(deriv_diff_squared) * delta
#       
#       # Calculer la distance D_1
#       D1_matrix[i, j] <- sqrt(integral)
#       D1_matrix[j, i] <- D1_matrix[i, j]  
#     }
#   }
#   
#   return(D1_matrix)
# }

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
  n <- length(fd_obj)
  Dp_matrix <- matrix(0, nrow = n, ncol = n)
  
  # Largeur des intervalles pour la somme de Riemann
  delta <- diff(range(fine_grid)) / (length(fine_grid) - 1)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      # Évaluer les courbes et leurs dérivées sur une grille fine
      fd_i_vals <- eval.fd(fine_grid, fd_obj[[i]])
      fd_j_vals <- eval.fd(fine_grid, fd_obj[[j]])
      
      fd_i_deriv_vals <- eval.fd(fine_grid, deriv.fd(fd_obj[[i]]))
      fd_j_deriv_vals <- eval.fd(fine_grid, deriv.fd(fd_obj[[j]]))
      
      # Calculer la différence au carré des courbes
      diff_squared <- (fd_i_vals - fd_j_vals)^2
      
      # Calculer la différence au carré des dérivées
      deriv_diff_squared <- (fd_i_deriv_vals - fd_j_deriv_vals)^2
      
      # Approcher l'intégrale par la somme de Riemann
      integral_values <-  (1-omega) * sum(diff_squared) * delta + omega * sum(deriv_diff_squared) * delta
      
      # Calculer la distance D_p
      Dp_matrix[i, j] <- sqrt(integral_values)
      Dp_matrix[j, i] <- Dp_matrix[i, j] 
    }
  }
  return(Dp_matrix)
}

