library(pracma)

################################################################################
#                       Fonction pour calculer D_0
################################################################################





D_0 <- function(X1, X2, T, n = 1000) {
  # Créer un vecteur de points dans l'intervalle T
  t_values <- seq(T[1], T[2], length.out = n)
  # Calculer l'intégrale du carré des différences
  integrand <- function(t) {
    (X1(t) - X2(t))^2
  }
  # Calculer la racine carrée de l'intégrale
  sqrt(integrate(integrand, lower = T[1], upper = T[2])$value)
}






################################################################################
#                       Fonction pour calculer D_1
################################################################################




D_1 <- function(X1, X2, T, n = 1000) {
  # Créer un vecteur de points dans l'intervalle T
  t_values <- seq(T[1], T[2], length.out = n)
  # Calculer les valeurs des fonctions sur les points
  X1_values <- X1(t_values)
  X2_values <- X2(t_values)
  # Calculer les dérivées numériques
  X1_deriv <- gradient(X1_values, t_values)
  X2_deriv <- gradient(X2_values, t_values)
  # Calculer l'intégrale du carré des différences des dérivées
  integrand <- function(t) {
    approx(t_values, (X1_deriv - X2_deriv)^2, t)$y
  }
  # Calculer la racine carrée de l'intégrale
  sqrt(integrate(integrand, lower = T[1], upper = T[2])$value)
}




################################################################################
#                       Fonction pour calculer D_p
################################################################################




D_p <- function(X1, X2, T, omega, n = 1000) {
  t_values <- seq(T[1], T[2], length.out = n)
  
  # Calculer les valeurs des fonctions sur les points
  X1_values <- X1(t_values)
  X2_values <- X2(t_values)
  
  # Calculer les dérivées numériques
  X1_deriv <- gradient(X1_values, t_values)
  X2_deriv <- gradient(X2_values, t_values)
  
  # Calculer les intégrales pondérées
  integrand_0 <- function(t) {
    (X1(t) - X2(t))^2
  }
  integrand_1 <- function(t) {
    approx(t_values, (X1_deriv - X2_deriv)^2, t)$y
  }
  
  integral_0 <- integrate(integrand_0, lower = T[1], upper = T[2])$value
  integral_1 <- integrate(integrand_1, lower = T[1], upper = T[2])$value
  
  # Calculer D_p
  sqrt((1 - omega) * integral_0 + omega * integral_1)
}





################################################################################
#                       Fonction pour calculer D_0 (matrice)
################################################################################



D_0_matrix <- function(df, T, n = 1000) {
  
  num_processes <- ncol(df) - 1  
  
  # Initialiser une matrice de distances
  distance_matrix <- matrix(0, nrow = num_processes, ncol = num_processes)
  
  # Calculer les distances pour chaque paire de processus
  for (i in 1:(num_processes - 1)) {
    for (j in (i + 1):num_processes) {
      X1 <- function(t) approx(df[, 1], df[, i+1], t)$y
      X2 <- function(t) approx(df[, 1], df[, j+1], t)$y
      
      # Calculer D_0, D_1, et D_p
      distance_matrix[i, j] <- D_0(X1, X2, T, n)
      distance_matrix[j, i] <- distance_matrix[i, j]
    }
  }
  
  return(distance_matrix)
}


################################################################################
#                       Fonction pour calculer D_1 (matrice)
################################################################################



D_1_matrix <- function(df, T, n = 1000) {
  
  num_processes <- ncol(df) - 1  
  
  # Initialiser une matrice de distances
  distance_matrix <- matrix(0, nrow = num_processes, ncol = num_processes)
  
  # Calculer les distances pour chaque paire de processus
  for (i in 1:(num_processes - 1)) {
    for (j in (i + 1):num_processes) {
      X1 <- function(t) approx(df[, 1], df[, i+1], t)$y
      X2 <- function(t) approx(df[, 1], df[, j+1], t)$y
      
      # Calculer D_0, D_1, et D_p
      distance_matrix[i, j] <- D_1(X1, X2, T, n)
      distance_matrix[j, i] <- distance_matrix[i, j]
    }
  }
  
  return(distance_matrix)
}





################################################################################
#                       Fonction pour calculer D_p (matrice)
################################################################################



D_p_matrix <- function(df, T, omega = 0.5, n = 1000) {

  num_processes <- ncol(df) - 1  
  
  # Initialiser une matrice de distances
  distance_matrix <- matrix(0, nrow = num_processes, ncol = num_processes)
  
  # Calculer les distances pour chaque paire de processus
  for (i in 1:(num_processes - 1)) {
    for (j in (i + 1):num_processes) {
      X1 <- function(t) approx(df[, 1], df[, i+1], t)$y
      X2 <- function(t) approx(df[, 1], df[, j+1], t)$y
      
      # Calculer D_0, D_1, et D_p
      distance_matrix[i, j] <- D_p(X1, X2, T, omega, n)
      distance_matrix[j, i] <- distance_matrix[i, j]
    }
  }
  
  return(distance_matrix)
}
