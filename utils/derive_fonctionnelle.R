library(fda)

derive_fd <- function(fd_obj) {
  basis <- fd_obj$basis  # Extraire la base fonctionnelle
  coefs <- fd_obj$coefs  # Extraire les coefficients
  
  # eval.basis(Lfdobj = 1) #evalue les valeurs de la dérivee d'ordre 1
  # basis_deriv <- ... # creer la bases derivee
  # 
  # ...
  # 
  # fd_obj_deriv <- ... # creer objet derive

  
  return(fd_obj_deriv)
}