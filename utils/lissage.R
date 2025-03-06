#           Fonctions de lissage par spline cubique avec une CV de lambda



################################################################################
#                           lissage mono-courbe
################################################################################

lissage_spline_cubique <- function(data, val, D){
  library(fda)
  
  q = quantile(as.numeric(val), na.rm = TRUE) 

  nb_breaks = D-2 
  
  nb_breaks_byQuartile = (nb_breaks%/%4)
  
  remainder = nb_breaks - (nb_breaks_byQuartile*4)
  
  breaks = c(seq(q[1], q[2], length.out = nb_breaks_byQuartile),
             seq(q[2] + 0.1, q[3], length.out = nb_breaks_byQuartile+ remainder),
             seq(q[3] + 0.1, q[4], length.out = nb_breaks_byQuartile),
             seq(q[4] + 0.1, q[5], length.out = nb_breaks_byQuartile))
  
  basis <- create.bspline.basis(rangeval = range(val, na.rm = TRUE), nbasis = D, norder = 4, breaks = breaks )
  
  lambda_gcv <- lambda2gcv(argvals = val, 
                           y = data, 
                           fdParobj = basis, 
                           log10lambda = seq(-4, 4, length.out = 100))
  
  fdobj <- smooth.basis(argvals = val,
                        y = data,
                        fdParobj = fdPar(basis, lambda = lambda_gcv))
  
  return(fdobj)
  
}





################################################################################
#                           lissage par bloc
################################################################################

lissage_spline_cubique_bloc <- function(data, val, D){
  library(fda)
  
  fdobj_list <- list()
  
  # Calcul des points de coupure (breaks) pour les splines
  q <- quantile(as.numeric(val))
  nb_breaks <- D - 2
  nb_breaks_byQuartile <- (nb_breaks %/% 4)
  remainder <- nb_breaks - (nb_breaks_byQuartile * 4)
  
  breaks <- c(
    seq(q[1], q[2], length.out = nb_breaks_byQuartile),
    seq(q[2] + 0.001, q[3], length.out = nb_breaks_byQuartile + remainder),
    seq(q[3] + 0.001, q[4], length.out = nb_breaks_byQuartile),
    seq(q[4] + 0.001, q[5], length.out = nb_breaks_byQuartile)
  )
  
  basis <- create.bspline.basis(rangeval = range(val), nbasis = D, norder = 4, breaks = breaks)
  
  l_grille <- 10^seq(-4, 4, length.out = 100)
  
  for (i in seq_len(nrow(data))) {
    
    y <- data[i,]
    
    # Calcul de lambda optimal via GCV
    #  lambda_gcv <- lambda2gcv(
    #  argvals = val,
    #  y = y,
    #  fdParobj = basis,
    #  log10lambda = 10^seq(-3, 5, length.out = 1000)
    #)
    
    gcv_vals <- sapply(l_grille, function(lambda) {
      gcv_MPC_bloc(y, lambda, basis)
    })
    
    
    fdobj <- smooth.basis(argvals = val,
                          y = y,
                          fdParobj =  fdPar(basis,
                                      Lfdobj = 2,
                                      lambda = lambda_gcv))$fd

    fdobj_list[[i]] <- fdobj
  }
  
  return(fdobj_list)
}


gcv_MPC_bloc = function(bloc, lambda, base){
  if(is.vector(bloc)){
    mat = smooth.basis(as.numeric(names(bloc)), bloc%>%unname(), fdPar(base, 2, lambda = lambda))$gcv
  } else{
    mat = smooth.basis(as.numeric(colnames(bloc)), bloc%>%t(), fdPar(base, 2, lambda = lambda))$gcv
  }
  
  return(mat)
}


