# ------- Fonctions de lissage par spline cubique avec une CV de lambda --------

gcv_MPC_bloc = function(bloc, lambda, base){
  if(is.vector(bloc)){
    mat = smooth.basis(as.numeric(names(bloc)), bloc%>%unname(), fdPar(base, 2, lambda = lambda))$gcv
  } else{
    mat = smooth.basis(as.numeric(colnames(bloc)), bloc%>%t(), fdPar(base, 2, lambda = lambda))$gcv
  }
  
  return(mat)
}
spline_lissage_bloc_quantile = function(bloc, l_grille, D, z){
  # bloc pas de valeurs manquantes
  # gestion du cas un bloc d'un vecteur
  # grille à pas variable
  q = quantile(z)
  d = (D-2)/4  # je sais que je prendrais D = 50 à adapter à l'avenir 

  lev = z

  breaks = c(seq(q[1], q[2], length.out = d), seq(q[2] + 0.1, q[3], length.out = d), seq(q[3] + 0.1, q[4], length.out = d), seq(q[4] + 0.1, q[5], length.out = d))
  base = create.bspline.basis(rangeval = range(lev), nbasis = D, breaks = breaks)
  
  if (is.vector(bloc)){
    mat = matrix(0, nrow = 1, ncol = length(l_grille))
    for (i in 1:length(l_grille)){
      mat[, i] = gcv_MPC_bloc(bloc, l_grille[i], base)
    }
    # lambda optimal pour chaque courbe
    lambda_optimaux = l_grille[which.min(mat%>%c())]
    liss_df = smooth.basis(lev, bloc%>%unname() , fdPar(base, 2, lambda = lambda_optimaux))$fd
    
  }else{
    # bloc --> matrice
    mat = matrix(0, nrow = dim(bloc)[1], ncol = length(l_grille))
    for (i in 1:length(l_grille)){
      mat[, i] = gcv_MPC_bloc(bloc, l_grille[i], base)
    }
    # lambda optimal pour chaque courbe
    lambda_optimaux = l_grille[apply(mat, MARGIN = 1, FUN = which.min)]
    liss_df = list()
    for (i in 1:dim(bloc)[1]){ 
      y = bloc[i,]
      liss = smooth.basis(lev, y%>%unname()%>%unlist()%>%as.numeric() , fdPar(base, 2, lambda = lambda_optimaux[i]))$fd 
      liss_df[[i]] = liss
      
    }
    
  }
  
  return(fd_obj = liss_df)
}

