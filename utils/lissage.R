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
    liss = smooth.basis(lev, bloc%>%unname() , fdPar(base, 2, lambda = lambda_optimaux))$fd
    bloc_hat = eval.fd(lev, liss)
    coefs_bloc = matrix(liss$coefs, ncol = D)
    rmse = sqrt(mean((bloc_hat - bloc)^2))
  }else{
    # bloc --> matrice
    mat = matrix(0, nrow = dim(bloc)[1], ncol = length(l_grille))
    for (i in 1:length(l_grille)){
      mat[, i] = gcv_MPC_bloc(bloc, l_grille[i], base)
    }
    # lambda optimal pour chaque courbe
    lambda_optimaux = l_grille[apply(mat, MARGIN = 1, FUN = which.min)]
    bloc_hat = matrix(0, ncol = dim(bloc)[2], nrow = dim(bloc)[1])
    coefs_bloc = matrix(0, nrow = dim(bloc)[1], ncol =  D)
    rmse = rep(0, dim(bloc)[1])
    liss_df = list()
    for (i in 1:dim(bloc)[1]){ 
      y = bloc[i,]
      liss = smooth.basis(lev, y%>%unname()%>%unlist()%>%as.numeric() , fdPar(base, 2, lambda = lambda_optimaux[i]))$fd 
      liss_df[[i]] = liss
      bloc_hat[i, ] = eval.fd(lev, liss)
      coefs_bloc[i, ] = liss$coefs
      rmse[i] = sqrt(mean((bloc_hat[i, ] - bloc[i, ])^2))
    }
    
  }
  
  return(list(fd_obj = liss_df, bloc = bloc_hat, coefs = coefs_bloc, basis= base, rmse = rmse))
}

