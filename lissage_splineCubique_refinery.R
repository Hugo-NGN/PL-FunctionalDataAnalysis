library(fda)

data(refinery)

#View(refinery)

D <- 10  #dimension de la base

yield_data <- refinery$Tray47

time_range <- range(refinery$Time)
################################################################################
# Créer une base pour les splines cubiques avec la plage correcte
basis <- create.bspline.basis(rangeval = time_range, nbasis = D, norder = 4)

#plot base de fonction
plot(basis, main="Fonctions de la base de B-spline", xlab="Temps (min)", ylab="Valeur de la base")


################################################################################
# Calculer les coefficients de spline cubique
lambda_gcv <- lambda2gcv(argvals = refinery$Time, 
                         y = yield_data, 
                         fdParobj = basis, 
                         log10lambda = seq(-4, 4, length.out = 100))

fdobj <- smooth.basis(argvals = refinery$Time, y = yield_data, fdParobj = fdPar(basis, lambda = lambda_gcv))
  
################################################################################
# Tracer le résultat
plot(refinery$Time, yield_data,
     main="Lissage de spline cubique", xlab="Temps (min)", ylab="Production pétrole",
     pch=19, col="gray")

lines(fdobj, col="blue", lwd=2)

legend("bottomright", legend=c("Données originales", "Lissage spline cubique"),
       col=c("gray", "blue"), lty=1:1, pch=c(19, NA), lwd=c(1, 2))





################################################################################
#plot base de fonctions ponderees
coefs <- coef(fdobj)
time_vals <- seq(time_range[1], time_range[2], length.out = 100)

# Créer une matrice pour stocker les valeurs des bases pondérées
weighted_basis <- matrix(0, nrow = length(time_vals), ncol = length(coefs))

for (i in 1:length(coefs)) {
  weighted_basis[, i] <- coefs[i] * eval.basis(time_vals, basis)[, i]
}

# Tracer chaque fonction de base pondérée
matplot(time_vals, weighted_basis, type = "l", lty = 1, col = rainbow(length(coefs)),
        main = "Fonctions de base pondérées par les coefficients",
        xlab = "Temps (min)", ylab = "Valeur pondérée")
legend("topright", legend = paste("Fonction", 1:length(coefs)), col = rainbow(length(coefs)), lty = 1, lwd = 2)

