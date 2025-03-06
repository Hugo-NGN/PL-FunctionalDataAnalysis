library(fda)

data("CanadianWeather")

temperature <- CanadianWeather$dailyAv[, , "Temperature.C"]
time_points <- 1:365  # Jours de l'année

# Définir une base B-spline avec 15 bases pour un lissage raisonnable
nbasis <- 15
basis <- create.bspline.basis(rangeval = c(1, 365), nbasis = nbasis, norder = 4)

# Choisir une séquence de paramètres de lissage pour GCV
log_lambda_seq <- seq(-4, 4, length.out = 50)

# Calculer le paramètre de lissage optimal par GCV
lambda_gcv <- lambda2gcv(argvals = time_points, 
                         y = temperature, 
                         fdParobj = basis, 
                         log10lambda = log_lambda_seq)

cat("Paramètre de lissage optimal (lambda) :", 10^lambda_gcv, "\n")

# Lissage des courbes avec le meilleur lambda
best_index <- which.min(lambda_gcv$gcv)  # Trouver l'indice avec le GCV minimal
lambda_value <- 10^lambda_gcv  # Convertir l'indice log10lambda en lambda

# Créer l'objet fdPar avec le bon lambda
fdParobj <- fdPar(basis, lambda = lambda_value)
smooth_temp <- smooth.basis(argvals = time_points, y = temperature, fdParobj = fdParobj)

# Visualiser les données brutes et lissées pour quelques villes
cities_to_plot <- c("Vancouver", "Toronto", "Montreal")
city_indices <- match(cities_to_plot, colnames(temperature))

# Tracer les données brutes et lissées
par(mfrow = c(1, length(cities_to_plot)))
for (i in city_indices) {
  plot(time_points, temperature[, i], type = "p", col = "gray",
       main = colnames(temperature)[i], xlab = "Jour", ylab = "Température (°C)")
  lines(time_points, eval.fd(time_points, smooth_temp$fd[, i]), col = "blue", lwd = 2)
  legend("topright", legend = c("Données brutes", "Lissage B-spline"),
         col = c("gray", "blue"), lty = 1, pch = c(1, NA), lwd = c(1, 2))
}

par(mfrow = c(1, 1))  # Réinitialiser la disposition des graphiques