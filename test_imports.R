################################################################################
#                                     LISSAGE
################################################################################
source("./utils/lissage.R")

data(refinery)

D <- 13  #dimension de la basE

data <- refinery$Tray47
val <- refinery$Time

res_lissage <- lissage_spline_cubique(data, val, D)


plot(val, data,
     main="Lissage de spline cubique", xlab="Temps (min)", ylab="Production pétrole",
     pch=19, col="gray")

lines(res_lissage, col="blue", lwd=2)

legend("bottomright", legend=c("Données originales", "Lissage spline cubique"),
       col=c("gray", "blue"), lty=1:1, pch=c(19, NA), lwd=c(1, 2))





######                  EXEMPLE VISUALISATION DES NOEUDS

# Remarque: ici les noeuds sont équidistant car les données simulées sont uniformément
#           "connues" sur l'intervalle de "val"

q = quantile(as.numeric(val)) 

nb_breaks = D-2 

nb_breaks_byQuartile = (nb_breaks%/%4)

remainder = nb_breaks - (nb_breaks_byQuartile*4)

breaks = c(seq(q[1], q[2], length.out = nb_breaks_byQuartile),
           seq(q[2] + 5, q[3], length.out = nb_breaks_byQuartile+ remainder),
           seq(q[3] + 5, q[4], length.out = nb_breaks_byQuartile),
           seq(q[4] + 5, q[5], length.out = nb_breaks_byQuartile))


basis <- create.bspline.basis(rangeval = range(val), nbasis = D, norder = 4, breaks = breaks )


# Tracer les noeuds
par(mar = c(4, 4, 2, 1))
plot(val, rep(0, length(val)), type = "n", xlab = "val", ylab = "Noeuds", main = "Placement des Noeuds")
abline(v = breaks, col = "red", lty = 2)
points(breaks, rep(0, length(breaks)), col = "red", pch = 16)




length(breaks)


################################################################################
#                               LISSAGE BLOC
################################################################################
source('./Code_R/utils/lissage.R')
library('fda')

data("CanadianWeather")

daily_avgtemp <- t(MontrealTemp[, 16:47]) #matrice 32*34

# Définir les valeurs pour 'val' et 'D'
jours <- ((16:47)+0.5)
D <- 20  
fdobj_list <- lissage_spline_cubique_bloc(t(daily_avgtemp), jours, D)


# Définir une grille de points pour l'évaluation des splines
x_vals <- seq(min(jours), max(jours), length.out = 1000)
colors <- rainbow(length(fdobj_list))

# Initialiser le graphique avec la première courbe
plot(x_vals, eval.fd(x_vals, fdobj_list[[1]]), type = "l",
     main = "Courbes lissées des températures à Montréal",
     ylab = "Température lissée (°C)",
     xlab = "Jours",
     col = colors[1],
     ylim = c(-30,5),
     lwd = 2)

# Ajouter les autres courbes
for (i in 2:length(fdobj_list)) {
  lines(x_vals, eval.fd(x_vals, fdobj_list[[i]]), col = colors[i], lwd = 2)
}

# Ajouter une légende
legend("topright", legend = paste("année", 1:32), col = colors, lwd = 2, cex = 0.6)

##############################################################################
#                                     DISTANCES
################################################################################

source("./utils/distances_fonctionnelles.R")

X1 <- function(t) {
  cos(t)
}
X2 <- function(t) {
  1+cos(t)
}


X11 <- function(t) {
  3*t +cos(t)
}


# Définir l'intervalle T et le paramètre omega
T <- c(0, 2 * pi)
omega <- 0.5

# Calculer D_0, D_1 et D_p
result_D0 <- D_0(X1, X2, T)
result_D1 <- D_1(X1, X2, T)
result_Dp <- D_p(X1, X2, T, omega)

# Afficher les résultats
cat("D_0:", result_D0, "\n")
cat("D_1:", result_D1, "\n")
cat("D_p:", result_Dp, "\n")




n= 1000
t_values <- seq(T[1], T[2], length.out = n)

# Créer un dataframe avec les valeurs des processus
df <- data.frame(
  t = t_values,
  X1 = X1(t_values),
  X2 = X2(t_values),
  X11 = X11(t_values)
)

# Calculer la matrice de distances
omega <- 0.5
distance_matrix <- D_p_matrix(df, T, omega)

# Afficher la matrice de distances
print(distance_matrix)


hc <- hclust(as.dist(distance_matrix), method = "ward.D2")

# Afficher le dendrogramme
plot(hc, main = "Dendrogramme du Clustering Hiérarchique", sub = "", xlab = "", ylab = "Distance", hang = -1)

# Découper en clusters
clusters <- cutree(hc, k = 2)  # Par exemple, pour 2 clusters
print(clusters)







