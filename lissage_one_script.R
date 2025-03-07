library(fda)

lissage_spline <- function(data, D, breaks = NULL) {
  # Calcul de la grille commune pour toutes les données
  n_max <- max(sapply(data, function(x) sum(!is.na(x))))
  grid <- seq(0, 1, length.out = n_max)
  
  # Base de spline cubique naturelle
  basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = D, norder = 4, breaks = breaks)
  
  # Stocker les résultats lissés
  lissage_resultats <- matrix(NA, nrow = nrow(data), ncol = n_max)
  lambdas_optimaux <- numeric(nrow(data))
  
  # Lissage de chaque profil
  for (i in 1:nrow(data)) {
    profil <- as.numeric(data[i, ])
    profil <- profil[!is.na(profil)] # Suppression des NA
    x <- seq(0, 1, length.out = length(profil))
    
    if (length(profil) < 4) {
      next  # Sauter les profils trop courts pour le lissage
    }
    
    # Sélection du lambda optimal via validation croisée
    gcv <- sapply(seq(0.001, 1, length.out = 100), function(lambda) {
      fdParObj <- fdPar(basis, Lfdobj = 2, lambda = lambda)
      smooth <- smooth.basis(x, profil, fdParObj)
      smooth$gcv
    })
    lambda_opt <- seq(0.001, 1, length.out = 100)[which.min(gcv)]
    lambdas_optimaux[i] <- lambda_opt
    
    # Lissage avec le lambda optimal
    fdParObj <- fdPar(basis, Lfdobj = 2, lambda = lambda_opt)
    smooth <- smooth.basis(x, profil, fdParObj)
    
    # Évaluation du lissage sur la grille commune
    lissage_resultats[i, ] <- eval.fd(grid, smooth$fd)
  }
  
  # Conversion en DataFrame avec gestion des noms de lignes
  lissage_df <- as.data.frame(lissage_resultats)
  rownames(lissage_df) <- rownames(data)
  
  return(list(lissage = lissage_df, lambdas = lambdas_optimaux))
}
# ------------------------------------------------------------------------------
#                           application
# ------------------------------------------------------------------------------


data <- read.csv("./data/var_SV_2018-01-01_00H_nan-forced_depth.csv", sep = ";")

rownames(data) <- data$pixel

data <- data[1:500, ]
data <- data[,-1]

row_lengths <- apply(data, 1, function(row) sum(!is.na(row)))

unique_lengths <- sort(unique(row_lengths))

D <- 20


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
data_lissse <- lissage_spline(data, D)

library(ggplot2)
library(reshape2)

# Récupération des données lissées
lissage_df <- data_lissse$lissage
grid <- seq(0, 1, length.out = ncol(lissage_df))

# Conversion en format long pour ggplot2 avec la colonne "Profil" correcte
lissage_df$Profil <- rownames(lissage_df)
lissage_melt <- melt(lissage_df, id.vars = "Profil", variable.name = "Grid", value.name = "Value")

# Ajout de la grille normalisée
lissage_melt$Grid <- as.numeric(lissage_melt$Grid)
lissage_melt$Grid <- grid[lissage_melt$Grid]

# Suppression des valeurs manquantes (si nécessaire)
lissage_melt <- na.omit(lissage_melt)

# Tracé des courbes lissées
ggplot(lissage_melt, aes(x = Value , y = Grid, color = factor(Profil), group = Profil)) +
  geom_line() +
  labs(title = "Courbes Lissées avec Splines Cubiques Naturelles",
       x = "Celerite",
       y = "profondeur",
       color = "Profil") +
  theme_minimal() +
  scale_y_reverse()+
  theme(legend.position = "none")  # Cacher la légende si trop de profils

