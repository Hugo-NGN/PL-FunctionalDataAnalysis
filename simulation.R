# Simulation name : fda - dataset : growth 
library(fda)
library(tidyverse)

source("./utils/clustering.R")
source("./utils/distances_fonctionnelles.R")

# -------------------------- Données Simulation --------------------------------


data_growth <- fda::growth


# ---------------------------- VISUALISATION -----------------------------------
## ----------------------------- Courbes M -------------------------------------
matplot(data_growth$age,
        data_growth$hgtm,
        type="l",
        lty=1,
        xlab="Age",
        ylab="Height")

title("Courbes de croissance (M)")

## ----------------------------- Courbes F -------------------------------------
matplot(data_growth$age,
        data_growth$hgtf,
        type="l",
        lty=1,
        xlab="Age",
        ylab="Height")

title("Courbes de croissance (F)")

## ---------------------------- Comparaison M/F --------------------------------

data_combined <- cbind(data_growth$hgtm, data_growth$hgtf)

colors <- c(rep("blue", ncol(data_growth$hgtm)), rep("red", ncol(data_growth$hgtf)))

matplot(data_growth$age, data_combined, type="l", lty=1, col=colors, xlab="Age", ylab="Height")

title("Courbes de croissance par genre")

legend("topright", legend=c("M", "F"), col=c("blue", "red"), lty=1)



# ---------------------------- LISSAGE B-SPLINE --------------------------------
# Base de spline
basis <- create.bspline.basis(rangeval=range(data_growth$age), nbasis=7)

# lissage par sexe
fd_m <- smooth.basis(data_growth$age, data_growth$hgtm, basis)$fd
fd_f <- smooth.basis(data_growth$age, data_growth$hgtf, basis)$fd


# fusion en une liste d'objet fd
fd_list <- list()
for (i in seq_len(ncol(data_growth$hgtm))) {
  fd_list[[i]] <- smooth.basis(data_growth$age, data_growth$hgtm[, i], basis)$fd
}
for (i in seq_len(ncol(data_growth$hgtf))) {
  fd_list[[length(fd_list) + 1]] <- smooth.basis(data_growth$age, data_growth$hgtf[, i], basis)$fd
}

# Visualisation 
plot(fd_m, col="blue", lwd=1, main="Courbes de croissance lissées (M)")
plot(fd_f, col="red", lwd=1, main="Courbes de croissance lissées (F)")

# Comparaison des courbes lissées
#plot(fd_m, col=colors, lwd=1, main="Courbes de croissance lissées par genre")
#plot(fd_f, add=TRUE, col=colors, lwd=2)
#legend("topright", legend=c("M", "F"), col=c("blue", "red"), lty=1)

# ---------------------------- DERIVEES DES COURBES ---------------------------

# Calcul des dérivées
fd_m_prime <- eval.fd(fdobj = fd_m,evalarg = data_growth$age,  Lfdobj=1)
fd_f_prime <- eval.fd(fdobj = fd_f, evalarg = data_growth$age, Lfdobj=1)

# Visualisation
fd_m_prime %>% matplot(type="l",
                       main = "Dérivées des courbes de croissance (M)",
                       xlab="Age",
                       ylab="Dérivée")

fd_f_prime %>% matplot(type="l",
                       main = "Dérivées des courbes de croissance (F)",
                       xlab="Age",
                       ylab="Dérivée")


# ---------------------------- MATRICE DES DISTANCES ---------------------------
fine_grid <- seq(min(data_growth$age), max(data_growth$age), length.out = 1000)

D0_matrix <- calculate_D0_matrix_parallel(fd_list, fine_grid=fine_grid)
D1_matrix <- calculate_D1_matrix_parallel(fd_list, fine_grid=fine_grid)
Dp_matrix <- calculate_Dp_matrix_parallel(fd_list, fine_grid=fine_grid, omega= 0.8)



# ---------------------------- CLUSTERING -------------------------------------
## ---------------------------- BASELINE ---------------------------------------

## ----------------------- K-MEANS fonctionnel ---------------------------------

## ----------------------- CAH fonctionnelle -----------------------------------

## ------------------ CAH + K-MEANS fonctionnel --------------------------------












