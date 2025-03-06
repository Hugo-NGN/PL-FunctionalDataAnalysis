# Lissage des profils de célérité
source("./Code_R/fonctionsR/lissage.R")

data <- read.csv("./data/var_SV_2018-01-01_00H_nan-forced_depth.csv", sep = ";")

rownames(data) <- data$pixel

data <- data[1:100, ]

row_lengths <- apply(data, 1, function(row) sum(!is.na(row)))

unique_lengths <- sort(unique(row_lengths))

D <- 20

rows_by_length <- list()
fdobj_list <- list()

val <- as.integer(sub("^X", "", colnames(data[-1])))
val <- val[-length(val)]

for (length in unique_lengths) {
  rows <- data[row_lengths == length, ]
  
  # Supprimer les valeurs NA en coupant les lignes à la bonne longueur
  trimmed_rows <- t(apply(rows, 1, function(row) {
    head(row[!is.na(row)], length) # Correction ici : trim NA correctement
  }))
  
  # Stocker les lignes nettoyées dans la liste
  rows_by_length[[as.character(length)]] <- trimmed_rows
  
  cat("Rows with length", length, ":\n")
  
  ##############################################################################
  #                               Lissage par bloc
  ##############################################################################
  
  # Application de la fonction de lissage par bloc
  fdobj_list_bloc <- lissage_spline_cubique_bloc(t(as.data.frame(trimmed_rows[-1])), val[1:length(trimmed_rows[-1])], D)
  
  for (i in seq_along(fdobj_list_bloc)) {
    fdobj_list[nrow(fdobj_list)+1] = fdobj_list_bloc[[i]]
  }
}