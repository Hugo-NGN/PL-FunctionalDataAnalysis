library(dplyr)
library(ggplot2)


show_localisation <- function(data, cluster, title = "Carte des clusters", data_path) {
  
  pixel <- get_pixel_bloc(data_path)
  
  # Récupère les coordonnées des pixels
  coords <- read.csv("data/coords.csv", sep = ";")
  geo <- get_geo_by_pixel(pixel, coords)
  
  geo[, 4] <- cluster
  
  # Ajustemnt pour la visualisation(aire des rectangles)
  geo <- geo %>%
    mutate(lon_min = lon - 0.01,
           lon_max = lon + 0.01,
           lat_min = lat - 0.01,
           lat_max = lat + 0.01)
  
  ggplot(geo, aes(xmin = lon_min, xmax = lon_max, ymin = lat_min, ymax = lat_max, fill = as.factor(V4))) +
    geom_rect() +
    labs(title = title,
         x = "Longitude",
         y = "Latitude",
         fill = "Cluster") +
    theme_minimal() +
    scale_fill_brewer(type = "qual", palette = "Set1")
}

get_pixel_bloc <- function(data_path, extract_n_data = 1000, bloc = 47){
  data <- read.csv(data_path, sep=";")
  taille_courbes <- apply(data, MARGIN = 1, function(x) which(is.na(x))[1]) - 1
  
  data <- data[taille_courbes == bloc+1, ]
  
  if (!is.null(extract_n_data)){
    data <- data[1:extract_n_data, ]
  }
  pixel <- data[,1]
  return(pixel)
}

library(dplyr)

get_geo_by_pixel <- function(pixels, coords) {
  if (!is.vector(pixels)) {
    stop("pixels doit être un vecteur d'indices.")
  }
  data_filtered <- coords %>% filter(.[[1]] %in% pixels)
  return(data_filtered)
}

