# preprocess les donnees pour le lissage

preprocess <- function(data, extract_n_data = NULL){
  taille_courbes <- apply(data, MARGIN = 1, function(x) which(is.na(x))[1]) - 1
  
  data <- data[taille_courbes == 48, ]
  
  if (!is.null(extract_n_data)){
    data <- data[1:extract_n_data, ]
  }
  
  data <- data[,-1]
  data <- data[, -48]
  
  data <- as.data.frame(lapply(data, as.numeric))

  
  return(data)
}
