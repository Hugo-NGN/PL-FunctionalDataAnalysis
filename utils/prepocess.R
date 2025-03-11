# preprocess les donnees pour le lissage

preprocess <- function(data){
  taille_courbes <- apply(data, MARGIN = 1, function(x) which(is.na(x))[1]) - 1
  
  data <- data[taille_courbes == 48, ]
  
  data <- data[1:100, ]
  
  data <- data[,-1]
  data <- data[, -48]
  
  data <- as.data.frame(lapply(data, as.numeric))

  
  return(data)
}
