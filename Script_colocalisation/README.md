If you just want the important fonction to read the colocalisation score .csv file :
```
Traitement_data_coloc <- function(tableau_data_souce){
  
  #Tableau de donné ou tout seras entreposer
  
  indice_coloc <- data.frame(Coef_pearson = rep(0,1),
                             
                             Coeaf_overlap = rep(0,1),
                             
                             k1  = rep(0,1),
                             
                             k2 = rep(0,1), 
                             
                             M1 = rep(0,1),
                             
                             M2 = rep(0,1),
                             
                             a_cytofluogram  = rep(0,1),
                             
                             b_cytofluogram = rep(0,1), 
                             
                             ICQ = rep(0,1),
                             
                             stringsAsFactors = FALSE)
  
  #indice d'ajout de ligne pour chaque occurence dans le fichier txt
  
  new_row <- rep(0, 9)
  
  #indice pour la boucle initial
  
  n <- 1
  
  
  for ( i in c(1:nrow(tableau_data_souce))){
    
    if (grepl("Pearson's Coefficient:", tableau_data_souce[i,1],fixed = TRUE)){
      
      new_row[1] <- tableau_data_souce[i+1,1]
      
    }
    
    if (grepl("Overlap Coefficient:", tableau_data_souce[i,1],fixed = TRUE)){
      
      new_row[2] <- tableau_data_souce[i+1,1]
      
    }
    
    if (grepl("r^2=k1xk2:", tableau_data_souce[i,1],fixed = TRUE)){
      
      new_row[3] <- tableau_data_souce[i+1,1]
      
      new_row[4] <- tableau_data_souce[i+2,1]
      
    }
    
    if (grepl("Manders' Coefficients (using threshold", tableau_data_souce[i,1], fixed = TRUE)){
      
      new_row[5] <- tableau_data_souce[i+1,1]
      
      new_row[6] <- tableau_data_souce[i+2,1]
      
    }
    
    if (grepl("Cytofluorogram's parameters:", tableau_data_souce[i,1],fixed = TRUE)){
      
      new_row[7] <- tableau_data_souce[i+1,1]
      
      new_row[8] <- tableau_data_souce[i+2,1]
      
    }
    
    if (grepl("ICQ: ",tableau_data_souce[i,1],fixed = TRUE)){
      
      new_row[9] <- tableau_data_souce[i,1]
      
      indice_coloc[n, ] <- new_row
      
      n <- n + 1
      new_row <- rep(0, 9)
      
    }
    
  }
  
  #Nettoyage de tout les caractère non numérique et transfromation en caratère numérique
  
  indice_coloc[] <- lapply(indice_coloc, function(x) {
    as.numeric(sub(".*[:=]\\s*([-0-9eE\\.]+).*", "\\1", x))
  })
    
  return(indice_coloc)
  
}
```
