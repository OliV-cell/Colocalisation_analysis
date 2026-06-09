#################################
################################
#Gaph des combinaison avec CTAC
###############################
###############################
#activation des packages

library(car)
library(multcompView)
library(dplyr)
library(agricolae)
library(multcomp)
library(rcompanion)
library(tidyr)
library(nlstools)
library(psych)
library(qpcR)
library(emmeans)
library(multcomp)
library(multcompView)

library(gamm4)
library(gvlma)
library(stargazer)
library(performance)
library(see)
#représentation graphique
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(pammtools)



decompte <- function(data_source, data_puit){
  for (i in c(3:12)){
    #null
    data_puit[[i-2,1]] <- length(which(data_source[,i] == "rien"))
    
    #faible
    data_puit[[i-2,2]] <- length(which(data_source[,i] == "faible"))
    
    #moyen
    data_puit[[i-2,3]] <- length(which(data_source[,i] == "moyen"))
    
    #elever
    data_puit[[i-2,4]] <- length(which(data_source[,i] == "elever"))
  }
  
  return(data_puit)
  
}

data_CTAC_aveccaac_antisens <- read.csv("CTAC_avec_caac.csv",
                                        sep=";", dec=",", header=FALSE)

data_CTAC_avecctat_antisens <- read.csv("CTAC_avec_ctat.csv",
                                        sep=";", dec=",", header=FALSE)

data_CTAC_avecctat_antisens_rep3 <- read.csv("CTAC_avec_ctat_rep3.csv",
                                             sep=";", dec=",", header=FALSE)

data_CTAC_avecttac_antisens <- read.csv("CTAC_avec_ttac.csv",
                                        sep=";", dec=",", header=FALSE)

data_CTAC_avecttac_antisens_rep2 <- read.csv("CTAC_avec_ttac_rep2.csv",
                                             sep=";", dec=",", header=FALSE)

data_CTAC_avectat_sens <- read.csv("CTAC_sens_avec_ctat_sens.csv",
                                   sep=";", dec=",", header=FALSE)

data_CTAC_avectat_sens_rep1 <- read.csv("CTAC_sens_avec_ctat_sens_rep1.csv",
                                        sep=";", dec=",", header=FALSE)

data_TTACavecCTAC_antisens <- read.csv("TTAC_avec_CTAC.csv",
                                       sep=";", dec=",", header=FALSE)

data_TTACavecCTAC_antisens_rep2 <- read.csv("TTAC_avec_CTAC_rep2.csv",
                                            sep=";", dec=",", header=FALSE)

data_CTATavecCTAC_antisens <- read.csv("CTAT_avec_CTAC.csv",
                                       sep=";", dec=",", header=FALSE)

data_CTATavecCTAC_antisens_rep3 <- read.csv("CTAT_avec_CTAC_rep3.csv",
                                            sep=";", dec=",", header=FALSE)

data_CTAT_avecCTAC_sens <- read.csv("CTAT_sens_avec_CTAC_sens.csv",
                                    sep=";", dec=",", header=FALSE)

data_CTAT_avecCTAC_sens_rep1 <- read.csv("CTAT_sens_avec_CTAC_sens_rep1.csv",
                                         sep=";", dec=",", header=FALSE)

data_CAACavecCTAC_antisens <- read.csv("CAAC_avec_CTAC.csv",
                                       sep=";", dec=",", header=FALSE)

#data_CAACavecCTAC_antisens <- read.csv("Analyse_photo_ovaire_R_ttac.csv",
#sep=";", dec=",", header=FALSE)

CTAC_aveccaac_antisens <- matrix( nrow = 10, ncol = 4)

CTAC_avecttac_antisens <- matrix( nrow = 10, ncol = 4)

CTAC_avecttac_antisens_rep2 <- matrix( nrow = 10, ncol = 4)

CTAC_avecctat_antisens <- matrix( nrow = 10, ncol = 4)

CTAC_avecctat_antisens_rep3 <- matrix( nrow = 10, ncol = 4)

TTACavecCTAC_antisens <- matrix( nrow = 10, ncol = 4)

TTACavecCTAC_antisens_rep2 <- matrix( nrow = 10, ncol = 4)

CTATavecCTAC_antisens <- matrix( nrow = 10, ncol = 4)

CTATavecCTAC_antisens_rep3 <- matrix( nrow = 10, ncol = 4)

CAACavecCTAC_antisens <- matrix( nrow = 10, ncol = 4)

CTAC_avecctat_sens <- matrix( nrow = 10, ncol = 4)

CTAC_avecctat_sens_rep1 <- matrix( nrow = 10, ncol = 4)

CTAT_avecCTAC_sens <- matrix( nrow = 10, ncol = 4)

CTAT_avecCTAC_sens_rep1 <- matrix( nrow = 10, ncol = 4)

#Décompte des occurences de différente intensitée

CTAC_aveccaac_antisens <- decompte(data_source = data_CTAC_aveccaac_antisens,data_puit = CTAC_aveccaac_antisens)

CTAC_avecttac_antisens <- decompte(data_source = data_CTAC_avecttac_antisens,data_puit = CTAC_avecttac_antisens)

CTAC_avecttac_antisens_rep2 <- decompte(data_source = data_CTAC_avecttac_antisens_rep2,data_puit = CTAC_avecttac_antisens_rep2)

CTAC_avecctat_antisens <- decompte(data_source = data_CTAC_avecctat_antisens,data_puit = CTAC_avecctat_antisens)

CTAC_avecctat_antisens_rep3 <- decompte(data_source = data_CTAC_avecctat_antisens_rep3,data_puit = CTAC_avecctat_antisens_rep3)

TTACavecCTAC_antisens <- decompte(data_source = data_TTACavecCTAC_antisens,data_puit = TTACavecCTAC_antisens )

TTACavecCTAC_antisens_rep2 <- decompte(data_source = data_TTACavecCTAC_antisens_rep2,data_puit = TTACavecCTAC_antisens_rep2 )

CTATavecCTAC_antisens <- decompte(data_source = data_CTATavecCTAC_antisens,data_puit = CTATavecCTAC_antisens)

CTATavecCTAC_antisens_rep3 <- decompte(data_source = data_CTATavecCTAC_antisens_rep3,data_puit = CTATavecCTAC_antisens_rep3)

CAACavecCTAC_antisens <- decompte(data_source = data_CAACavecCTAC_antisens,data_puit = CAACavecCTAC_antisens)

CTAC_avecctat_sens <- decompte(data_source = data_CTAC_avectat_sens,data_puit = CTAC_avecctat_sens)

CTAC_avecctat_sens_rep1 <- decompte(data_source = data_CTAC_avectat_sens_rep1,data_puit = CTAC_avecctat_sens_rep1)

CTAT_avecCTAC_sens <- decompte(data_source = data_CTAT_avecCTAC_sens ,data_puit = CTAT_avecCTAC_sens)

CTAT_avecCTAC_sens_rep1 <- decompte(data_source = data_CTAT_avecCTAC_sens_rep1 ,data_puit = CTAT_avecCTAC_sens_rep1)

noms_matrices <- c("CTAC_aveccaac_antisens", "CTAC_avecttac_antisens", "CTAC_avecttac_antisens_rep2",
                   "CTAC_avecctat_antisens", "CTAC_avecctat_antisens_rep3")

total_observations <- sum(sapply(mget(noms_matrices), sum, na.rm = TRUE))

print(total_observations)

matrice_fusionnee <- Reduce("+", mget(noms_matrices))

comptage_par_ligne <- rowSums(matrice_fusionnee, na.rm = TRUE)

noms_stades <- paste("Stade", 1:(1 + length(comptage_par_ligne) - 1))
names(comptage_par_ligne) <- noms_stades


print(comptage_par_ligne)


#renome les colone

colnames(CTAC_aveccaac_antisens) <- c("null", "faible", "moyen", "élevé")
colnames(CTAC_avecttac_antisens) <- c("null", "faible", "moyen", "élevé")
colnames(CTAC_avecttac_antisens_rep2) <- c("null", "faible", "moyen", "élevé")
colnames(CTAC_avecctat_antisens) <- c("null", "faible", "moyen", "élevé")
colnames(CTAC_avecctat_antisens_rep3) <- c("null", "faible", "moyen", "élevé")
colnames(TTACavecCTAC_antisens) <- c("null", "faible", "moyen", "élevé")
colnames(TTACavecCTAC_antisens_rep2) <- c("null", "faible", "moyen", "élevé")
colnames(CTATavecCTAC_antisens) <- c("null", "faible", "moyen", "élevé")
colnames(CTATavecCTAC_antisens_rep3) <- c("null", "faible", "moyen", "élevé")
colnames(CAACavecCTAC_antisens) <- c("null", "faible", "moyen", "élevé")
colnames(CTAC_avecctat_sens) <- c("null", "faible", "moyen", "élevé")
colnames(CTAT_avecCTAC_sens) <- c("null", "faible", "moyen", "élevé")
colnames(CTAC_avecctat_sens_rep1) <- c("null", "faible", "moyen", "élevé")
colnames(CTAT_avecCTAC_sens_rep1) <- c("null", "faible", "moyen", "élevé")

#Boucle de calcule du score pour chaque stade et chaque satellite

score_occurence_data_combinaison <- data.frame(score_CTAC_aveccaac_antisens = rep(0,10),
                                               
                                               score_CTAC_avecttac_antisens = rep(0,10),
                                               
                                               score_CTAC_avecctat_antisens = rep(0,10),
                                               
                                               score_TTACavecCTAC_antisens = rep(0,10),
                                               
                                               score_CTATavecCTAC_antisens = rep(0,10),
                                               
                                               score_CAACavecCTAC_antisens = rep(0,10), 
                                               
                                               score_CTAC_avecctat_sens = rep(0,10), 
                                               
                                               score_CTAT_avecCTAC_sens = rep(0,10),
                                               
                                               score_CTATavecCTAC_antisens_rep3 = rep(0,10), 
                                               
                                               score_CTAC_avecctat_antisens_rep3 = rep(0,10),
                                               
                                               score_CTAC_avecttat_antisens_rep2 = rep(0,10),
                                               
                                               score_TTACavecCTAC_antisens_rep2 = rep(0,10), 
                                               
                                               score_CTAC_avecctat_sens_rep1 = rep(0,10), 
                                               
                                               score_CTAT_avecCTAC_sens_rep1 = rep(0,10),
                                               
                                               stade = c(1:10),
                                               
                                               row.names = c("stade 1","stade 2","stade 3",
                                                             "stade 4","stade 5","stade 6",
                                                             "stade 7","stade 8","stade 9",
                                                             "stade 10"))


for ( i in c(1:10)){
  
  score_occurence_data_combinaison[[i,1]] <- (CTAC_aveccaac_antisens[i,2] + 2*CTAC_aveccaac_antisens[i,3] + 3*CTAC_aveccaac_antisens[i,4])
  
  score_occurence_data_combinaison[[i,2]] <- (CTAC_avecttac_antisens[i,2] + 2*CTAC_avecttac_antisens[i,3] + 3*CTAC_avecttac_antisens[i,4]) 
  
  score_occurence_data_combinaison[[i,3]] <- (CTAC_avecctat_antisens[i,2] + 2*CTAC_avecctat_antisens[i,3] + 3*CTAC_avecctat_antisens[i,4])
  
  score_occurence_data_combinaison[[i,4]] <- (TTACavecCTAC_antisens[i,2] + 2*TTACavecCTAC_antisens[i,3] + 3*TTACavecCTAC_antisens[i,4])
  
  score_occurence_data_combinaison[[i,5]] <- (CTATavecCTAC_antisens[i,2] + 2*CTATavecCTAC_antisens[i,3] + 3*CTATavecCTAC_antisens[i,4])
  
  score_occurence_data_combinaison[[i,6]] <- (CAACavecCTAC_antisens[i,2] + 2*CAACavecCTAC_antisens[i,3] + 3*CAACavecCTAC_antisens[i,4])
  
  score_occurence_data_combinaison[[i,7]] <- (CTAC_avecctat_sens[i,2] + 2*CTAC_avecctat_sens[i,3] + 3*CTAC_avecctat_sens[i,4])
  
  score_occurence_data_combinaison[[i,8]] <- (CTAT_avecCTAC_sens[i,2] + 2*CTAT_avecCTAC_sens[i,3] + 3*CTAT_avecCTAC_sens[i,4])
  
  score_occurence_data_combinaison[[i,9]] <- (CTATavecCTAC_antisens_rep3[i,2] + 2*CTATavecCTAC_antisens_rep3[i,3] + 3*CTATavecCTAC_antisens_rep3[i,4])
  
  score_occurence_data_combinaison[[i,10]] <- (CTAC_avecctat_antisens_rep3[i,2] + 2*CTAC_avecctat_antisens_rep3[i,3] + 3*CTAC_avecctat_antisens_rep3[i,4])
  
  score_occurence_data_combinaison[[i,11]] <- (CTAC_avecttac_antisens_rep2[i,2] + 2*CTAC_avecttac_antisens_rep2[i,3] + 3*CTAC_avecttac_antisens_rep2[i,4]) 
  
  score_occurence_data_combinaison[[i,12]] <- (TTACavecCTAC_antisens_rep2[i,2] + 2*TTACavecCTAC_antisens_rep2[i,3] + 3*TTACavecCTAC_antisens_rep2[i,4])
  
  score_occurence_data_combinaison[[i,13]] <- (CTAC_avecctat_sens[i,2] + 2*CTAC_avecctat_sens[i,3] + 3*CTAC_avecctat_sens[i,4])
  
  score_occurence_data_combinaison[[i,14]] <- (CTAT_avecCTAC_sens[i,2] + 2*CTAT_avecCTAC_sens[i,3] + 3*CTAT_avecCTAC_sens[i,4])
  
  print(paste("Boucle",i,"fait"))
}

score_occurence_data_normaliser_combinaison <- score_occurence_data_combinaison

for (i in c(1:10)){
  
  score_occurence_data_normaliser_combinaison[[i,1]] <- score_occurence_data_normaliser_combinaison[i,1]/(sum(CTAC_aveccaac_antisens[i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[[i,2]] <- score_occurence_data_normaliser_combinaison[i,2]/(sum(CTAC_avecttac_antisens[i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[[i,3]] <- score_occurence_data_normaliser_combinaison[i,3]/(sum(CTAC_avecctat_antisens[i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[[i,4]] <- score_occurence_data_normaliser_combinaison[i,4]/(sum(TTACavecCTAC_antisens[i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[[i,5]] <- score_occurence_data_normaliser_combinaison[i,5]/(sum(CTATavecCTAC_antisens[i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[[i,6]] <- score_occurence_data_normaliser_combinaison[i,6]/(sum(CAACavecCTAC_antisens[i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[[i,7]] <- score_occurence_data_normaliser_combinaison[i,7]/(sum(CTAC_avecctat_sens[i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[[i,8]] <- score_occurence_data_normaliser_combinaison[i,8]/(sum(CTAT_avecCTAC_sens[i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[[i,9]] <- score_occurence_data_normaliser_combinaison[i,9]/(sum(CTATavecCTAC_antisens_rep3[i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[[i,10]] <- score_occurence_data_normaliser_combinaison[i,10]/(sum(CTAC_avecctat_antisens_rep3[i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[[i,11]] <- score_occurence_data_normaliser_combinaison[i,11]/(sum(CTAC_avecttac_antisens_rep2[i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[[i,12]] <- score_occurence_data_normaliser_combinaison[i,12]/(sum(TTACavecCTAC_antisens_rep2[i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[[i,13]] <- score_occurence_data_normaliser_combinaison[i,13]/(sum(CTAC_avecctat_sens_rep1[i,1:4]) +1)
  
  score_occurence_data_normaliser_combinaison[[i,14]] <- score_occurence_data_normaliser_combinaison[i,14]/(sum(CTAT_avecCTAC_sens_rep1[i,1:4]) +1)
  
}

#############################
#compilationde tout les CTAC
#CTAT sens et antisens
############################

CTAC_antisens_compiler <- data.frame(Null = rep(0,10),
                                     
                                     Faible = rep(0,10),
                                     
                                     Moyen  = rep(0,10),
                                     
                                     Élever = rep(0,10), 
                                     
                                     row.names = c("stade 1","stade 2","stade 3",
                                                   "stade 4","stade 5","stade 6",
                                                   "stade 7","stade 8","stade 9",
                                                   "stade 10"))

CTAT_antisens_compiler <- data.frame(Null = rep(0,10),
                                     
                                     Faible = rep(0,10),
                                     
                                     Moyen  = rep(0,10),
                                     
                                     Élever = rep(0,10), 
                                     
                                     row.names = c("stade 1","stade 2","stade 3",
                                                   "stade 4","stade 5","stade 6",
                                                   "stade 7","stade 8","stade 9",
                                                   "stade 10"))

CTAC_sens_compiler <- data.frame(Null = rep(0,10),
                                 
                                 Faible = rep(0,10),
                                 
                                 Moyen  = rep(0,10),
                                 
                                 Élever = rep(0,10), 
                                 
                                 row.names = c("stade 1","stade 2","stade 3",
                                               "stade 4","stade 5","stade 6",
                                               "stade 7","stade 8","stade 9",
                                               "stade 10"))

CTAT_sens_compiler <- data.frame(Null = rep(0,10),
                                 
                                 Faible = rep(0,10),
                                 
                                 Moyen  = rep(0,10),
                                 
                                 Élever = rep(0,10), 
                                 
                                 row.names = c("stade 1","stade 2","stade 3",
                                               "stade 4","stade 5","stade 6",
                                               "stade 7","stade 8","stade 9",
                                               "stade 10"))

score_occurence_data_CTAC_antisens_compiler <- data.frame(Score_CTAC_antisens_compiler = rep(0,10), 
                                                          
                                                          stade = c(1:10),
                                                          
                                                          row.names = c("stade 1","stade 2","stade 3",
                                                                        "stade 4","stade 5","stade 6",
                                                                        "stade 7","stade 8","stade 9",
                                                                        "stade 10"))

score_occurence_data_CTAT_antisens_compiler <- data.frame(Score_CTAT_antisens_compiler = rep(0,10), 
                                                          
                                                          stade = c(1:10),
                                                          
                                                          row.names = c("stade 1","stade 2","stade 3",
                                                                        "stade 4","stade 5","stade 6",
                                                                        "stade 7","stade 8","stade 9",
                                                                        "stade 10"))

score_occurence_data_CTAT_sens_compiler <- data.frame(Score_CTAT_sens_compiler = rep(0,10), 
                                                      
                                                      stade = c(1:10),
                                                      
                                                      row.names = c("stade 1","stade 2","stade 3",
                                                                    "stade 4","stade 5","stade 6",
                                                                    "stade 7","stade 8","stade 9",
                                                                    "stade 10"))

score_occurence_data_CTAC_sens_compiler <- data.frame(Score_CTAC_sens_compiler = rep(0,10), 
                                                      
                                                      stade = c(1:10),
                                                      
                                                      row.names = c("stade 1","stade 2","stade 3",
                                                                    "stade 4","stade 5","stade 6",
                                                                    "stade 7","stade 8","stade 9",
                                                                    "stade 10"))

for (i in c(1:4)) {
  for (j in c(1:10)) {
    
    CTAC_antisens_compiler[[j,i]] <- mean(c(CTAC_aveccaac_antisens[j,i], CTAC_avecctat_antisens[j,i], 
                                         CTAC_avecttac_antisens[j,i],CTAC_avecctat_antisens_rep3[j,i], 
                                         CTAC_avecttac_antisens_rep2[j,i]))
    
    CTAT_antisens_compiler[[j,i]] <- mean(c(CTATavecCTAC_antisens[j,i], CTATavecCTAC_antisens_rep3[j,i]))
    
    CTAC_sens_compiler[[j,i]] <- mean(c(CTAC_avecctat_sens[j,i], CTAC_avecctat_sens_rep1[j,i]))
    
    CTAT_sens_compiler[[j,i]] <- mean(c(CTAT_avecCTAC_sens[j,i], CTAT_avecCTAC_sens_rep1[j,i]))
    
  }
}

for (i in c(1:10)){
  
  score_occurence_data_CTAC_antisens_compiler[[i,1]] <- (CTAC_antisens_compiler[i,2] + 2*CTAC_antisens_compiler[i,3] + 3*CTAC_antisens_compiler[i,4])
  
  score_occurence_data_CTAT_antisens_compiler[[i,1]] <- (CTAT_antisens_compiler[i,2] + 2*CTAT_antisens_compiler[i,3] + 3*CTAT_antisens_compiler[i,4])
  
  score_occurence_data_CTAC_sens_compiler[[i,1]] <- (CTAC_sens_compiler[i,2] + 2*CTAC_sens_compiler[i,3] + 3*CTAC_sens_compiler[i,4])
  
  score_occurence_data_CTAT_sens_compiler[[i,1]] <- (CTAT_sens_compiler[i,2] + 2*CTAT_sens_compiler[i,3] + 3*CTAT_sens_compiler[i,4])
  
}

score_occurence_data_CTAC_antisens_compiler_normaliser <- score_occurence_data_CTAC_antisens_compiler

score_occurence_data_CTAT_antisens_compiler_normaliser <- score_occurence_data_CTAT_antisens_compiler

score_occurence_data_CTAC_sens_compiler_normaliser <- score_occurence_data_CTAC_sens_compiler

score_occurence_data_CTAT_sens_compiler_normaliser <- score_occurence_data_CTAT_sens_compiler

for (i in c(1:10)) {
  
  score_occurence_data_CTAC_antisens_compiler_normaliser[[i,1]] <- score_occurence_data_CTAC_antisens_compiler_normaliser[i,1]/(sum(CTAC_antisens_compiler[i,1:4]) +1)
  
  score_occurence_data_CTAT_antisens_compiler_normaliser[[i,1]] <- score_occurence_data_CTAT_antisens_compiler_normaliser[i,1]/(sum(CTAT_antisens_compiler[i,1:4]) +1)
  
  score_occurence_data_CTAC_sens_compiler_normaliser[[i,1]] <- score_occurence_data_CTAC_sens_compiler_normaliser[i,1]/(sum(CTAC_sens_compiler[i,1:4]) +1)
  
  score_occurence_data_CTAT_sens_compiler_normaliser[[i,1]] <- score_occurence_data_CTAT_sens_compiler_normaliser[i,1]/(sum(CTAT_sens_compiler[i,1:4]) +1)
  
}


CTAC_sens_compiler_moyenne <- data.frame(moyenne = rep(0,10),
                                      
                                      SD = rep(0,10),
                                 
                                      stade = c(1:10),
                                 
                                      row.names = c("stade 1","stade 2","stade 3",
                                               "stade 4","stade 5","stade 6",
                                               "stade 7","stade 8","stade 9",
                                               "stade 10"))

CTAT_sens_compiler_moyenne <- data.frame(moyenne = rep(0,10),
                                      
                                      SD = rep(0,10),
                                      
                                      stade = c(1:10),
                                      
                                      row.names = c("stade 1","stade 2","stade 3",
                                                    "stade 4","stade 5","stade 6",
                                                    "stade 7","stade 8","stade 9",
                                                    "stade 10"))

CTAT_antisens_compiler_moyenne <- data.frame(moyenne = rep(0,10),
                                          
                                          SD = rep(0,10),
                                      
                                          stade = c(1:10),
                                      
                                          row.names = c("stade 1","stade 2","stade 3",
                                                    "stade 4","stade 5","stade 6",
                                                    "stade 7","stade 8","stade 9",
                                                    "stade 10"))

CTAC_antisens_compiler_moyenne <- data.frame(moyenne = rep(0,10),
                                          
                                          SD = rep(0,10),
                                          
                                          stade = c(1:10),
                                          
                                          row.names = c("stade 1","stade 2","stade 3",
                                                        "stade 4","stade 5","stade 6",
                                                        "stade 7","stade 8","stade 9",
                                                        "stade 10"))



for (j in c(1:10)) {
    
    CTAC_sens_compiler_moyenne[[j,1]] <- mean(c(score_occurence_data_normaliser_combinaison$score_CTAC_avecctat_sens[j], 
                                             score_occurence_data_normaliser_combinaison$score_CTAC_avecctat_sens_rep1[j]))
    
    CTAC_sens_compiler_moyenne[[j,2]] <- sd(c(score_occurence_data_normaliser_combinaison$score_CTAC_avecctat_sens[j], 
                                             score_occurence_data_normaliser_combinaison$score_CTAC_avecctat_sens_rep1[j]))
    
    CTAT_sens_compiler_moyenne[[j,1]] <- mean(c(score_occurence_data_normaliser_combinaison$score_CTAT_avecCTAC_sens[j], 
                                             score_occurence_data_normaliser_combinaison$score_CTAT_avecCTAC_sens_rep1[j]))
    
    CTAT_sens_compiler_moyenne[[j,2]] <- sd(c(score_occurence_data_normaliser_combinaison$score_CTAT_avecCTAC_sens[j], 
                                             score_occurence_data_normaliser_combinaison$score_CTAT_avecCTAC_sens_rep1[j]))
    
    CTAT_antisens_compiler_moyenne[[j,1]] <- mean(c(score_occurence_data_normaliser_combinaison$score_CTATavecCTAC_antisens[j], 
                                             score_occurence_data_normaliser_combinaison$score_CTATavecCTAC_antisens_rep3[j]))
    
    CTAT_antisens_compiler_moyenne[[j,2]] <- sd(c(score_occurence_data_normaliser_combinaison$score_CTATavecCTAC_antisens[j], 
                                                 score_occurence_data_normaliser_combinaison$score_CTATavecCTAC_antisens_rep3[j]))
    
    CTAC_antisens_compiler_moyenne[[j,1]] <- mean(c(score_occurence_data_normaliser_combinaison$score_CTAC_aveccaac_antisens[j],
                                                 score_occurence_data_normaliser_combinaison$score_CTAC_avecttac_antisens[j],
                                                 score_occurence_data_normaliser_combinaison$score_CTAC_avecctat_antisens[j],
                                                 score_occurence_data_normaliser_combinaison$score_CTAC_avecctat_antisens_rep3[j],
                                                 score_occurence_data_normaliser_combinaison$score_CTAC_avecttat_antisens_rep2[j]))
    
    CTAC_antisens_compiler_moyenne[[j,2]] <- sd(c(score_occurence_data_normaliser_combinaison$score_CTAC_aveccaac_antisens[j],
                                                 score_occurence_data_normaliser_combinaison$score_CTAC_avecttac_antisens[j],
                                                 score_occurence_data_normaliser_combinaison$score_CTAC_avecctat_antisens[j],
                                                 score_occurence_data_normaliser_combinaison$score_CTAC_avecctat_antisens_rep3[j],
                                                 score_occurence_data_normaliser_combinaison$score_CTAC_avecttat_antisens_rep2[j]))
    
    
}




########################
#Graphique avec GGplot
######################


transparence = 0.2

taile_ligne = 0.8


windowsFonts( A = windowsFont("Arial"))

graph <- ggplot(score_occurence_data_normaliser_combinaison,
                aes(x = stade)) +
  
  #ligne des CTAC
  
  geom_line(aes(y = score_CTAC_aveccaac_antisens,
                color = "Forward AAACTAC"), size = taile_ligne, alpha = transparence,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTAC_avecttac_antisens,
                color = "Forward AAACTAC"), size = taile_ligne,alpha = transparence,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTAC_avecttat_antisens_rep2,
                color = "Forward AAACTAC"), size = taile_ligne,alpha = transparence,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTAC_avecctat_antisens,
                color = "Forward AAACTAC"), size = taile_ligne,alpha = transparence,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTAC_avecctat_antisens_rep3,
                color = "Forward AAACTAC"), size = taile_ligne,alpha = transparence,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTAC_avecctat_sens,
                color = "Reverse AAACTAC"), size = taile_ligne, alpha = transparence,
            linetype = "dotdash") +
  
  geom_line(aes(y = score_CTAC_avecctat_sens_rep1,
                color = "Reverse AAACTAC"), size = taile_ligne, alpha = transparence,
            linetype = "dotdash") +
  
  #ligne des CTAT
  
  geom_line(aes(y = score_CTATavecCTAC_antisens,
                color = "Forward AAACTAT"), size = taile_ligne, alpha = transparence,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTATavecCTAC_antisens_rep3,
                color = "Forward AAACTAT"), size = taile_ligne, alpha = transparence,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTAT_avecCTAC_sens,
                color = "Reverse AAACTAT"), size = taile_ligne, alpha = transparence,
            linetype = "dotdash") +
  
  geom_line(aes(y = score_CTAT_avecCTAC_sens_rep1,
                color = "Reverse AAACTAT"), size = taile_ligne, alpha = transparence,
            linetype = "dotdash") +
  
  #Ligne des compilation
  
  geom_line(data = CTAC_antisens_compiler_moyenne,
            aes(y = moyenne,
                color = "Compilation foward AAACTAC"),
            linetype = "solid", size = 1.5) +
  
  geom_line(data = CTAT_antisens_compiler_moyenne,
            aes(y = moyenne,
                color = "Compilation foward AAACTAT"),
            linetype = "solid", size = 1.5) +
  
  geom_line(data = CTAC_sens_compiler_moyenne,
            aes(y = moyenne,
                color = "Compilation reverse AAACTAC"),
            linetype = "dotdash", size = 1.5) +
  
  geom_line(data = CTAT_sens_compiler_moyenne,
            aes(y = moyenne,
                color = "Compilation reverse AAACTAT"),
            linetype = "dotdash", size = 1.5) +
  
  labs(
    x = "Egg chamber stages",
    y = "Normalized expression score",
    title = "The lncRNA expression levels of satellite DNA sequences \n across egg chamber development in Drosophila virllis",
    subtitle = "The experiments multiplexed the AAACTAC probe with the other related probes",
    color = "LncRNA strand"
  ) +
  
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  ylim(0, 4) +
  
  scale_color_manual(values = c(
    "Forward AAACTAC" = "limegreen",
    "Reverse AAACTAC" = "limegreen",
    "Forward AAATTAC"  = "yellow",
    "Forward AAACTAT"  = "magenta",
    "Reverse AAACTAT"  = "magenta",
    "Forward AAACAAC"  = "blue",
    "Compilation foward AAACTAC" = "limegreen",
    "Compilation foward AAACTAT" = "magenta",
    "Compilation reverse AAACTAC" = "limegreen",
    "Compilation reverse AAACTAT" = "magenta"
  )) +
  
  theme_bw() +
  theme(
    legend.position = "right",
    legend.background = element_rect(fill = FALSE, color = FALSE), 
    legend.title = element_text(size = 12, face = "bold",family = "A"),
    legend.text = element_text(size = 9),
    text = element_text(size = 12, family = "A"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "A"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain", family = "A"),
    plot.margin = margin(20, 20, 20, 20)
  )

print(graph)


#############################################


#################################
graph_personalisable <- ggplot(score_occurence_data_normaliser_combinaison,
                               aes(x = stade)) +
  
  geom_line(aes(y = score_CTAC_avecctat_sens,
                color = "Reverse AAACTAC"), size = 0.75,
            linetype = "dotdash") +
  
  geom_line(aes(y = score_CTAT_avecCTAC_sens,
                color = "Reverse AAACTAT"), size = 0.75,
            linetype = "dotdash") +
  
  geom_line(data = score_occurence_data_CTAC_antisens_compiler_normaliser,
            aes(y = Score_CTAC_antisens_compiler,
                color = "Compilation AAACTAC"),
            linetype = "solid", size = 1.3) +
  
  geom_line(data = score_occurence_data_CTAT_antisens_compiler_normaliser,
            aes(y = Score_CTAT_antisens_compiler,
                color = "Compilation AAACTAT"),
            linetype = "solid", size = 1.3) +
  
  labs(
    x = "Egg chamber stages",
    y = "Normalized expression score",
    title = "The lncRNA expression levels of satellite DNA sequences \n across egg chamber development in Drosophila virllis",
    subtitle = "The experiments multiplexed the AAACTAC probe with the other related probes",
    color = "LncRNA strand"
  ) +
  
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  ylim(0, 4) +
  
  scale_color_manual(values = c(
    "Forward AAACTAC" = "limegreen",
    "Reverse AAACTAC" = "limegreen",
    "Forward AAATTAC"  = "yellow",
    "Forward AAACTAT"  = "magenta",
    "Reverse AAACTAT"  = "magenta",
    "Forward AAACAAC"  = "blue",
    "Compilation AAACTAC" = "limegreen",
    "Compilation AAACTAT" = "magenta"
  )) +
  
  theme_bw() +
  theme(
    legend.position = "right",
    legend.background = element_rect(fill = FALSE, color = FALSE), 
    legend.title = element_text(size = 10, face = "bold",family = "A"),
    legend.text = element_text(size = 9),
    text = element_text(size = 12, family = "A"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "A"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain", family = "A"),
    plot.margin = margin(20, 20, 20, 20)
  )

print(graph_personalisable)


############################################
#Exploration data
###########################################

par(mfrow=c(4,2), oma = c(0, 0, 3, 0))

for ( i in c(3:10)){
  
    plot(CTAC_aveccaac_antisens[i,],
         xlab= " null = 1, faible = 2, moyen = 3, élever = 4",
         ylab = "occurence des cathégorie",
         pch = 19,
         type = "o",
         main = paste0("Stade",i))
  
  tot <- sum(CTAC_aveccaac_antisens[i,])
  
  text( x = 2.5, y = (0.75*max(CTAC_aveccaac_antisens[i,])), label = paste0("N = ", tot))
  
}

mtext("CTAC_aveccaac_antisens nombre d'obersvation", outer = TRUE, cex = 1.5, font = 2)

par(mfrow=c(4,2), oma = c(0, 0, 3, 0))

for ( i in c(3:10)){
  
  plot(CTAC_avecctat_antisens_rep3[i,],
       xlab= " null = 1, faible = 2, moyen = 3, élever = 4",
       ylab = "occurence des cathégorie",
       pch = 19,
       type = "o",
       main = paste0("Stade",i))
  
  tot <- sum(CTAC_avecctat_antisens_rep3[i,])
  
  text( x = 2.5, y = (0.75*max(CTAC_avecctat_antisens_rep3[i,])), label = paste0("N = ", tot))
  
}

mtext("CTAC_avecctat_antisens_rep3 nombre d'obersvation", outer = TRUE, cex = 1.5, font = 2)

par(mfrow=c(4,2), oma = c(0, 0, 3, 0))

for ( i in c(3:10)){
  
  plot(CTAC_avecctat_antisens[i,],
       xlab= " null = 1, faible = 2, moyen = 3, élever = 4",
       ylab = "occurence des cathégorie",
       pch = 19,
       type = "o",
       main = paste0("Stade",i))
  
  tot <- sum(CTAC_avecctat_antisens[i,])
  
  text( x = 2.5, y = (0.75*max(CTAC_avecctat_antisens[i,])), label = paste0("N = ", tot))
  
}

mtext("CTAC_avecctat_antisens nombre d'obersvation", outer = TRUE, cex = 1.5, font = 2)

#GLM distribution Gamma test

#confection des data_frame mutate de sens et antisens CTAC et CTAT

score_CTAC_antisens <- score_occurence_data_normaliser_combinaison[,c(1,2,3,10,11,15)] %>%
  pivot_longer(
    cols = !stade,
    names_to = "brin",
    values_to = "score"
  )

score_CTAC_sens <- score_occurence_data_normaliser_combinaison[,c(7,13,15)] %>%
  pivot_longer(
    cols = !stade,
    names_to = "brin",
    values_to = "score"
  )

score_CTAT_antisens <- score_occurence_data_normaliser_combinaison[,c(5,9,15)] %>%
  pivot_longer(
    cols = !stade,
    names_to = "brin",
    values_to = "score"
  )

score_CTAT_sens <- score_occurence_data_normaliser_combinaison[,c(8,14,15)] %>%
  pivot_longer(
    cols = !stade,
    names_to = "brin",
    values_to = "score"
  )

hist(score_CTAC_antisens$score)
hist(score_CTAC_sens$score)
hist(score_CTAT_antisens$score)
hist(score_CTAT_sens$score)

boxplot(score_CTAC_antisens$score ~ score_CTAC_antisens$stade)
boxplot(score_CTAC_sens$score ~ score_CTAC_sens$stade)
boxplot(score_CTAT_antisens$score ~ score_CTAT_antisens$stade)
boxplot(score_CTAT_sens$score ~ score_CTAT_sens$stade)

concatene <- qpcR:::cbind.na(score_CTAC_antisens[,3],score_CTAC_sens[,3],
                              score_CTAT_antisens[,3],score_CTAT_sens[,3])

pairs.panels(concatene)



model_CTAC_antisens <- lm(data = score_CTAC_antisens, score ~ stade)
model_CTAC_sens <- lm(data = score_CTAC_sens, score ~ stade)
model_CTAT_antisens <- lm(data = score_CTAT_antisens, score ~ stade)
model_CTAT_sens <- lm(data = score_CTAT_sens, score ~ stade)


gvlma(model_CTAC_antisens)
check_model(model_CTAC_antisens)

gvlma(model_CTAC_sens)
check_model(model_CTAC_sens)


gvlma(model_CTAT_antisens)
check_model(model_CTAT_antisens)


gvlma(model_CTAT_sens)
check_model(model_CTAT_sens)


#vérification de modele gam ainsi que de leurs condition application
#####################################
#Faire test de comparaions des courbes

score_CTAC_antisens$brin <- "CTAC_antisens" 

score_CTAC_sens$brin <- "CTAC_sens" 

score_CTAT_antisens$brin <- "CTAT_antisens" 

score_CTAT_sens$brin <- "CTAT_sens" 

data_combier <- rbind(score_CTAC_antisens, score_CTAC_sens,score_CTAT_antisens,score_CTAT_sens)

#exportation des donner

data_combier$brin <- as.factor(data_combier$brin)

modele_combiné_tw <- gam(score ~ brin + s(stade, by = as.factor(brin), bs = "cr", k = 10), 
                      data = data_combier, 
                      method = "REML", family = tw(link = "log"),control = list(
                        maxit = 1000))

modele_combiné <- gam(list(score ~ brin + s(stade, by = brin, bs = "cr", k = 10), 
                           ~ s(stade)), 
                      data = data_combier,
                      optimizer = c("outer", "newton"),
                      method = "REML", 
                      family = gaulss(),
                      control = list(maxit = 1000)
)

summary(modele_combiné)
par(mfrow=c(2,2))
gam.check(modele_combiné)
k.check(modele_combiné)

concurvity(modele_combiné, full = FALSE)


model_performance(modele_combiné)
plot(modele_combiné)

anova <- anova.gam(modele_combiné)


stade_dense <- seq(min(data_combier$stade), 
                   max(data_combier$stade), 
                   length.out = 200)

df_predict_combine <- expand.grid(
  stade = stade_dense,
  brin = unique(data_combier$brin) 
)

preds_smooth <- predict(modele_combiné, newdata = df_predict_combine, 
                        se.fit = TRUE, type = "link")

ilink <- family(modele_combiné)$linkinv

df_predict_combine <- df_predict_combine %>%
  mutate(
    fit   = preds_smooth$fit[,1],  # Moyenne
    upper = preds_smooth$fit[,1] + (1.96 * preds_smooth$se.fit[,1]),
    lower = preds_smooth$fit[,1] - (1.96 * preds_smooth$se.fit[,1])
  )

# calcule des R^2 par sp

data_combier$preds <- predict(modele_combiné, type = "response")[,1]

tableau_R2 <- data_combier %>%
  group_by(brin) %>%
  summarize(
    RSS = sum((score - preds)^2),                  
    TSS = sum((score - mean(score))^2),           
    R2  = round(1 - (RSS / TSS), 4)               
  )

print(tableau_R2)

noms_legendes <- c(
  "CTAC_antisens" = "Forward AAACTAC",
  "CTAC_sens"     = "Reverse AAACTAC",
  "CTAT_antisens" = "Foward AAACTAT",
  "CTAT_sens"     = "Reverse AAACTAT"
)

p_val_brin <- anova$pTerms.pv[1]

#####################################################

ploT_GAM <- ggplot() + 
  
  geom_point(data = data_combier, 
             aes(x = stade, y = score, colour = brin),
             position = position_jitter(width = 0.05)) + 
  
  geom_ribbon(data = df_predict_combine, 
              aes(x = stade, ymin = lower, ymax = upper, fill = brin), 
              alpha = 0.09) +
  
  geom_line(data = df_predict_combine, 
            aes(x = stade, y = fit, colour = brin, linetype = brin),
            linewidth = 1.2) +
  
  labs(
    x = "Ovariole developmental stages",
    y = "Normalize expression score",
    title = "Transcription evolution of forward and reverse LncRNA of AAACTAC and AAACTAT for D. virilis",
    subtitle = "Used probe were forward and reverse",
    colour = "LncRNA strand",
    linetype = "LncRNA strand" 
  ) +
  scale_color_manual(labels = noms_legendes,
                     values = c("CTAC_antisens" = "limegreen", 
                                "CTAC_sens"     = "limegreen", 
                                "CTAT_antisens" = "magenta", 
                                "CTAT_sens"     = "magenta")) +
  
  scale_fill_manual(labels = noms_legendes,
                    values = c("CTAC_antisens" = "limegreen", 
                               "CTAC_sens"     = "limegreen", 
                               "CTAT_antisens" = "magenta", 
                               "CTAT_sens"     = "magenta")) +
  
  scale_linetype_manual(labels = noms_legendes,
                        values = c("CTAC_antisens" = "solid", 
                                   "CTAC_sens"     = "dotted", 
                                   "CTAT_antisens" = "solid", 
                                   "CTAT_sens"     = "dotted")) +
  
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  ylim(c(-0.5, 5)) + 
  
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.background = element_rect(fill = "transparent", color = NA), 
    legend.title = element_text(size = 12, face = "bold", family = "A"),
    legend.text = element_text(size = 9),
    text = element_text(size = 12, family = "A"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "A"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain", family = "A"),
    plot.margin = margin(20, 20, 20, 20),
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold")
  ) + 
  
  guides(fill = "none") +
  
  annotate("rect", xmin = 6.5, xmax = Inf, ymin = 4.2, ymax = 5,
           alpha = 1, fill = "white") +
  
  #annotate("text", x = Inf, y = Inf, label = paste0("Global ANOVA: p = ", format.pval(p_val_brin, digits = 2)), 
           #hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
  
  annotate("text", x = 8.5, y = 4.9, label = paste0("R² S.CTAC = ", round(tableau_R2[2,4], digits = 3), "***"), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic",
           color = "limegreen") + 
  
  annotate("text", x = 8.5, y = 4.6, label = paste0("R² A.CTAC = ", round(tableau_R2[1,4], digits = 3), "***" ), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic",
           color = "limegreen") +
  
  annotate("text", x = Inf, y = 4.9, label = paste0("R² S.CTAT = ", round(tableau_R2[4,4], digits = 3), "***"), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic",
           color = "magenta") + 
  
  annotate("text", x = Inf, y = 4.6, label = paste0("R² A.CTAT = ", round(tableau_R2[3,4], digits = 3), "***" ), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic",
           color = "magenta")
  

print(ploT_GAM)

data_combiner_CTAC <- rbind(data_combier[data_combier$brin == "CTAC_antisens", ],
                              data_combier[data_combier$brin == "CTAC_sens", ])

df_predict_combiner_CTAC <- rbind(df_predict_combine[df_predict_combine$brin == "CTAC_antisens", ],
                                  df_predict_combine[df_predict_combine$brin == "CTAC_sens", ])


plot_GAM_CTAC_seul <- ggplot() + 
  
  geom_point(data = data_combiner_CTAC, 
             aes(x = stade, y = score, colour = brin),
             position = position_jitter(width = 0.05)) + 
  
  geom_ribbon(data = df_predict_combiner_CTAC, 
              aes(x = stade, ymin = lower, ymax = upper, fill = brin), 
              alpha = 0.09) +
  
  geom_line(data = df_predict_combiner_CTAC, 
            aes(x = stade, y = fit, colour = brin, linetype = brin),
            linewidth = 1.2) +
  
  labs(
    x = "Ovariole developmental stages",
    y = "Normalize expression score",
    title = "Transcription evolution of forward and reverse LncRNA of AAACTAC and AAACTAT for D. virilis",
    subtitle = "Used probe were forward and reverse",
    colour = "LncRNA strand",
    linetype = "LncRNA strand" 
  ) +
  scale_color_manual(labels = noms_legendes,
                     values = c("CTAC_antisens" = "darkgreen", 
                                "CTAC_sens"     = "limegreen")) +
  
  scale_fill_manual(labels = noms_legendes,
                    values = c("CTAC_antisens" = "darkgreen", 
                               "CTAC_sens"     = "limegreen")) +
  
  scale_linetype_manual(labels = noms_legendes,
                        values = c("CTAC_antisens" = "solid", 
                                   "CTAC_sens"     = "dotted")) +
  
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  ylim(c(-0.5, 5)) + 
  
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.background = element_rect(fill = "transparent", color = NA), 
    legend.title = element_text(size = 12, face = "bold", family = "A"),
    legend.text = element_text(size = 12),
    text = element_text(size = 12, family = "A"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "A"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain", family = "A"),
    plot.margin = margin(20, 20, 20, 20),
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold")
  ) + 
  
  guides(fill = "none") +
  
  annotate("rect", xmin = 8.5, xmax = Inf, ymin = 4.2, ymax = 5,
           alpha = 1, fill = "white") +
  
  #annotate("text", x = Inf, y = Inf, label = paste0("Global ANOVA: p = ", format.pval(p_val_brin, digits = 2)), 
  #hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
  
  annotate("text", x = Inf, y = 4.9, label = paste0("R² F.CTAC = ", round(tableau_R2[2,4], digits = 3), "***"), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic",
           color = "darkgreen") + 
  
  annotate("text", x = Inf, y = 4.6, label = paste0("R² R.CTAC = ", round(tableau_R2[1,4], digits = 3), "***" ), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic",
           color = "limegreen") 

print(plot_GAM_CTAC_seul)

############################################################


tableau_pvalues <- data.frame()

for (s in 3:10) {
  donnees_stade <- subset(data_combier, stade == s)
  
  if(length(unique(donnees_stade$brin)) > 1) {
    
    fit <- aov(score ~ brin, data = donnees_stade)
    
    tukey <- TukeyHSD(fit)$brin
    
    temp <- data.frame(
      Stade = s,
      Comparaison = rownames(tukey),
      p_value = round(tukey[, "p adj"], 4)
    )
    tableau_pvalues <- rbind(tableau_pvalues, temp)
  }
}

print(tableau_pvalues)

signif_data <- tableau_pvalues %>%
  filter(p_value < 0.05) %>%
  mutate(
    label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

















