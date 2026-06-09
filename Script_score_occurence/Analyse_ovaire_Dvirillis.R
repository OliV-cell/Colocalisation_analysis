##########################
#Analyse donnée ovaire
#########################

library(ggplot2)




data_TTAC_antisens <- read.csv("Analyse_photo_ovaire_R_ttac.csv",
                               sep=";", dec=",", header=FALSE)

data_TTAC_sens <- read.csv("Analyse_photo_ovaire_TTAC_sens.csv",
                               sep=";", dec=",", header=FALSE)

data_CTAT_antisens <- read.csv("Analyse_photo_ovaire_CTAT.csv",
                               sep=";", dec=",", header=FALSE)

data_CTAT_sens <- read.csv("Analyse_photo_ovaire_CTAT_sens.csv",
                           sep=";", dec=",", header=FALSE)
  
data_CTAC_sens <- read.csv("Analyse_photo_ovaire_CTAC_sens.csv",
                              sep=";", dec=",", header=FALSE)

data_CTAC_antisens <-read.csv("Analyse_photo_ovaire_CTAC_antisens.csv",
                              sep=";", dec=",", header=FALSE)
  
data_CAAC_antisens <- read.csv("Analyse_photo_ovaire_CAAC_antisens.csv",
                               sep=";", dec=",", header=FALSE)
  
#Calcule du nombre d'occurence par stade pour chaque brin   

TTAC_antisens <- matrix( nrow = 10, ncol = 4)

TTAC_sens <- matrix( nrow = 10, ncol = 4)

CTAC_antisens <-matrix( nrow = 10, ncol = 4)

CTAC_sens <-matrix( nrow = 10, ncol = 4)

CTAC_antisens <-matrix( nrow = 10, ncol = 4)

CTAT_antisens <-matrix( nrow = 10, ncol = 4)

CTAT_sens <-matrix( nrow = 10, ncol = 4)

CAAC_antisens <-matrix( nrow = 10, ncol = 4)

#boucle for afin de faire le décompte

#for (i in c(3:12)){
    #null
   #TTAC_antisens[[i-2,1]] <- length(which(data_TTAC_antisens[,i] == "rien"))
    
    #faible
   #TTAC_antisens[[i-2,2]] <- length(which(data_TTAC_antisens[,i] == "faible"))
    
    #moyen
    #TTAC_antisens[[i-2,3]] <- length(which(data_TTAC_antisens[,i] == "moyen"))
    
    #elever
    #TTAC_antisens[[i-2,4]] <- length(which(data_TTAC_antisens[,i] == "elever"))
#}

#fonction

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

TTAC_antisens <- decompte(data_source = data_TTAC_antisens,data_puit = TTAC_antisens)
TTAC_sens <- decompte(data_source = data_TTAC_sens, data_puit = TTAC_sens)
CTAC_antisens <- decompte(data_source = data_CTAC_antisens, data_puit = CTAC_antisens)
CTAC_sens <- decompte(data_source = data_CTAC_sens, data_puit = CTAC_sens)
CTAT_antisens <- decompte(data_source = data_CTAT_antisens, data_puit = CTAT_antisens)
CTAT_sens <- decompte(data_source = data_CTAT_sens, data_puit = CTAT_sens)
CAAC_antisens <- decompte(data_source = data_CAAC_antisens, data_puit = CAAC_antisens)

colnames(TTAC_antisens) <- c("null", "faible", "moyen", "élevé")
colnames(TTAC_sens) <- c("null", "faible", "moyen", "élevé")
colnames(CTAC_antisens) <- c("null", "faible", "moyen", "élevé")
colnames(CTAC_sens) <- c("null", "faible", "moyen", "élevé")
colnames(CTAC_sens) <- c("null", "faible", "moyen", "élevé")
colnames(CTAT_antisens) <- c("null", "faible", "moyen", "élevé")
colnames(CTAT_sens) <- c("null", "faible", "moyen", "élevé")
colnames(CAAC_antisens) <- c("null", "faible", "moyen", "élevé")

#Boucle de calcule du score pour chaque stade et chaque satellite

score_occurence_data <- data.frame(score_CTAC_antisens = rep(0,10),
                                   score_CTAC_sens = rep(0,10),
                                   score_CTAT_antisens = rep(0,10),
                                   score_CTAT_sens = rep(0,10),
                                   score_TTAC_antisens = rep(0,10),
                                   score_TTAC_sens = rep(0,10),
                                   score_CAAC_antisens = rep(0,10),
                                   stade = c(1:10),
                                   row.names = c("stade 1","stade 2","stade 3",
                                                 "stade 4","stade 5","stade 6",
                                                 "stade 7","stade 8","stade 9",
                                                 "stade 10"))

for ( i in c(1:10)){
  
  score_occurence_data[[i,1]] <- (CTAC_antisens[i,2] + 2*CTAC_antisens[i,3] + 3*CTAC_antisens[i,4])
  
  score_occurence_data[[i,2]] <- (CTAC_sens[i,2] + 2*CTAC_sens[i,3] + 3*CTAC_sens[i,4]) 
  
  score_occurence_data[[i,3]] <- (CTAT_antisens[i,2] + 2*CTAT_antisens[i,3] + 3*CTAT_antisens[i,4])
  
  score_occurence_data[[i,4]] <- (CTAT_sens[i,2] + 2*CTAT_sens[i,3] + 3*CTAT_sens[i,4])
  
  score_occurence_data[[i,5]] <- (TTAC_antisens[i,2] + 2*TTAC_antisens[i,3] + 3*TTAC_antisens[i,4])
  
  score_occurence_data[[i,6]] <- (TTAC_sens[i,2] + 2*TTAC_sens[i,3] + 3*TTAC_sens[i,4])
  
  score_occurence_data[[i,7]] <- (CAAC_antisens[i,2] + 2*CAAC_antisens[i,3] + 3*CAAC_antisens[i,4])
  
  print(paste("Boucle",i,"fait"))
}



#Couleur
#CTAC -> "green"
#TTAC -> "jaune"
#CTAT -> "magenta
#Antisens -> signe pleine
#ligne pointiller -> Sens

par(mar = c(4, 5, 4, 4))

par(xpd = TRUE)

plot(score_occurence_data$score_CTAC_antisens,
     type = "l",
     col = "limegreen",
     ylim = c(0,100),
     lwd = 3, 
     ylab = "Score d'expression",
     xlab = "Stade du dévelopement des ovarioles")

matlines(score_occurence_data$score_CTAC_sens,
         type = "l",
         lty = 4,
         col = "limegreen",
         lwd = 3)

matlines(score_occurence_data$score_TTAC_antisens,
         type = "l",
         col = "yellow",
         lwd = 3)

matlines(score_occurence_data$score_TTAC_sens,
         type = "l",
         lty = 4,
         col = "yellow",
         lwd = 3)

matlines(score_occurence_data$score_CTAT_antisens,
         type = "l",
         col = "magenta",
         lwd = 3)

matlines(score_occurence_data$score_CTAT_sens,
         type = "l",
         lty = 4,
         col = "magenta",
         lwd = 3)

matlines(score_occurence_data$score_CAAC_antisens,
         type = "l",
         col = "darkblue",
         lwd = 3)


axis(1, at=seq(1, 10, 1))

axis(2, at=seq(0, 100, 10))

legend("topright", 
       legend = c("Sens AAACTAC","antisens AAACTAC",
                  "Sens AAATTAC","Antisens AAATTAC","Sens AAACTAT",
                  "Antisens AAACTAT","Sens AAACAAC"),
       col = c("limegreen","limegreen",
               "yellow","yellow",
               "magenta","magenta",
               "darkblue"),
       lty = c(1,3,1,3,1,3,1),
       lwd = 3,
       xpd = TRUE,
       horiz = FALSE)

##############################################
#nomalisation par division du nombre total d'observation des stades 
#+1 afin d'évité les divisions 0 si il y a lieux

score_occurence_data_normaliser <- score_occurence_data

for (i in c(1:10)){
  
  score_occurence_data_normaliser[[i,1]] <- score_occurence_data_normaliser[i,1]/(sum(CTAC_antisens[i,1:4]) +1)
  
  score_occurence_data_normaliser[[i,2]] <- score_occurence_data_normaliser[i,2]/(sum(CTAC_sens[i,1:4]) +1)
  
  score_occurence_data_normaliser[[i,3]] <- score_occurence_data_normaliser[i,3]/(sum(CTAT_antisens[i,1:4]) +1)
  
  score_occurence_data_normaliser[[i,4]] <- score_occurence_data_normaliser[i,4]/(sum(CTAT_sens[i,1:4]) +1)
  
  score_occurence_data_normaliser[[i,5]] <- score_occurence_data_normaliser[i,5]/(sum(TTAC_antisens[i,1:4]) +1)
  
  score_occurence_data_normaliser[[i,6]] <- score_occurence_data_normaliser[i,6]/(sum(TTAC_sens[i,1:4]) +1)
  
  score_occurence_data_normaliser[[i,7]] <- score_occurence_data_normaliser[i,7]/(sum(CAAC_antisens[i,1:4]) +1)
}

par(mar = c(4, 5, 4, 4))

par(xpd = TRUE)

plot(score_occurence_data_normaliser$score_CTAC_antisens,
     type = "l",
     col = "limegreen",
     lwd = 3, 
     ylim = c(0,5),
     ylab = "Score d'expression normalisé",
     xlab = "Stade du dévelopement des ovarioles")

matlines(score_occurence_data_normaliser$score_CTAC_sens,
         type = "l",
         lty = 4,
         col = "limegreen",
         lwd = 3)

matlines(score_occurence_data_normaliser$score_TTAC_antisens,
         type = "l",
         col = "yellow",
         lwd = 3)

matlines(score_occurence_data_normaliser$score_TTAC_sens,
         type = "l",
         lty = 4,
         col = "yellow",
         lwd = 3)

matlines(score_occurence_data_normaliser$score_CTAT_antisens,
         type = "l",
         col = "magenta",
         lwd = 3)

matlines(score_occurence_data_normaliser$score_CTAT_sens,
         type = "l",
         lty = 4,
         col = "magenta",
         lwd = 3)

matlines(score_occurence_data_normaliser$score_CAAC_antisens,
         type = "l",
         col = "darkblue",
         lwd = 3)


axis(1, at=seq(1, 10, 1))

axis(2, at=seq(0, 5, 1))

legend("topright", 
       legend = c("Sens AAACTAC","antisens AAACTAC",
                  "Sens AAATTAC","Antisens AAATTAC","Sens AAACTAT",
                  "Antisens AAACTAT","Sens AAACAAC"),
       col = c("limegreen","limegreen",
               "yellow","yellow",
               "magenta","magenta",
               "darkblue"),
       lty = c(1,3,1,3,1,3,1),
       lwd = 3,
       xpd = TRUE,
       horiz = FALSE)


##################
#Graph avec ggplot
##################

windowsFonts( A = windowsFont("baskerville old face"))

ggplot(score_occurence_data_normaliser,
       aes(x = stade)) +
  
  geom_line(aes(y = score_CTAC_antisens,
                color = "Forward AAACTAC"), size = 0.75,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTAC_sens,
                color = "Reverse AAACTAC"), size = 0.75,
            linetype = "dotdash") +
  
  geom_line(aes(y = score_TTAC_antisens,
                color = "Forward AAATTAC"), size = 0.75,
            linetype = "solid") +
  
  geom_line(aes(y = score_TTAC_sens,
                color = "Reverse AAATTAC"), size = 0.75,
            linetype = "dotdash") +
  
  geom_line(aes(y = score_CTAT_antisens,
                color = "Forward AAACTAT"), size = 0.75,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTAT_sens,
                color = "Reverse AAACTAT"), size = 0.75,
            linetype = "dotdash") +
  
  geom_line(aes(y = score_CAAC_antisens,
                color = "Reverse AAACAAC"),
            linetype = "solid", size = 0.75) +
  
  labs(
    x = "Egg chamber stages",
    y = "Normalized expression score",
    title = "The lncRNA expression levels of satellite DNA sequences across egg chamber development",
    subtitle = "The experiments multiplexed the AAACTAC probe with the other related probes",
    color = "LncRNA strand"
  ) +
  
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  ylim(0, 4) +
  
  scale_color_manual(values = c(
    "Forward AAACTAC" = "limegreen",
    "Reverse AAACTAC" = "limegreen",
    "Forward AAATTAC"  = "yellow",
    "Reverse AAATTAC"  = "yellow",
    "Forward AAACTAT"  = "magenta",
    "Reverse AAACTAT"  = "magenta",
    "Forward AAACAAC"  = "blue"
  )) +
  
  theme_bw() +
  theme(
    legend.position = c(0.845, 0.805),
    legend.background = element_rect(fill = FALSE, color = FALSE), 
    legend.title = element_text(size = 10, face = "bold",family = "A"),
    legend.text = element_text(size = 9),
    text = element_text(size = 12, family = "A"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "A"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain", family = "A"),
    plot.margin = margin(20, 20, 20, 20)
  )


############################
#Graph CTAC et CTAT antisens
############################


data_CTAT_CTAC_antisens <- read.csv("Analyse_photo_ovaire_CTAC_CTAT_antisens.csv",
                               sep=";", dec=",", header=FALSE)

CTAC <- matrix( nrow = 10, ncol = 4)
  
CTAT <- matrix( nrow = 10, ncol = 4)

for (i in c(3:12)){
  #null
  CTAC[[i-2,1]] <- length(which(data_CTAT_CTAC_antisens[,i] == "rien"))
  
  #faible
  CTAC[[i-2,2]] <- length(which(data_CTAT_CTAC_antisens[,i] == "faible"))
  
  #moyen
  CTAC[[i-2,3]] <- length(which(data_CTAT_CTAC_antisens[,i] == "moyen"))
  
  #elever
  CTAC[[i-2,4]] <- length(which(data_CTAT_CTAC_antisens[,i] == "elever"))
}

colnames(CTAC) <- c("null", "faible", "moyen", "élevé")

for (i in c(17:26)){
  #null
  CTAT[[i-16,1]] <- length(which(data_CTAT_CTAC_antisens[,i] == "rien"))
  
  #faible
  CTAT[[i-16,2]] <- length(which(data_CTAT_CTAC_antisens[,i] == "faible"))
  
  #moyen
  CTAT[[i-16,3]] <- length(which(data_CTAT_CTAC_antisens[,i] == "moyen"))
  
  #elever
  CTAT[[i-16,4]] <- length(which(data_CTAT_CTAC_antisens[,i] == "elever"))
}

colnames(CTAT) <- c("null", "faible", "moyen", "élevé")

score_occurence_data_CTAT_CTAC <- data.frame(score_CTAT_antisens_both = rep(0,10),
                                   score_CTAC_antisens_both = rep(0,10),
                                   row.names = c("stade 1","stade 2","stade 3",
                                                 "stade 4","stade 5","stade 6",
                                                 "stade 7","stade 8","stade 9",
                                                 "stade 10"))

for ( i in c(1:10)){
  
  score_occurence_data_CTAT_CTAC[[i,1]] <- (CTAT[i,2] + 2*CTAT[i,3] + 3*CTAT[i,4])
  
  score_occurence_data_CTAT_CTAC[[i,2]] <- (CTAC[i,2] + 2*CTAC[i,3] + 3*CTAC[i,4])
}

score_occurence_data_CTAT_CTAC_normaliser <- score_occurence_data_CTAT_CTAC

for (i in c(1:10)){
  
  score_occurence_data_CTAT_CTAC_normaliser[[i,1]] <- score_occurence_data_CTAT_CTAC_normaliser[i,1]/(sum(CTAC[i,1:4]) +1)
  
  score_occurence_data_CTAT_CTAC_normaliser[[i,2]] <- score_occurence_data_CTAT_CTAC_normaliser[i,2]/(sum(CTAC[i,1:4]) +1)
  
}

par(mar = c(4, 5, 4, 4))

par(xpd = TRUE)

plot(score_occurence_data_CTAT_CTAC$score_CTAT_antisens_both,
     type = "l",
     col = "magenta",
     lwd = 3, 
     ylim = c(0,150),
     ylab = "Score d'expression",
     xlab = "Stade du dévelopement des ovarioles")

matlines(score_occurence_data_CTAT_CTAC$score_CTAC_antisens_both,
         type = "l",
         col = "limegreen",
         lwd = 3)

legend("topright", 
       legend = c("Antisens AAACTAC","Antisens AAACTAT"),
       col = c("limegreen","magenta"),
       lty = c(1,1,1,1),
       lwd = 3,
       xpd = TRUE,
       horiz = FALSE)


plot(score_occurence_data_CTAT_CTAC_normaliser$score_CTAT_antisens_both,
     type = "l",
     col = "magenta",
     lwd = 3, 
     ylim = c(0,5),
     ylab = "Score d'expression normalisé",
     xlab = "Stade du dévelopement des ovarioles")

matlines(score_occurence_data_CTAT_CTAC_normaliser$score_CTAC_antisens_both,
         type = "l",
         col = "limegreen",
         lwd = 3)

legend("topright", 
       legend = c("Antisens AAACTAC","Antisens AAACTAT"),
       col = c("limegreen","magenta"),
       lty = c(1,1,1,1),
       lwd = 3,
       xpd = TRUE,
       horiz = FALSE)

plot(as.numeric(data_CTAT_CTAC_antisens[4:143,13]))

matlines(data_CTAT_CTAC_antisens[4:143,27],
        type = "p",
        pch = 19)

#Faire des test de distribution de poissons car expression de gène suit
#la lois de poissons donc, l'a distribution de l'aire et donc les obs de l'aire
#devrais suivre une lois de poissons. 

#################################
################################
#Gaph des combinaison avec CTAC
###############################
###############################

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
    
    CTAC_antisens_compiler[[j,i]] <- sum(CTAC_aveccaac_antisens[j,i], CTAC_avecctat_antisens[j,i], 
                                         CTAC_avecttac_antisens[j,i],CTAC_avecctat_antisens_rep3[j,i], 
                                         CTAC_avecttac_antisens_rep2[j,i])
    
    CTAT_antisens_compiler[[j,i]] <- sum(CTATavecCTAC_antisens[j,i], CTATavecCTAC_antisens_rep3[j,i])
    
    CTAC_sens_compiler[[j,i]] <- sum(CTAC_avecctat_sens[j,i], CTAC_avecctat_sens_rep1[j,i])
    
    CTAT_sens_compiler[[j,i]] <- sum(CTAT_avecCTAC_sens[j,i], CTAT_avecCTAC_sens_rep1[j,i])
    
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




########################
#Graphique avec GGplot
######################


transparence = 0.3


windowsFonts( A = windowsFont("baskerville old face"))

graph <- ggplot(score_occurence_data_normaliser_combinaison,
       aes(x = stade)) +
  
  geom_line(aes(y = score_CTAC_aveccaac_antisens,
                color = "Forward AAACTAC"), size = 0.75, alpha = transparence,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTAC_avecttac_antisens,
                color = "Forward AAACTAC"), size = 0.75,alpha = transparence,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTAC_avecttat_antisens_rep2,
                color = "Forward AAACTAC"), size = 0.75,alpha = transparence,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTAC_avecctat_antisens,
                color = "Forward AAACTAC"), size = 0.75,alpha = transparence,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTAC_avecctat_antisens_rep3,
                color = "Forward AAACTAC"), size = 0.75,alpha = transparence,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTAC_avecctat_sens,
                color = "Reverse AAACTAC"), size = 0.75, alpha = transparence,
            linetype = "dotdash") +
  
  geom_line(aes(y = score_CTAC_avecctat_sens_rep1,
                color = "Reverse AAACTAC"), size = 0.75, alpha = transparence,
            linetype = "dotdash") +
  
  #geom_line(aes(y = score_TTACavecCTAC_antisens,
                #color = "Forward AAATTAC"), size = 0.75,
            #linetype = "solid") +
  
  geom_line(aes(y = score_CTATavecCTAC_antisens,
                color = "Forward AAACTAT"), size = 0.75, alpha = transparence,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTATavecCTAC_antisens_rep3,
                color = "Forward AAACTAT"), size = 0.75, alpha = transparence,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTAT_avecCTAC_sens,
                color = "Reverse AAACTAT"), size = 0.75, alpha = transparence,
            linetype = "dotdash") +
  
  geom_line(aes(y = score_CTAT_avecCTAC_sens_rep1,
                color = "Reverse AAACTAT"), size = 0.75, alpha = transparence,
            linetype = "dotdash") +
  
  #geom_line(aes(y = score_CAACavecCTAC_antisens,
                #color = "Forward AAACAAC"), size = 0.75,
            #linetype = "solid") +
  
  geom_line(data = score_occurence_data_CTAC_antisens_compiler_normaliser,
            aes(y = Score_CTAC_antisens_compiler,
                color = "Compilation foward AAACTAC"),
            linetype = "solid", size = 1.3) +
  
  geom_line(data = score_occurence_data_CTAT_antisens_compiler_normaliser,
            aes(y = Score_CTAT_antisens_compiler,
                color = "Compilation foward AAACTAT"),
            linetype = "solid", size = 1.3) +
  
  geom_line(data = score_occurence_data_CTAC_sens_compiler_normaliser,
            aes(y = Score_CTAC_sens_compiler,
                color = "Compilation reverse AAACTAC"),
            linetype = "dotdash", size = 1.3) +
  
  geom_line(data = score_occurence_data_CTAT_sens_compiler_normaliser,
            aes(y = Score_CTAT_sens_compiler,
                color = "Compilation reverse AAACTAT"),
            linetype = "dotdash", size = 1.3) +
  
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
    legend.title = element_text(size = 10, face = "bold",family = "A"),
    legend.text = element_text(size = 9),
    text = element_text(size = 12, family = "A"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "A"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain", family = "A"),
    plot.margin = margin(20, 20, 20, 20)
  )

print(graph)

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












