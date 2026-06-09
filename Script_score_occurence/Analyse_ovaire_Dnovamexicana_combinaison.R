##################################
#Analyse donnée ovaire Dnovamexicana
##################################

library(ggplot2)

data_CAAC_antisens_avec_ctac <- read.csv("Analyse_photo_ovaire_Dnovamexicana_CAAC_avec_ctac.csv",
                                         sep=";", dec=",", header=FALSE)

data_TTAC_antisens_avec_ctac <- read.csv("Analyse_photo_ovaire_Dnovamexicana_TTAC_avec_ctac.csv",
                                         sep=";", dec=",", header=FALSE)

data_CTAC_antisens_avec_caac <- read.csv("Analyse_photo_ovaire_Dnovamexicana_CTAC_avec_caac.csv",
                                         sep=";", dec=",", header=FALSE)

data_CTAC_antisens_avec_ttac <- read.csv("Analyse_photo_ovaire_Dnovamexicana_CTAC_avec_ttac.csv",
                                         sep=";", dec=",", header=FALSE)


####################

#CAAC_antisens <-matrix( nrow = 10, ncol = 4)
#TTAC_antisens <-matrix( nrow = 10, ncol = 4)
#CTAC_antisens <-matrix( nrow = 10, ncol = 4)

CAAC_antisens_ctac <-matrix( nrow = 10, ncol = 4)
TTAC_antisens_ctac <-matrix( nrow = 10, ncol = 4)
CTAC_antisens_caac <-matrix( nrow = 10, ncol = 4)
CTAC_antisens_ttac <-matrix( nrow = 10, ncol = 4)

###################

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

#CAAC_antisens <- decompte(data_source = data_CAAC_antisens, data_puit = CAAC_antisens)
#TTAC_antisens <- decompte(data_source = data_TTAC_antisens, data_puit = TTAC_antisens)
#CTAC_antisens <- decompte(data_source = data_CTAC_antisens, data_puit = CTAC_antisens)

CAAC_antisens_ctac <- decompte(data_source = data_CAAC_antisens_avec_ctac, data_puit = CAAC_antisens_ctac)
TTAC_antisens_ctac <- decompte(data_source = data_TTAC_antisens_avec_ctac, data_puit = TTAC_antisens_ctac)
CTAC_antisens_caac <- decompte(data_source = data_CTAC_antisens_avec_caac, data_puit = CTAC_antisens_caac)
CTAC_antisens_ttac <- decompte(data_source = data_CTAC_antisens_avec_ttac, data_puit = CTAC_antisens_ttac)

#colnames(CAAC_antisens) <- c("null", "faible", "moyen", "élevé")
#colnames(TTAC_antisens) <- c("null", "faible", "moyen", "élevé")
#colnames(CTAC_antisens) <- c("null", "faible", "moyen", "élevé")

colnames(CAAC_antisens_ctac) <- c("null", "faible", "moyen", "élevé")
colnames(TTAC_antisens_ctac) <- c("null", "faible", "moyen", "élevé")
colnames(CTAC_antisens_caac) <- c("null", "faible", "moyen", "élevé")
colnames(CTAC_antisens_ttac) <- c("null", "faible", "moyen", "élevé")

noms_matrices <- c("CTAC_antisens_ttac", "CTAC_antisens_caac"
)

total_observations <- sum(sapply(mget(noms_matrices), sum, na.rm = TRUE))

print(total_observations)

matrice_fusionnee <- Reduce("+", mget(noms_matrices))

comptage_par_ligne <- rowSums(matrice_fusionnee, na.rm = TRUE)

noms_stades <- paste("Stade", 1:(1 + length(comptage_par_ligne) - 1))
names(comptage_par_ligne) <- noms_stades


print(comptage_par_ligne)

#Boucle de calcule du score pour chaque stade et chaque satellite

score_occurence_data <- data.frame(score_CAAC_antisens = rep(0,10),
                                   score_TTAC_antisens = rep(0,10),
                                   score_CTAC_antisens = rep(0,10),
                                   stade = rep(1:10),
                                   row.names = c("stade 1","stade 2","stade 3",
                                                 "stade 4","stade 5","stade 6",
                                                 "stade 7","stade 8","stade 9",
                                                 "stade 10"))

score_occurence_data_combinaison <- data.frame(score_CTAC_antisens_ttac = rep(0,10),
                                               score_CTAC_antisens_caac = rep(0,10),
                                               score_TTAC_antisens_ctac = rep(0,10),
                                               score_CAAC_antisens_ctac = rep(0,10),
                                               stade = rep(1:10),
                                               row.names = c("stade 1","stade 2","stade 3",
                                                             "stade 4","stade 5","stade 6",
                                                             "stade 7","stade 8","stade 9",
                                                             "stade 10"))


for ( i in c(1:10)){
  
  #score_occurence_data[[i,3]] <- (CTAC_antisens[i,2] + 2*CTAC_antisens[i,3] + 3*CTAC_antisens[i,4])
  
  #score_occurence_data[[i,2]] <- (TTAC_antisens[i,2] + 2*TTAC_antisens[i,3] + 3*TTAC_antisens[i,4])
  
  #score_occurence_data[[i,1]] <- (CAAC_antisens[i,2] + 2*CAAC_antisens[i,3] + 3*CAAC_antisens[i,4])
  
  score_occurence_data_combinaison[[i,1]] <- (CTAC_antisens_ttac[i,2] + 2*CTAC_antisens_ttac[i,3] + 3*CTAC_antisens_ttac[i,4])
  
  score_occurence_data_combinaison[[i,2]] <- (CTAC_antisens_caac[i,2] + 2*CTAC_antisens_caac[i,3] + 3*CTAC_antisens_caac[i,4])                                            
  
  score_occurence_data_combinaison[[i,3]] <- (TTAC_antisens_ctac[i,2] + 2*TTAC_antisens_ctac[i,3] + 3*TTAC_antisens_ctac[i,4])
  
  score_occurence_data_combinaison[[i,4]] <- (CAAC_antisens_ctac[i,2] + 2*CAAC_antisens_ctac[i,3] + 3*CAAC_antisens_ctac[i,4])
  
  
  
}

#####################
#Plot
####################

par(mar = c(4, 5, 4, 4))

par(xpd = TRUE)

plot(score_occurence_data$score_CAAC_antisens,
     type = "l",
     col = "darkblue",
     ylim = c(0,150),
     lwd = 3, 
     ylab = "Score d'expression",
     xlab = "Stade du dévelopement des ovarioles")

matlines(score_occurence_data$score_TTAC_antisens,
         type = "l",
         col = "yellow",
         lwd = 3)

matlines(score_occurence_data$score_CTAC_antisens,
         type = "l",
         col = "green",
         lwd = 3)

axis(1, at=seq(1, 10, 1))

legend("topright", 
       legend = c("Sens AAACTAC",
                  #"antisens AAACTAC",
                  "Sens AAATTAC",#"Antisens AAATTAC","Sens AAACTAT",
                  #"Antisens AAACTAT",
                  "Sens AAACAAC"),
       col = c("limegreen",#"limegreen",
               "yellow",#"yellow",
               #"magenta","magenta",
               "darkblue"),
       lty = c(1,#3,
               1,
               #3,
               1,
               #3,
               1),
       lwd = 3,
       xpd = TRUE,
       horiz = FALSE)

#################################
#normalisation
################################

#score_occurence_data_normaliser <- score_occurence_data

score_occurence_data_combiniason_normaliser <- score_occurence_data_combinaison

for (i in c(1:10)){
  
  #score_occurence_data_normaliser[[i,1]] <- score_occurence_data_normaliser[i,1]/(sum(CAAC_antisens[i,1:4]) +1)
  
  #score_occurence_data_normaliser[[i,2]] <- score_occurence_data_normaliser[i,2]/(sum(TTAC_antisens[i,1:4]) +1)
  
  #score_occurence_data_normaliser[[i,3]] <- score_occurence_data_normaliser[i,3]/(sum(CTAC_antisens[i,1:4]) +1)
  
  score_occurence_data_combiniason_normaliser[[i,1]] <- score_occurence_data_combiniason_normaliser[i,1]/(sum(CTAC_antisens_ttac[i,1:4]) +1)
  
  score_occurence_data_combiniason_normaliser[[i,2]] <- score_occurence_data_combiniason_normaliser[i,2]/(sum(CTAC_antisens_caac[i,1:4]) +1)
  
  score_occurence_data_combiniason_normaliser[[i,3]] <- score_occurence_data_combiniason_normaliser[i,3]/(sum(CAAC_antisens_ctac[i,1:4]) +1)
  
  score_occurence_data_combiniason_normaliser[[i,4]] <- score_occurence_data_combiniason_normaliser[i,4]/(sum(TTAC_antisens_ctac[i,1:4]) +1)
  
  
}

############################
#compilationde tout les CTAC
############################

score_occurence_data_CTAC_moyenne <- data.frame(moyenne = rep(0,10), 
                                                
                                                SD = rep(0,10),
                                                
                                                stade = c(1:10),
                                                
                                                row.names = c("stade 1","stade 2","stade 3",
                                                              "stade 4","stade 5","stade 6",
                                                              "stade 7","stade 8","stade 9",
                                                              "stade 10"))



for (j in c(1:10)) {
  
  score_occurence_data_CTAC_moyenne[[j,1]] <- mean(c(score_occurence_data_combiniason_normaliser$score_CTAC_antisens_ttac[j],
                                                     score_occurence_data_combiniason_normaliser$score_CTAC_antisens_caac[j]))
  
  score_occurence_data_CTAC_moyenne[[j,2]] <- sd(c(score_occurence_data_combiniason_normaliser$score_CTAC_antisens_ttac[j],
                                                   score_occurence_data_combiniason_normaliser$score_CTAC_antisens_caac[j]))
  
}


########################
#Graphique avec GGplot
######################


transparence = 0.2


windowsFonts( A = windowsFont("baskerville old face"))

graph <- ggplot(score_occurence_data_combiniason_normaliser,
       aes(x = stade)) +
  
  geom_line(aes(y = score_CTAC_antisens_ttac,
                color = "Forward AAACTAC"), size = 0.8, alpha = transparence,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTAC_antisens_caac,
                color = "Forward AAACTAC"), size = 0.8,alpha = transparence,
            linetype = "solid") +
  
  geom_line(aes(y = score_CAAC_antisens_ctac,
                color = "Forward AAACAAC"), size = 0.8,
            linetype = "solid") +
  
  geom_line(data = score_occurence_data_CTAC_moyenne,
            aes(y = moyenne,
                color = "Compilation AAACTAC"),
            linetype = "solid", size = 1.3) +
  
  labs(
    x = "Egg chamber stages",
    y = "Normalized expression score",
    title = "The lncRNA expression levels of satellite DNA sequences \n across egg chamber development in Drosophila novamexicana",
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
    "Compilation AAACTAC" = "limegreen"
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

#################
#exportation des donner

data_combier_Dnova_CTAC_antisens <- score_occurence_data_combiniason_normaliser[,c(1,2,5)] %>%
  pivot_longer(
    cols = !stade,
    names_to = "brin",
    values_to = "score"
  )

data_combier_Dnova_CTAC_antisens$sp <- "Dnova"

write.csv2(data_combier_Dnova_CTAC_antisens, "Donner_Dnova_CTAC_antisens.csv", row.names = FALSE)

















