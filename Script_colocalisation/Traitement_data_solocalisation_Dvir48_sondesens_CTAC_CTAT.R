###########################################
#Script analyse résultats colocalisation
#CTAC+CTAT sonde sens (brin antisens)
##########################################

#test satistique divers


library(car)
library(multcompView)
library(dplyr)
library(agricolae)
library(coin)
library(FSA)
library(multcomp)
library(rcompanion)
library(psych)
library(nlstools)
library(lmtest)
library(tidyr)
library(SuppDists)

#regression
library(DHARMa)
library(splines)
library(gamm4)
library(mgcv)
library(gvlma)
library(stargazer)
library(performance)
library(see)
library(glmmTMB)
library(MASS)
library(robustbase)

#representation graphique
library(ggeffects)
library(ggplot2)
library(ggpubr)

#citation 
library(grateful)

#########################################################################
#ouverture des fichier contenant les indices, est sous forme txt dans csv

#Vecteur des nom des futur fichier utiliser

nom_des_fichier <- vector()
nom_des_fichier_rep1 <- vector()

for ( i in (3:10)){
  
  nom_des_fichier[[i-2]] <- paste0("score_coloc_stade",i)
  nom_des_fichier_rep1[[i-2]] <- paste0("score_coloc_stade",i,"_rep1")
  
}

#ouverture de la répétition 0

for ( i in c(3:10)) {
  
  ID <- read.csv(paste0("Coef_coloc_Dvir48_CTAC_CTAT_sensboth_stade",i,".csv"),
                              sep = ",", dec = ".", header =TRUE)
  assign(nom_des_fichier[i-2],ID)
  
}

#ouverture de la répétition 1

for ( i in c(3:10)) {
  
  ID <- read.csv(paste0("Coef_coloc_Dvir48_CTAC_CTAT_sensboth_stade",i,"_rep1.csv"),
                                            sep = ",", dec = ".", header =TRUE)
  assign(nom_des_fichier_rep1[i-2],ID)
  
}

#en suivant la même logique, on peut faire ensuite une boulce j à l'intérieure de la boucle i pour i stade
#et J répétiton

###############################################################################
#Nettotage des donnée afin de ne récupérer que les données importante au analyse

#Fonction de nettoyage

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

###################################################################
#Nettoyage des données
indice_colocalisation_stade3 <- Traitement_data_coloc(score_coloc_stade3)
indice_colocalisation_stade4 <- Traitement_data_coloc(score_coloc_stade4)
indice_colocalisation_stade5 <- Traitement_data_coloc(score_coloc_stade5)
indice_colocalisation_stade6 <- Traitement_data_coloc(score_coloc_stade6)
indice_colocalisation_stade7 <- Traitement_data_coloc(score_coloc_stade7)
indice_colocalisation_stade8 <- Traitement_data_coloc(score_coloc_stade8)
indice_colocalisation_stade9 <- Traitement_data_coloc(score_coloc_stade9)
indice_colocalisation_stade10 <- Traitement_data_coloc(score_coloc_stade10)

indice_colocalisation_stade3_rep1 <- Traitement_data_coloc(score_coloc_stade3_rep1)
indice_colocalisation_stade4_rep1 <- Traitement_data_coloc(score_coloc_stade4_rep1)
indice_colocalisation_stade5_rep1 <- Traitement_data_coloc(score_coloc_stade5_rep1)
indice_colocalisation_stade6_rep1 <- Traitement_data_coloc(score_coloc_stade6_rep1)
indice_colocalisation_stade7_rep1 <- Traitement_data_coloc(score_coloc_stade7_rep1)
indice_colocalisation_stade8_rep1 <- Traitement_data_coloc(score_coloc_stade8_rep1)
indice_colocalisation_stade9_rep1 <- Traitement_data_coloc(score_coloc_stade9_rep1)
indice_colocalisation_stade10_rep1 <- Traitement_data_coloc(score_coloc_stade10_rep1)

####################################################################################

#Fusion des row de tout les essaies

indice_colocalisation_stade3_totale <- rbind(indice_colocalisation_stade3_rep1,
                                      indice_colocalisation_stade3)

indice_colocalisation_stade4_totale <- rbind(indice_colocalisation_stade4_rep1,
                                      indice_colocalisation_stade4)

indice_colocalisation_stade5_totale <- rbind(indice_colocalisation_stade5_rep1,
                                      indice_colocalisation_stade5)

indice_colocalisation_stade6_totale <- rbind(indice_colocalisation_stade6_rep1,
                                      indice_colocalisation_stade6)

indice_colocalisation_stade7_totale <- rbind(indice_colocalisation_stade7_rep1,
                                      indice_colocalisation_stade7)

indice_colocalisation_stade8_totale <- rbind(indice_colocalisation_stade8_rep1,
                                      indice_colocalisation_stade8)

indice_colocalisation_stade9_totale <- rbind(indice_colocalisation_stade9_rep1,
                                      indice_colocalisation_stade9)

indice_colocalisation_stade10_totale <- rbind(indice_colocalisation_stade10_rep1,
                                       indice_colocalisation_stade10)

##################################################################################

#Calcule des moyennes de tout les scores

#Fonction de calcule des moyenne pour chaque score de chaque stade

Calcule_moyene <- function(tableau_data_souce){
  
  #Tableau de donné ou tout seras entreposer
  
  indice_coloc_moyenne <- data.frame(Coef_pearson = rep(0,1),
                                     
                                     Coeaf_overlap = rep(0,1),
                                     
                                     k1  = rep(0,1),
                                     
                                     k2 = rep(0,1), 
                                     
                                     M1 = rep(0,1),
                                     
                                     M2 = rep(0,1),
                                     
                                     a_cytofluogram  = rep(0,1),
                                     
                                     b_cytofluogram = rep(0,1), 
                                     
                                     ICQ = rep(0,1))
  
  for (i in c(1:ncol(tableau_data_souce))) {
    
    indice_coloc_moyenne[,i] <- mean(tableau_data_souce[,i])
    
    
  }
  
  return(indice_coloc_moyenne)
}  

moyenne_coloc_coef_stade3 <- Calcule_moyene(indice_colocalisation_stade3_totale)
moyenne_coloc_coef_stade4 <- Calcule_moyene(indice_colocalisation_stade4_totale)
moyenne_coloc_coef_stade5 <- Calcule_moyene(indice_colocalisation_stade5_totale)
moyenne_coloc_coef_stade6 <- Calcule_moyene(indice_colocalisation_stade6_totale)
moyenne_coloc_coef_stade7 <- Calcule_moyene(indice_colocalisation_stade7_totale)
moyenne_coloc_coef_stade8 <- Calcule_moyene(indice_colocalisation_stade8_totale)
moyenne_coloc_coef_stade9 <- Calcule_moyene(indice_colocalisation_stade9_totale)
moyenne_coloc_coef_stade10 <- Calcule_moyene(indice_colocalisation_stade10_totale)


moyenne_coloc_coef <- rbind(moyenne_coloc_coef_stade3,moyenne_coloc_coef_stade4,moyenne_coloc_coef_stade5,
                            moyenne_coloc_coef_stade6, moyenne_coloc_coef_stade7, moyenne_coloc_coef_stade8,
                            moyenne_coloc_coef_stade9, moyenne_coloc_coef_stade10)

rownames(moyenne_coloc_coef) <- c(
  "Stade3",
  "Stade4",
  "Stade5",
  "Stade6",
  "Stade7",
  "Stade8",
  "Stade9",
  "Stade10"
)

Stade <- c(
  "Stade3",
  "Stade4",
  "Stade5",
  "Stade6",
  "Stade7",
  "Stade8",
  "Stade9",
  "Stade10")

moyenne_coloc_coef <- cbind(moyenne_coloc_coef,Stade)


######################################
#data long des score combinés ensemble

indice_colocalisation_stade3_totale$stade <- rep(3,length(indice_colocalisation_stade3_totale[,1]))

indice_colocalisation_stade4_totale$stade <- rep(4,length(indice_colocalisation_stade4_totale[,1]))

indice_colocalisation_stade5_totale$stade <- rep(5,length(indice_colocalisation_stade5_totale[,1]))

indice_colocalisation_stade6_totale$stade <- rep(6,length(indice_colocalisation_stade6_totale[,1]))

indice_colocalisation_stade7_totale$stade <- rep(7,length(indice_colocalisation_stade7_totale[,1]))

indice_colocalisation_stade8_totale$stade <- rep(8,length(indice_colocalisation_stade8_totale[,1]))

indice_colocalisation_stade9_totale$stade <- rep(9,length(indice_colocalisation_stade9_totale[,1]))

indice_colocalisation_stade10_totale$stade <- rep(10,length(indice_colocalisation_stade10_totale[,1]))

Indice_colocalisation_tout_stade <- rbind(indice_colocalisation_stade3_totale,
                                          indice_colocalisation_stade4_totale,
                                          indice_colocalisation_stade5_totale,
                                          indice_colocalisation_stade6_totale,
                                          indice_colocalisation_stade7_totale,
                                          indice_colocalisation_stade8_totale,
                                          indice_colocalisation_stade9_totale,
                                          indice_colocalisation_stade10_totale)

Indice_colocalisation_tout_stade_lon <- Indice_colocalisation_tout_stade[,c(1,2,5,6,9,10)] %>%
  pivot_longer(
    cols = colnames(Indice_colocalisation_tout_stade[,c(1,2,5,6,9)]), 
    names_to = "ID_score",             
    values_to = "valeurs"
  )

####################################################################################
#Représentation graphique de l'évolution des scores en fonction des stades
#Ainsi que leurs variances interstade
####################################################################################


#####################################
#Test modele LM
###################################

#confection de model préléminaire pour tester les conditions d'application

model_pearson <- lm(data = Indice_colocalisation_tout_stade, Coef_pearson ~ bs(stade,7))

#on test les conditions
check_model(model_pearson)
check_autocorrelation(model_pearson)
check_heteroscedasticity(model_pearson)

shapiro.test(model_pearson$residuals)
mean(model_pearson$residuals)
cor.test(Indice_colocalisation_tout_stade$stade, model_pearson$residuals)
model_performance(model_pearson)

#vérification linéarité 

newdat <- Indice_colocalisation_tout_stade_lon %>% 
  filter(ID_score == "Coef_pearson")
newdat$fit <- fitted(model_pearson)
newdat$res <- resid(model_pearson)
newdat$weights <- model_pearson$w

plot_1 <- ggplot(newdat, aes(fit, res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('linéarité pearson sens')

#sur les valeurs x exacte

plot_2 <- ggplot(newdat, aes(x = Indice_colocalisation_tout_stade$stade, y = res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('linéarité pearson sens sur stade')

ggarrange(plot_1,plot_2)

######################

for (i in c(1:length(Indice_colocalisation_tout_stade$Coeaf_overlap))) {
  
  if (Indice_colocalisation_tout_stade$Coeaf_overlap[i] == 0){
    
    Indice_colocalisation_tout_stade$Coeaf_overlap[i] <- NA
    
  }
  
}

model_Overlap <- rlm(data = Indice_colocalisation_tout_stade, Coeaf_overlap ~ bs(stade,5))

shapiro.test(model_Overlap$residuals)

check_heteroscedasticity(model_Overlap)

qqPlot(resid(model_Overlap),
       main = "Overlap")
#on test les conditions

newdat <- Indice_colocalisation_tout_stade_lon %>% 
  filter(ID_score == "Coeaf_overlap",
         valeurs != 0)
newdat$fit <- fitted(model_Overlap)
newdat$res <- resid(model_Overlap)
newdat$weights <- model_Overlap$w


ggplot(newdat, aes(x = fit, y = res)) +
  geom_point(aes(color = weights, alpha = weights),
             size = 4) +
  scale_color_gradient(low = "cyan", high = "blue") +
  geom_smooth(color = "black", se = TRUE, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'red4', size = 1) +
  labs(
    title = 'Overlap',
    x = 'Valeurs ajustées (Fitted)',
    y = 'Résidus',
    color = 'Poids'
  ) +
  theme_minimal()


plot_2 <- ggplot(newdat, aes(x = stade, y = res)) +
  geom_point(aes(color = weights, alpha = weights),
             size = 4) +
  scale_color_gradient(low = "cyan", high = "blue") +
  geom_smooth(color = "black", se = TRUE, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'red4', size = 1) +
  labs(
    title = 'Overlap sur stade',
    x = 'Valeurs ajustées (Fitted)',
    y = 'Résidus',
    color = 'Poids'
  ) +
  theme_minimal()

ggarrange(plot_1,plot_2)

mean(model_Overlap$residuals)
dwtest(model_Overlap)
lmtest::bptest(model_Overlap)
cor.test(newdat$stade, model_Overlap$residuals)


#on replace les 0 

for (i in c(1:length(Indice_colocalisation_tout_stade$Coeaf_overlap))) {
  
  if (is.na(Indice_colocalisation_tout_stade$Coeaf_overlap[i])){
    
    Indice_colocalisation_tout_stade$Coeaf_overlap[i] <- 0
    
  }
  
}


######################
model_ICQ <- rlm(data = Indice_colocalisation_tout_stade, ICQ ~ bs(stade,4))

shapiro.test(model_ICQ$residuals)

check_heteroscedasticity(model_ICQ)

mean(model_ICQ$residuals)

dwtest(model_ICQ)

cor.test(Indice_colocalisation_tout_stade$stade, model_ICQ$residuals)

qqPlot(resid(model_ICQ),
       main = "ICQ")

newdat <- Indice_colocalisation_tout_stade_lon %>% 
  filter(ID_score == "ICQ")
newdat$fit <- fitted(model_ICQ)
newdat$res <- resid(model_ICQ)
newdat$weights <- model_ICQ$w


plot_1 <- ggplot(newdat, aes(x = fit, y = res)) +
  geom_point(aes(color = weights, alpha = weights),
             size = 4) +
  scale_color_gradient(low = "cyan", high = "blue") +
  geom_smooth(color = "black", se = TRUE, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'red4', size = 1) +
  labs(
    title = 'ICQ',
    x = 'Valeurs ajustées (Fitted)',
    y = 'Résidus',
    color = 'Poids'
  ) +
  theme_minimal()

plot_2 <- ggplot(newdat, aes(x = Indice_colocalisation_tout_stade$stade, y = res)) +
  geom_point(aes(color = weights, alpha = weights),
             size = 4) +
  scale_color_gradient(low = "cyan", high = "blue") +
  geom_smooth(color = "black", se = TRUE, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'red4', size = 1) +
  labs(
    title = 'ICQ sur stade',
    x = 'Valeurs ajustées (Fitted)',
    y = 'Résidus',
    color = 'Poids'
  ) +
  theme_minimal()

ggarrange(plot_1,plot_2)

######################
model_M1 <- rlm(data = Indice_colocalisation_tout_stade, M1 ~  bs(stade,5))

#on test les conditions
shapiro.test(model_M1$residuals)
mean(model_M1$residuals)
check_heteroscedasticity(model_M1)
cor.test(Indice_colocalisation_tout_stade$stade, model_M1$residuals)

qqPlot(resid(model_M1),
       main = "M1")

newdat <- Indice_colocalisation_tout_stade_lon %>% 
  filter(ID_score == "M1")
newdat$fit <- fitted(model_M1)
newdat$res <- resid(model_M1)
newdat$weights <- model_M1$w


plot_1 <- ggplot(newdat, aes(x = fit, y = res)) +
  geom_point(aes(color = weights, alpha = weights),
             size = 4) +
  scale_color_gradient(low = "cyan", high = "blue") +
  geom_smooth(color = "black", se = TRUE, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'red4', size = 1) +
  labs(
    title = 'M1',
    x = 'Valeurs ajustées (Fitted)',
    y = 'Résidus',
    color = 'Poids'
  ) +
  theme_minimal()

plot_2 <- ggplot(newdat, aes(x = Indice_colocalisation_tout_stade$stade, y = res)) +
  geom_point(aes(color = weights, alpha = weights),
             size = 4) +
  scale_color_gradient(low = "cyan", high = "blue") +
  geom_smooth(color = "black", se = TRUE, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'red4', size = 1) +
  labs(
    title = 'M1 sur stade',
    x = 'Valeurs ajustées (Fitted)',
    y = 'Résidus',
    color = 'Poids'
  ) +
  theme_minimal()

ggarrange(plot_1,plot_2)

######################
model_M2 <- rlm(data = Indice_colocalisation_tout_stade,M2 ~ bs(stade,6))

#on test les conditions
mean(model_M2$residuals)
shapiro.test(model_M2$residuals)
check_heteroscedasticity(model_M2)
dwtest(model_M2)
cor.test(Indice_colocalisation_tout_stade$stade, model_M2$residuals)

qqPlot(resid(model_M2),
       main = "M2")

newdat <- Indice_colocalisation_tout_stade_lon %>% 
  filter(ID_score == "M2")
newdat$fit <- fitted(model_M2)
newdat$res <- resid(model_M2)
newdat$weights <- model_M2$w


plot_1 <- ggplot(newdat, aes(x = fit, y = res)) +
  geom_point(aes(color = weights, alpha = weights),
             size = 4) +
  scale_color_gradient(low = "cyan", high = "blue") +
  geom_smooth(color = "black", se = TRUE, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'red4', size = 1) +
  labs(
    title = 'M2',
    x = 'Valeurs ajustées (Fitted)',
    y = 'Résidus',
    color = 'Poids'
  ) +
  theme_minimal()

plot_2 <- ggplot(newdat, aes(x = Indice_colocalisation_tout_stade$stade, y = res)) +
  geom_point(aes(color = weights, alpha = weights),
             size = 4) +
  scale_color_gradient(low = "cyan", high = "blue") +
  geom_smooth(color = "black", se = TRUE, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'red4', size = 1) +
  labs(
    title = 'M2 sur stade',
    x = 'Valeurs ajustées (Fitted)',
    y = 'Résidus',
    color = 'Poids'
  ) +
  theme_minimal()

ggarrange(plot_1,plot_2)

#on test les conditions
mean(model_M2$residuals)
shapiro.test(model_M2$residuals)
dwtest(model_M2)
lmtest::bptest(model_M2)
cor.test(Indice_colocalisation_tout_stade$stade, model_M2$residuals)


#mean de tout le monde
mean(model_pearson$residuals)
mean(model_Overlap$residuals)
mean(model_ICQ$residuals)
mean(model_M1$residuals)
mean(model_M2$residuals)


####################################
#Prédiction de chaque model
####################################

summary(model_pearson)
summary(model_Overlap)
summary(model_ICQ)
summary(model_M1)
summary(model_M2)

#calcule des R2 et speudo R2

data_r2 <- data.frame("Pearson" = c(0,0),
                      "Overlap" = c(0,0),
                      "ICQ" = c(0,0),
                      "M1" = c(0,0),
                      "M2" = c(0,0),
                      row.names = c("R^2","P.value"))

data_r2[1,1] <- summary(model_pearson)$r.squared
data_r2[1,2] <- (cor(model_Overlap$model$Coeaf_overlap, predict(model_Overlap)))^2
data_r2[1,3] <- (cor(model_ICQ$model$ICQ, predict(model_ICQ)))^2
data_r2[1,4] <- (cor(model_M1$model$M1, predict(model_M1)))^2
data_r2[1,5] <- (cor(model_M2$model$M2, predict(model_M2)))^2

pred_pearson <- predict_response(model_pearson, terms = "stade [3:10]")
pred_overlap <- predict_response(model_Overlap, terms = "stade [3:10]")
pred_ICQ <- predict_response(model_ICQ, terms = "stade [3:10]")
pred_M1 <- predict_response(model_M1, terms = "stade [3:10]")
pred_M2 <- predict_response(model_M2, terms = "stade [3:10]")

plot_pred_pearson <- plot(pred_pearson,
                          show_data = TRUE,
                          jitter = 0.2)

plot_pred_overlap <- plot(pred_overlap,
                          show_residuals = TRUE,
                          jitter = 0.2)

plot_pred_ICQ <- plot(pred_ICQ,
                      show_data = TRUE,
                      jitter = 0.2)

plot_pred_M1 <- plot(pred_M1,
                     show_data = TRUE,
                     jitter = 0.2)

plot_pred_M2 <- plot(pred_M2,
                     show_data = TRUE,
                     jitter = 0.2)

plot_fusion <- ggarrange(plot_pred_pearson,plot_pred_overlap,
                         plot_pred_ICQ,plot_pred_M1,plot_pred_M2,
                         nrow = 3,
                         ncol = 2
)

annotate_figure(plot_fusion, top = text_grob("Prediction score Sens", 
                                             color = "red", face = "bold", size = 14))


#on augmente le détail de la prédiction
#afin de comblé les cathégorie discrète

stade_dense <- seq(min(Indice_colocalisation_tout_stade$stade), 
                   max(Indice_colocalisation_tout_stade$stade), 
                   length.out = 200)

###########################
#Pearson

df_predict_pearson <- expand.grid(
  stade = stade_dense,
  Score = "Coef_pearson" 
)

pred_pearson <- predict(model_pearson, newdata = df_predict_pearson , 
                        se.fit = TRUE, interval = "confidence", level = 0.95)

df_predict_pearson <- df_predict_pearson %>%
  bind_cols(as_tibble(pred_pearson))


#Overlap

df_predict_Overlap <- expand.grid(
  stade = stade_dense,
  Score = "Coef_Overlap" 
)

pred_Overlap <- predict(model_Overlap, newdata = df_predict_pearson , 
                        se.fit = TRUE, interval = "confidence", level = 0.95)

df_predict_Overlap <- df_predict_Overlap %>%
  bind_cols(as_tibble(pred_Overlap))


#ICQ

df_predict_ICQ <- expand.grid(
  stade = stade_dense,
  Score = "Coef_ICQ" 
)

pred_ICQ <- predict(model_ICQ, newdata = df_predict_pearson , 
                    se.fit = TRUE, interval = "confidence", level = 0.95)

df_predict_ICQ <- df_predict_ICQ %>%
  bind_cols(as_tibble(pred_ICQ))

#M1

df_predict_M1 <- expand.grid(
  stade = stade_dense,
  Score = "Coef_M1" 
)

pred_M1 <- predict(model_M1, newdata = df_predict_pearson , 
                   se.fit = TRUE, interval = "confidence", level = 0.95)

df_predict_M1 <- df_predict_M1 %>%
  bind_cols(as_tibble(pred_M1))


#M2

df_predict_M2 <- expand.grid(
  stade = stade_dense,
  Score = "Coef_M2" 
)

pred_M2 <- predict(model_M2, newdata = df_predict_pearson , 
                   se.fit = TRUE, interval = "confidence", level = 0.95)


df_predict_M2 <- df_predict_M2 %>%
  bind_cols(as_tibble(pred_M2))

#####################
#On plot tout le monde pour observer le comportement 
#général

par(mfrow=c(3,2), oma = c(0, 0, 4, 0))

plot(pred_pearson[["fit"]][,1],
     main = "Pearson")
plot(pred_Overlap[["fit"]][,1],
     main = "Overlap")
plot(pred_ICQ[["fit"]][,1],
     main = "ICQ")
plot(pred_M1[["fit"]][,1],
     main = "M1")
plot(pred_M2[["fit"]][,1],
     main = "M2")

##############################
#Graphique ggplot pearson coef
##############################
#on le mets avec les donnée

#pearson 

windowsFonts( A = windowsFont("baskerville old face"))

graph_pearson <- ggplot(data = Indice_colocalisation_tout_stade)+ 
  geom_boxplot(aes(x = as.factor(stade), y = Coef_pearson), color = "blue", alpha = 0.5, fill = NA) + 
  geom_point(aes(x = stade - 2, y = Coef_pearson), color = "blue", position = position_jitter(width = 0.1),
             size = 1) + 
  geom_ribbon(data = df_predict_pearson, 
              color = "blue",
              aes(x = stade -2 , ymin = fit[,"lwr"], ymax = fit[,"upr"]), 
              alpha = 0.1) +
  geom_line(data = df_predict_pearson,
            aes(stade -2, fit[,"lwr"]), color = "grey30", size = 0.1) + 
  geom_line(data = df_predict_pearson,
            aes(stade -2, fit[,"upr"]), color = "grey30", size = 0.1) + 
  geom_line(data = df_predict_pearson, 
            aes(x = stade -2 , y = fit[,"fit"]),
            linewidth = 1.2)   +
  
  scale_x_discrete("Stade",labels = c("Stade3","Stade4","Stade5","Stade6",
                                       "Stade7","Stade8","Stade9","Stade10"))+
  ylim(0, 1)+
  
  labs(
    x = "Egg chamber stages",
    y = "coefficient de Pearson",
    title = "Évolution du score de Pearson (r) des LncARN antisens de AAACTAT et AAACTAC",
    subtitle = "Les sondes utilisées étaient sens"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "right",
    text = element_text(size = 12, family = "A"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "A"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain", family = "A"),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  
  annotate("text", x = Inf, y = Inf, label = paste0("R²= ", round(data_r2[1,1], digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic")

#Overlap

graph_overlap <- ggplot(data = Indice_colocalisation_tout_stade)+ 
  geom_boxplot(aes(x = as.factor(stade), y = Coeaf_overlap), color = "orange", alpha = 0.5,
               fill = NA) + 
  geom_point(aes(x = stade - 2, y = Coeaf_overlap), color = "orange", 
             position = position_jitter(width = 0.1), size  = 1) + 
  geom_ribbon(data = df_predict_Overlap, 
              color = "orange",
              aes(x = stade -2 , ymin = fit[,"lwr"], ymax = fit[,"upr"]), 
              alpha = 0.1) +
  geom_line(data = df_predict_Overlap,
            aes(stade -2, fit[,"lwr"]), color = "grey30", size = 0.1) + 
  geom_line(data = df_predict_Overlap,
            aes(stade -2, fit[,"upr"]), color = "grey30", size = 0.1) +
  geom_line(data = df_predict_Overlap, 
            aes(x = stade -2 , y = fit),
            linewidth = 1.2)   +
  
  scale_x_discrete("Stade",labels = c("Stade3","Stade4","Stade5","Stade6",
                                       "Stade7","Stade8","Stade9","Stade10"))+
  ylim(0, 1)+
  
  labs(
    x = "Stade dévellopemental",
    y = "Coeficient de MOC",
    title = "Évolution du score de overlap (MOC) des LncARN antisens de AAACTAT et AAACTAC",
    subtitle = "Les sondes utilisées étaient sens"
  )+
  
  theme_bw() +
  theme(
    legend.position = "right",
    text = element_text(size = 12, family = "A"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "A"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain", family = "A"),
    plot.margin = margin(20, 20, 20, 20)
  )+
  
  annotate("text", x = Inf, y = Inf, label = paste0("R²= ", round(data_r2[1,2],digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic")

#ICQ

graph_ICQ <- ggplot(data = Indice_colocalisation_tout_stade)+ 
  geom_boxplot(aes(x = as.factor(stade), y = ICQ), color = "red", alpha = 0.5, fill = NA) + 
  geom_point(aes(x = stade - 2, y = ICQ),color = "red",
             position = position_jitter(width = 0.1), size = 1) + 
  geom_ribbon(data = df_predict_ICQ,
              color = "red",
              aes(x = stade -2 , ymin = fit[,"lwr"], ymax = fit[,"upr"]), 
              alpha = 0.1) +
  geom_line(data = df_predict_ICQ,
            aes(stade -2, fit[,"lwr"]), color = "grey30", size = 0.1) + 
  geom_line(data = df_predict_ICQ,
            aes(stade -2, fit[,"upr"]), color = "grey30", size = 0.1) +
  geom_line(data = df_predict_ICQ, 
            aes(x = stade -2 , y = fit[,"fit"]),
            linewidth = 1.2)    +
  
  annotate("text", x = Inf, y = Inf, label = paste0("R²= ", round(data_r2[1,3], digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") + 
  
  scale_x_discrete("Stages",labels = c("Stade3","Stade4","Stade5","Stade6",
                                       "Stade7","Stade8","Stade9","Stade10"))+
  ylim(-0.5, 0.5)+
  
  labs(
    x = "Egg chamber stages",
    y = "coefficient ICQ",
    title = "Évolution du score de ICQ des LncARN antisens de AAACTAT et AAACTAC",
    subtitle = "Les sondes utilisées étaient sens"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "right",
    text = element_text(size = 12, family = "A"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "A"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain", family = "A"),
    plot.margin = margin(20, 20, 20, 20)
  )


#menders

#M1 = CTAC, M2 = CTAT

data_M1_M2 <- Indice_colocalisation_tout_stade[,c(5,6,10)] %>%
  pivot_longer(
    cols = colnames(Indice_colocalisation_tout_stade[,c(5,6)]), 
    names_to = "ID_score",             
    values_to = "valeurs"
  )

transparence = 1

windowsFonts( A = windowsFont("baskerville old face"))

graph_menders <- ggplot(data_M1_M2, aes(x = factor(stade), y = valeurs,  color  = ID_score)) +
  
  geom_boxplot(aes(),
               fill = NA,
               outlier.shape = NA, 
               position = position_dodge(width = 0.8)) +
  
  geom_point(aes(),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.1),
             alpha = 1, size = 1) +
  
  scale_fill_manual(values = c("M2" = "magenta", "M1" = "green")) +
  scale_color_manual(values = c("M2" = "magenta", "M1" = "green")) +
  
  geom_ribbon(data = df_predict_M2, 
              aes(x = stade - 2, ymin = fit[,"lwr"], ymax = fit[,"upr"], y = NULL, fill = NULL), 
              alpha = 0.15, fill = "magenta", inherit.aes = FALSE) +
  geom_line(data = df_predict_M2,
            aes(stade -2, fit[,"lwr"]), 
            color = "darkmagenta", size = 0.1, inherit.aes = FALSE) + 
  
  geom_line(data = df_predict_M2,
            aes(stade -2, fit[,"upr"]), 
            color = "darkmagenta", size = 0.1, inherit.aes = FALSE) +
  
  geom_line(data = df_predict_M2, 
            aes(x = stade - 2, y = fit), 
            color = "magenta", linewidth = 1.2, inherit.aes = FALSE) +
  
  geom_ribbon(data = df_predict_M1, 
              aes(x = stade - 2, ymin = fit[,"lwr"], ymax = fit[,"upr"], y = NULL, fill = NULL), 
              alpha = 0.15, fill = "green", inherit.aes = FALSE) +
  
  geom_line(data = df_predict_M1,
            aes(stade -2, fit[,"lwr"]), color = "darkgreen", size = 0.1, inherit.aes = FALSE) + 
  
  geom_line(data = df_predict_M1,
            aes(stade -2, fit[,"upr"]), color = "darkgreen", size = 0.1, inherit.aes = FALSE) +
  
  geom_line(data = df_predict_M1, 
            aes(x = stade - 2, y = fit), 
            color = "green", linewidth = 1.2, inherit.aes = FALSE) +
  
  annotate("text", x = Inf, y = Inf, label = paste0("R²= ", round(data_r2[1,4], digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic", color = "green")+
  
  annotate("text", x = Inf, y = 0.95, label = paste0("R²= ", round(data_r2[1,5], digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic", color = "magenta") + 
  
  ylim(0, 1)+
  
  labs(
    x = "Egg chamber stages",
    y = "coefficient M1 & M2",
    title = "Évolution du score M1 & M2 des LncARN antisens de AAACTAT et AAACTAC",
    subtitle = "Les sondes utilisées étaient sens"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "right",
    text = element_text(size = 12, family = "A"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "A"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain", family = "A"),
    plot.margin = margin(20, 20, 20, 20)
  )



####################################################
#Impression des graphiques
###################################################

print(graph_pearson)
print(graph_overlap)
print(graph_ICQ)
print(graph_menders)

ggarrange(graph_pearson,graph_overlap,graph_ICQ,graph_menders)

#################################################
#################################################
#test de variance entre Y

data_concatener <- rbind(indice_colocalisation_stade10, indice_colocalisation_stade3,
                         indice_colocalisation_stade4, indice_colocalisation_stade5,
                         indice_colocalisation_stade6, indice_colocalisation_stade7, 
                         indice_colocalisation_stade8, indice_colocalisation_stade9)
data_concatener <- data_concatener[,c(-3,-4,-7,-8)]

pairs.panels(data_concatener[,c(1:5)])


par(mfrow=c(3,2), oma = c(0, 0, 3, 0))

hist(data_concatener$Coef_pearson,
     main = "Pearson sens")

hist(data_concatener$Coeaf_overlap,
     main = "Overlap sens")

hist(data_concatener$M1,
     main = "M1 sens")

hist(data_concatener$M2,
     main = "M2 sens")

hist(data_concatener$ICQ,
     main = "ICQ sens")

mtext("Sens", outer = TRUE, cex = 1.5, font = 2)

boxplot(data_concatener,
        main = " Sonde sens")

for ( i in c(1:5)){
  
  name = colnames(data_concatener)[i]
  
  nombre_0 <- (length(which(data_concatener[,i] == 0)))
  
  print(paste("Nombre de 0 pour",name , nombre_0))
  
}

############################

#on le mets avec les donnée

#pearson 


data_Per_ICQ<- Indice_colocalisation_tout_stade[,c(1,9,10)] %>%
  pivot_longer(
    cols = colnames(Indice_colocalisation_tout_stade[,c(1,9)]), 
    names_to = "ID_score",             
    values_to = "valeurs"
  )

plot_ICQ_per <- ggplot(data_Per_ICQ, aes(x = factor(stade), y = valeurs, fill = ID_score))+ 
  
  geom_boxplot(alpha = 0.5, 
               outlier.shape = NA, 
               position = position_dodge(width = 0.8)) +
  
  geom_point(color = "black",
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.1),
             alpha = 1, size = 1,
             shape = 1) +
  
  scale_fill_manual(values = c("Coef_pearson" = "blue", "ICQ" = "red")) +
  scale_color_manual(values = c("Coef_pearson" = "blue", "ICQ" = "red")) +
  
  geom_ribbon(data = df_predict_pearson, 
              aes(x = stade -2 , ymin = fit[,"lwr"], ymax = fit[,"upr"]), 
              alpha = 0.09, fill = "blue", inherit.aes = FALSE) +
  geom_line(data = df_predict_pearson,
            aes(stade -2, fit[,"lwr"]), color = "darkblue", size = 0.1,
            inherit.aes = FALSE) + 
  geom_line(data = df_predict_pearson,
            aes(stade -2, fit[,"upr"]), color = "darkblue", size = 0.1,
            inherit.aes = FALSE) + 
  geom_line(data = df_predict_pearson, 
            aes(x = stade -2 , y = fit[,"fit"]),
            linewidth = 1.2, color = "blue",inherit.aes = FALSE)+
  
  geom_ribbon(data = df_predict_ICQ, 
              aes(x = stade -2 , ymin = fit[,"lwr"], ymax = fit[,"upr"]), 
              alpha = 0.09, fill = "red",inherit.aes = FALSE) +
  geom_line(data = df_predict_ICQ,
            aes(stade -2, fit[,"lwr"]), color = "darkred", size = 0.1,
            inherit.aes = FALSE) + 
  geom_line(data = df_predict_ICQ,
            aes(stade -2, fit[,"upr"]), color = "darkred", size = 0.1,
            inherit.aes = FALSE) +
  geom_line(data = df_predict_ICQ, 
            aes(x = stade -2 , y = fit[,"fit"]),
            linewidth = 1.2, color = "red",inherit.aes = FALSE)    +
  
  scale_x_discrete("Stages",labels = c("Stade3","Stade4","Stade5","Stade6",
                                       "Stade7","Stade8","Stade9","Stade10"))+
  ylim(0, 1)+
  
  labs(
    x = "Stade de dévelopement des oeufs",
    y = "coeficient Pearson & ICQ",
    title = "Évolution de la colocalisaton des LncARN antisens de AAACTAT et AAACTAC",
    subtitle = "Les sondes utilisée étais sens"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "right",
    text = element_text(size = 12, family = "A"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "A"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain", family = "A"),
    plot.margin = margin(20, 20, 20, 20)
  )+
  
  annotate("text", x = Inf, y = Inf, label = paste0("R²= ", round(data_r2[1,3], digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic", color = "red") +
  annotate("text", x = Inf, y = 0.98, label = paste0("R²= ", round(data_r2[1,1], digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic", color = "blue")


#Grahique poster

transparence = 1

windowsFonts( A = windowsFont("baskerville old face"))

plot_mender_pres <- ggplot(data_M1_M2, aes(x = factor(stade), y = valeurs, color = ID_score, fill = NA)) +
  
  geom_boxplot(aes(),
               alpha = 0.7, 
               outlier.shape = NA, 
               position = position_dodge(width = 0.8)) +
  
  geom_point(aes(),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.1),
             alpha = 1, size = 1) +
  
  scale_fill_manual(values = c("M1" = "green", "M2" = "magenta")) +
  scale_color_manual(values = c("M1" = "green", "M2" = "magenta")) +
  
  geom_ribbon(data = df_predict_M1, 
              aes(x = stade - 2, ymin = fit[,"lwr"], ymax = fit[,"upr"], y = NULL, fill = NULL), 
              alpha = 0.15, fill = "green", inherit.aes = FALSE) +
  geom_line(data = df_predict_M1,
            aes(stade -2, fit[,"lwr"]), 
            color = "darkgreen", size = 0.1, inherit.aes = FALSE) + 
  
  geom_line(data = df_predict_M1,
            aes(stade -2, fit[,"upr"]), 
            color = "darkgreen", size = 0.1, inherit.aes = FALSE) +
  
  geom_line(data = df_predict_M1, 
            aes(x = stade - 2, y = fit[,"fit"]), 
            color = "green", linewidth = 1.2, inherit.aes = FALSE) +
  
  geom_ribbon(data = df_predict_M2, 
              aes(x = stade - 2, ymin = fit[,"lwr"], ymax = fit[,"upr"], y = NULL, fill = NULL), 
              alpha = 0.15, fill = "magenta", inherit.aes = FALSE) +
  
  geom_line(data = df_predict_M2,
            aes(stade -2, fit[,"lwr"]), color = "darkmagenta", size = 0.1, inherit.aes = FALSE) + 
  
  geom_line(data = df_predict_M2,
            aes(stade -2, fit[,"upr"]), color = "darkmagenta", size = 0.1, inherit.aes = FALSE) +
  
  geom_line(data = df_predict_M2, 
            aes(x = stade - 2, y = fit[,"fit"]), 
            color = "magenta", linewidth = 1.2, inherit.aes = FALSE)  +
  
  scale_x_discrete("Stages",labels = c("Stade3","Stade4","Stade5","Stade6",
                                       "Stade7","Stade8","Stade9","Stade10"))+
  ylim(0, 1)+
  
  labs(
    x = "Stade de dévelopement des oeufs",
    y = "coeficient M1 & M2",
    title = "Évolution de la colocalisaton des LncARN antisens de AAACTAT et AAACTAC",
    subtitle = "Les sondes utilisée étais sens"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "right",
    text = element_text(size = 12, family = "A"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "A"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain", family = "A"),
    plot.margin = margin(20, 20, 20, 20)
  )+
  
  annotate("text", x = Inf, y = Inf, label = paste0("R²= ", round(data_r2[1,4], digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic", color = "magenta")+
  
  annotate("text", x = Inf, y = 0.95, label = paste0("R²= ", round(data_r2[1,5], digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic", color = "green")

print(plot_ICQ_per)
print(plot_mender_pres)

#bocplot pour chaque score

par(mfrow=c(3,2), oma = c(0, 0, 3, 0))

boxplot(data = Indice_colocalisation_tout_stade, Coef_pearson ~ stade)
boxplot(data = Indice_colocalisation_tout_stade, Coeaf_overlap ~ stade)
boxplot(data = Indice_colocalisation_tout_stade, ICQ ~ stade)
boxplot(data = Indice_colocalisation_tout_stade, M1 ~ stade)
boxplot(data = Indice_colocalisation_tout_stade, M2 ~ stade)

mtext("Distribution de la variance en fonction du stade \n pour les sondes sens"
      , outer = TRUE, cex = 1.2, font = 1.5)

################################
#anova pour spline
###############################


#pearson
model_pearson_1 <- lm(data = Indice_colocalisation_tout_stade, Coef_pearson ~ stade)
model_pearson_2 <- lm(data = Indice_colocalisation_tout_stade, Coef_pearson ~ poly(stade,2))
model_pearson_3 <- lm(data = Indice_colocalisation_tout_stade, Coef_pearson ~ bs(stade,3))
model_pearson_4 <- lm(data = Indice_colocalisation_tout_stade, Coef_pearson ~ bs(stade,4))
model_pearson_5 <- lm(data = Indice_colocalisation_tout_stade, Coef_pearson ~ bs(stade,5))
model_pearson_6 <- lm(data = Indice_colocalisation_tout_stade, Coef_pearson ~ bs(stade,6))
model_pearson_7 <- lm(data = Indice_colocalisation_tout_stade, Coef_pearson ~ bs(stade,7))

anova(model_pearson_1,model_pearson_2, model_pearson_3, model_pearson_4, model_pearson_5, model_pearson_6, model_pearson_7)

#overlap


model_Overlap_1 <- rlm(data = Indice_colocalisation_tout_stade, Coeaf_overlap ~ stade)
model_Overlap_2 <- rlm(data = Indice_colocalisation_tout_stade, Coeaf_overlap ~ poly(stade,2))
model_Overlap_3 <- rlm(data = Indice_colocalisation_tout_stade, Coeaf_overlap ~ bs(stade,3))
model_Overlap_4 <- rlm(data = Indice_colocalisation_tout_stade, Coeaf_overlap ~ bs(stade,4))
model_Overlap_5 <- rlm(data = Indice_colocalisation_tout_stade, Coeaf_overlap ~ bs(stade,5))
model_Overlap_6 <- rlm(data = Indice_colocalisation_tout_stade, Coeaf_overlap ~ bs(stade,6))
model_Overlap_7 <- rlm(data = Indice_colocalisation_tout_stade, Coeaf_overlap ~ bs(stade,7))

anova(model_Overlap_1,model_Overlap_2, model_Overlap_3, model_Overlap_4, model_Overlap_5, model_Overlap_6, model_Overlap_7)


#ICQ

model_ICQ_1 <- rlm(data = Indice_colocalisation_tout_stade, ICQ ~ stade)
model_ICQ_2 <- rlm(data = Indice_colocalisation_tout_stade, ICQ ~ poly(stade,2))
model_ICQ_3 <- rlm(data = Indice_colocalisation_tout_stade, ICQ ~ bs(stade,3))
model_ICQ_4 <- rlm(data = Indice_colocalisation_tout_stade, ICQ ~ bs(stade,4))
model_ICQ_5 <- rlm(data = Indice_colocalisation_tout_stade, ICQ ~ bs(stade,5))
model_ICQ_6 <- rlm(data = Indice_colocalisation_tout_stade, ICQ ~ bs(stade,6))
model_ICQ_7 <- rlm(data = Indice_colocalisation_tout_stade, ICQ ~ bs(stade,7))

anova(model_ICQ_1,model_ICQ_2, model_ICQ_3, model_ICQ_4, model_ICQ_5, model_ICQ_6, model_ICQ_7)

#M1
model_M1_1 <- rlm(data = Indice_colocalisation_tout_stade, M1 ~  stade)
model_M1_2 <- rlm(data = Indice_colocalisation_tout_stade, M1 ~  poly(stade,2))
model_M1_3 <- rlm(data = Indice_colocalisation_tout_stade, M1 ~  bs(stade,3))
model_M1_4 <- rlm(data = Indice_colocalisation_tout_stade, M1 ~  bs(stade,4))
model_M1_5 <- rlm(data = Indice_colocalisation_tout_stade, M1 ~  bs(stade,5))
model_M1_6 <- rlm(data = Indice_colocalisation_tout_stade, M1 ~  bs(stade,6))
model_M1_7 <- rlm(data = Indice_colocalisation_tout_stade, M1 ~  bs(stade,7))

anova(model_M1_1,model_M1_2, model_M1_3, model_M1_4, model_M1_5, model_M1_6, model_M1_7)

#M2

model_M2_1 <- rlm(data = Indice_colocalisation_tout_stade, M1 ~  stade)
model_M2_2 <- rlm(data = Indice_colocalisation_tout_stade, M1 ~  poly(stade,2))
model_M2_3 <- rlm(data = Indice_colocalisation_tout_stade, M1 ~  bs(stade,3))
model_M2_4 <- rlm(data = Indice_colocalisation_tout_stade, M1 ~  bs(stade,4))
model_M2_5 <- rlm(data = Indice_colocalisation_tout_stade, M1 ~  bs(stade,5))
model_M2_6 <- rlm(data = Indice_colocalisation_tout_stade, M1 ~  bs(stade,6))
model_M2_7 <- rlm(data = Indice_colocalisation_tout_stade, M1 ~  bs(stade,7))

anova(model_M2_1,model_M2_2, model_M2_3, model_M2_4, model_M2_5, model_M2_6, model_M2_7)













