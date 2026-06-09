###########################################
#Script analyse résultats colocalisation
#CTAC+CTAT sonde antisens (brin sens)
##########################################

library(feasts)
library(tsibble)
library(tidyverse)
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
#representation graphique
library(ggeffects)
library(ggplot2)
library(ggpubr)

#citation 
library(grateful)

mes_packages <- c(
  "feasts", "tsibble", "tidyverse", "car",  "dplyr", 
  "agricolae", "coin", "multcomp", "rcompanion", "psych", 
  "nlstools", "lmtest", "tidyr", "SuppDists", "DHARMa", "splines", 
  "gamm4", "mgcv", "gvlma", "performance", "see", 
  "glmmTMB", "MASS", "ggeffects", "ggplot2", "ggpubr", "grateful"
)

# 2. Appel de la fonction de citation en utilisant ce vecteur
library(grateful)
cite_packages(
  pkgs = mes_packages,
  output = "file", 
  out.format = "docx", 
  out.dir = ".", 
  cite.tidyverse = TRUE,
  dependencies = FALSE
)

#################################
#importation des données
score_coloc_stade3_rep1 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_rep1_stade3.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stade4_rep1 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antosensboth_stade4_rep1.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stade5_rep1 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antosensboth_stade5_rep1.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stade6_rep1 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade6_rep1.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stade7_rep1 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade7_rep1.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stade8_rep1 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade8_rep1.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stade9_rep1 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade9_rep1.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stade10_rep1 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade10_rep1.csv",
                                    sep = ",", dec = ".", header =TRUE)

score_coloc_stade3_rep2 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade3_rep2.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stade4_rep2 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade4_rep2.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stade5_rep2 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade5_rep2.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stade6_rep2 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade6_rep2.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stade7_rep2 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade7_rep2.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stade8_rep2 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade8_rep2.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stade9_rep2 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade9_rep2.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stade10_rep2 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade10_rep2.csv",
                                     sep = ",", dec = ".", header =TRUE)




#Fonction de traitement et des donnée, ou, il serons classer en ordre dans un beau
#tableau clair pour faire des graphiques de comparaison de moyenne.

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



indice_colocalisation_stade3_rep1 <- Traitement_data_coloc(score_coloc_stade3_rep1)
indice_colocalisation_stade4_rep1 <- Traitement_data_coloc(score_coloc_stade4_rep1)
indice_colocalisation_stade5_rep1 <- Traitement_data_coloc(score_coloc_stade5_rep1)
indice_colocalisation_stade6_rep1 <- Traitement_data_coloc(score_coloc_stade6_rep1)
indice_colocalisation_stade7_rep1 <- Traitement_data_coloc(score_coloc_stade7_rep1)
indice_colocalisation_stade8_rep1 <- Traitement_data_coloc(score_coloc_stade8_rep1)
indice_colocalisation_stade9_rep1 <- Traitement_data_coloc(score_coloc_stade9_rep1)
indice_colocalisation_stade10_rep1 <- Traitement_data_coloc(score_coloc_stade10_rep1)

indice_colocalisation_stade3_rep2 <- Traitement_data_coloc(score_coloc_stade3_rep2)
indice_colocalisation_stade4_rep2 <- Traitement_data_coloc(score_coloc_stade4_rep2)
indice_colocalisation_stade5_rep2 <- Traitement_data_coloc(score_coloc_stade5_rep2)
indice_colocalisation_stade6_rep2 <- Traitement_data_coloc(score_coloc_stade6_rep2)
indice_colocalisation_stade7_rep2 <- Traitement_data_coloc(score_coloc_stade7_rep2)
indice_colocalisation_stade8_rep2 <- Traitement_data_coloc(score_coloc_stade8_rep2)
indice_colocalisation_stade9_rep2 <- Traitement_data_coloc(score_coloc_stade9_rep2)
indice_colocalisation_stade10_rep2 <- Traitement_data_coloc(score_coloc_stade10_rep2)

indice_colocalisation_stade3 <- rbind(indice_colocalisation_stade3_rep1,
                                      indice_colocalisation_stade3_rep2)

indice_colocalisation_stade4 <- rbind(indice_colocalisation_stade4_rep1,
                                      indice_colocalisation_stade4_rep2)

indice_colocalisation_stade5 <- rbind(indice_colocalisation_stade5_rep1,
                                      indice_colocalisation_stade5_rep2)

indice_colocalisation_stade6 <- rbind(indice_colocalisation_stade6_rep1,
                                      indice_colocalisation_stade6_rep2)

indice_colocalisation_stade7 <- rbind(indice_colocalisation_stade7_rep1,
                                      indice_colocalisation_stade7_rep2)

indice_colocalisation_stade8 <- rbind(indice_colocalisation_stade8_rep1,
                                      indice_colocalisation_stade8_rep2)

indice_colocalisation_stade9 <- rbind(indice_colocalisation_stade9_rep1,
                                      indice_colocalisation_stade9_rep2)

indice_colocalisation_stade10 <- rbind(indice_colocalisation_stade10_rep1,
                                      indice_colocalisation_stade10_rep2)

######################################
#data long des score combinés ensemble

indice_colocalisation_stade3$stade <- rep(3,length(indice_colocalisation_stade3[,1]))

indice_colocalisation_stade4$stade <- rep(4,length(indice_colocalisation_stade4[,1]))

indice_colocalisation_stade5$stade <- rep(5,length(indice_colocalisation_stade5[,1]))

indice_colocalisation_stade6$stade <- rep(6,length(indice_colocalisation_stade6[,1]))

indice_colocalisation_stade7$stade <- rep(7,length(indice_colocalisation_stade7[,1]))

indice_colocalisation_stade8$stade <- rep(8,length(indice_colocalisation_stade8[,1]))

indice_colocalisation_stade9$stade <- rep(9,length(indice_colocalisation_stade9[,1]))

indice_colocalisation_stade10$stade <- rep(10,length(indice_colocalisation_stade10[,1]))


Indice_colocalisation_tout_stade <- rbind(indice_colocalisation_stade3,
                                          indice_colocalisation_stade4,
                                          indice_colocalisation_stade5,
                                          indice_colocalisation_stade6,
                                          indice_colocalisation_stade7,
                                          indice_colocalisation_stade8,
                                          indice_colocalisation_stade9,
                                          indice_colocalisation_stade10)

Indice_colocalisation_tout_stade_lon <- Indice_colocalisation_tout_stade[,c(1,2,5,6,9,10)] %>%
  pivot_longer(
    cols = colnames(Indice_colocalisation_tout_stade[,c(1,2,5,6,9)]), 
    names_to = "ID_score",             
    values_to = "valeurs"
  )

str(Indice_colocalisation_tout_stade_lon)

#####################################
#Test modele LM
###################################

#confection de model préléminaire pour tester les conditions d'application
model_pearson_1 <- rlm(data = Indice_colocalisation_tout_stade, Coef_pearson ~ stade)
model_pearson_2 <- rlm(data = Indice_colocalisation_tout_stade, Coef_pearson ~ poly(stade,6))
model_pearson_3 <- rlm(data = Indice_colocalisation_tout_stade, Coef_pearson ~ bs(stade,3))
model_pearson_4 <- rlm(data = Indice_colocalisation_tout_stade, Coef_pearson ~ bs(stade,4))
model_pearson_5 <- rlm(data = Indice_colocalisation_tout_stade, Coef_pearson ~ bs(stade,5))
model_pearson_6 <- rlm(data = Indice_colocalisation_tout_stade, Coef_pearson ~ bs(stade,6))
model_pearson_7 <- rlm(data = Indice_colocalisation_tout_stade, Coef_pearson ~ bs(stade,7))

anova(model_pearson_1,model_pearson_2,model_pearson_3,model_pearson_4,model_pearson_5,model_pearson_6,model_pearson_7)

summary(model_pearson_6)

mean(model_pearson_6$residuals)

dwtest(model_pearson_6)

check_heteroscedasticity(model_pearson_6)

cor.test(Indice_colocalisation_tout_stade$stade, model_pearson_6$residuals)

model_performance(model_pearson_6)

shapiro.test(model_pearson_6$residuals)

qqPlot(resid(model_pearson_6),
       main = "Pearson")

newdat <- Indice_colocalisation_tout_stade_lon %>% 
  filter(ID_score == "Coef_pearson")
newdat$fit <- fitted(model_pearson_6)
newdat$res <- resid(model_pearson_6)
newdat$weights <- model_pearson_6$w

plot_1 <- ggplot(newdat, aes(x = fit, y = res)) +
  geom_point(aes(color = weights, alpha = weights),
             size = 4) +
  scale_color_gradient(low = "cyan", high = "blue") +
  geom_smooth(color = "black", se = TRUE, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'red4', size = 1) +
  labs(
    title = 'Pearson',
    x = 'Valeurs ajustées (Fitted)',
    y = 'Résidus',
    color = 'Poids'
  ) +
  theme_minimal()

#sur les valeurs x exactes
plot_2 <- ggplot(newdat, aes(x = Indice_colocalisation_tout_stade$stade, y = res)) +
  geom_point(aes(color = weights, alpha = weights),
             size = 4) +
  scale_color_gradient(low = "cyan", high = "blue") +
  geom_smooth(color = "black", se = TRUE, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'red4', size = 1) +
  labs(
    title = 'Pearson sur x stade',
    x = 'Valeurs ajustées (Fitted)',
    y = 'Résidus',
    color = 'Poids'
  ) +
  theme_minimal()

ggarrange(plot_1,plot_2)

#pearson ok
######################

#test on coupe les 0 

for (i in c(1:length(Indice_colocalisation_tout_stade$Coeaf_overlap))) {
  
  if (Indice_colocalisation_tout_stade$Coeaf_overlap[i] == 0){
    
    Indice_colocalisation_tout_stade$Coeaf_overlap[i] <- NA
    
  }
  
}

model_Overlap <- glmmTMB(Coeaf_overlap ~ bs(stade,5), 
                         family = ordbeta(link = "logit"), 
                         data = Indice_colocalisation_tout_stade)
#on test les conditions
summary(model_Overlap)

mod_sim_overlap <- simulateResiduals(fittedModel = model_Overlap, n = 250)

par(mfrow=c(2,2), oma = c(0, 0, 4, 0))

plot(mod_sim_overlap ,rank = F)

plotResiduals(mod_sim_overlap)

testOutliers(mod_sim_overlap)

testDispersion(mod_sim_overlap) 

testZeroInflation(mod_sim_overlap) 

#on test la linéarité

newdat <- Indice_colocalisation_tout_stade_lon %>% 
  filter(ID_score == "Coeaf_overlap",
         valeurs != 0)

newdat$fit = fitted(model_Overlap)
newdat$res = resid(model_Overlap)


plot_1 <- ggplot(newdat, aes(fit, res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('linéarité overlap antisens')

#sur les valeurs x exacte

plot_2 <- ggplot(newdat, aes(x = stade, y = res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('linéarité overlap antisens sur stade')

ggarrange(plot_1,plot_2)

dwtest(model_Overlap)
lmtest::bptest(model_Overlap)

######################

model_ICQ_1 <- lm(data = Indice_colocalisation_tout_stade, ICQ ~ stade)
model_ICQ_2 <- lm(data = Indice_colocalisation_tout_stade, ICQ ~ poly(stade, 2))
model_ICQ_3 <- lm(data = Indice_colocalisation_tout_stade, ICQ ~ bs(stade, 3))
model_ICQ_4 <- lm(data = Indice_colocalisation_tout_stade, ICQ ~ bs(stade, 4))
model_ICQ_5 <- lm(data = Indice_colocalisation_tout_stade, ICQ ~ bs(stade, 5))
model_ICQ_6 <- lm(data = Indice_colocalisation_tout_stade, ICQ ~ bs(stade, 6))
model_ICQ_7 <- lm(data = Indice_colocalisation_tout_stade, ICQ ~ bs(stade, 7))

anova(model_ICQ_1,model_ICQ_2,model_ICQ_3,model_ICQ_4,model_ICQ_5,model_ICQ_6,model_ICQ_7)

model_ICQ <- lm(data = Indice_colocalisation_tout_stade, ICQ ~ bs(stade, 7))

#on test les conditions

check_model(model_ICQ)

shapiro.test(model_ICQ$residuals)

check_heteroscedasticity(model_ICQ)

shapiro.test(model_ICQ$residuals)

qqPlot(resid(model_ICQ))

mean(model_ICQ$residuals)

dwtest(model_ICQ)

cor.test(Indice_colocalisation_tout_stade$stade, model_ICQ$residuals)


newdat <- Indice_colocalisation_tout_stade_lon %>% 
  filter(ID_score == "ICQ")
newdat$fit = fitted(model_ICQ)
newdat$res = resid(model_ICQ )

plot_1 <- ggplot(newdat, aes(fit, res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('linéarité ICQ antisens')

#sur les vrai x

plot_2 <- ggplot(newdat, aes(x = Indice_colocalisation_tout_stade$stade, res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('linéarité ICQ antisens sur stade')

ggarrange(plot_1,plot_2)

######################

model_M1_1 <- glmmTMB(M1 ~ stade, 
                    family = ordbeta(link = "logit"), 
                    data = Indice_colocalisation_tout_stade)

model_M1_2 <- glmmTMB(M1 ~ poly(stade,2), 
                    family = ordbeta(link = "logit"), 
                    data = Indice_colocalisation_tout_stade)

model_M1_3<- glmmTMB(M1 ~ bs(stade,3), 
                    family = ordbeta(link = "logit"), 
                    data = Indice_colocalisation_tout_stade)

model_M1_4 <- glmmTMB(M1 ~ bs(stade,4), 
                    family = ordbeta(link = "logit"), 
                    data = Indice_colocalisation_tout_stade)

model_M1_5 <- glmmTMB(M1 ~ bs(stade,5), 
                    family = ordbeta(link = "logit"), 
                    data = Indice_colocalisation_tout_stade)

model_M1_6 <- glmmTMB(M1 ~ bs(stade,6), 
                    family = ordbeta(link = "logit"), 
                    data = Indice_colocalisation_tout_stade)

model_M1_7 <- glmmTMB(M1 ~ bs(stade,7), 
                    family = ordbeta(link = "logit"), 
                    data = Indice_colocalisation_tout_stade)

anova(model_M1_1,model_M1_2,model_M1_3,model_M1_4,model_M1_5,model_M1_6,model_M1_7)

model_M1 <- glmmTMB(M1 ~ bs(stade,5), 
                    family = ordbeta(link = "logit"), 
                    data = Indice_colocalisation_tout_stade)

summary(model_M1)

mod_sim_M1 <- simulateResiduals(fittedModel = model_M1, n = 250)

plot(mod_sim_M1 ,rank = F)

par(mfrow=c(2,2), oma = c(0, 0, 4, 0))

plotResiduals(mod_sim_M1)

testOutliers(mod_sim_M1)

testDispersion(mod_sim_M1) 

testZeroInflation(mod_sim_M1) 

#on test la linéarité

newdat <- Indice_colocalisation_tout_stade_lon %>% 
  filter(ID_score == "M1")
newdat$fit = fitted(model_M1)
newdat$res = resid(model_M1)

plot_1 <- ggplot(newdat, aes(fit, res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('linéarité M1 antisens')

plot_2 <- ggplot(newdat, aes(x = Indice_colocalisation_tout_stade$stade, res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('linéarité M1 antisens sur stade')

ggarrange(plot_1,plot_2)

######################

model_M2_1 <- glmmTMB(M2 ~ stade, 
                      family = ordbeta(link = "logit"), 
                      data = Indice_colocalisation_tout_stade)

model_M2_2 <- glmmTMB(M2 ~ poly(stade,2), 
                      family = ordbeta(link = "logit"), 
                      data = Indice_colocalisation_tout_stade)

model_M2_3<- glmmTMB(M2 ~ bs(stade,3), 
                     family = ordbeta(link = "logit"), 
                     data = Indice_colocalisation_tout_stade)

model_M2_4 <- glmmTMB(M2 ~ bs(stade,4), 
                      family = ordbeta(link = "logit"), 
                      data = Indice_colocalisation_tout_stade)

model_M2_5 <- glmmTMB(M2 ~ bs(stade,5), 
                      family = ordbeta(link = "logit"), 
                      data = Indice_colocalisation_tout_stade)

model_M2_6 <- glmmTMB(M2 ~ bs(stade,6), 
                      family = ordbeta(link = "logit"), 
                      data = Indice_colocalisation_tout_stade)

model_M2_7 <- glmmTMB(M2 ~ bs(stade,7), 
                      family = ordbeta(link = "logit"), 
                      data = Indice_colocalisation_tout_stade)

anova(model_M2_1,model_M2_2,model_M2_3,model_M2_4,model_M2_5,model_M2_6,model_M2_7)


model_beta_M2 <- glmmTMB(M2 ~ bs(stade, 7), 
                      family = ordbeta(link = "logit"), 
                      data = Indice_colocalisation_tout_stade)

summary(model_beta_M2)

mod_sim_M2 <- simulateResiduals(fittedModel = model_beta_M2, n = 250)

plot(mod_sim_M2 ,rank = F)

par(mfrow=c(2,2), oma = c(0, 0, 4, 0))

plotResiduals(mod_sim_M2)

testOutliers(mod_sim_M2)

testDispersion(mod_sim_M2) 

testZeroInflation(mod_sim_M2) 

#on test la linéarité

newdat <- Indice_colocalisation_tout_stade_lon %>% 
  filter(ID_score == "M2")
newdat$fit = fitted(model_beta_M2)
newdat$res = resid(model_beta_M2)

plot_1 <- ggplot(newdat, aes(fit, res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('linéarité M2 antisens')

plot_2 <- ggplot(newdat, aes(x = Indice_colocalisation_tout_stade$stade, res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('linéarité M2 antisens sur stade')


ggarrange(plot_1,plot_2)

#on replace les 0 

for (i in c(1:length(Indice_colocalisation_tout_stade$Coeaf_overlap))) {
  
  if (is.na(Indice_colocalisation_tout_stade$Coeaf_overlap[i])){
    
    Indice_colocalisation_tout_stade$Coeaf_overlap[i] <- 0
    
  }
  
}

####################################
#Prédiction de chaque model
#ainsi que la démonstation des effets
####################################

summary(model_pearson_6)

summary(model_Overlap)

summary(model_ICQ)

summary(model_M1)

summary(model_beta_M2)


#calcule des R2 et speudo R2

data_r2 <- data.frame("Pearson" = c(0,0),
                      "Overlap" = c(0,0),
                      "ICQ" = c(0,0),
                      "M1" = c(0,0),
                      "M2" = c(0,0),
                      row.names = c("R^2","P.value"))


data_r2[1,1] <- (cor(model_pearson_6$model$Coef_pearson, predict(model_pearson_6)))^2
data_r2[1,2] <- efronRSquared(model_Overlap)
data_r2[1,3] <- summary(model_ICQ)$r.squared
data_r2[1,4] <- efronRSquared(model_M1)
data_r2[1,5] <- efronRSquared(model_beta_M2)

pred_pearson <- predict_response(model_pearson_6, terms = "stade")
pred_overlap <- predict_response(model_Overlap, terms = "stade [3:10]")
pred_ICQ <- predict_response(model_ICQ, terms = "stade")
pred_M1 <- predict_response(model_M1, terms = "stade [all]")
pred_M2 <- predict_response(model_beta_M2, terms = "stade [all]")

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

annotate_figure(plot_fusion, top = text_grob("Prediction score Antisens", 
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

pred_pearson <- predict(model_pearson_6, newdata = df_predict_pearson , 
                        se.fit = TRUE, interval = "confidence", level = 0.95)

df_predict_pearson <- df_predict_pearson %>%
  bind_cols(as_tibble(pred_pearson))


#Overlap

df_predict_Overlap <- expand.grid(
  stade = stade_dense,
  Score = "Coef_Overlap" 
)

pred_Overlap <- predict(model_Overlap, newdata = df_predict_Overlap , 
                        se.fit = TRUE , type = "link")

df_predict_Overlap <- df_predict_Overlap %>%
  mutate(
  fit   = model_Overlap$modelInfo$family$linkinv(pred_Overlap$fit), 
  upper = model_Overlap$modelInfo$family$linkinv(pred_Overlap$fit + (1.96 * pred_Overlap$se.fit)),
  lower = model_Overlap$modelInfo$family$linkinv(pred_Overlap$fit - (1.96 * pred_Overlap$se.fit))
)

#ICQ

df_predict_ICQ <- expand.grid(
  stade = stade_dense,
  Score = "Coef_ICQ" 
)

pred_ICQ <- predict(model_ICQ, newdata = df_predict_ICQ , 
                        se.fit = TRUE, interval = "confidence", level = 0.95)

df_predict_ICQ <- df_predict_ICQ %>%
  bind_cols(as_tibble(pred_ICQ))


#M1

df_predict_M1 <- expand.grid(
  stade = stade_dense,
  Score = "Coef_M1" 
)

pred_M1 <- predict(model_M1, newdata = df_predict_M1 , 
                        se.fit = TRUE, type = "link")

df_predict_M1 <- df_predict_M1 %>%
  mutate(
    fit   = model_M1$modelInfo$family$linkinv(pred_M1$fit), 
    upper = model_M1$modelInfo$family$linkinv(pred_M1$fit + (1.96 * pred_M1$se.fit)),
    lower = model_M1$modelInfo$family$linkinv(pred_M1$fit - (1.96 * pred_M1$se.fit))
  )


#M2

df_predict_M2 <- expand.grid(
  stade = stade_dense,
  Score = "Coef_M2" 
)

pred_M2 <- predict(model_beta_M2, newdata = df_predict_M2 , 
                        se.fit = TRUE, type = "link")

df_predict_M2 <- df_predict_M2 %>%
  mutate(
    fit   = model_beta_M2$modelInfo$family$linkinv(pred_M2$fit), 
    upper = model_beta_M2$modelInfo$family$linkinv(pred_M2$fit + (1.96 * pred_M2$se.fit)),
    lower = model_beta_M2$modelInfo$family$linkinv(pred_M2$fit - (1.96 * pred_M2$se.fit))
  )

#####################
#On plot tout le monde pour observer le comportement 
#général

par(mfrow=c(3,2), oma = c(0, 0, 4, 0))

plot(pred_pearson[["fit"]][,1],
     main = "Pearson")
plot(pred_Overlap$fit,
     main = "Overlap")
plot(pred_ICQ[["fit"]][,1],
     main = "ICQ")
plot(pred_M1$fit,
     main = "M1")
plot(pred_M2$fit,
     main = "M2")

#on le mets avec les donnée

#pearson 

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
    y = "Coeficient de Pearson (r)",
    title = "Évolution do score de Pearson (r) des LncARN sens de AAACTAT et AAACTAC",
    subtitle = "Les sondes utilisées étaient antisens"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "right",
    text = element_text(size = 12, family = "A"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "A"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain", family = "A"),
    plot.margin = margin(20, 20, 20, 20)
  )+
  
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
              aes(x = stade -2 , ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_line(data = df_predict_Overlap,
            aes(stade -2, lower), color = "grey30", size = 0.1) + 
  geom_line(data = df_predict_Overlap,
            aes(stade -2, upper), color = "grey30", size = 0.1) +
  geom_line(data = df_predict_Overlap, 
            aes(x = stade -2 , y = fit),
            linewidth = 1.2)   +
  
  scale_x_discrete("Stade",labels = c("Stade3","Stade4","Stade5","Stade6",
                                       "Stade7","Stade8","Stade9","Stade10"))+
  ylim(0,1)+
  
  labs(
    x = "Stade dévellopemental",
    y = "Coeficient de MOC",
    title = "Évolution du score de overlap (MOC) des LncARN sens de AAACTAT et AAACTAC",
    subtitle = "Les sondes utilisées étaient antisens"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "right",
    text = element_text(size = 12, family = "A"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "A"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain", family = "A"),
    plot.margin = margin(20, 20, 20, 20)
  )+
  
  annotate("text", x = Inf, y = Inf, label = paste0("R²= ", round(data_r2[1,2], digits = 3)), 
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
  
  scale_x_discrete("Stade",labels = c("Stade3","Stade4","Stade5","Stade6",
                                       "Stade7","Stade8","Stade9","Stade10"))+
  ylim(-0.5, 0.5)+
  
  labs(
    x = "Egg chamber stages",
    y = "coeficient de ICQ",
    title = "Évolution du score de ICQ  des LncARN sens de AAACTAT et AAACTAC",
    subtitle = "Les sondes utilisées étaient antisens"
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
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic")



######################
#Tableau contenant les moyennes de chaque valeurs

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

moyenne_coloc_coef_stade3 <- Calcule_moyene(indice_colocalisation_stade3)
moyenne_coloc_coef_stade4 <- Calcule_moyene(indice_colocalisation_stade4)
moyenne_coloc_coef_stade5 <- Calcule_moyene(indice_colocalisation_stade5)
moyenne_coloc_coef_stade6 <- Calcule_moyene(indice_colocalisation_stade6)
moyenne_coloc_coef_stade7 <- Calcule_moyene(indice_colocalisation_stade7)
moyenne_coloc_coef_stade8 <- Calcule_moyene(indice_colocalisation_stade8)
moyenne_coloc_coef_stade9 <- Calcule_moyene(indice_colocalisation_stade9)
moyenne_coloc_coef_stade10 <- Calcule_moyene(indice_colocalisation_stade10)


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


##############################
#Graphique ggplot Menders coef
##############################

#on refet un dataframe avec les 2 valeurs d'intéret

#M1 = CTAT, M2 = CTAC

data_M1_M2 <- Indice_colocalisation_tout_stade[,c(5,6,10)] %>%
  pivot_longer(
    cols = colnames(Indice_colocalisation_tout_stade[,c(5,6)]), 
    names_to = "ID_score",             
    values_to = "valeurs"
  )

transparence = 1

windowsFonts( A = windowsFont("baskerville old face"))

graph_menders <- ggplot(data_M1_M2, aes(x = factor(stade), y = valeurs, color = ID_score, fill - NA)) +

  geom_boxplot(alpha = 0.7, 
               outlier.shape = NA, 
               position = position_dodge(width = 0.8)) +
  
  geom_point(aes(),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.1),
             alpha = 1, size = 1) +
  
  scale_fill_manual(values = c("M1" = "magenta", "M2" = "green")) +
  scale_color_manual(values = c("M1" = "magenta", "M2" = "green")) +

  geom_ribbon(data = df_predict_M2, 
              aes(x = stade - 2, ymin = lower, ymax = upper, y = NULL, fill = NULL), 
              alpha = 0.15, fill = "green", inherit.aes = FALSE) +
  geom_line(data = df_predict_M2,
            aes(stade -2, lower), 
            color = "darkgreen", size = 0.1, inherit.aes = FALSE) + 
  
  geom_line(data = df_predict_M2,
            aes(stade -2, upper), 
            color = "darkgreen", size = 0.1, inherit.aes = FALSE) +
  
  geom_line(data = df_predict_M2, 
            aes(x = stade - 2, y = fit), 
            color = "green", linewidth = 1.2, inherit.aes = FALSE) +

  geom_ribbon(data = df_predict_M1, 
              aes(x = stade - 2, ymin = lower, ymax = upper, y = NULL, fill = NULL), 
              alpha = 0.15, fill = "magenta", inherit.aes = FALSE) +
  
  geom_line(data = df_predict_M1,
            aes(stade -2, lower), color = "darkmagenta", size = 0.1, inherit.aes = FALSE) + 
  
  geom_line(data = df_predict_M1,
            aes(stade -2, upper), color = "darkmagenta", size = 0.1, inherit.aes = FALSE) +
  
  geom_line(data = df_predict_M1, 
            aes(x = stade - 2, y = fit), 
            color = "magenta", linewidth = 1.2, inherit.aes = FALSE)  +
  
  scale_x_discrete("Stages",labels = c("Stade3","Stade4","Stade5","Stade6",
                                       "Stade7","Stade8","Stade9","Stade10"))+
  ylim(0, 1)+
  
  labs(
    x = "Egg chamber stages",
    y = "Coefficient de M1 & M2",
    title = "Évolution du score M1 & M2 des LncARN sens de AAACTAT et AAACTAC",
    subtitle = "Les sondes utilisées étaient antisens"
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
  
####################################################
#Impression des graphiques boxplot
###################################################

print(graph_pearson)
print(graph_overlap)
print(graph_ICQ)
print(graph_menders)

ggarrange(graph_pearson,graph_overlap,graph_ICQ,graph_menders)

#################################################
#test de variance entre Y

par(mfrow=c(4,2), oma = c(0, 0, 2, 0))

data_concatener <- rbind(indice_colocalisation_stade10, indice_colocalisation_stade3,
                         indice_colocalisation_stade4, indice_colocalisation_stade5,
                         indice_colocalisation_stade6, indice_colocalisation_stade7, 
                         indice_colocalisation_stade8, indice_colocalisation_stade9)
data_concatener <- data_concatener[,c(-3,-4,-7,-8)]

pairs.panels(data_concatener[,c(1:5)])

par(mfrow=c(3,2), oma = c(0, 0, 3, 0))

hist(data_concatener$Coef_pearson,
     main = "Pearson antisens")

hist(data_concatener$Coeaf_overlap,
     main = "Overlap antisens")

hist(data_concatener$M1,
     main = "M1 antisens")

hist(data_concatener$M2,
     main = "M2 antisens")

hist(data_concatener$ICQ,
     main = "ICQ antisens")

mtext("Antisens", outer = TRUE, cex = 1.5, font = 2)


library(reshape2)

boxplot(data_concatener,
        main = "sonde antisens")

for ( i in c(1:5)){
  
  name = colnames(data_concatener)[i]
  
  nombre_0 <- (length(which(data_concatener[,i] == 0)))
  
  print(paste("Nombre de 0 pour",name , nombre_0))
  
}

#bocplot pour chaque score

par(mfrow=c(3,2), oma = c(0, 0, 3, 0))

boxplot(data = Indice_colocalisation_tout_stade, Coef_pearson ~ stade)
boxplot(data = Indice_colocalisation_tout_stade, Coeaf_overlap ~ stade)
boxplot(data = Indice_colocalisation_tout_stade, ICQ ~ stade)
boxplot(data = Indice_colocalisation_tout_stade, M1 ~ stade)
boxplot(data = Indice_colocalisation_tout_stade, M2 ~ stade)

mtext("Distribution de la variance en fonction du stade \n pour les sondes antisens"
      , outer = TRUE, cex = 1.2, font = 1.5)

##############################################
#vérification ANOVA spline
#############################################

#pearson
model_pearson_1 <- rlm(data = Indice_colocalisation_tout_stade, Coef_pearson ~ stade)
model_pearson_2 <- rlm(data = Indice_colocalisation_tout_stade, Coef_pearson ~ poly(stade,6))
model_pearson_3 <- rlm(data = Indice_colocalisation_tout_stade, Coef_pearson ~ bs(stade,3))
model_pearson_4 <- rlm(data = Indice_colocalisation_tout_stade, Coef_pearson ~ bs(stade,4))
model_pearson_5 <- rlm(data = Indice_colocalisation_tout_stade, Coef_pearson ~ bs(stade,5))
model_pearson_6 <- rlm(data = Indice_colocalisation_tout_stade, Coef_pearson ~ bs(stade,6))
model_pearson_7 <- rlm(data = Indice_colocalisation_tout_stade, Coef_pearson ~ bs(stade,7))

anova(model_pearson_1,model_pearson_2,model_pearson_3,model_pearson_4,model_pearson_5,model_pearson_6,model_pearson_7)


#overlap
model_Overlap_1 <- glmmTMB(Coeaf_overlap ~ stade, 
                         family = ordbeta(link = "logit"), 
                         data = Indice_colocalisation_tout_stade)
model_Overlap_2 <- glmmTMB(Coeaf_overlap ~ poly(stade,2), 
                         family = ordbeta(link = "logit"), 
                         data = Indice_colocalisation_tout_stade)
model_Overlap_3 <- glmmTMB(Coeaf_overlap ~ bs(stade,3), 
                         family = ordbeta(link = "logit"), 
                         data = Indice_colocalisation_tout_stade)
model_Overlap_4 <- glmmTMB(Coeaf_overlap ~ bs(stade,4), 
                         family = ordbeta(link = "logit"), 
                         data = Indice_colocalisation_tout_stade)
model_Overlap_5 <- glmmTMB(Coeaf_overlap ~ bs(stade,5), 
                         family = ordbeta(link = "logit"), 
                        data = Indice_colocalisation_tout_stade)
model_Overlap_6 <- glmmTMB(Coeaf_overlap ~ bs(stade,6), 
                           family = ordbeta(link = "logit"), 
                           data = Indice_colocalisation_tout_stade)
model_Overlap_7 <- glmmTMB(Coeaf_overlap ~ bs(stade,7), 
                           family = ordbeta(link = "logit"), 
                           data = Indice_colocalisation_tout_stade)

anova(model_Overlap_1,model_Overlap_2,model_Overlap_3,model_Overlap_4,model_Overlap_5,model_Overlap_6,model_Overlap_7)



#ICQ
model_ICQ_1 <- lm(data = Indice_colocalisation_tout_stade, ICQ ~ stade)
model_ICQ_2 <- lm(data = Indice_colocalisation_tout_stade, ICQ ~ poly(stade, 2))
model_ICQ_3 <- lm(data = Indice_colocalisation_tout_stade, ICQ ~ bs(stade, 3))
model_ICQ_4 <- lm(data = Indice_colocalisation_tout_stade, ICQ ~ bs(stade, 4))
model_ICQ_5 <- lm(data = Indice_colocalisation_tout_stade, ICQ ~ bs(stade, 5))
model_ICQ_6 <- lm(data = Indice_colocalisation_tout_stade, ICQ ~ bs(stade, 6))
model_ICQ_7 <- lm(data = Indice_colocalisation_tout_stade, ICQ ~ bs(stade, 7))

anova(model_ICQ_1,model_ICQ_2,model_ICQ_3,model_ICQ_4,model_ICQ_5,model_ICQ_6,model_ICQ_7)


#M1

model_M1_1 <- glmmTMB(M1 ~ stade, 
                      family = ordbeta(link = "logit"), 
                      data = Indice_colocalisation_tout_stade)

model_M1_2 <- glmmTMB(M1 ~ poly(stade,2), 
                      family = ordbeta(link = "logit"), 
                      data = Indice_colocalisation_tout_stade)

model_M1_3<- glmmTMB(M1 ~ bs(stade,3), 
                     family = ordbeta(link = "logit"), 
                     data = Indice_colocalisation_tout_stade)

model_M1_4 <- glmmTMB(M1 ~ bs(stade,4), 
                      family = ordbeta(link = "logit"), 
                      data = Indice_colocalisation_tout_stade)

model_M1_5 <- glmmTMB(M1 ~ bs(stade,5), 
                      family = ordbeta(link = "logit"), 
                      data = Indice_colocalisation_tout_stade)

model_M1_6 <- glmmTMB(M1 ~ bs(stade,6), 
                      family = ordbeta(link = "logit"), 
                      data = Indice_colocalisation_tout_stade)

model_M1_7 <- glmmTMB(M1 ~ bs(stade,7), 
                      family = ordbeta(link = "logit"), 
                      data = Indice_colocalisation_tout_stade)

anova(model_M1_1,model_M1_2,model_M1_3,model_M1_4,model_M1_5,model_M1_6,model_M1_7)


#M2

model_M2_1 <- glmmTMB(M2 ~ stade, 
                      family = ordbeta(link = "logit"), 
                      data = Indice_colocalisation_tout_stade)

model_M2_2 <- glmmTMB(M2 ~ poly(stade,2), 
                      family = ordbeta(link = "logit"), 
                      data = Indice_colocalisation_tout_stade)

model_M2_3<- glmmTMB(M2 ~ bs(stade,3), 
                     family = ordbeta(link = "logit"), 
                     data = Indice_colocalisation_tout_stade)

model_M2_4 <- glmmTMB(M2 ~ bs(stade,4), 
                      family = ordbeta(link = "logit"), 
                      data = Indice_colocalisation_tout_stade)

model_M2_5 <- glmmTMB(M2 ~ bs(stade,5), 
                      family = ordbeta(link = "logit"), 
                      data = Indice_colocalisation_tout_stade)

model_M2_6 <- glmmTMB(M2 ~ bs(stade,6), 
                      family = ordbeta(link = "logit"), 
                      data = Indice_colocalisation_tout_stade)

model_M2_7 <- glmmTMB(M2 ~ bs(stade,7), 
                      family = ordbeta(link = "logit"), 
                      data = Indice_colocalisation_tout_stade)

anova(model_M2_1,model_M2_2,model_M2_3,model_M2_4,model_M2_5,model_M2_6,model_M2_7)























