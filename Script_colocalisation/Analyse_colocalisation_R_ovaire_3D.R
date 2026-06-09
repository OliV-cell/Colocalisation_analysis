###################################
#Analyse coloc en 3D
#Dvir48 sonde antisens CTAC et CTAT
###################################

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
library(ggsignif)
#citation 
library(grateful)

#########################################################################
#ouverture des fichier contenant les indices, est sous forme txt dans csv

#Vecteur des nom des futur fichier utiliser

stade_3 = c(1, 4, 6, 10, 15, 22, 30, 33, 35, 48, 59)
stade_4 = c(2, 3, 17, 28, 44, 52, 58)
stade_5 = c(1, 3, 4, 7, 9, 11, 13, 15, 26, 29, 32, 35, 46, 47, 57)
stade_6 = c(2, 14,19, 25, 27, 38, 50, 53)
stade_7 = c(5, 8, 12, 20, 31, 34, 45)
stade_8 = c(21, 24, 23, 37, 39, 42, 49)
stade_9 = c(36, 40, 41, 51)
stade_10 = c(54, 55, 56)



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

#fonction ouverture 

ouverture <- function(n_serie,stade,nom_ficher){
  
  donne_concatene <- data.frame(Coef_pearson = rep(0,1),
                                
                                Coeaf_overlap = rep(0,1),
                                
                                k1  = rep(0,1),
                                
                                k2 = rep(0,1), 
                                
                                M1 = rep(0,1),
                                
                                M2 = rep(0,1),
                                
                                a_cytofluogram  = rep(0,1),
                                
                                b_cytofluogram = rep(0,1), 
                                
                                ICQ = rep(0,1),
                                
                                stringsAsFactors = FALSE)
  
  nom_des_fichier <- vector()
  
  for ( i in (1:length(n_serie))){
    
    nom_des_fichier[[i]] <- paste0("Coef_coloc_zstack_Series_",n_serie[i],"stade",stade)
    
  }
  
  
  for ( i in (1:length(n_serie))){
    
    ID <- read.csv(paste0("Coef_coloc_zstack_Series_",n_serie[i],"Dvir48_CTAC_CTAT_antisens_stade",stade,".csv"),
                   sep = ",", dec = ".", header =TRUE)
    
    ID <- Traitement_data_coloc(ID)
    
    donne_concatene[i,] <- ID
    
  }
  
  assign(nom_ficher, donne_concatene, .GlobalEnv)
  
}


for (i in c(3:10)){
  
  x <- get(paste0("stade_",i))
  
  ouverture(x,i,paste0("Coloc_stade_",i))
  
}

######################################
#data long des score combinés ensemble

Coloc_stade_3$stade <- rep(3,length(Coloc_stade_3[,1]))

Coloc_stade_4$stade <- rep(4,length(Coloc_stade_4[,1]))

Coloc_stade_5$stade <- rep(5,length(Coloc_stade_5[,1]))

Coloc_stade_6$stade <- rep(6,length(Coloc_stade_6[,1]))

Coloc_stade_7$stade <- rep(7,length(Coloc_stade_7[,1]))

Coloc_stade_8$stade <- rep(8,length(Coloc_stade_8[,1]))

Coloc_stade_9$stade <- rep(9,length(Coloc_stade_9[,1]))

Coloc_stade_10$stade <- rep(10,length(Coloc_stade_10[,1]))

Indice_colocalisation_tout_stade <- rbind(Coloc_stade_3,
                                          Coloc_stade_4,
                                          Coloc_stade_5,
                                          Coloc_stade_6,
                                          Coloc_stade_7,
                                          Coloc_stade_8,
                                          Coloc_stade_9,
                                          Coloc_stade_10)

Indice_colocalisation_tout_stade_lon <- Indice_colocalisation_tout_stade[,c(1,2,5,6,9,10)] %>%
  pivot_longer(
    cols = colnames(Indice_colocalisation_tout_stade[,c(1,2,5,6,9)]), 
    names_to = "ID_score",             
    values_to = "valeurs"
  )

###################################################################
#orbeta
hist(Indice_colocalisation_tout_stade$Coef_pearson)

#on essaie gausienne
hist(Indice_colocalisation_tout_stade$Coeaf_overlap)

#gausienne aussi
hist(Indice_colocalisation_tout_stade$ICQ)

#ordbeta
hist(Indice_colocalisation_tout_stade$M1)

#ordbeta
hist(Indice_colocalisation_tout_stade$M2)

boxplot(Indice_colocalisation_tout_stade[,c(1,2,5,6,9)])

data_concatener <- Indice_colocalisation_tout_stade[,c(-3,-4,-7,-8)]

for ( i in c(1:5)){
  
  name = colnames(data_concatener)[i]
  
  nombre_0 <- (length(which(data_concatener[,i] == 0)))
  
  print(paste("Nombre de 0 pour",name , nombre_0))
  
}

###############################
#on essaie des modèles et des fitting

boxplot(Indice_colocalisation_tout_stade$Coef_pearson ~ Indice_colocalisation_tout_stade$stade)


#################################
#pearson Gamma
#a voir, rlm meilleur fit en général sur les donnée.
model_pearson <- glmmTMB(data = Indice_colocalisation_tout_stade,
                         Coef_pearson ~ bs(stade,5),
                         family = ziGamma(link = "inverse"))

qqPlot(Indice_colocalisation_tout_stade$Coef_pearson)

summary(model_pearson)

mod_sim_M1 <- simulateResiduals(fittedModel = model_pearson, n = 250)

plot(mod_sim_M1 ,rank = F)

par(mfrow=c(2,2), oma = c(0, 0, 4, 0))

plotResiduals(mod_sim_M1)

testOutliers(mod_sim_M1)

testDispersion(mod_sim_M1) 

testZeroInflation(mod_sim_M1) 

#on test la linéarité
#####
newdat <- Indice_colocalisation_tout_stade_lon %>% 
  filter(ID_score == "Coef_pearson")
newdat$fit = fitted(model_pearson)
newdat$res = resid(model_pearson)

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


#tout est ok :)

#########################
#Overlap

model_overlap <- lm(Coeaf_overlap ~ bs(stade,3), data = Indice_colocalisation_tout_stade)

check_model(model_overlap)

shapiro.test(model_overlap$residuals)

check_heteroscedasticity(model_overlap)

mean(model_overlap$residuals)


#parfait, lm avec polynome de degré 6 :)

##########################
#ICQ

model_ICQ <- lm(ICQ ~ bs(stade,4), data = Indice_colocalisation_tout_stade)

check_model(model_ICQ)

shapiro.test(model_ICQ$residuals)

check_heteroscedasticity(model_ICQ)

mean(model_ICQ$residuals)


#tout est ok :)


########################
#M1
model_M1 <- glmmTMB(data = Indice_colocalisation_tout_stade,
                         M1 ~ bs(stade,5),
                         family = ordbeta(link = "logit"))

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
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('linéarité M1  antisens sur stade')

ggarrange(plot_1,plot_2)

#tout semble en ordre :)

#################################
#M2
model_M2 <- glmmTMB(data = Indice_colocalisation_tout_stade,
                    M2 ~ bs(stade,5),
                    family = ordbeta(link = "logit"))

summary(model_M2)

mod_sim_M1 <- simulateResiduals(fittedModel = model_M2, n = 250)

plot(mod_sim_M1 ,rank = F)

par(mfrow=c(2,2), oma = c(0, 0, 4, 0))

plotResiduals(mod_sim_M1)

testOutliers(mod_sim_M1)

testDispersion(mod_sim_M1) 

testZeroInflation(mod_sim_M1) 

#on test la linéarité

newdat <- Indice_colocalisation_tout_stade_lon %>% 
  filter(ID_score == "M2")
newdat$fit = fitted(model_M2)
newdat$res = resid(model_M2)

plot_1 <- ggplot(newdat, aes(fit, res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('linéarité M1 antisens')

plot_2 <- ggplot(newdat, aes(x = Indice_colocalisation_tout_stade$stade, res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('linéarité M1  antisens sur stade')

ggarrange(plot_1,plot_2)

#la linéarité de la réponse des résidue est étrange mais! ces chill

############################
#On lance les projections sur les 
#vrai donnée.


summary(model_pearson)

summary(model_overlap)

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


data_r2[1,1] <- r2_coxsnell(model_pearson)
data_r2[1,2] <- summary(model_overlap)$r.squared
data_r2[1,3] <- summary(model_ICQ)$r.squared
data_r2[1,4] <- r2_ferrari(model_M1)
data_r2[1,5] <- r2(model_M2)

predict_response(model_pearson, terms = "stade")
predict_response(model_overlap, terms = "stade")
predict_response(model_ICQ, terms = "stade")
predict_response(model_M1, terms = "stade [all]")
predict_response(model_M2, terms = "stade [all]")


#on augmente le détail de la prédiction
#afin de comblé les cathégorie discrète

stade_dense <- seq(min(Indice_colocalisation_tout_stade$stade), 
                   max(Indice_colocalisation_tout_stade$stade), 
                   length.out = 200)

stade_dense_pour_pearson <- seq(min(Indice_colocalisation_tout_stade$stade), 
                   max(Indice_colocalisation_tout_stade$stade), 
                   length.out = 200)

###########################
#Pearson

df_predict_pearson <- expand.grid(
  stade = stade_dense_pour_pearson,
  Score = "Coef_pearson" 
)

pred_pearson <- predict(model_pearson, newdata = df_predict_pearson , 
                   se.fit = TRUE, type = "link")

df_predict_pearson <- df_predict_pearson %>%
  mutate(
    fit   = model_pearson$modelInfo$family$linkinv(pred_pearson$fit), 
    upper = model_pearson$modelInfo$family$linkinv(pred_pearson$fit + (1.96 * pred_pearson$se.fit)),
    lower = model_pearson$modelInfo$family$linkinv(pred_pearson$fit - (1.96 * pred_pearson$se.fit))
  )


#Overlap

df_predict_Overlap <- expand.grid(
  stade = stade_dense,
  Score = "Coef_Overlap" 
)

pred_Overlap <- predict(model_overlap, newdata = df_predict_Overlap , 
                        se.fit = TRUE , interval = "confidence", level = 0.95)

df_predict_Overlap <- df_predict_Overlap %>%
  bind_cols(as_tibble(pred_Overlap))

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

pred_M2 <- predict(model_M2, newdata = df_predict_M2 , 
                   se.fit = TRUE, type = "link")

df_predict_M2 <- df_predict_M2 %>%
  mutate(
    fit   = model_M2$modelInfo$family$linkinv(pred_M2$fit), 
    upper = model_M2$modelInfo$family$linkinv(pred_M2$fit + (1.96 * pred_M2$se.fit)),
    lower = model_M2$modelInfo$family$linkinv(pred_M2$fit - (1.96 * pred_M2$se.fit))
  )

#####################################
#Ploting des graphiques


#####################
#On plot tout le monde pour observer le comportement 
#général

par(mfrow=c(3,2), oma = c(0, 0, 4, 0))

plot(pred_pearson$fit,
     main = "Pearson")
plot(pred_Overlap[["fit"]][,1],
     main = "Overlap")
plot(pred_ICQ[["fit"]][,1],
     main = "ICQ")
plot(pred_M1$fit,
     main = "M1")
plot(pred_M2$fit,
     main = "M2")

##############################
#Graphique ggplot pearson coef
##############################
#on le mets avec les donnée

#pearson 

windowsFonts( A = windowsFont("Arial"))

graph_pearson <- ggplot(data = Indice_colocalisation_tout_stade)+ 
  geom_boxplot(aes(x = as.factor(stade), y = Coef_pearson), color = "blue", alpha = 0.2, fill = NA) + 
  geom_point(aes(x = stade - 2, y = Coef_pearson), color = "blue", position = position_jitter(width = 0.1)) + 
  geom_ribbon(data = df_predict_pearson, 
              aes(x = stade -2 , ymin = lower, ymax = upper), 
              alpha = 0.2, inherit.aes = FALSE, fill = "blue") +
  geom_line(data = df_predict_pearson, 
            aes(x = stade -2 , y = fit),
            linewidth = 1.2, color = "black")   +
  geom_line(data = df_predict_pearson, 
            aes(x = stade -2 , y = upper),
            linewidth = 0.5, color = "darkblue")   +
  geom_line(data = df_predict_pearson, 
            aes(x = stade -2 , y = lower),
            linewidth = 0.5, color = "darkblue")   +
  
  scale_x_discrete("Stade",labels = c("Stade3","Stade4","Stade5","Stade6",
                                       "Stade7","Stade8","Stade9","Stade10"))+
  ylim(0, 1)+
  
  labs(
    x = "Egg chamber stages",
    y = "coeficient de Pearson (r)",
    title = "Évolution du coefficient de Pearson des LncARN sens de AAACTAT et AAACTAC",
    subtitle = "Les sondes utilisées étaient antisens"
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
  geom_boxplot(aes(x = as.factor(stade), y = Coeaf_overlap), color = "orange", alpha = 0.5, fill = NA) + 
  geom_point(aes(x = stade - 2, y = Coeaf_overlap), color = "orange", 
             position = position_jitter(width = 0.1)) + 
  geom_ribbon(data = df_predict_Overlap, 
              aes(x = stade -2 , ymin = fit[,"lwr"], ymax = fit[,"upr"]), 
              alpha = 0.2, inherit.aes = FALSE, fill = "orange") +
  geom_line(data = df_predict_Overlap, 
            aes(x = stade -2 , y = fit[,"fit"]),
            linewidth = 1.2, color = "black")   +
  geom_line(data = df_predict_Overlap, 
            aes(x = stade -2 , y = fit[,"lwr"]),
            linewidth = 0.5, color = "darkorange",
            inherit.aes = FALSE)   +
  geom_line(data = df_predict_Overlap, 
            aes(x = stade -2 , y = fit[,"upr"]),
            linewidth = 0.5, color = "darkorange",
            inherit.aes = FALSE)   +
  
  scale_x_discrete("Stade",labels = c("Stade3","Stade4","Stade5","Stade6",
                                       "Stade7","Stade8","Stade9","Stade10"))+
  ylim(0, 1)+
  
  labs(
    x = "Egg chamber stages",
    y = "coeficient de MOC",
    title = "Évolution du coefficient de MOC des LncARN sens de AAACTAT et AAACTAC",
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
  
  annotate("text", x = Inf, y = Inf, label = paste0("R²= ", round(data_r2[1,2],digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic")

#ICQ

graph_ICQ <- ggplot(data = Indice_colocalisation_tout_stade)+ 
  geom_boxplot(aes(x = as.factor(stade), y = ICQ), color = "red", alpha = 0.5,fill = NA) + 
  geom_point(aes(x = stade - 2, y = ICQ),color = "red",
             position = position_jitter(width = 0.1)) + 
  geom_ribbon(data = df_predict_ICQ, 
              aes(x = stade -2 , ymin = fit[,"lwr"], ymax = fit[,"upr"]), 
              alpha = 0.2, fill = "red",inherit.aes = FALSE) +
  geom_line(data = df_predict_ICQ, 
            aes(x = stade -2 , y = fit[,"fit"]),
            linewidth = 1.2, color = "black")    +
  geom_line(data = df_predict_ICQ, 
            aes(x = stade -2 , y = fit[,"lwr"]),
            linewidth = 0.5, color = "darkred",
            inherit.aes = FALSE)   +
  geom_line(data = df_predict_ICQ, 
            aes(x = stade -2 , y = fit[,"upr"]),
            linewidth = 0.5, color = "darkred",
            inherit.aes = FALSE)   +
  geom_hline(yintercept = 0, linetype = "dotted",
             size = 1) + 
  
  annotate("text", x = Inf, y = Inf, label = paste0("R²= ", round(data_r2[1,3], digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") + 
  
  scale_x_discrete("Egg chamber stages",labels = c("Stade3","Stade4","Stade5","Stade6",
                                       "Stade7","Stade8","Stade9","Stade10"))+
  ylim(-0.5, 0.5)+
  
  labs(
    x = "Egg chamber stages",
    y = "ICQ score",
    title = "Evolution of the ICQ score through development \n of oocyte for the forward transcript AAACTAT et AAACTAC",
    subtitle = "Both probe were reversed"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "right",
    text = element_text(size = 12, family = "A"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "A"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain", family = "A"),
    plot.margin = margin(20, 20, 20, 20),,
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold")
  )


#menders

#M1 = CTAT, M2 = CTAC

data_M1_M2 <- Indice_colocalisation_tout_stade[,c(5,6,10)] %>%
  pivot_longer(
    cols = colnames(Indice_colocalisation_tout_stade[,c(5,6)]), 
    names_to = "ID_score",             
    values_to = "valeurs"
  )

transparence = 1

windowsFonts( A = windowsFont("Arial"))

graph_menders <- ggplot(data_M1_M2, aes(x = factor(stade), y = valeurs, color = ID_score)) +
  
  geom_boxplot(aes(),
               fill = NA,
               alpha = 0.4, 
               outlier.shape = NA, 
               position = position_dodge(width = 0.8)) +
  
  geom_point(aes(),
             fill = NA,
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.1),
             alpha = 1, size = 1) +
  
  scale_fill_manual(values = c("M1" = "magenta", "M2" = "green")) +
  
  scale_color_manual(values = c("M1" = "magenta", "M2" = "green")) +
  
  geom_ribbon(data = df_predict_M2, 
              aes(x = stade - 2, ymin = lower, ymax = upper, y = NULL, fill = NULL), 
              alpha = 0.15, fill = "green", inherit.aes = FALSE) +
  
  geom_line(data = df_predict_M2, 
            aes(x = stade - 2, y = fit), 
            color = "green", linewidth = 1.2, inherit.aes = FALSE) +
  geom_line(data = df_predict_M2, 
            aes(x = stade - 2, y = upper), 
            color = "darkgreen", linewidth = 0.5, inherit.aes = FALSE) +
  geom_line(data = df_predict_M2, 
            aes(x = stade - 2, y = lower), 
            color = "darkgreen", linewidth = 0.5, inherit.aes = FALSE) +
  
  geom_ribbon(data = df_predict_M1, 
              aes(x = stade - 2, ymin = lower, ymax = upper, y = NULL, fill = NULL), 
              alpha = 0.15, fill = "magenta", inherit.aes = FALSE) +
  geom_line(data = df_predict_M1, 
            aes(x = stade - 2, y = fit), 
            color = "magenta", linewidth = 1.2, inherit.aes = FALSE)  +
  geom_line(data = df_predict_M1, 
            aes(x = stade - 2, y = upper), 
            color = "darkmagenta", linewidth = 0.5, inherit.aes = FALSE)  +
  geom_line(data = df_predict_M1, 
            aes(x = stade - 2, y = lower), 
            color = "darkmagenta", linewidth = 0.5, inherit.aes = FALSE)  +
  
  annotate("text", x = Inf, y = Inf, label = paste0("R²= ", round(data_r2[1,4],digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic", 
           color = "magenta") + 
  
  annotate("text", x = Inf, y = 0.95, label = paste0("R²= ", round(data_r2[1,5],digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic", 
           color = "green") + 
  
  scale_x_discrete("Stages",labels = c("Stade3","Stade4","Stade5","Stade6",
                                       "Stade7","Stade8","Stade9","Stade10"))+
  ylim(0, 1)+
  
  labs(
    x = "Egg chamber stages",
    y = "coefficient de M1 & M2",
    title = "Évolution du coefficient de M1 & M2 des LncARN sens de AAACTAT et AAACTAC",
    subtitle = "Les sondes utilisées étaient antisens"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "right",
    text = element_text(size = 12, family = "A"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "A"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain", family = "A"),
    plot.margin = margin(20, 20, 20, 20)
  )


print(graph_pearson)
print(graph_overlap)
print(graph_ICQ)
print(graph_menders)

ggarrange(graph_pearson,graph_overlap,graph_ICQ,graph_menders)


#########################################################



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
  
  geom_ribbon(data = df_predict_pearson, 
              aes(x = stade -2 , ymin = lower, ymax = upper), 
              alpha = 0.09, fill = "blue", inherit.aes = FALSE) +
  geom_line(data = df_predict_pearson,
            aes(stade -2, lower), color = "darkblue", size = 0.1,
            inherit.aes = FALSE) + 
  geom_line(data = df_predict_pearson,
            aes(stade -2, upper), color = "darkblue", size = 0.1,
            inherit.aes = FALSE) + 
  geom_line(data = df_predict_pearson, 
            aes(x = stade -2 , y = fit),
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
    y = "coeficient p & ICQ",
    title = "Évolution de la colocalisaton 3D des LncARN sens de AAACTAT et AAACTAC",
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
  
  scale_fill_manual(values = c("Coef_pearson" = "blue", "ICQ" = "red")) +
  scale_color_manual(values = c("Coef_pearson" = "blue", "ICQ" = "red")) +
  
  annotate("text", x = Inf, y = Inf, label = paste0("R²= ", round(data_r2[1,3], digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic", color = "red") +
  annotate("text", x = Inf, y = 0.98, label = paste0("R²= ", round(data_r2[1,1], digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic", color = "blue")

print(plot_ICQ_per)


#########################################
#Donnée de Volume DIANA

#stade 9 et 10

#Vecteur des nom des futur fichier utiliser

setwd("C:/Users/verme/Desktop/Data_ovaire/Dvir48/Analyse_coloc_3D_sonde_antisens_Dvir48_CTAC_CTAT/Data_graphique_3D")

nom_des_data_9 <- vector()
nom_des_data_10 <- vector()

for ( i in seq_along(stade_9)){
  
  nom_des_data_9[[i]] <- paste0("Volume_stade_9_",i)
  
}

for ( i in seq_along(stade_10)){
  
  nom_des_data_10[[i]] <- paste0("Volume_stade_10_",i)
  
}

#ouverture de 9

for ( i in seq_along(stade_9)) {
  
  n_image <- stade_9[i]
  ID <- read.csv(paste0("Resultat_DIANA_Series_",n_image,"Dvir48_CTAC_CTAT_antisens_stade9.csv"),
                              sep = ",", dec = ".", header =TRUE)
  assign(nom_des_data_9[i],ID)
  
}

for ( i in seq_along(stade_10)) {
  
  n_image <- stade_10[i]
  
  ID <- read.csv(paste0("Resultat_DIANA_Series_",n_image,"Dvir48_CTAC_CTAT_antisens_stade10.csv"),
                 sep = ",", dec = ".", header =TRUE)
  assign(nom_des_data_10[i],ID)
  
}

Volume_stade_9 <- rbind(Volume_stade_9_1,Volume_stade_9_2,Volume_stade_9_3,Volume_stade_9_4)
Volume_stade_9$Stade <- rep(9,4)

Volume_stade_10 <- rbind(Volume_stade_10_1,Volume_stade_10_2,Volume_stade_10_3)
Volume_stade_10$Stade <- rep(10,3)

Volume_9_10 <- rbind(Volume_stade_10,Volume_stade_9)

data_Volume_9_10<- Volume_9_10[,c(3,4,5,9)] %>%
  pivot_longer(
    cols = colnames(Volume_9_10[,c(3,4,5)]), 
    names_to = "Mesure",             
    values_to = "Valeurs"
  )

##################
#graphique

ggplot(data = data_Volume_9_10, aes(x = as.factor(Stade), y = Valeurs, fill = Mesure)) + 
  
  scale_fill_manual(values = c("ColocFromAvolume" = "magenta", 
                               "ColocFromBvolume" = "green",
                               "ColocFromABvolume" = "red")) +
  
  geom_boxplot(alpha = 0.6) + 
  geom_point(color = "black",
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.1),
             alpha = 1, size = 2) + 
  labs(
    y = "Pourcentage",
    x = "Stade"
  ) +
  
  scale_y_continuous(breaks = seq(0, 110, by = 10),
                     limits = c(0, 110))+
  
  theme_bw() +
  
  theme(
    legend.position = "right",
    text = element_text(size = 12, family = "A"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "A"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain", family = "A"),
    plot.margin = margin(20, 20, 20, 20)
  )+
  
  annotate("rect", xmin = as.factor(10), xmax = Inf, ymin = 100, ymax = 110,
           alpha = 1, fill = "white") + 
  
annotate("text", x = as.factor(9), y = 105, label = "**",size = 4, fontface = "italic")+

annotate("text", x = as.factor(10), y = 105, label = "**",size = 4, fontface = "italic")


anov_test <- aov(data = data_Volume_9_10,anov_test <- aov(data = dStadeata_Volume_9_10,
                 Valeurs ~ Mesure*Stade))

shapiro.test(anov_test$residuals)

summary(anov_test)

#test de comparaison paramétrique

tukey_volume <- TukeyHSD(anov_test)

print(anov_test)

print(tukey_volume)



