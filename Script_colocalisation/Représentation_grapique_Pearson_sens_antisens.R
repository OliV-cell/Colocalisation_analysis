###########################################
#Script graphique pearson
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

############################################
load("data_model_sens.RData")
load("data_model_antisens.RData")

load("data_antisens_score.RData")
df_predict_pearson_antisens <- df_predict_pearson
Indice_colocalisation_tout_stade_antisens <- Indice_colocalisation_tout_stade
load("data_sens_score.RData")

Indice_colocalisation_tout_stade_antisens$Sens <- rep("antisens", length(Indice_colocalisation_tout_stade_antisens$stade))

Indice_colocalisation_tout_stade$Sens <- rep("sens", length(Indice_colocalisation_tout_stade$stade))

Indice_colocalisation_tout_stade_tout <- rbind(Indice_colocalisation_tout_stade_antisens, Indice_colocalisation_tout_stade)

data_Per_sens_antisens<- Indice_colocalisation_tout_stade_tout[,c(1,10,11)]

ggplot(Indice_colocalisation_tout_stade_antisens, aes(x = factor(stade), y = Coef_pearson, fill = Sens))+ 
  
  geom_boxplot(alpha = 0.2, 
               outlier.shape = NA, 
               position = position_dodge(width = 0.8)) +
  
  geom_point(color = "black",
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.1),
             alpha = 1, size = 1,
             shape = 1) +
  
  geom_ribbon(data = df_predict_pearson_antisens, 
              aes(x = stade -2 , ymin = fit[,"lwr"], ymax = fit[,"upr"]), 
              alpha = 0.3, fill = "blue", inherit.aes = FALSE) +
  geom_line(data = df_predict_pearson_antisens,
            aes(stade -2, fit[,"lwr"]), color = "darkblue", size = 0.1,
            inherit.aes = FALSE) + 
  geom_line(data = df_predict_pearson_antisens,
            aes(stade -2, fit[,"upr"]), color = "darkblue", size = 0.1,
            inherit.aes = FALSE) + 
  geom_line(data = df_predict_pearson_antisens, 
            aes(x = stade -2 , y = fit[,"fit"]),
            linewidth = 1.2, color = "blue",inherit.aes = FALSE)+
  
  scale_x_discrete("Stade de dévelopement des oeufs",labels = c("Stade3","Stade4","Stade5","Stade6",
                                       "Stade7","Stade8","Stade9","Stade10"))+
  ylim(0, 1)+
  
  labs(
    x = "Stade de dévelopement des oeufs",
    y = "coeficient de Pearson (r)",
    title = "Évolution de la colocalisaton des LncARN sens de AAACTAT et AAACTAC",
    subtitle = "Se sont le sens des sondes"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "right",
    text = element_text(size = 12, family = "A"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "A"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain", family = "A"),
    plot.margin = margin(20, 20, 20, 20)
  )+
  
  scale_fill_manual(values = c("antisens" = "blue")) +
  scale_color_manual(values = c("antisens" = "blue")) +
  
  annotate("text", x = Inf, y = 0.98, label = "R² = 0.039", 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic", color = "blue")


















