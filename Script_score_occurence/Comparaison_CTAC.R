###########################
#Script analyse comparaison
#CTAC avec Virlis, américana
#novamexicana et lumeil
################################

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
library(marginaleffects)
library(FSA)

library(gvlma)
library(stargazer)
library(performance)
library(see)
#représentation graphique
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(pammtools)

library(Distance)
library(dsm)

########################################################
#ouverture des fichier

data_dlumei <- read.csv("Donner_Dlumei_CTAC_antisens.csv",
                        sep = ";", dec = ",", header =TRUE)

data_dvir  <- read.csv("Donner_Dvir_CTAC_antisens.csv",
                        sep = ";", dec = ",", header =TRUE)

data_dnova  <- read.csv("Donner_Dnova_CTAC_antisens.csv",
                       sep = ";", dec = ",", header =TRUE)

data_damer  <- read.csv("Donner_Damer_CTAC_antisens.csv",
                       sep = ";", dec = ",", header =TRUE)


#traitement Dlumei
decompte <- function(data_source, data_puit){
  for (i in c(3:10)){
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

Dlumei_rep1 <- data_dlumei[1:37,]
  
Dlumei_rep2 <- data_dlumei[38:98,]

Dlumei_1 <- matrix( nrow = 8, ncol = 4)

Dlumei_2 <- matrix( nrow = 8, ncol = 4)

Dlumei_1 <- decompte(data_source = data_dlumei,data_puit = Dlumei_1)
Dlumei_2 <- decompte(data_source = data_dlumei,data_puit = Dlumei_2)

noms_matrices <- c("Dlumei_1","Dlumei_2")

total_observations <- sum(sapply(mget(noms_matrices), sum, na.rm = TRUE))

print(total_observations)

matrice_fusionnee <- Reduce("+", mget(noms_matrices))

comptage_par_ligne <- rowSums(matrice_fusionnee, na.rm = TRUE)

noms_stades <- paste("Stade", 3:(3 + length(comptage_par_ligne) - 1))
names(comptage_par_ligne) <- noms_stades


print(comptage_par_ligne)

data_score_lumei <- data.frame(score_1 = rep(0,8),
                               
                               score_2 = rep(0,8),
                               
                               stade = c(3:10),
                               
                               row.names = c("stade 3",
                                             "stade 4","stade 5","stade 6",
                                             "stade 7","stade 8","stade 9",
                                             "stade 10"))

for ( i in c(1:8)){
  
  data_score_lumei[[i,1]] <- (Dlumei_1[i,2] + 2*Dlumei_1[i,3] + 3*Dlumei_1[i,4])
  
  data_score_lumei[[i,2]] <- (Dlumei_2[i,2] + 2*Dlumei_2[i,3] + 3*Dlumei_2[i,4])
  
  print(paste("Boucle",i,"fait"))
}

data_score_lumei_normalisé <- data_score_lumei

for (i in c(1:8)){
  
  data_score_lumei_normalisé[[i,1]] <- data_score_lumei_normalisé[i,1]/(sum(Dlumei_1[i,1:4]) +1)
  data_score_lumei_normalisé[[i,2]] <- data_score_lumei_normalisé[i,2]/(sum(Dlumei_2[i,1:4]) +1)
  
}

data_dlumei <- data_score_lumei_normalisé

data_dlumei[nrow(data_dlumei) + 1, ] <- list(1,0)

data_dlumei[nrow(data_dlumei) + 1, ] <- list(1,0)

data_dlumei[9:10, 1] <- 0

data_dlumei$sp <- rep("Dlumei", length(data_dlumei$score_1))

data_dlumei <- data_dlumei %>%
  pivot_longer(
    cols = !stade & !sp,
    names_to = "brin",
    values_to = "score"
  )

#explo data

plot(data_dvir$score)

plot(data_dnova$score)

plot(data_damer$score)

plot(data_dlumei$score)

hist(data_dvir$score)

hist(data_dnova$score)

hist(data_damer$score)

hist(data_dlumei$score)

data_concaténé <- rbind(data_dvir, data_dnova, data_damer, data_dlumei)

boxplot(data_concaténé$score ~ data_concaténé$sp)
leveneTest(data_concaténé$score ~ data_concaténé$sp)

#essaie avec GAM

data_concaténé$sp <- as.factor(data_concaténé$sp)

modele_combiné <- gam(list(score ~ sp + s(stade, by = sp, k = 10), 
    ~ s(stade)), 
  data = data_concaténé,
  optimizer = c("outer", "newton"),
  method = "REML", 
  family = gaulss(),
  control = list(maxit = 1000)
)

#Fonction choix du modèle

par(mfrow=c(2,2))

plot(modele_combiné, 
     pages = 1, 
     scheme = 1, 
     all.terms = TRUE)

summary.gam(modele_combiné)

par(mfrow=c(2,2))

gam.check(modele_combiné)

k.check(modele_combiné)

performance(modele_combiné)

concurvity(modele_combiné, full = FALSE)

vis_concurvity(modele_combiné, type = "estimate")

model_performance(modele_combiné)

stade_dense <- seq(min(data_concaténé$stade), 
                   max(data_concaténé$stade), 
                   length.out = 200)

df_predict_combine <- expand.grid(
  stade = stade_dense,
  sp = unique(data_concaténé$sp) 
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

data_concaténé$preds <- predict(modele_combiné, type = "response")[,1]

tableau_R2 <- data_concaténé %>%
  group_by(sp) %>%
  summarize(
    RSS = sum((score - preds)^2),                  
    TSS = sum((score - mean(score))^2),           
    R2  = round(1 - (RSS / TSS), 4)               
  )

print(tableau_R2)

######################
#statistique de différence général, par stade et pente d'expression

anova <- anova(modele_combiné)
anova.gam(modele_combiné)

p_val_brin <- anova$pTerms.pv[1]

noms_legendes <- c(
  "Dvir48" = "D.Virilis",
  "Dnova"  = "D.novamexicana",
  "Damer" = "D.americana",
  "Dlumei" = "D.lummei"
)

#########################

windowsFonts( A = windowsFont("Arial"))

#représentation graphique
plot_GAM <- ggplot() + 
  
  geom_point(data = data_concaténé, 
             aes(x = stade, y = score, color = sp),
             alpha = 0.6, position = position_jitter(width = 0.1)) + 
  geom_ribbon(data = df_predict_combine, 
              aes(x = stade, ymin = lower, ymax = upper, fill = sp), 
              alpha = 0.09) +
  geom_line(data = df_predict_combine, 
            aes(x = stade, y = fit, color = sp),
            linewidth = 1.2) +
  
  scale_color_manual(labels = noms_legendes,
                     values = c("Dvir48" = "red", 
                                "Dnova" = "blue", 
                                "Damer" = "orange",
                                "Dlumei" = "darkgreen")) +
  
  scale_fill_manual(labels = noms_legendes,
                    values = c("Dvir48" = "red", 
                               "Dnova" = "blue", 
                               "Damer" = "orange",
                               "Dlumei" = "darkgreen")) +
  
  labs(
    x = "Stade de dévelopement des oeufs",
    y = "Score d'expression normalisé",
    title = "Évolution de la transcription des LncARN sens de AAACTAC pour \n D.Virilis, D.americana, D.nocamexicana & D.lummei",
    subtitle = "Les sondes utilisé était antisens",
    color = "Espèces",
    fill = "Espèces"
  ) +
    
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  ylim(c(-1,4.5)) + 
  
  guides(fill = "none") +
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
  
  annotate("rect", xmin = 9.2, xmax = Inf, ymin = 2.8, ymax = 4.5,
           alpha = 1, fill = "white") +
  
  #annotate("text", x = Inf, y = Inf, label = paste0("Global ANOVA: p = ", format.pval(p_val_brin, digits = 2)), 
           #hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") + 
  
  annotate("text", x = Inf, y = 4.4, label = paste0("R² = ", round(tableau_R2[4,4], digits = 3), "***"), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic",
           color = "red") + 
  
  annotate("text", x = Inf, y = 4, label = paste0("R² = ", round(tableau_R2[3,4], digits = 3), "***" ), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic",
           color = "blue") +
  
  annotate("text", x = Inf, y = 3.6, label = paste0("R² = ", round(tableau_R2[1,4], digits = 3) , "***"), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic",
           color = "orange") + 
  
  annotate("text", x = Inf, y = 3.2, label = paste0("R² = ", round(tableau_R2[2,4], digits = 3) , "***"), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic",
           color = "darkgreen")
  

print(plot_GAM)

##########################################
#comparaison des courbe avec tets non-paramétrique 
#manque de donnée

tableau_pvalues <- data.frame()

for (s in 3:10) {
  donnees_stade <- subset(data_concaténé, stade == s)
  
  if(length(unique(donnees_stade$sp)) > 1) {
    
    fit <- aov(score ~ sp, data = donnees_stade)
    
    tukey <- TukeyHSD(fit)$sp

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

h_crochet <- 3

plot_GAM +
  
  #stade3-4-5 
  
  annotate("rect", xmin = 2.5, xmax = 5.2, ymin = 0, ymax = 3,
           linetype = "dashed", color = "black", fill = NA) +
  
  annotate("text", x = 4, y = 3 + 0.3, 
           label = "***", size = 5) +
  
  
  #stade7
  
  annotate("segment", x = 6.8, xend = 7.2, y = h_crochet, yend = h_crochet, 
           linetype = "dashed", color = "black") +
  
  annotate("segment", x = 6.8, xend = 6.8, y = 0.49, yend = h_crochet, 
           linetype = "dashed", color = "black") +
  
  annotate("segment", x = 7.2, xend = 7.2, y = 1.025, yend = h_crochet, 
           linetype = "dashed", color = "black") +
  
  annotate("text", x = 7, y = h_crochet + 0.3, 
           label = signif_data$label[4], size = 5) 
s






