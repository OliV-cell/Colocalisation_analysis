###########################################
# Colocalization results analysis script
# CTAC+CTAT antisense probe (sense strand)
##########################################

library(tidyverse)
# various statistical tests

library(car)
library(dplyr)
library(coin)
library(rcompanion)
library(psych)
library(nlstools)
library(tidyr)

# Regression
library(DHARMa)
library(splines)
library(gamm4)
library(mgcv)
library(gvlma)
library(performance)
library(see)
library(MASS)
library(lmtest)
# graphical representation
library(ggplot2)
library(ggpubr)
library(ggeffects)

# Citation 
library(grateful)

mes_packages <- c(
  "tidyverse", "car",  "dplyr", 
  "coin", "rcompanion", "psych", 
  "nlstools", "tidyr", "splines", 
  "gamm4", "mgcv", "gvlma", "performance", "see", 
  "MASS", "ggplot2", "ggpubr", "grateful"
)

# Call citation function using this vector
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
# Import data
score_coloc_stage3_rep1 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_rep1_stade3.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stage4_rep1 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antosensboth_stade4_rep1.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stage5_rep1 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antosensboth_stade5_rep1.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stage6_rep1 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade6_rep1.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stage7_rep1 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade7_rep1.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stage8_rep1 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade8_rep1.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stage9_rep1 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade9_rep1.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stage10_rep1 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade10_rep1.csv",
                                    sep = ",", dec = ".", header =TRUE)

score_coloc_stage3_rep2 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade3_rep2.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stage4_rep2 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade4_rep2.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stage5_rep2 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade5_rep2.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stage6_rep2 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade6_rep2.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stage7_rep2 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade7_rep2.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stage8_rep2 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade8_rep2.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stage9_rep2 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade9_rep2.csv",
                                    sep = ",", dec = ".", header =TRUE)
score_coloc_stage10_rep2 <- read.csv("Coef_coloc_Dvir48_CTAC_CTAT_antisensboth_stade10_rep2.csv",
                                     sep = ",", dec = ".", header =TRUE)




# Data processing function, where they will be ordered in a nice
# clear table to make mean comparison graphics.

Process_coloc_data <- function(source_data_table){
  
  # Data table where everything will be stored
  
  coloc_index <- data.frame(Coef_pearson = rep(0,1),
                             
                             Coef_overlap = rep(0,1),
                             
                             k1  = rep(0,1),
                             
                             k2 = rep(0,1), 
                             
                             M1 = rep(0,1),
                             
                             M2 = rep(0,1),
                             
                             a_cytofluorogram  = rep(0,1),
                             
                             b_cytofluorogram = rep(0,1), 
                             
                             ICQ = rep(0,1),
                             
                             stringsAsFactors = FALSE)
  
  # index to add line for each occurrence in the txt file
  
  new_row <- rep(0, 9)
  
  # index for the initial loop
  
  n <- 1
  
  
  for ( i in c(1:nrow(source_data_table))){
    
    if (grepl("Pearson's Coefficient:", source_data_table[i,1],fixed = TRUE)){
      
      new_row[1] <- source_data_table[i+1,1]
      
    }
    
    if (grepl("Overlap Coefficient:", source_data_table[i,1],fixed = TRUE)){
      
      new_row[2] <- source_data_table[i+1,1]
      
    }
    
    if (grepl("r^2=k1xk2:", source_data_table[i,1],fixed = TRUE)){
      
      new_row[3] <- source_data_table[i+1,1]
      
      new_row[4] <- source_data_table[i+2,1]
      
    }
    
    if (grepl("Manders' Coefficients (using threshold", source_data_table[i,1], fixed = TRUE)){
      
      new_row[5] <- source_data_table[i+1,1]
      
      new_row[6] <- source_data_table[i+2,1]
      
    }
    
    if (grepl("Cytofluorogram's parameters:", source_data_table[i,1],fixed = TRUE)){
      
      new_row[7] <- source_data_table[i+1,1]
      
      new_row[8] <- source_data_table[i+2,1]
      
    }
    
    if (grepl("ICQ: ",source_data_table[i,1],fixed = TRUE)){
      
      new_row[9] <- source_data_table[i,1]
      
      coloc_index[n, ] <- new_row
      
      n <- n + 1
      new_row <- rep(0, 9)
      
    }
    
  }
  
  # Clean all non-numeric characters and transform to numeric
  
  coloc_index[] <- lapply(coloc_index, function(x) {
    as.numeric(sub(".*[:=]\\s*([-0-9eE\\.]+).*", "\\1", x))
  })
    
  return(coloc_index)
  
}



colocalization_index_stage3_rep1 <- Process_coloc_data(score_coloc_stage3_rep1)
colocalization_index_stage4_rep1 <- Process_coloc_data(score_coloc_stage4_rep1)
colocalization_index_stage5_rep1 <- Process_coloc_data(score_coloc_stage5_rep1)
colocalization_index_stage6_rep1 <- Process_coloc_data(score_coloc_stage6_rep1)
colocalization_index_stage7_rep1 <- Process_coloc_data(score_coloc_stage7_rep1)
colocalization_index_stage8_rep1 <- Process_coloc_data(score_coloc_stage8_rep1)
colocalization_index_stage9_rep1 <- Process_coloc_data(score_coloc_stage9_rep1)
colocalization_index_stage10_rep1 <- Process_coloc_data(score_coloc_stage10_rep1)

colocalization_index_stage3_rep2 <- Process_coloc_data(score_coloc_stage3_rep2)
colocalization_index_stage4_rep2 <- Process_coloc_data(score_coloc_stage4_rep2)
colocalization_index_stage5_rep2 <- Process_coloc_data(score_coloc_stage5_rep2)
colocalization_index_stage6_rep2 <- Process_coloc_data(score_coloc_stage6_rep2)
colocalization_index_stage7_rep2 <- Process_coloc_data(score_coloc_stage7_rep2)
colocalization_index_stage8_rep2 <- Process_coloc_data(score_coloc_stage8_rep2)
colocalization_index_stage9_rep2 <- Process_coloc_data(score_coloc_stage9_rep2)
colocalization_index_stage10_rep2 <- Process_coloc_data(score_coloc_stage10_rep2)

colocalization_index_stage3 <- rbind(colocalization_index_stage3_rep1,
                                      colocalization_index_stage3_rep2)

colocalization_index_stage4 <- rbind(colocalization_index_stage4_rep1,
                                      colocalization_index_stage4_rep2)

colocalization_index_stage5 <- rbind(colocalization_index_stage5_rep1,
                                      colocalization_index_stage5_rep2)

colocalization_index_stage6 <- rbind(colocalization_index_stage6_rep1,
                                      colocalization_index_stage6_rep2)

colocalization_index_stage7 <- rbind(colocalization_index_stage7_rep1,
                                      colocalization_index_stage7_rep2)

colocalization_index_stage8 <- rbind(colocalization_index_stage8_rep1,
                                      colocalization_index_stage8_rep2)

colocalization_index_stage9 <- rbind(colocalization_index_stage9_rep1,
                                      colocalization_index_stage9_rep2)

colocalization_index_stage10 <- rbind(colocalization_index_stage10_rep1,
                                      colocalization_index_stage10_rep2)

######################################
# long data of combined scores together

colocalization_index_stage3$stage <- rep(3,length(colocalization_index_stage3[,1]))

colocalization_index_stage4$stage <- rep(4,length(colocalization_index_stage4[,1]))

colocalization_index_stage5$stage <- rep(5,length(colocalization_index_stage5[,1]))

colocalization_index_stage6$stage <- rep(6,length(colocalization_index_stage6[,1]))

colocalization_index_stage7$stage <- rep(7,length(colocalization_index_stage7[,1]))

colocalization_index_stage8$stage <- rep(8,length(colocalization_index_stage8[,1]))

colocalization_index_stage9$stage <- rep(9,length(colocalization_index_stage9[,1]))

colocalization_index_stage10$stage <- rep(10,length(colocalization_index_stage10[,1]))


Colocalization_index_all_stages <- rbind(colocalization_index_stage3,
                                          colocalization_index_stage4,
                                          colocalization_index_stage5,
                                          colocalization_index_stage6,
                                          colocalization_index_stage7,
                                          colocalization_index_stage8,
                                          colocalization_index_stage9,
                                          colocalization_index_stage10)

Colocalization_index_all_stages_lon <- Colocalization_index_all_stages[,c(1,2,5,6,9,10)] %>%
  pivot_longer(
    cols = colnames(Colocalization_index_all_stages[,c(1,2,5,6,9)]), 
    names_to = "score_ID",             
    values_to = "values"
  )

str(Colocalization_index_all_stages_lon)

#####################################
# Test LM model
###################################

# preliminary model creation to test application conditions
model_pearson_1 <- rlm(data = Colocalization_index_all_stages, Coef_pearson ~ stage)
model_pearson_2 <- rlm(data = Colocalization_index_all_stages, Coef_pearson ~ poly(stage,6))
model_pearson_3 <- rlm(data = Colocalization_index_all_stages, Coef_pearson ~ bs(stage,3))
model_pearson_4 <- rlm(data = Colocalization_index_all_stages, Coef_pearson ~ bs(stage,4))
model_pearson_5 <- rlm(data = Colocalization_index_all_stages, Coef_pearson ~ bs(stage,5))
model_pearson_6 <- rlm(data = Colocalization_index_all_stages, Coef_pearson ~ bs(stage,6))
model_pearson_7 <- rlm(data = Colocalization_index_all_stages, Coef_pearson ~ bs(stage,7))

anova(model_pearson_1,model_pearson_2,model_pearson_3,model_pearson_4,model_pearson_5,model_pearson_6,model_pearson_7)

summary(model_pearson_6)

mean(model_pearson_6$residuals)

dwtest(model_pearson_6)

check_heteroscedasticity(model_pearson_6)

cor.test(Colocalization_index_all_stages$stage, model_pearson_6$residuals)

model_performance(model_pearson_6)

shapiro.test(model_pearson_6$residuals)

qqPlot(resid(model_pearson_6),
       main = "Pearson")

newdat <- Colocalization_index_all_stages_lon %>% 
  filter(score_ID == "Coef_pearson")
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
    x = 'Fitted values',
    y = 'Residuals',
    color = 'Weights'
  ) +
  theme_minimal()

# on exact x values
plot_2 <- ggplot(newdat, aes(x = Colocalization_index_all_stages$stage, y = res)) +
  geom_point(aes(color = weights, alpha = weights),
             size = 4) +
  scale_color_gradient(low = "cyan", high = "blue") +
  geom_smooth(color = "black", se = TRUE, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'red4', size = 1) +
  labs(
    title = 'Pearson on stage',
    x = 'Fitted values',
    y = 'Residuals',
    color = 'Weights'
  ) +
  theme_minimal()

ggarrange(plot_1,plot_2)

# pearson ok
######################

# test cutting out 0s 

for (i in c(1:length(Colocalization_index_all_stages$Coef_overlap))) {
  
  if (Colocalization_index_all_stages$Coef_overlap[i] == 0){
    
    Colocalization_index_all_stages$Coef_overlap[i] <- NA
    
  }
  
}

model_Overlap <- glmmTMB(Coef_overlap ~ bs(stage,5), 
                         family = ordbeta(link = "logit"), 
                         data = Colocalization_index_all_stages)
# test conditions
summary(model_Overlap)

mod_sim_overlap <- simulateResiduals(fittedModel = model_Overlap, n = 250)

par(mfrow=c(2,2), oma = c(0, 0, 4, 0))

plot(mod_sim_overlap ,rank = F)

plotResiduals(mod_sim_overlap)

testOutliers(mod_sim_overlap)

testDispersion(mod_sim_overlap) 

testZeroInflation(mod_sim_overlap) 

# test linearity

newdat <- Colocalization_index_all_stages_lon %>% 
  filter(score_ID == "Coef_overlap",
         values != 0)

newdat$fit = fitted(model_Overlap)
newdat$res = resid(model_Overlap)


plot_1 <- ggplot(newdat, aes(fit, res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('overlap linearity antisense')

# on exact x values

plot_2 <- ggplot(newdat, aes(x = stage, y = res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('overlap linearity antisense on stage')

ggarrange(plot_1,plot_2)

dwtest(model_Overlap)
lmtest::bptest(model_Overlap)

######################

model_ICQ_1 <- lm(data = Colocalization_index_all_stages, ICQ ~ stage)
model_ICQ_2 <- lm(data = Colocalization_index_all_stages, ICQ ~ poly(stage, 2))
model_ICQ_3 <- lm(data = Colocalization_index_all_stages, ICQ ~ bs(stage, 3))
model_ICQ_4 <- lm(data = Colocalization_index_all_stages, ICQ ~ bs(stage, 4))
model_ICQ_5 <- lm(data = Colocalization_index_all_stages, ICQ ~ bs(stage, 5))
model_ICQ_6 <- lm(data = Colocalization_index_all_stages, ICQ ~ bs(stage, 6))
model_ICQ_7 <- lm(data = Colocalization_index_all_stages, ICQ ~ bs(stage, 7))

anova(model_ICQ_1,model_ICQ_2,model_ICQ_3,model_ICQ_4,model_ICQ_5,model_ICQ_6,model_ICQ_7)

model_ICQ <- lm(data = Colocalization_index_all_stages, ICQ ~ bs(stage, 7))

# test conditions

check_model(model_ICQ)

shapiro.test(model_ICQ$residuals)

check_heteroscedasticity(model_ICQ)

shapiro.test(model_ICQ$residuals)

qqPlot(resid(model_ICQ))

mean(model_ICQ$residuals)

dwtest(model_ICQ)

cor.test(Colocalization_index_all_stages$stage, model_ICQ$residuals)


newdat <- Colocalization_index_all_stages_lon %>% 
  filter(score_ID == "ICQ")
newdat$fit = fitted(model_ICQ)
newdat$res = resid(model_ICQ )

plot_1 <- ggplot(newdat, aes(fit, res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('ICQ linearity antisense')

# on true x

plot_2 <- ggplot(newdat, aes(x = Colocalization_index_all_stages$stage, res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('ICQ linearity antisense on stage')

ggarrange(plot_1,plot_2)

######################

model_M1_1 <- glmmTMB(M1 ~ stage, 
                    family = ordbeta(link = "logit"), 
                    data = Colocalization_index_all_stages)

model_M1_2 <- glmmTMB(M1 ~ poly(stage,2), 
                    family = ordbeta(link = "logit"), 
                    data = Colocalization_index_all_stages)

model_M1_3<- glmmTMB(M1 ~ bs(stage,3), 
                    family = ordbeta(link = "logit"), 
                    data = Colocalization_index_all_stages)

model_M1_4 <- glmmTMB(M1 ~ bs(stage,4), 
                    family = ordbeta(link = "logit"), 
                    data = Colocalization_index_all_stages)

model_M1_5 <- glmmTMB(M1 ~ bs(stage,5), 
                    family = ordbeta(link = "logit"), 
                    data = Colocalization_index_all_stages)

model_M1_6 <- glmmTMB(M1 ~ bs(stage,6), 
                    family = ordbeta(link = "logit"), 
                    data = Colocalization_index_all_stages)

model_M1_7 <- glmmTMB(M1 ~ bs(stage,7), 
                    family = ordbeta(link = "logit"), 
                    data = Colocalization_index_all_stages)

anova(model_M1_1,model_M1_2,model_M1_3,model_M1_4,model_M1_5,model_M1_6,model_M1_7)

model_M1 <- glmmTMB(M1 ~ bs(stage,5), 
                    family = ordbeta(link = "logit"), 
                    data = Colocalization_index_all_stages)

summary(model_M1)

mod_sim_M1 <- simulateResiduals(fittedModel = model_M1, n = 250)

plot(mod_sim_M1 ,rank = F)

par(mfrow=c(2,2), oma = c(0, 0, 4, 0))

plotResiduals(mod_sim_M1)

testOutliers(mod_sim_M1)

testDispersion(mod_sim_M1) 

testZeroInflation(mod_sim_M1) 

# test linearity

newdat <- Colocalization_index_all_stages_lon %>% 
  filter(score_ID == "M1")
newdat$fit = fitted(model_M1)
newdat$res = resid(model_M1)

plot_1 <- ggplot(newdat, aes(fit, res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('M1 linearity antisense')

plot_2 <- ggplot(newdat, aes(x = Colocalization_index_all_stages$stage, res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('M1 linearity antisense on stage')

ggarrange(plot_1,plot_2)

######################

model_M2_1 <- glmmTMB(M2 ~ stage, 
                      family = ordbeta(link = "logit"), 
                      data = Colocalization_index_all_stages)

model_M2_2 <- glmmTMB(M2 ~ poly(stage,2), 
                      family = ordbeta(link = "logit"), 
                      data = Colocalization_index_all_stages)

model_M2_3<- glmmTMB(M2 ~ bs(stage,3), 
                     family = ordbeta(link = "logit"), 
                     data = Colocalization_index_all_stages)

model_M2_4 <- glmmTMB(M2 ~ bs(stage,4), 
                      family = ordbeta(link = "logit"), 
                      data = Colocalization_index_all_stages)

model_M2_5 <- glmmTMB(M2 ~ bs(stage,5), 
                      family = ordbeta(link = "logit"), 
                      data = Colocalization_index_all_stages)

model_M2_6 <- glmmTMB(M2 ~ bs(stage,6), 
                      family = ordbeta(link = "logit"), 
                      data = Colocalization_index_all_stages)

model_M2_7 <- glmmTMB(M2 ~ bs(stage,7), 
                      family = ordbeta(link = "logit"), 
                      data = Colocalization_index_all_stages)

anova(model_M2_1,model_M2_2,model_M2_3,model_M2_4,model_M2_5,model_M2_6,model_M2_7)


model_beta_M2 <- glmmTMB(M2 ~ bs(stage, 7), 
                      family = ordbeta(link = "logit"), 
                      data = Colocalization_index_all_stages)

summary(model_beta_M2)

mod_sim_M2 <- simulateResiduals(fittedModel = model_beta_M2, n = 250)

plot(mod_sim_M2 ,rank = F)

par(mfrow=c(2,2), oma = c(0, 0, 4, 0))

plotResiduals(mod_sim_M2)

testOutliers(mod_sim_M2)

testDispersion(mod_sim_M2) 

testZeroInflation(mod_sim_M2) 

# test linearity

newdat <- Colocalization_index_all_stages_lon %>% 
  filter(score_ID == "M2")
newdat$fit = fitted(model_beta_M2)
newdat$res = resid(model_beta_M2)

plot_1 <- ggplot(newdat, aes(fit, res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('M2 linearity antisense')

plot_2 <- ggplot(newdat, aes(x = Colocalization_index_all_stages$stage, res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('M2 linearity antisense on stage')


ggarrange(plot_1,plot_2)

# replace the 0s 

for (i in c(1:length(Colocalization_index_all_stages$Coef_overlap))) {
  
  if (is.na(Colocalization_index_all_stages$Coef_overlap[i])){
    
    Colocalization_index_all_stages$Coef_overlap[i] <- 0
    
  }
  
}

####################################
# Prediction of each model
# as well as the demonstration of effects
####################################

summary(model_pearson_6)

summary(model_Overlap)

summary(model_ICQ)

summary(model_M1)

summary(model_beta_M2)


# calculation of R2 and pseudo R2

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

pred_pearson <- predict_response(model_pearson_6, terms = "stage")
pred_overlap <- predict_response(model_Overlap, terms = "stage [3:10]")
pred_ICQ <- predict_response(model_ICQ, terms = "stage")
pred_M1 <- predict_response(model_M1, terms = "stage [all]")
pred_M2 <- predict_response(model_beta_M2, terms = "stage [all]")

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

# increase the detail of the prediction
# in order to fill the discrete categories

stage_dense <- seq(min(Colocalization_index_all_stages$stage), 
                   max(Colocalization_index_all_stages$stage), 
                   length.out = 200)

###########################
# Pearson

df_predict_pearson <- expand.grid(
  stage = stage_dense,
  Score = "Coef_pearson" 
)

pred_pearson <- predict(model_pearson_6, newdata = df_predict_pearson , 
                        se.fit = TRUE, interval = "confidence", level = 0.95)

df_predict_pearson <- df_predict_pearson %>%
  bind_cols(as_tibble(pred_pearson))


# Overlap

df_predict_Overlap <- expand.grid(
  stage = stage_dense,
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

# ICQ

df_predict_ICQ <- expand.grid(
  stage = stage_dense,
  Score = "Coef_ICQ" 
)

pred_ICQ <- predict(model_ICQ, newdata = df_predict_ICQ , 
                        se.fit = TRUE, interval = "confidence", level = 0.95)

df_predict_ICQ <- df_predict_ICQ %>%
  bind_cols(as_tibble(pred_ICQ))


# M1

df_predict_M1 <- expand.grid(
  stage = stage_dense,
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


# M2

df_predict_M2 <- expand.grid(
  stage = stage_dense,
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
# Plot everyone to observe behavior
# general

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

# put it with the data

# pearson 

graph_pearson <- ggplot(data = Colocalization_index_all_stages)+ 
  geom_boxplot(aes(x = as.factor(stage), y = Coef_pearson), color = "blue", alpha = 0.5, fill = NA) + 
  geom_point(aes(x = stage - 2, y = Coef_pearson), color = "blue", position = position_jitter(width = 0.1),
             size = 1) + 
  geom_ribbon(data = df_predict_pearson, 
              color = "blue",
              aes(x = stage -2 , ymin = fit[,"lwr"], ymax = fit[,"upr"]), 
              alpha = 0.1) +
  geom_line(data = df_predict_pearson,
            aes(stage -2, fit[,"lwr"]), color = "grey30", size = 0.1) + 
  geom_line(data = df_predict_pearson,
            aes(stage -2, fit[,"upr"]), color = "grey30", size = 0.1) + 
  geom_line(data = df_predict_pearson, 
            aes(x = stage -2 , y = fit[,"fit"]),
            linewidth = 1.2)   +
  
  scale_x_discrete("Stage",labels = c("Stage3","Stage4","Stage5","Stage6",
                                       "Stage7","Stage8","Stage9","Stage10"))+
  ylim(0, 1)+
  
  labs(
    x = "Egg chamber stages",
    y = "Pearson coefficient (r)",
    title = "Evolution of the Pearson score (r) of sense LncRNA AAACTAT and AAACTAC",
    subtitle = "The probes used were antisense"
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

# Overlap

graph_overlap <- ggplot(data = Colocalization_index_all_stages)+ 
  geom_boxplot(aes(x = as.factor(stage), y = Coef_overlap), color = "orange", alpha = 0.5,
               fill = NA) + 
  geom_point(aes(x = stage - 2, y = Coef_overlap), color = "orange", 
             position = position_jitter(width = 0.1), size  = 1) + 
  geom_ribbon(data = df_predict_Overlap, 
              color = "orange",
              aes(x = stage -2 , ymin = lower, ymax = upper), 
              alpha = 0.1) +
  geom_line(data = df_predict_Overlap,
            aes(stage -2, lower), color = "grey30", size = 0.1) + 
  geom_line(data = df_predict_Overlap,
            aes(stage -2, upper), color = "grey30", size = 0.1) +
  geom_line(data = df_predict_Overlap, 
            aes(x = stage -2 , y = fit),
            linewidth = 1.2)   +
  
  scale_x_discrete("Stage",labels = c("Stage3","Stage4","Stage5","Stage6",
                                       "Stage7","Stage8","Stage9","Stage10"))+
  ylim(0,1)+
  
  labs(
    x = "Developmental stage",
    y = "MOC coefficient",
    title = "Evolution of the overlap score (MOC) of sense LncRNA AAACTAT and AAACTAC",
    subtitle = "The probes used were antisense"
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

# ICQ

graph_ICQ <- ggplot(data = Colocalization_index_all_stages)+ 
  geom_boxplot(aes(x = as.factor(stage), y = ICQ), color = "red", alpha = 0.5, fill = NA) + 
  geom_point(aes(x = stage - 2, y = ICQ),color = "red",
             position = position_jitter(width = 0.1), size = 1) + 
  geom_ribbon(data = df_predict_ICQ,
              color = "red",
              aes(x = stage -2 , ymin = fit[,"lwr"], ymax = fit[,"upr"]), 
              alpha = 0.1) +
  geom_line(data = df_predict_ICQ,
            aes(stage -2, fit[,"lwr"]), color = "grey30", size = 0.1) + 
  geom_line(data = df_predict_ICQ,
            aes(stage -2, fit[,"upr"]), color = "grey30", size = 0.1) +
  geom_line(data = df_predict_ICQ, 
            aes(x = stage -2 , y = fit[,"fit"]),
            linewidth = 1.2)    +
  
  scale_x_discrete("Stage",labels = c("Stage3","Stage4","Stage5","Stage6",
                                       "Stage7","Stage8","Stage9","Stage10"))+
  ylim(-0.5, 0.5)+
  
  labs(
    x = "Egg chamber stages",
    y = "ICQ coefficient",
    title = "Evolution of the ICQ score of sense LncRNA AAACTAT and AAACTAC",
    subtitle = "The probes used were antisense"
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
# Table containing the means of each value

Calculate_mean <- function(source_data_table){
  
  # Data table where everything will be stored
  
  mean_coloc_index <- data.frame(Coef_pearson = rep(0,1),
                             
                             Coef_overlap = rep(0,1),
                             
                             k1  = rep(0,1),
                             
                             k2 = rep(0,1), 
                             
                             M1 = rep(0,1),
                             
                             M2 = rep(0,1),
                             
                             a_cytofluorogram  = rep(0,1),
                             
                             b_cytofluorogram = rep(0,1), 
                             
                             ICQ = rep(0,1))
  
  for (i in c(1:ncol(source_data_table))) {
    
    mean_coloc_index[,i] <- mean(source_data_table[,i])
    
    
  }
  
  return(mean_coloc_index)
}  

mean_coloc_coef_stage3 <- Calculate_mean(colocalization_index_stage3)
mean_coloc_coef_stage4 <- Calculate_mean(colocalization_index_stage4)
mean_coloc_coef_stage5 <- Calculate_mean(colocalization_index_stage5)
mean_coloc_coef_stage6 <- Calculate_mean(colocalization_index_stage6)
mean_coloc_coef_stage7 <- Calculate_mean(colocalization_index_stage7)
mean_coloc_coef_stage8 <- Calculate_mean(colocalization_index_stage8)
mean_coloc_coef_stage9 <- Calculate_mean(colocalization_index_stage9)
mean_coloc_coef_stage10 <- Calculate_mean(colocalization_index_stage10)


mean_coloc_coef <- rbind(mean_coloc_coef_stage3,mean_coloc_coef_stage4,mean_coloc_coef_stage5,
                            mean_coloc_coef_stage6, mean_coloc_coef_stage7, mean_coloc_coef_stage8,
                            mean_coloc_coef_stage9, mean_coloc_coef_stage10)

rownames(mean_coloc_coef) <- c(
  "Stage3",
  "Stage4",
  "Stage5",
  "Stage6",
  "Stage7",
  "Stage8",
  "Stage9",
  "Stage10"
)

Stage <- c(
  "Stage3",
  "Stage4",
  "Stage5",
  "Stage6",
  "Stage7",
  "Stage8",
  "Stage9",
  "Stage10")

mean_coloc_coef <- cbind(mean_coloc_coef,Stage)


##############################
# Manders coef ggplot Graphic
##############################

# remake a dataframe with the 2 values of interest

# M1 = CTAT, M2 = CTAC

data_M1_M2 <- Colocalization_index_all_stages[,c(5,6,10)] %>%
  pivot_longer(
    cols = colnames(Colocalization_index_all_stages[,c(5,6)]), 
    names_to = "score_ID",             
    values_to = "values"
  )

transparency = 1

windowsFonts( A = windowsFont("baskerville old face"))

graph_manders <- ggplot(data_M1_M2, aes(x = factor(stage), y = values, color = score_ID, fill - NA)) +

  geom_boxplot(alpha = 0.7, 
               outlier.shape = NA, 
               position = position_dodge(width = 0.8)) +
  
  geom_point(aes(),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.1),
             alpha = 1, size = 1) +
  
  scale_fill_manual(values = c("M1" = "magenta", "M2" = "green")) +
  scale_color_manual(values = c("M1" = "magenta", "M2" = "green")) +

  geom_ribbon(data = df_predict_M2, 
              aes(x = stage - 2, ymin = lower, ymax = upper, y = NULL, fill = NULL), 
              alpha = 0.15, fill = "green", inherit.aes = FALSE) +
  geom_line(data = df_predict_M2,
            aes(stage -2, lower), 
            color = "darkgreen", size = 0.1, inherit.aes = FALSE) + 
  
  geom_line(data = df_predict_M2,
            aes(stage -2, upper), 
            color = "darkgreen", size = 0.1, inherit.aes = FALSE) +
  
  geom_line(data = df_predict_M2, 
            aes(x = stage - 2, y = fit), 
            color = "green", linewidth = 1.2, inherit.aes = FALSE) +

  geom_ribbon(data = df_predict_M1, 
              aes(x = stage - 2, ymin = lower, ymax = upper, y = NULL, fill = NULL), 
              alpha = 0.15, fill = "magenta", inherit.aes = FALSE) +
  
  geom_line(data = df_predict_M1,
            aes(stage -2, lower), color = "darkmagenta", size = 0.1, inherit.aes = FALSE) + 
  
  geom_line(data = df_predict_M1,
            aes(stage -2, upper), color = "darkmagenta", size = 0.1, inherit.aes = FALSE) +
  
  geom_line(data = df_predict_M1, 
            aes(x = stage - 2, y = fit), 
            color = "magenta", linewidth = 1.2, inherit.aes = FALSE)  +
  
  scale_x_discrete("Stages",labels = c("Stage3","Stage4","Stage5","Stage6",
                                       "Stage7","Stage8","Stage9","Stage10"))+
  ylim(0, 1)+
  
  labs(
    x = "Egg chamber stages",
    y = "M1 & M2 coefficient",
    title = "Evolution of the M1 & M2 score of sense LncRNA AAACTAT and AAACTAC",
    subtitle = "The probes used were antisense"
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
# Print boxplot graphics
###################################################

print(graph_pearson)
print(graph_overlap)
print(graph_ICQ)
print(graph_manders)

ggarrange(graph_pearson,graph_overlap,graph_ICQ,graph_manders)

#################################################
# variance test between Y

par(mfrow=c(4,2), oma = c(0, 0, 2, 0))

concatenated_data <- rbind(colocalization_index_stage10, colocalization_index_stage3,
                         colocalization_index_stage4, colocalization_index_stage5,
                         colocalization_index_stage6, colocalization_index_stage7, 
                         colocalization_index_stage8, colocalization_index_stage9)
concatenated_data <- concatenated_data[,c(-3,-4,-7,-8)]

pairs.panels(concatenated_data[,c(1:5)])

par(mfrow=c(3,2), oma = c(0, 0, 3, 0))

hist(concatenated_data$Coef_pearson,
     main = "Pearson antisense")

hist(concatenated_data$Coef_overlap,
     main = "Overlap antisense")

hist(concatenated_data$M1,
     main = "M1 antisense")

hist(concatenated_data$M2,
     main = "M2 antisense")

hist(concatenated_data$ICQ,
     main = "ICQ antisense")

mtext("Antisense", outer = TRUE, cex = 1.5, font = 2)


library(reshape2)

boxplot(concatenated_data,
        main = "antisense probe")

for ( i in c(1:5)){
  
  name = colnames(concatenated_data)[i]
  
  number_0 <- (length(which(concatenated_data[,i] == 0)))
  
  print(paste("Number of 0 for",name , number_0))
  
}

# boxplot for each score

par(mfrow=c(3,2), oma = c(0, 0, 3, 0))

boxplot(data = Colocalization_index_all_stages, Coef_pearson ~ stage)
boxplot(data = Colocalization_index_all_stages, Coef_overlap ~ stage)
boxplot(data = Colocalization_index_all_stages, ICQ ~ stage)
boxplot(data = Colocalization_index_all_stages, M1 ~ stage)
boxplot(data = Colocalization_index_all_stages, M2 ~ stage)

mtext("Distribution of variance according to stage \n for antisense probes"
      , outer = TRUE, cex = 1.2, font = 1.5)

##############################################
# spline ANOVA verification
#############################################

# pearson
model_pearson_1 <- rlm(data = Colocalization_index_all_stages, Coef_pearson ~ stage)
model_pearson_2 <- rlm(data = Colocalization_index_all_stages, Coef_pearson ~ poly(stage,6))
model_pearson_3 <- rlm(data = Colocalization_index_all_stages, Coef_pearson ~ bs(stage,3))
model_pearson_4 <- rlm(data = Colocalization_index_all_stages, Coef_pearson ~ bs(stage,4))
model_pearson_5 <- rlm(data = Colocalization_index_all_stages, Coef_pearson ~ bs(stage,5))
model_pearson_6 <- rlm(data = Colocalization_index_all_stages, Coef_pearson ~ bs(stage,6))
model_pearson_7 <- rlm(data = Colocalization_index_all_stages, Coef_pearson ~ bs(stage,7))

anova(model_pearson_1,model_pearson_2,model_pearson_3,model_pearson_4,model_pearson_5,model_pearson_6,model_pearson_7)


# overlap
model_Overlap_1 <- glmmTMB(Coef_overlap ~ stage, 
                         family = ordbeta(link = "logit"), 
                         data = Colocalization_index_all_stages)
model_Overlap_2 <- glmmTMB(Coef_overlap ~ poly(stage,2), 
                         family = ordbeta(link = "logit"), 
                         data = Colocalization_index_all_stages)
model_Overlap_3 <- glmmTMB(Coef_overlap ~ bs(stage,3), 
                         family = ordbeta(link = "logit"), 
                         data = Colocalization_index_all_stages)
model_Overlap_4 <- glmmTMB(Coef_overlap ~ bs(stage,4), 
                         family = ordbeta(link = "logit"), 
                         data = Colocalization_index_all_stages)
model_Overlap_5 <- glmmTMB(Coef_overlap ~ bs(stage,5), 
                         family = ordbeta(link = "logit"), 
                        data = Colocalization_index_all_stages)
model_Overlap_6 <- glmmTMB(Coef_overlap ~ bs(stage,6), 
                           family = ordbeta(link = "logit"), 
                           data = Colocalization_index_all_stages)
model_Overlap_7 <- glmmTMB(Coef_overlap ~ bs(stage,7), 
                           family = ordbeta(link = "logit"), 
                           data = Colocalization_index_all_stages)

anova(model_Overlap_1,model_Overlap_2,model_Overlap_3,model_Overlap_4,model_Overlap_5,model_Overlap_6,model_Overlap_7)



# ICQ
model_ICQ_1 <- lm(data = Colocalization_index_all_stages, ICQ ~ stage)
model_ICQ_2 <- lm(data = Colocalization_index_all_stages, ICQ ~ poly(stage, 2))
model_ICQ_3 <- lm(data = Colocalization_index_all_stages, ICQ ~ bs(stage, 3))
model_ICQ_4 <- lm(data = Colocalization_index_all_stages, ICQ ~ bs(stage, 4))
model_ICQ_5 <- lm(data = Colocalization_index_all_stages, ICQ ~ bs(stage, 5))
model_ICQ_6 <- lm(data = Colocalization_index_all_stages, ICQ ~ bs(stage, 6))
model_ICQ_7 <- lm(data = Colocalization_index_all_stages, ICQ ~ bs(stage, 7))

anova(model_ICQ_1,model_ICQ_2,model_ICQ_3,model_ICQ_4,model_ICQ_5,model_ICQ_6,model_ICQ_7)


# M1

model_M1_1 <- glmmTMB(M1 ~ stage, 
                      family = ordbeta(link = "logit"), 
                      data = Colocalization_index_all_stages)

model_M1_2 <- glmmTMB(M1 ~ poly(stage,2), 
                      family = ordbeta(link = "logit"), 
                      data = Colocalization_index_all_stages)

model_M1_3<- glmmTMB(M1 ~ bs(stage,3), 
                     family = ordbeta(link = "logit"), 
                     data = Colocalization_index_all_stages)

model_M1_4 <- glmmTMB(M1 ~ bs(stage,4), 
                      family = ordbeta(link = "logit"), 
                      data = Colocalization_index_all_stages)

model_M1_5 <- glmmTMB(M1 ~ bs(stage,5), 
                      family = ordbeta(link = "logit"), 
                      data = Colocalization_index_all_stages)

model_M1_6 <- glmmTMB(M1 ~ bs(stage,6), 
                      family = ordbeta(link = "logit"), 
                      data = Colocalization_index_all_stages)

model_M1_7 <- glmmTMB(M1 ~ bs(stage,7), 
                      family = ordbeta(link = "logit"), 
                      data = Colocalization_index_all_stages)

anova(model_M1_1,model_M1_2,model_M1_3,model_M1_4,model_M1_5,model_M1_6,model_M1_7)


# M2

model_M2_1 <- glmmTMB(M2 ~ stage, 
                      family = ordbeta(link = "logit"), 
                      data = Colocalization_index_all_stages)

model_M2_2 <- glmmTMB(M2 ~ poly(stage,2), 
                      family = ordbeta(link = "logit"), 
                      data = Colocalization_index_all_stages)

model_M2_3<- glmmTMB(M2 ~ bs(stage,3), 
                     family = ordbeta(link = "logit"), 
                     data = Colocalization_index_all_stages)

model_M2_4 <- glmmTMB(M2 ~ bs(stage,4), 
                      family = ordbeta(link = "logit"), 
                      data = Colocalization_index_all_stages)

model_M2_5 <- glmmTMB(M2 ~ bs(stage,5), 
                      family = ordbeta(link = "logit"), 
                      data = Colocalization_index_all_stages)

model_M2_6 <- glmmTMB(M2 ~ bs(stage,6), 
                      family = ordbeta(link = "logit"), 
                      data = Colocalization_index_all_stages)

model_M2_7 <- glmmTMB(M2 ~ bs(stage,7), 
                      family = ordbeta(link = "logit"), 
                      data = Colocalization_index_all_stages)

anova(model_M2_1,model_M2_2,model_M2_3,model_M2_4,model_M2_5,model_M2_6,model_M2_7)

