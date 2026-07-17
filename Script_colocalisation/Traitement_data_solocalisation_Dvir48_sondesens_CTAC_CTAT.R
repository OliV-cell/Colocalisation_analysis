###########################################
# Colocalization results analysis script
# CTAC+CTAT sense probe (antisense strand)
##########################################

# various statistical tests


library(car)
library(dplyr)
library(coin)
library(rcompanion)
library(psych)
library(nlstools)
library(tidyr)

# regression
library(DHARMa)
library(splines)
library(gamm4)
library(mgcv)
library(gvlma)
library(performance)
library(see)
library(MASS)
library(robustbase)

# graphical representation
library(ggplot2)
library(ggpubr)
library(ggeffects)

# citation 
library(grateful)

#########################################################################
# opening of files containing the indices, is in txt format inside csv

# Vector of future file names to be used

file_names <- vector()
file_names_rep1 <- vector()

for ( i in (3:10)){
  
  file_names[[i-2]] <- paste0("score_coloc_stage",i)
  file_names_rep1[[i-2]] <- paste0("score_coloc_stage",i,"_rep1")
  
}

# opening of repetition 0

for ( i in c(3:10)) {
  
  ID <- read.csv(paste0("Coef_coloc_Dvir48_CTAC_CTAT_sensboth_stade",i,".csv"),
                              sep = ",", dec = ".", header =TRUE)
  assign(file_names[i-2],ID)
  
}

# opening of repetition 1

for ( i in c(3:10)) {
  
  ID <- read.csv(paste0("Coef_coloc_Dvir48_CTAC_CTAT_sensboth_stade",i,"_rep1.csv"),
                                            sep = ",", dec = ".", header =TRUE)
  assign(file_names_rep1[i-2],ID)
  
}

# following the same logic, we can then make a loop j inside the loop i for i stage
# and J repetition

###############################################################################
# Data cleaning to recover only the data important to the analysis

# Cleaning function

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

###################################################################
# Data cleaning
colocalization_index_stage3 <- Process_coloc_data(score_coloc_stage3)
colocalization_index_stage4 <- Process_coloc_data(score_coloc_stage4)
colocalization_index_stage5 <- Process_coloc_data(score_coloc_stage5)
colocalization_index_stage6 <- Process_coloc_data(score_coloc_stage6)
colocalization_index_stage7 <- Process_coloc_data(score_coloc_stage7)
colocalization_index_stage8 <- Process_coloc_data(score_coloc_stage8)
colocalization_index_stage9 <- Process_coloc_data(score_coloc_stage9)
colocalization_index_stage10 <- Process_coloc_data(score_coloc_stage10)

colocalization_index_stage3_rep1 <- Process_coloc_data(score_coloc_stage3_rep1)
colocalization_index_stage4_rep1 <- Process_coloc_data(score_coloc_stage4_rep1)
colocalization_index_stage5_rep1 <- Process_coloc_data(score_coloc_stage5_rep1)
colocalization_index_stage6_rep1 <- Process_coloc_data(score_coloc_stage6_rep1)
colocalization_index_stage7_rep1 <- Process_coloc_data(score_coloc_stage7_rep1)
colocalization_index_stage8_rep1 <- Process_coloc_data(score_coloc_stage8_rep1)
colocalization_index_stage9_rep1 <- Process_coloc_data(score_coloc_stage9_rep1)
colocalization_index_stage10_rep1 <- Process_coloc_data(score_coloc_stage10_rep1)

####################################################################################

# Merging rows of all trials

colocalization_index_stage3_total <- rbind(colocalization_index_stage3_rep1,
                                      colocalization_index_stage3)

colocalization_index_stage4_total <- rbind(colocalization_index_stage4_rep1,
                                      colocalization_index_stage4)

colocalization_index_stage5_total <- rbind(colocalization_index_stage5_rep1,
                                      colocalization_index_stage5)

colocalization_index_stage6_total <- rbind(colocalization_index_stage6_rep1,
                                      colocalization_index_stage6)

colocalization_index_stage7_total <- rbind(colocalization_index_stage7_rep1,
                                      colocalization_index_stage7)

colocalization_index_stage8_total <- rbind(colocalization_index_stage8_rep1,
                                      colocalization_index_stage8)

colocalization_index_stage9_total <- rbind(colocalization_index_stage9_rep1,
                                      colocalization_index_stage9)

colocalization_index_stage10_total <- rbind(colocalization_index_stage10_rep1,
                                       colocalization_index_stage10)

##################################################################################

# Calculate means of all scores

# Function to calculate means for each score of each stage

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

mean_coloc_coef_stage3 <- Calculate_mean(colocalization_index_stage3_total)
mean_coloc_coef_stage4 <- Calculate_mean(colocalization_index_stage4_total)
mean_coloc_coef_stage5 <- Calculate_mean(colocalization_index_stage5_total)
mean_coloc_coef_stage6 <- Calculate_mean(colocalization_index_stage6_total)
mean_coloc_coef_stage7 <- Calculate_mean(colocalization_index_stage7_total)
mean_coloc_coef_stage8 <- Calculate_mean(colocalization_index_stage8_total)
mean_coloc_coef_stage9 <- Calculate_mean(colocalization_index_stage9_total)
mean_coloc_coef_stage10 <- Calculate_mean(colocalization_index_stage10_total)


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


######################################
# long data of combined scores together

colocalization_index_stage3_total$stage <- rep(3,length(colocalization_index_stage3_total[,1]))

colocalization_index_stage4_total$stage <- rep(4,length(colocalization_index_stage4_total[,1]))

colocalization_index_stage5_total$stage <- rep(5,length(colocalization_index_stage5_total[,1]))

colocalization_index_stage6_total$stage <- rep(6,length(colocalization_index_stage6_total[,1]))

colocalization_index_stage7_total$stage <- rep(7,length(colocalization_index_stage7_total[,1]))

colocalization_index_stage8_total$stage <- rep(8,length(colocalization_index_stage8_total[,1]))

colocalization_index_stage9_total$stage <- rep(9,length(colocalization_index_stage9_total[,1]))

colocalization_index_stage10_total$stage <- rep(10,length(colocalization_index_stage10_total[,1]))

Colocalization_index_all_stages <- rbind(colocalization_index_stage3_total,
                                          colocalization_index_stage4_total,
                                          colocalization_index_stage5_total,
                                          colocalization_index_stage6_total,
                                          colocalization_index_stage7_total,
                                          colocalization_index_stage8_total,
                                          colocalization_index_stage9_total,
                                          colocalization_index_stage10_total)

Colocalization_index_all_stages_lon <- Colocalization_index_all_stages[,c(1,2,5,6,9,10)] %>%
  pivot_longer(
    cols = colnames(Colocalization_index_all_stages[,c(1,2,5,6,9)]), 
    names_to = "score_ID",             
    values_to = "values"
  )

####################################################################################
# Graphical representation of the evolution of scores according to stages
# As well as their inter-stage variances
####################################################################################


#####################################
# Test LM model
###################################

# preliminary model creation to test application conditions

model_pearson <- lm(data = Colocalization_index_all_stages, Coef_pearson ~ bs(stage,7))

# test conditions
check_model(model_pearson)
check_autocorrelation(model_pearson)
check_heteroscedasticity(model_pearson)

shapiro.test(model_pearson$residuals)
mean(model_pearson$residuals)
cor.test(Colocalization_index_all_stages$stage, model_pearson$residuals)
model_performance(model_pearson)

# linearity verification 

newdat <- Colocalization_index_all_stages_lon %>% 
  filter(score_ID == "Coef_pearson")
newdat$fit <- fitted(model_pearson)
newdat$res <- resid(model_pearson)
newdat$weights <- model_pearson$w

plot_1 <- ggplot(newdat, aes(fit, res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('pearson linearity sense')

# on exact x values

plot_2 <- ggplot(newdat, aes(x = Colocalization_index_all_stages$stage, y = res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('pearson linearity sense on stage')

ggarrange(plot_1,plot_2)

######################

for (i in c(1:length(Colocalization_index_all_stages$Coef_overlap))) {
  
  if (Colocalization_index_all_stages$Coef_overlap[i] == 0){
    
    Colocalization_index_all_stages$Coef_overlap[i] <- NA
    
  }
  
}

model_Overlap <- rlm(data = Colocalization_index_all_stages, Coef_overlap ~ bs(stage,5))

shapiro.test(model_Overlap$residuals)

check_heteroscedasticity(model_Overlap)

qqPlot(resid(model_Overlap),
       main = "Overlap")
# test conditions

newdat <- Colocalization_index_all_stages_lon %>% 
  filter(score_ID == "Coef_overlap",
         values != 0)
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
    x = 'Fitted values',
    y = 'Residuals',
    color = 'Weights'
  ) +
  theme_minimal()


plot_2 <- ggplot(newdat, aes(x = stage, y = res)) +
  geom_point(aes(color = weights, alpha = weights),
             size = 4) +
  scale_color_gradient(low = "cyan", high = "blue") +
  geom_smooth(color = "black", se = TRUE, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'red4', size = 1) +
  labs(
    title = 'Overlap on stage',
    x = 'Fitted values',
    y = 'Residuals',
    color = 'Weights'
  ) +
  theme_minimal()

ggarrange(plot_1,plot_2)

mean(model_Overlap$residuals)
dwtest(model_Overlap)
lmtest::bptest(model_Overlap)
cor.test(newdat$stage, model_Overlap$residuals)


# replace the 0s 

for (i in c(1:length(Colocalization_index_all_stages$Coef_overlap))) {
  
  if (is.na(Colocalization_index_all_stages$Coef_overlap[i])){
    
    Colocalization_index_all_stages$Coef_overlap[i] <- 0
    
  }
  
}


######################
model_ICQ <- rlm(data = Colocalization_index_all_stages, ICQ ~ bs(stage,4))

shapiro.test(model_ICQ$residuals)

check_heteroscedasticity(model_ICQ)

mean(model_ICQ$residuals)

dwtest(model_ICQ)

cor.test(Colocalization_index_all_stages$stage, model_ICQ$residuals)

qqPlot(resid(model_ICQ),
       main = "ICQ")

newdat <- Colocalization_index_all_stages_lon %>% 
  filter(score_ID == "ICQ")
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
    x = 'Fitted values',
    y = 'Residuals',
    color = 'Weights'
  ) +
  theme_minimal()

plot_2 <- ggplot(newdat, aes(x = Colocalization_index_all_stages$stage, y = res)) +
  geom_point(aes(color = weights, alpha = weights),
             size = 4) +
  scale_color_gradient(low = "cyan", high = "blue") +
  geom_smooth(color = "black", se = TRUE, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'red4', size = 1) +
  labs(
    title = 'ICQ on stage',
    x = 'Fitted values',
    y = 'Residuals',
    color = 'Weights'
  ) +
  theme_minimal()

ggarrange(plot_1,plot_2)

######################
model_M1 <- rlm(data = Colocalization_index_all_stages, M1 ~  bs(stage,5))

# test conditions
shapiro.test(model_M1$residuals)
mean(model_M1$residuals)
check_heteroscedasticity(model_M1)
cor.test(Colocalization_index_all_stages$stage, model_M1$residuals)

qqPlot(resid(model_M1),
       main = "M1")

newdat <- Colocalization_index_all_stages_lon %>% 
  filter(score_ID == "M1")
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
    x = 'Fitted values',
    y = 'Residuals',
    color = 'Weights'
  ) +
  theme_minimal()

plot_2 <- ggplot(newdat, aes(x = Colocalization_index_all_stages$stage, y = res)) +
  geom_point(aes(color = weights, alpha = weights),
             size = 4) +
  scale_color_gradient(low = "cyan", high = "blue") +
  geom_smooth(color = "black", se = TRUE, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'red4', size = 1) +
  labs(
    title = 'M1 on stage',
    x = 'Fitted values',
    y = 'Residuals',
    color = 'Weights'
  ) +
  theme_minimal()

ggarrange(plot_1,plot_2)

######################
model_M2 <- rlm(data = Colocalization_index_all_stages,M2 ~ bs(stage,6))

# test conditions
mean(model_M2$residuals)
shapiro.test(model_M2$residuals)
check_heteroscedasticity(model_M2)
dwtest(model_M2)
cor.test(Colocalization_index_all_stages$stage, model_M2$residuals)

qqPlot(resid(model_M2),
       main = "M2")

newdat <- Colocalization_index_all_stages_lon %>% 
  filter(score_ID == "M2")
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
    x = 'Fitted values',
    y = 'Residuals',
    color = 'Weights'
  ) +
  theme_minimal()

plot_2 <- ggplot(newdat, aes(x = Colocalization_index_all_stages$stage, y = res)) +
  geom_point(aes(color = weights, alpha = weights),
             size = 4) +
  scale_color_gradient(low = "cyan", high = "blue") +
  geom_smooth(color = "black", se = TRUE, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'red4', size = 1) +
  labs(
    title = 'M2 on stage',
    x = 'Fitted values',
    y = 'Residuals',
    color = 'Weights'
  ) +
  theme_minimal()

ggarrange(plot_1,plot_2)

# test conditions
mean(model_M2$residuals)
shapiro.test(model_M2$residuals)
dwtest(model_M2)
lmtest::bptest(model_M2)
cor.test(Colocalization_index_all_stages$stage, model_M2$residuals)


# mean of everyone
mean(model_pearson$residuals)
mean(model_Overlap$residuals)
mean(model_ICQ$residuals)
mean(model_M1$residuals)
mean(model_M2$residuals)


####################################
# Prediction of each model
####################################

summary(model_pearson)
summary(model_Overlap)
summary(model_ICQ)
summary(model_M1)
summary(model_M2)

# calculation of R2 and pseudo R2

data_r2 <- data.frame("Pearson" = c(0,0),
                      "Overlap" = c(0,0),
                      "ICQ" = c(0,0),
                      "M1" = c(0,0),
                      "M2" = c(0,0),
                      row.names = c("R^2","P.value"))

data_r2[1,1] <- summary(model_pearson)$r.squared
data_r2[1,2] <- (cor(model_Overlap$model$Coef_overlap, predict(model_Overlap)))^2
data_r2[1,3] <- (cor(model_ICQ$model$ICQ, predict(model_ICQ)))^2
data_r2[1,4] <- (cor(model_M1$model$M1, predict(model_M1)))^2
data_r2[1,5] <- (cor(model_M2$model$M2, predict(model_M2)))^2

pred_pearson <- predict_response(model_pearson, terms = "stage [3:10]")
pred_overlap <- predict_response(model_Overlap, terms = "stage [3:10]")
pred_ICQ <- predict_response(model_ICQ, terms = "stage [3:10]")
pred_M1 <- predict_response(model_M1, terms = "stage [3:10]")
pred_M2 <- predict_response(model_M2, terms = "stage [3:10]")

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

annotate_figure(plot_fusion, top = text_grob("Prediction score Sense", 
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

pred_pearson <- predict(model_pearson, newdata = df_predict_pearson , 
                        se.fit = TRUE, interval = "confidence", level = 0.95)

df_predict_pearson <- df_predict_pearson %>%
  bind_cols(as_tibble(pred_pearson))


# Overlap

df_predict_Overlap <- expand.grid(
  stage = stage_dense,
  Score = "Coef_Overlap" 
)

pred_Overlap <- predict(model_Overlap, newdata = df_predict_pearson , 
                        se.fit = TRUE, interval = "confidence", level = 0.95)

df_predict_Overlap <- df_predict_Overlap %>%
  bind_cols(as_tibble(pred_Overlap))


# ICQ

df_predict_ICQ <- expand.grid(
  stage = stage_dense,
  Score = "Coef_ICQ" 
)

pred_ICQ <- predict(model_ICQ, newdata = df_predict_pearson , 
                    se.fit = TRUE, interval = "confidence", level = 0.95)

df_predict_ICQ <- df_predict_ICQ %>%
  bind_cols(as_tibble(pred_ICQ))

# M1

df_predict_M1 <- expand.grid(
  stage = stage_dense,
  Score = "Coef_M1" 
)

pred_M1 <- predict(model_M1, newdata = df_predict_pearson , 
                   se.fit = TRUE, interval = "confidence", level = 0.95)

df_predict_M1 <- df_predict_M1 %>%
  bind_cols(as_tibble(pred_M1))


# M2

df_predict_M2 <- expand.grid(
  stage = stage_dense,
  Score = "Coef_M2" 
)

pred_M2 <- predict(model_M2, newdata = df_predict_pearson , 
                   se.fit = TRUE, interval = "confidence", level = 0.95)


df_predict_M2 <- df_predict_M2 %>%
  bind_cols(as_tibble(pred_M2))

#####################
# Plot everyone to observe behavior
# general

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
# pearson coef ggplot Graphic
##############################
# put it with the data

# pearson 

windowsFonts( A = windowsFont("baskerville old face"))

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
    y = "Pearson coefficient",
    title = "Evolution of the Pearson score (r) of antisense LncRNA AAACTAT and AAACTAC",
    subtitle = "The probes used were sense"
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

# Overlap

graph_overlap <- ggplot(data = Colocalization_index_all_stages)+ 
  geom_boxplot(aes(x = as.factor(stage), y = Coef_overlap), color = "orange", alpha = 0.5,
               fill = NA) + 
  geom_point(aes(x = stage - 2, y = Coef_overlap), color = "orange", 
             position = position_jitter(width = 0.1), size  = 1) + 
  geom_ribbon(data = df_predict_Overlap, 
              color = "orange",
              aes(x = stage -2 , ymin = fit[,"lwr"], ymax = fit[,"upr"]), 
              alpha = 0.1) +
  geom_line(data = df_predict_Overlap,
            aes(stage -2, fit[,"lwr"]), color = "grey30", size = 0.1) + 
  geom_line(data = df_predict_Overlap,
            aes(stage -2, fit[,"upr"]), color = "grey30", size = 0.1) +
  geom_line(data = df_predict_Overlap, 
            aes(x = stage -2 , y = fit),
            linewidth = 1.2)   +
  
  scale_x_discrete("Stage",labels = c("Stage3","Stage4","Stage5","Stage6",
                                       "Stage7","Stage8","Stage9","Stage10"))+
  ylim(0, 1)+
  
  labs(
    x = "Developmental stage",
    y = "MOC coefficient",
    title = "Evolution of the overlap score (MOC) of antisense LncRNA AAACTAT and AAACTAC",
    subtitle = "The probes used were sense"
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
  
  annotate("text", x = Inf, y = Inf, label = paste0("R²= ", round(data_r2[1,3], digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") + 
  
  scale_x_discrete("Stages",labels = c("Stage3","Stage4","Stage5","Stage6",
                                       "Stage7","Stage8","Stage9","Stage10"))+
  ylim(-0.5, 0.5)+
  
  labs(
    x = "Egg chamber stages",
    y = "ICQ coefficient",
    title = "Evolution of the ICQ score of antisense LncRNA AAACTAT and AAACTAC",
    subtitle = "The probes used were sense"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "right",
    text = element_text(size = 12, family = "A"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "A"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain", family = "A"),
    plot.margin = margin(20, 20, 20, 20)
  )


# manders

# M1 = CTAC, M2 = CTAT

data_M1_M2 <- Colocalization_index_all_stages[,c(5,6,10)] %>%
  pivot_longer(
    cols = colnames(Colocalization_index_all_stages[,c(5,6)]), 
    names_to = "score_ID",             
    values_to = "values"
  )

transparency = 1

windowsFonts( A = windowsFont("baskerville old face"))

graph_manders <- ggplot(data_M1_M2, aes(x = factor(stage), y = values,  color  = score_ID)) +
  
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
              aes(x = stage - 2, ymin = fit[,"lwr"], ymax = fit[,"upr"], y = NULL, fill = NULL), 
              alpha = 0.15, fill = "magenta", inherit.aes = FALSE) +
  geom_line(data = df_predict_M2,
            aes(stage -2, fit[,"lwr"]), 
            color = "darkmagenta", size = 0.1, inherit.aes = FALSE) + 
  
  geom_line(data = df_predict_M2,
            aes(stage -2, fit[,"upr"]), 
            color = "darkmagenta", size = 0.1, inherit.aes = FALSE) +
  
  geom_line(data = df_predict_M2, 
            aes(x = stage - 2, y = fit), 
            color = "magenta", linewidth = 1.2, inherit.aes = FALSE) +
  
  geom_ribbon(data = df_predict_M1, 
              aes(x = stage - 2, ymin = fit[,"lwr"], ymax = fit[,"upr"], y = NULL, fill = NULL), 
              alpha = 0.15, fill = "green", inherit.aes = FALSE) +
  
  geom_line(data = df_predict_M1,
            aes(stage -2, fit[,"lwr"]), color = "darkgreen", size = 0.1, inherit.aes = FALSE) + 
  
  geom_line(data = df_predict_M1,
            aes(stage -2, fit[,"upr"]), color = "darkgreen", size = 0.1, inherit.aes = FALSE) +
  
  geom_line(data = df_predict_M1, 
            aes(x = stage - 2, y = fit), 
            color = "green", linewidth = 1.2, inherit.aes = FALSE) +
  
  annotate("text", x = Inf, y = Inf, label = paste0("R²= ", round(data_r2[1,4], digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic", color = "green")+
  
  annotate("text", x = Inf, y = 0.95, label = paste0("R²= ", round(data_r2[1,5], digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic", color = "magenta") + 
  
  ylim(0, 1)+
  
  labs(
    x = "Egg chamber stages",
    y = "M1 & M2 coefficient",
    title = "Evolution of the M1 & M2 score of antisense LncRNA AAACTAT and AAACTAC",
    subtitle = "The probes used were sense"
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
# Print graphics
###################################################

print(graph_pearson)
print(graph_overlap)
print(graph_ICQ)
print(graph_manders)

ggarrange(graph_pearson,graph_overlap,graph_ICQ,graph_manders)

#################################################
#################################################
# variance test between Y

concatenated_data <- rbind(colocalization_index_stage10, colocalization_index_stage3,
                         colocalization_index_stage4, colocalization_index_stage5,
                         colocalization_index_stage6, colocalization_index_stage7, 
                         colocalization_index_stage8, colocalization_index_stage9)
concatenated_data <- concatenated_data[,c(-3,-4,-7,-8)]

pairs.panels(concatenated_data[,c(1:5)])


par(mfrow=c(3,2), oma = c(0, 0, 3, 0))

hist(concatenated_data$Coef_pearson,
     main = "Pearson sense")

hist(concatenated_data$Coef_overlap,
     main = "Overlap sense")

hist(concatenated_data$M1,
     main = "M1 sense")

hist(concatenated_data$M2,
     main = "M2 sense")

hist(concatenated_data$ICQ,
     main = "ICQ sense")

mtext("Sense", outer = TRUE, cex = 1.5, font = 2)

boxplot(concatenated_data,
        main = " Sense probe")

for ( i in c(1:5)){
  
  name = colnames(concatenated_data)[i]
  
  number_0 <- (length(which(concatenated_data[,i] == 0)))
  
  print(paste("Number of 0 for",name , number_0))
  
}

############################

# put it with the data

# pearson 


data_Per_ICQ<- Colocalization_index_all_stages[,c(1,9,10)] %>%
  pivot_longer(
    cols = colnames(Colocalization_index_all_stages[,c(1,9)]), 
    names_to = "score_ID",             
    values_to = "values"
  )

plot_ICQ_per <- ggplot(data_Per_ICQ, aes(x = factor(stage), y = values, fill = score_ID))+ 
  
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
              aes(x = stage -2 , ymin = fit[,"lwr"], ymax = fit[,"upr"]), 
              alpha = 0.09, fill = "blue", inherit.aes = FALSE) +
  geom_line(data = df_predict_pearson,
            aes(stage -2, fit[,"lwr"]), color = "darkblue", size = 0.1,
            inherit.aes = FALSE) + 
  geom_line(data = df_predict_pearson,
            aes(stage -2, fit[,"upr"]), color = "darkblue", size = 0.1,
            inherit.aes = FALSE) + 
  geom_line(data = df_predict_pearson, 
            aes(x = stage -2 , y = fit[,"fit"]),
            linewidth = 1.2, color = "blue",inherit.aes = FALSE)+
  
  geom_ribbon(data = df_predict_ICQ, 
              aes(x = stage -2 , ymin = fit[,"lwr"], ymax = fit[,"upr"]), 
              alpha = 0.09, fill = "red",inherit.aes = FALSE) +
  geom_line(data = df_predict_ICQ,
            aes(stage -2, fit[,"lwr"]), color = "darkred", size = 0.1,
            inherit.aes = FALSE) + 
  geom_line(data = df_predict_ICQ,
            aes(stage -2, fit[,"upr"]), color = "darkred", size = 0.1,
            inherit.aes = FALSE) +
  geom_line(data = df_predict_ICQ, 
            aes(x = stage -2 , y = fit[,"fit"]),
            linewidth = 1.2, color = "red",inherit.aes = FALSE)    +
  
  scale_x_discrete("Stages",labels = c("Stage3","Stage4","Stage5","Stage6",
                                       "Stage7","Stage8","Stage9","Stage10"))+
  ylim(0, 1)+
  
  labs(
    x = "Egg development stage",
    y = "Pearson & ICQ coefficient",
    title = "Evolution of colocalization of antisense LncRNA AAACTAT and AAACTAC",
    subtitle = "The probes used were sense"
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


# Poster graphic

transparency = 1

windowsFonts( A = windowsFont("baskerville old face"))

plot_manders_pres <- ggplot(data_M1_M2, aes(x = factor(stage), y = values, color = score_ID)) +
  
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
              aes(x = stage - 2, ymin = fit[,"lwr"], ymax = fit[,"upr"], y = NULL, fill = NULL), 
              alpha = 0.15, fill = "green", inherit.aes = FALSE) +
  geom_line(data = df_predict_M1,
            aes(stage -2, fit[,"lwr"]), 
            color = "darkgreen", size = 0.1, inherit.aes = FALSE) + 
  
  geom_line(data = df_predict_M1,
            aes(stage -2, fit[,"upr"]), 
            color = "darkgreen", size = 0.1, inherit.aes = FALSE) +
  
  geom_line(data = df_predict_M1, 
            aes(x = stage - 2, y = fit[,"fit"]), 
            color = "green", linewidth = 1.2, inherit.aes = FALSE) +
  
  geom_ribbon(data = df_predict_M2, 
              aes(x = stage - 2, ymin = fit[,"lwr"], ymax = fit[,"upr"], y = NULL, fill = NULL), 
              alpha = 0.15, fill = "magenta", inherit.aes = FALSE) +
  
  geom_line(data = df_predict_M2,
            aes(stage -2, fit[,"lwr"]), color = "darkmagenta", size = 0.1, inherit.aes = FALSE) + 
  
  geom_line(data = df_predict_M2,
            aes(stage -2, fit[,"upr"]), color = "darkmagenta", size = 0.1, inherit.aes = FALSE) +
  
  geom_line(data = df_predict_M2, 
            aes(x = stage - 2, y = fit[,"fit"]), 
            color = "magenta", linewidth = 1.2, inherit.aes = FALSE)  +
  
  scale_x_discrete("Stages",labels = c("Stage3","Stage4","Stage5","Stage6",
                                       "Stage7","Stage8","Stage9","Stage10"))+
  ylim(0, 1)+
  
  labs(
    x = "Egg development stage",
    y = "M1 & M2 coefficient",
    title = "Evolution of colocalization of antisense LncRNA AAACTAT and AAACTAC",
    subtitle = "The probes used were sense"
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
print(plot_manders_pres)

# boxplot for each score

par(mfrow=c(3,2), oma = c(0, 0, 3, 0))

boxplot(data = Colocalization_index_all_stages, Coef_pearson ~ stage)
boxplot(data = Colocalization_index_all_stages, Coef_overlap ~ stage)
boxplot(data = Colocalization_index_all_stages, ICQ ~ stage)
boxplot(data = Colocalization_index_all_stages, M1 ~ stage)
boxplot(data = Colocalization_index_all_stages, M2 ~ stage)

mtext("Distribution of variance according to stage \n for sense probes"
      , outer = TRUE, cex = 1.2, font = 1.5)

################################
# anova for spline
###############################


# pearson
model_pearson_1 <- lm(data = Colocalization_index_all_stages, Coef_pearson ~ stage)
model_pearson_2 <- lm(data = Colocalization_index_all_stages, Coef_pearson ~ poly(stage,2))
model_pearson_3 <- lm(data = Colocalization_index_all_stages, Coef_pearson ~ bs(stage,3))
model_pearson_4 <- lm(data = Colocalization_index_all_stages, Coef_pearson ~ bs(stage,4))
model_pearson_5 <- lm(data = Colocalization_index_all_stages, Coef_pearson ~ bs(stage,5))
model_pearson_6 <- lm(data = Colocalization_index_all_stages, Coef_pearson ~ bs(stage,6))
model_pearson_7 <- lm(data = Colocalization_index_all_stages, Coef_pearson ~ bs(stage,7))

anova(model_pearson_1,model_pearson_2, model_pearson_3, model_pearson_4, model_pearson_5, model_pearson_6, model_pearson_7)

# overlap


model_Overlap_1 <- rlm(data = Colocalization_index_all_stages, Coef_overlap ~ stage)
model_Overlap_2 <- rlm(data = Colocalization_index_all_stages, Coef_overlap ~ poly(stage,2))
model_Overlap_3 <- rlm(data = Colocalization_index_all_stages, Coef_overlap ~ bs(stage,3))
model_Overlap_4 <- rlm(data = Colocalization_index_all_stages, Coef_overlap ~ bs(stage,4))
model_Overlap_5 <- rlm(data = Colocalization_index_all_stages, Coef_overlap ~ bs(stage,5))
model_Overlap_6 <- rlm(data = Colocalization_index_all_stages, Coef_overlap ~ bs(stage,6))
model_Overlap_7 <- rlm(data = Colocalization_index_all_stages, Coef_overlap ~ bs(stage,7))

anova(model_Overlap_1,model_Overlap_2, model_Overlap_3, model_Overlap_4, model_Overlap_5, model_Overlap_6, model_Overlap_7)


# ICQ

model_ICQ_1 <- rlm(data = Colocalization_index_all_stages, ICQ ~ stage)
model_ICQ_2 <- rlm(data = Colocalization_index_all_stages, ICQ ~ poly(stage,2))
model_ICQ_3 <- rlm(data = Colocalization_index_all_stages, ICQ ~ bs(stage,3))
model_ICQ_4 <- rlm(data = Colocalization_index_all_stages, ICQ ~ bs(stage,4))
model_ICQ_5 <- rlm(data = Colocalization_index_all_stages, ICQ ~ bs(stage,5))
model_ICQ_6 <- rlm(data = Colocalization_index_all_stages, ICQ ~ bs(stage,6))
model_ICQ_7 <- rlm(data = Colocalization_index_all_stages, ICQ ~ bs(stage,7))

anova(model_ICQ_1,model_ICQ_2, model_ICQ_3, model_ICQ_4, model_ICQ_5, model_ICQ_6, model_ICQ_7)

# M1
model_M1_1 <- rlm(data = Colocalization_index_all_stages, M1 ~  stage)
model_M1_2 <- rlm(data = Colocalization_index_all_stages, M1 ~  poly(stage,2))
model_M1_3 <- rlm(data = Colocalization_index_all_stages, M1 ~  bs(stage,3))
model_M1_4 <- rlm(data = Colocalization_index_all_stages, M1 ~  bs(stage,4))
model_M1_5 <- rlm(data = Colocalization_index_all_stages, M1 ~  bs(stage,5))
model_M1_6 <- rlm(data = Colocalization_index_all_stages, M1 ~  bs(stage,6))
model_M1_7 <- rlm(data = Colocalization_index_all_stages, M1 ~  bs(stage,7))

anova(model_M1_1,model_M1_2, model_M1_3, model_M1_4, model_M1_5, model_M1_6, model_M1_7)

# M2

model_M2_1 <- rlm(data = Colocalization_index_all_stages, M1 ~  stage)
model_M2_2 <- rlm(data = Colocalization_index_all_stages, M1 ~  poly(stage,2))
model_M2_3 <- rlm(data = Colocalization_index_all_stages, M1 ~  bs(stage,3))
model_M2_4 <- rlm(data = Colocalization_index_all_stages, M1 ~  bs(stage,4))
model_M2_5 <- rlm(data = Colocalization_index_all_stages, M1 ~  bs(stage,5))
model_M2_6 <- rlm(data = Colocalization_index_all_stages, M1 ~  bs(stage,6))
model_M2_7 <- rlm(data = Colocalization_index_all_stages, M1 ~  bs(stage,7))

anova(model_M2_1,model_M2_2, model_M2_3, model_M2_4, model_M2_5, model_M2_6, model_M2_7)
