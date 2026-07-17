###################################
# 3D colocalization analysis
# Dvir48 antisense probe CTAC and CTAT
###################################

# Various statistical tests
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
library(robustbase)
library(glmmTMB)

# Graphical representation
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggeffects)

# Citation 
library(grateful)

#########################################################################
# Open files containing the indices, it is in txt format within csv

# Vector of the image numbers for each stage

stage_3 = c(1, 4, 6, 10, 15, 22, 30, 33, 35, 48, 59)
stage_4 = c(2, 3, 17, 28, 44, 52, 58)
stage_5 = c(1, 3, 4, 7, 9, 11, 13, 15, 26, 29, 32, 35, 46, 47, 57)
stage_6 = c(2, 14, 19, 25, 27, 38, 50, 53)
stage_7 = c(5, 8, 12, 20, 31, 34, 45)
stage_8 = c(21, 24, 23, 37, 39, 42, 49)
stage_9 = c(36, 40, 41, 51)
stage_10 = c(54, 55, 56)



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
  
  # Index to add a row for each occurrence in the txt file
  
  new_row <- rep(0, 9)
  
  # Initial loop index
  
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
  
  # Clean all non-numeric characters and transform into numeric characters
  
  coloc_index[] <- lapply(coloc_index, function(x) {
    as.numeric(sub(".*[:=]\\s*([-0-9eE\\.]+).*", "\\1", x))
  })
  
  return(coloc_index)
  
}

# Opening function 

opening <- function(series_n,stage,file_name){
  
  concatenated_data <- data.frame(Coef_pearson = rep(0,1),
                                
                                Coef_overlap = rep(0,1),
                                
                                k1  = rep(0,1),
                                
                                k2 = rep(0,1), 
                                
                                M1 = rep(0,1),
                                
                                M2 = rep(0,1),
                                
                                a_cytofluorogram  = rep(0,1),
                                
                                b_cytofluorogram = rep(0,1), 
                                
                                ICQ = rep(0,1),
                                
                                stringsAsFactors = FALSE)
  
  file_names <- vector()
  
  for ( i in (1:length(series_n))){
    
    file_names[[i]] <- paste0("Coef_coloc_zstack_Series_",series_n[i],"stage",stage)
    
  }
  
  
  for ( i in (1:length(series_n))){
    
    ID <- read.csv(paste0("Coef_coloc_zstack_Series_",series_n[i],"Dvir48_CTAC_CTAT_antisens_stade",stage,".csv"),
                   sep = ",", dec = ".", header =TRUE)
    
    ID <- Process_coloc_data(ID)
    
    concatenated_data[i,] <- ID
    
  }
  
  assign(file_name, concatenated_data, .GlobalEnv)
  
}


for (i in c(3:10)){
  
  x <- get(paste0("stage_",i))
  
  opening(x,i,paste0("Coloc_stage_",i))
  
}

######################################
# Long data format of scores combined together

Coloc_stage_3$stage <- rep(3,length(Coloc_stage_3[,1]))

Coloc_stage_4$stage <- rep(4,length(Coloc_stage_4[,1]))

Coloc_stage_5$stage <- rep(5,length(Coloc_stage_5[,1]))

Coloc_stage_6$stage <- rep(6,length(Coloc_stage_6[,1]))

Coloc_stage_7$stage <- rep(7,length(Coloc_stage_7[,1]))

Coloc_stage_8$stage <- rep(8,length(Coloc_stage_8[,1]))

Coloc_stage_9$stage <- rep(9,length(Coloc_stage_9[,1]))

Coloc_stage_10$stage <- rep(10,length(Coloc_stage_10[,1]))

Colocalization_index_all_stages <- rbind(Coloc_stage_3,
                                          Coloc_stage_4,
                                          Coloc_stage_5,
                                          Coloc_stage_6,
                                          Coloc_stage_7,
                                          Coloc_stage_8,
                                          Coloc_stage_9,
                                          Coloc_stage_10)

Colocalization_index_all_stages_lon <- Colocalization_index_all_stages[,c(1,2,5,6,9,10)] %>%
  pivot_longer(
    cols = colnames(Colocalization_index_all_stages[,c(1,2,5,6,9)]), 
    names_to = "score_ID",             
    values_to = "values"
  )

###################################################################
# ordbeta
hist(Colocalization_index_all_stages$Coef_pearson)

# try gaussian
hist(Colocalization_index_all_stages$Coef_overlap)

# gaussian as well
hist(Colocalization_index_all_stages$ICQ)

# ordbeta
hist(Colocalization_index_all_stages$M1)

# ordbeta
hist(Colocalization_index_all_stages$M2)

boxplot(Colocalization_index_all_stages[,c(1,2,5,6,9)])

concatenated_data <- Colocalization_index_all_stages[,c(-3,-4,-7,-8)]

for ( i in c(1:5)){
  
  name = colnames(concatenated_data)[i]
  
  number_0 <- (length(which(concatenated_data[,i] == 0)))
  
  print(paste("Number of 0 for",name , number_0))
  
}

###############################
# Try models and fitting

boxplot(Colocalization_index_all_stages$Coef_pearson ~ Colocalization_index_all_stages$stage)


#################################
# Pearson Gamma
# To see, rlm better fit in general on the data.
model_pearson <- glmmTMB(data = Colocalization_index_all_stages,
                         Coef_pearson ~ bs(stage,5),
                         family = ziGamma(link = "inverse"))

qqPlot(Colocalization_index_all_stages$Coef_pearson)

summary(model_pearson)

mod_sim_M1 <- simulateResiduals(fittedModel = model_pearson, n = 250)

plot(mod_sim_M1 ,rank = F)

par(mfrow=c(2,2), oma = c(0, 0, 4, 0))

plotResiduals(mod_sim_M1)

testOutliers(mod_sim_M1)

testDispersion(mod_sim_M1) 

testZeroInflation(mod_sim_M1) 

# Test linearity
#####
newdat <- Colocalization_index_all_stages_lon %>% 
  filter(score_ID == "Coef_pearson")
newdat$fit = fitted(model_pearson)
newdat$res = resid(model_pearson)

plot_1 <- ggplot(newdat, aes(fit, res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('Pearson linearity sense')

# On exact x values

plot_2 <- ggplot(newdat, aes(x = Colocalization_index_all_stages$stage, y = res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('Pearson linearity sense on stage')

ggarrange(plot_1,plot_2)


# Everything is ok :)

#########################
# Overlap

model_overlap <- lm(Coef_overlap ~ bs(stage,3), data = Colocalization_index_all_stages)

check_model(model_overlap)

shapiro.test(model_overlap$residuals)

check_heteroscedasticity(model_overlap)

mean(model_overlap$residuals)


# Perfect, lm with degree 6 polynomial :)

##########################
# ICQ

model_ICQ <- lm(ICQ ~ bs(stage,4), data = Colocalization_index_all_stages)

check_model(model_ICQ)

shapiro.test(model_ICQ$residuals)

check_heteroscedasticity(model_ICQ)

mean(model_ICQ$residuals)


# Everything is ok :)


########################
# M1
model_M1 <- glmmTMB(data = Colocalization_index_all_stages,
                         M1 ~ bs(stage,5),
                         family = ordbeta(link = "logit"))

summary(model_M1)

mod_sim_M1 <- simulateResiduals(fittedModel = model_M1, n = 250)

plot(mod_sim_M1 ,rank = F)

par(mfrow=c(2,2), oma = c(0, 0, 4, 0))

plotResiduals(mod_sim_M1)

testOutliers(mod_sim_M1)

testDispersion(mod_sim_M1) 

testZeroInflation(mod_sim_M1) 

# Test linearity

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

# Everything seems in order :)

#################################
# M2
model_M2 <- glmmTMB(data = Colocalization_index_all_stages,
                    M2 ~ bs(stage,5),
                    family = ordbeta(link = "logit"))

summary(model_M2)

mod_sim_M1 <- simulateResiduals(fittedModel = model_M2, n = 250)

plot(mod_sim_M1 ,rank = F)

par(mfrow=c(2,2), oma = c(0, 0, 4, 0))

plotResiduals(mod_sim_M1)

testOutliers(mod_sim_M1)

testDispersion(mod_sim_M1) 

testZeroInflation(mod_sim_M1) 

# Test linearity

newdat <- Colocalization_index_all_stages_lon %>% 
  filter(score_ID == "M2")
newdat$fit = fitted(model_M2)
newdat$res = resid(model_M2)

plot_1 <- ggplot(newdat, aes(fit, res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('M1 linearity antisense')

plot_2 <- ggplot(newdat, aes(x = Colocalization_index_all_stages$stage, res))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0, linetype='dashed', col='red4', size=1.2)+ ggtitle('M1 linearity antisense on stage')

ggarrange(plot_1,plot_2)

# The linearity of the residual response is strange but! it's chill

############################
# Launch projections on 
# real data.


summary(model_pearson)

summary(model_overlap)

summary(model_ICQ)

summary(model_M1)

summary(model_M2)

# Calculate R2 and pseudo R2

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

predict_response(model_pearson, terms = "stage")
predict_response(model_overlap, terms = "stage")
predict_response(model_ICQ, terms = "stage")
predict_response(model_M1, terms = "stage [all]")
predict_response(model_M2, terms = "stage [all]")


# Increase the detail of the prediction
# to fill the discrete categories

stage_dense <- seq(min(Colocalization_index_all_stages$stage), 
                   max(Colocalization_index_all_stages$stage), 
                   length.out = 200)

stage_dense_for_pearson <- seq(min(Colocalization_index_all_stages$stage), 
                   max(Colocalization_index_all_stages$stage), 
                   length.out = 200)

###########################
# Pearson

df_predict_pearson <- expand.grid(
  stage = stage_dense_for_pearson,
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


# Overlap

df_predict_Overlap <- expand.grid(
  stage = stage_dense,
  Score = "Coef_Overlap" 
)

pred_Overlap <- predict(model_overlap, newdata = df_predict_Overlap , 
                        se.fit = TRUE , interval = "confidence", level = 0.95)

df_predict_Overlap <- df_predict_Overlap %>%
  bind_cols(as_tibble(pred_Overlap))

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

pred_M2 <- predict(model_M2, newdata = df_predict_M2 , 
                   se.fit = TRUE, type = "link")

df_predict_M2 <- df_predict_M2 %>%
  mutate(
    fit   = model_M2$modelInfo$family$linkinv(pred_M2$fit), 
    upper = model_M2$modelInfo$family$linkinv(pred_M2$fit + (1.96 * pred_M2$se.fit)),
    lower = model_M2$modelInfo$family$linkinv(pred_M2$fit - (1.96 * pred_M2$se.fit))
  )

#####################################
# Plotting the graphs


#####################
# Plot everyone to observe general behavior

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
# Graph ggplot pearson coef
##############################
# Put it with the data

# Pearson 

windowsFonts( A = windowsFont("Arial"))

graph_pearson <- ggplot(data = Colocalization_index_all_stages)+ 
  geom_boxplot(aes(x = as.factor(stage), y = Coef_pearson), color = "blue", alpha = 0.2, fill = NA) + 
  geom_point(aes(x = stage - 2, y = Coef_pearson), color = "blue", position = position_jitter(width = 0.1)) + 
  geom_ribbon(data = df_predict_pearson, 
              aes(x = stage -2 , ymin = lower, ymax = upper), 
              alpha = 0.2, inherit.aes = FALSE, fill = "blue") +
  geom_line(data = df_predict_pearson, 
            aes(x = stage -2 , y = fit),
            linewidth = 1.2, color = "black")   +
  geom_line(data = df_predict_pearson, 
            aes(x = stage -2 , y = upper),
            linewidth = 0.5, color = "darkblue")   +
  geom_line(data = df_predict_pearson, 
            aes(x = stage -2 , y = lower),
            linewidth = 0.5, color = "darkblue")   +
  
  scale_x_discrete("Stage",labels = c("Stage3","Stage4","Stage5","Stage6",
                                       "Stage7","Stage8","Stage9","Stage10"))+
  ylim(0, 1)+
  
  labs(
    x = "Egg chamber stages",
    y = "Pearson coefficient (r)",
    title = "Evolution of the Pearson coefficient of sense LncRNA AAACTAT and AAACTAC",
    subtitle = "The probes used were antisense"
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
  geom_boxplot(aes(x = as.factor(stage), y = Coef_overlap), color = "orange", alpha = 0.5, fill = NA) + 
  geom_point(aes(x = stage - 2, y = Coef_overlap), color = "orange", 
             position = position_jitter(width = 0.1)) + 
  geom_ribbon(data = df_predict_Overlap, 
              aes(x = stage -2 , ymin = fit[,"lwr"], ymax = fit[,"upr"]), 
              alpha = 0.2, inherit.aes = FALSE, fill = "orange") +
  geom_line(data = df_predict_Overlap, 
            aes(x = stage -2 , y = fit[,"fit"]),
            linewidth = 1.2, color = "black")   +
  geom_line(data = df_predict_Overlap, 
            aes(x = stage -2 , y = fit[,"lwr"]),
            linewidth = 0.5, color = "darkorange",
            inherit.aes = FALSE)   +
  geom_line(data = df_predict_Overlap, 
            aes(x = stage -2 , y = fit[,"upr"]),
            linewidth = 0.5, color = "darkorange",
            inherit.aes = FALSE)   +
  
  scale_x_discrete("Stage",labels = c("Stage3","Stage4","Stage5","Stage6",
                                       "Stage7","Stage8","Stage9","Stage10"))+
  ylim(0, 1)+
  
  labs(
    x = "Egg chamber stages",
    y = "MOC coefficient",
    title = "Evolution of the MOC coefficient of sense LncRNA AAACTAT and AAACTAC",
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
  
  annotate("text", x = Inf, y = Inf, label = paste0("R²= ", round(data_r2[1,2],digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic")

# ICQ

graph_ICQ <- ggplot(data = Colocalization_index_all_stages)+ 
  geom_boxplot(aes(x = as.factor(stage), y = ICQ), color = "red", alpha = 0.5,fill = NA) + 
  geom_point(aes(x = stage - 2, y = ICQ),color = "red",
             position = position_jitter(width = 0.1)) + 
  geom_ribbon(data = df_predict_ICQ, 
              aes(x = stage -2 , ymin = fit[,"lwr"], ymax = fit[,"upr"]), 
              alpha = 0.2, fill = "red",inherit.aes = FALSE) +
  geom_line(data = df_predict_ICQ, 
            aes(x = stage -2 , y = fit[,"fit"]),
            linewidth = 1.2, color = "black")    +
  geom_line(data = df_predict_ICQ, 
            aes(x = stage -2 , y = fit[,"lwr"]),
            linewidth = 0.5, color = "darkred",
            inherit.aes = FALSE)   +
  geom_line(data = df_predict_ICQ, 
            aes(x = stage -2 , y = fit[,"upr"]),
            linewidth = 0.5, color = "darkred",
            inherit.aes = FALSE)   +
  geom_hline(yintercept = 0, linetype = "dotted",
             size = 1) + 
  
  annotate("text", x = Inf, y = Inf, label = paste0("R²= ", round(data_r2[1,3], digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") + 
  
  scale_x_discrete("Egg chamber stages",labels = c("Stage3","Stage4","Stage5","Stage6",
                                       "Stage7","Stage8","Stage9","Stage10"))+
  ylim(-0.5, 0.5)+
  
  labs(
    x = "Egg chamber stages",
    y = "ICQ score",
    title = "Evolution of the ICQ score through development \n of oocyte for the forward transcript AAACTAT and AAACTAC",
    subtitle = "Both probes were reversed"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "right",
    text = element_text(size = 12, family = "A"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "A"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain", family = "A"),
    plot.margin = margin(20, 20, 20, 20),
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold")
  )


# Manders

# M1 = CTAT, M2 = CTAC

data_M1_M2 <- Colocalization_index_all_stages[,c(5,6,10)] %>%
  pivot_longer(
    cols = colnames(Colocalization_index_all_stages[,c(5,6)]), 
    names_to = "score_ID",             
    values_to = "values"
  )

transparency = 1

windowsFonts( A = windowsFont("Arial"))

graph_manders <- ggplot(data_M1_M2, aes(x = factor(stage), y = values, color = score_ID)) +
  
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
              aes(x = stage - 2, ymin = lower, ymax = upper, y = NULL, fill = NULL), 
              alpha = 0.15, fill = "green", inherit.aes = FALSE) +
  
  geom_line(data = df_predict_M2, 
            aes(x = stage - 2, y = fit), 
            color = "green", linewidth = 1.2, inherit.aes = FALSE) +
  geom_line(data = df_predict_M2, 
            aes(x = stage - 2, y = upper), 
            color = "darkgreen", linewidth = 0.5, inherit.aes = FALSE) +
  geom_line(data = df_predict_M2, 
            aes(x = stage - 2, y = lower), 
            color = "darkgreen", linewidth = 0.5, inherit.aes = FALSE) +
  
  geom_ribbon(data = df_predict_M1, 
              aes(x = stage - 2, ymin = lower, ymax = upper, y = NULL, fill = NULL), 
              alpha = 0.15, fill = "magenta", inherit.aes = FALSE) +
  geom_line(data = df_predict_M1, 
            aes(x = stage - 2, y = fit), 
            color = "magenta", linewidth = 1.2, inherit.aes = FALSE)  +
  geom_line(data = df_predict_M1, 
            aes(x = stage - 2, y = upper), 
            color = "darkmagenta", linewidth = 0.5, inherit.aes = FALSE)  +
  geom_line(data = df_predict_M1, 
            aes(x = stage - 2, y = lower), 
            color = "darkmagenta", linewidth = 0.5, inherit.aes = FALSE)  +
  
  annotate("text", x = Inf, y = Inf, label = paste0("R²= ", round(data_r2[1,4],digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic", 
           color = "magenta") + 
  
  annotate("text", x = Inf, y = 0.95, label = paste0("R²= ", round(data_r2[1,5],digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic", 
           color = "green") + 
  
  scale_x_discrete("Stages",labels = c("Stage3","Stage4","Stage5","Stage6",
                                       "Stage7","Stage8","Stage9","Stage10"))+
  ylim(0, 1)+
  
  labs(
    x = "Egg chamber stages",
    y = "coefficient of M1 & M2",
    title = "Evolution of the coefficient of M1 & M2 for sense LncRNA AAACTAT and AAACTAC",
    subtitle = "The probes used were antisense"
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
print(graph_manders)

ggarrange(graph_pearson,graph_overlap,graph_ICQ,graph_manders)


#########################################################



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
  
  geom_ribbon(data = df_predict_pearson, 
              aes(x = stage -2 , ymin = lower, ymax = upper), 
              alpha = 0.09, fill = "blue", inherit.aes = FALSE) +
  geom_line(data = df_predict_pearson,
            aes(stage -2, lower), color = "darkblue", size = 0.1,
            inherit.aes = FALSE) + 
  geom_line(data = df_predict_pearson,
            aes(stage -2, upper), color = "darkblue", size = 0.1,
            inherit.aes = FALSE) + 
  geom_line(data = df_predict_pearson, 
            aes(x = stage -2 , y = fit),
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
    y = "coefficient p & ICQ",
    title = "Evolution of the 3D colocalization of sense LncRNA AAACTAT and AAACTAC",
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
  
  scale_fill_manual(values = c("Coef_pearson" = "blue", "ICQ" = "red")) +
  scale_color_manual(values = c("Coef_pearson" = "blue", "ICQ" = "red")) +
  
  annotate("text", x = Inf, y = Inf, label = paste0("R²= ", round(data_r2[1,3], digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic", color = "red") +
  annotate("text", x = Inf, y = 0.98, label = paste0("R²= ", round(data_r2[1,1], digits = 3)), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic", color = "blue")

print(plot_ICQ_per)


#########################################

