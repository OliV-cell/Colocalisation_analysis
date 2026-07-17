###########################
#Comparison analysis script
#CTAC with Virlis, americana
#novamexicana and lummei
################################

#activate packages

library(car)
library(dplyr)
library(tidyr)
library(mgcv)
library(performance)
library(see)
library(dsm)

#graphical representation
library(ggplot2)

########################################################
#opening files

data_dlumei <- read.csv("Donner_Dlumei_CTAC_antisens.csv",
                        sep = ";", dec = ",", header =TRUE)

data_dvir  <- read.csv("Donner_Dvir_CTAC_antisens.csv",
                        sep = ";", dec = ",", header =TRUE)

data_dnova  <- read.csv("Donner_Dnova_CTAC_antisens.csv",
                       sep = ";", dec = ",", header =TRUE)

data_damer  <- read.csv("Donner_Damer_CTAC_antisens.csv",
                       sep = ";", dec = ",", header =TRUE)


#Dlumei processing
count_occurrences <- function(data_source, data_well){
  for (i in c(3:10)){
    #null
    data_well[[i-2,1]] <- length(which(data_source[,i] == "rien"))
    
    #low
    data_well[[i-2,2]] <- length(which(data_source[,i] == "faible"))
    
    #medium
    data_well[[i-2,3]] <- length(which(data_source[,i] == "moyen"))
    
    #high
    data_well[[i-2,4]] <- length(which(data_source[,i] == "elever"))
  }
  
  return(data_well)
  
}

Dlumei_rep1 <- data_dlumei[1:37,]
  
Dlumei_rep2 <- data_dlumei[38:98,]

Dlumei_1 <- matrix( nrow = 8, ncol = 4)

Dlumei_2 <- matrix( nrow = 8, ncol = 4)

Dlumei_1 <- count_occurrences(data_source = data_dlumei,data_well = Dlumei_1)
Dlumei_2 <- count_occurrences(data_source = data_dlumei,data_well = Dlumei_2)

matrix_names <- c("Dlumei_1","Dlumei_2")

total_observations <- sum(sapply(mget(matrix_names), sum, na.rm = TRUE))

print(total_observations)

merged_matrix <- Reduce("+", mget(matrix_names))

count_per_row <- rowSums(merged_matrix, na.rm = TRUE)

stage_names <- paste("Stage", 3:(3 + length(count_per_row) - 1))
names(count_per_row) <- stage_names


print(count_per_row)

data_score_lumei <- data.frame(score_1 = rep(0,8),
                               
                               score_2 = rep(0,8),
                               
                               stade = c(3:10),
                               
                               row.names = c("stage 3",
                                             "stage 4","stage 5","stage 6",
                                             "stage 7","stage 8","stage 9",
                                             "stage 10"))

for ( i in c(1:8)){
  
  data_score_lumei[[i,1]] <- (Dlumei_1[i,2] + 2*Dlumei_1[i,3] + 3*Dlumei_1[i,4])
  
  data_score_lumei[[i,2]] <- (Dlumei_2[i,2] + 2*Dlumei_2[i,3] + 3*Dlumei_2[i,4])
  
  print(paste("Loop",i,"done"))
}

data_score_lumei_normalized <- data_score_lumei

for (i in c(1:8)){
  
  data_score_lumei_normalized[[i,1]] <- data_score_lumei_normalized[i,1]/(sum(Dlumei_1[i,1:4]) +1)
  data_score_lumei_normalized[[i,2]] <- data_score_lumei_normalized[i,2]/(sum(Dlumei_2[i,1:4]) +1)
  
}

data_dlumei <- data_score_lumei_normalized

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

#data exploration

plot(data_dvir$score)

plot(data_dnova$score)

plot(data_damer$score)

plot(data_dlumei$score)

hist(data_dvir$score)

hist(data_dnova$score)

hist(data_damer$score)

hist(data_dlumei$score)

concatenated_data <- rbind(data_dvir, data_dnova, data_damer, data_dlumei)

concatenated_data <- rename(concatenated_data, stage = stade)
concatenated_data <- rename(concatenated_data, strand = brin)

boxplot(concatenated_data$score ~ concatenated_data$sp)
leveneTest(concatenated_data$score ~ concatenated_data$sp)

#try with GAM

concatenated_data$sp <- as.factor(concatenated_data$sp)

combined_model <- gam(list(score ~ sp + s(stage, by = sp, k = 10), 
    ~ s(stage)), 
  data = concatenated_data,
  optimizer = c("outer", "newton"),
  method = "REML", 
  family = gaulss(),
  control = list(maxit = 1000)
)

#Model selection function

par(mfrow=c(2,2))

plot(combined_model, 
     pages = 1, 
     scheme = 1, 
     all.terms = TRUE)

summary.gam(combined_model)

par(mfrow=c(2,2))

gam.check(combined_model)

k.check(combined_model)

performance(combined_model)

concurvity(combined_model, full = FALSE)

vis_concurvity(combined_model, type = "estimate")

model_performance(combined_model)

dense_stage <- seq(min(concatenated_data$stage), 
                   max(concatenated_data$stage), 
                   length.out = 200)

df_predict_combine <- expand.grid(
  stage = dense_stage,
  sp = unique(concatenated_data$sp) 
)

preds_smooth <- predict(combined_model, newdata = df_predict_combine, 
                        se.fit = TRUE, type = "link")

ilink <- family(combined_model)$linkinv

df_predict_combine <- df_predict_combine %>%
  mutate(
    fit   = preds_smooth$fit[,1],  # Mean
    upper = preds_smooth$fit[,1] + (1.96 * preds_smooth$se.fit[,1]),
    lower = preds_smooth$fit[,1] - (1.96 * preds_smooth$se.fit[,1])
  )

# R^2 calculation per sp

concatenated_data$preds <- predict(combined_model, type = "response")[,1]

R2_table <- concatenated_data %>%
  group_by(sp) %>%
  summarize(
    RSS = sum((score - preds)^2),                  
    TSS = sum((score - mean(score))^2),           
    R2  = round(1 - (RSS / TSS), 4)               
  )

print(R2_table)

######################
#general difference statistics, by stage and expression slope

anova <- anova(combined_model)
anova.gam(combined_model)

p_val_strand <- anova$pTerms.pv[1]

legend_names <- c(
  "Dvir48" = "D.Virilis",
  "Dnova"  = "D.novamexicana",
  "Damer" = "D.americana",
  "Dlumei" = "D.lummei"
)

#########################

windowsFonts( A = windowsFont("Arial"))

#graphical representation
plot_GAM <- ggplot() + 
  
  geom_point(data = concatenated_data, 
             aes(x = stage, y = score, color = sp),
             alpha = 0.6, position = position_jitter(width = 0.1)) + 
  geom_ribbon(data = df_predict_combine, 
              aes(x = stage, ymin = lower, ymax = upper, fill = sp), 
              alpha = 0.09) +
  geom_line(data = df_predict_combine, 
            aes(x = stage, y = fit, color = sp),
            linewidth = 1.2) +
  
  scale_color_manual(labels = legend_names,
                     values = c("Dvir48" = "red", 
                                "Dnova" = "blue", 
                                "Damer" = "orange",
                                "Dlumei" = "darkgreen")) +
  
  scale_fill_manual(labels = legend_names,
                    values = c("Dvir48" = "red", 
                               "Dnova" = "blue", 
                               "Damer" = "orange",
                               "Dlumei" = "darkgreen")) +
  
  labs(
    x = "Egg developmental stage",
    y = "Normalized expression score",
    title = "Evolution of forward LncRNA transcription of AAACTAC for \n D.Virilis, D.americana, D.novamexicana & D.lummei",
    subtitle = "The probes used were antisense",
    color = "Species",
    fill = "Species"
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
  
  #annotate("text", x = Inf, y = Inf, label = paste0("Global ANOVA: p = ", format.pval(p_val_strand, digits = 2)), 
           #hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") + 
  
  annotate("text", x = Inf, y = 4.4, label = paste0("R² = ", round(R2_table[4,4], digits = 3), "***"), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic",
           color = "red") + 
  
  annotate("text", x = Inf, y = 4, label = paste0("R² = ", round(R2_table[3,4], digits = 3), "***" ), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic",
           color = "blue") +
  
  annotate("text", x = Inf, y = 3.6, label = paste0("R² = ", round(R2_table[1,4], digits = 3) , "***"), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic",
           color = "orange") + 
  
  annotate("text", x = Inf, y = 3.2, label = paste0("R² = ", round(R2_table[2,4], digits = 3) , "***"), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic",
           color = "darkgreen")
  

print(plot_GAM)

##########################################
#comparison of curves with non-parametric tests 
#lack of data

pvalues_table <- data.frame()

for (s in 3:10) {
  stage_data <- subset(concatenated_data, stage == s)
  
  if(length(unique(stage_data$sp)) > 1) {
    
    fit <- aov(score ~ sp, data = stage_data)
    
    tukey <- TukeyHSD(fit)$sp

    temp <- data.frame(
      Stage = s,
      Comparison = rownames(tukey),
      p_value = round(tukey[, "p adj"], 4)
    )
    pvalues_table <- rbind(pvalues_table, temp)
  }
}

print(pvalues_table)

signif_data <- pvalues_table %>%
  filter(p_value < 0.05) %>%
  mutate(
    label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

h_bracket <- 3

plot_GAM +
  
  #stage3-4-5 
  
  annotate("rect", xmin = 2.5, xmax = 5.2, ymin = 0, ymax = 3,
           linetype = "dashed", color = "black", fill = NA) +
  
  annotate("text", x = 4, y = 3 + 0.3, 
           label = "***", size = 5) +
  
  
  #stage7
  
  annotate("segment", x = 6.8, xend = 7.2, y = h_bracket, yend = h_bracket, 
           linetype = "dashed", color = "black") +
  
  annotate("segment", x = 6.8, xend = 6.8, y = 0.49, yend = h_bracket, 
           linetype = "dashed", color = "black") +
  
  annotate("segment", x = 7.2, xend = 7.2, y = 1.025, yend = h_bracket, 
           linetype = "dashed", color = "black") +
  
  annotate("text", x = 7, y = h_bracket + 0.3, 
           label = signif_data$label[4], size = 5) 
