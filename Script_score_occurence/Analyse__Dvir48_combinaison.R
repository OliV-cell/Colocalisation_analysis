#################################
################################
#Graph of combinations with CTAC
###############################
###############################
#package activation

library(car)
library(dplyr)
library(tidyr)
library(psych)
library(qpcR)
library(mgcv)
library(gvlma)
library(performance)
library(see)
#graphical representation
library(ggplot2)



count_occurrences <- function(data_source, data_well){
  for (i in c(3:12)){
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

data_CTAC_withcaac_antisense <- read.csv("CTAC_avec_caac.csv",
                                        sep=";", dec=",", header=FALSE)

data_CTAC_withctat_antisense <- read.csv("CTAC_avec_ctat.csv",
                                        sep=";", dec=",", header=FALSE)

data_CTAC_withctat_antisense_rep3 <- read.csv("CTAC_avec_ctat_rep3.csv",
                                             sep=";", dec=",", header=FALSE)

data_CTAC_withttac_antisense <- read.csv("CTAC_avec_ttac.csv",
                                        sep=";", dec=",", header=FALSE)

data_CTAC_withttac_antisense_rep2 <- read.csv("CTAC_avec_ttac_rep2.csv",
                                             sep=";", dec=",", header=FALSE)

data_CTAC_withtat_sense <- read.csv("CTAC_sens_avec_ctat_sens.csv",
                                   sep=";", dec=",", header=FALSE)

data_CTAC_withtat_sense_rep1 <- read.csv("CTAC_sens_avec_ctat_sens_rep1.csv",
                                        sep=";", dec=",", header=FALSE)

data_TTACwithCTAC_antisense <- read.csv("TTAC_avec_CTAC.csv",
                                       sep=";", dec=",", header=FALSE)

data_TTACwithCTAC_antisense_rep2 <- read.csv("TTAC_avec_CTAC_rep2.csv",
                                            sep=";", dec=",", header=FALSE)

data_CTATwithCTAC_antisense <- read.csv("CTAT_avec_CTAC.csv",
                                       sep=";", dec=",", header=FALSE)

data_CTATwithCTAC_antisense_rep3 <- read.csv("CTAT_avec_CTAC_rep3.csv",
                                            sep=";", dec=",", header=FALSE)

data_CTAT_withCTAC_sense <- read.csv("CTAT_sens_avec_CTAC_sens.csv",
                                    sep=";", dec=",", header=FALSE)

data_CTAT_withCTAC_sense_rep1 <- read.csv("CTAT_sens_avec_CTAC_sens_rep1.csv",
                                         sep=";", dec=",", header=FALSE)

data_CAACwithCTAC_antisense <- read.csv("CAAC_avec_CTAC.csv",
                                       sep=";", dec=",", header=FALSE)

#data_CAACwithCTAC_antisense <- read.csv("Analyse_photo_ovaire_R_ttac.csv",
#sep=";", dec=",", header=FALSE)

CTAC_withcaac_antisense <- matrix( nrow = 10, ncol = 4)

CTAC_withttac_antisense <- matrix( nrow = 10, ncol = 4)

CTAC_withttac_antisense_rep2 <- matrix( nrow = 10, ncol = 4)

CTAC_withctat_antisense <- matrix( nrow = 10, ncol = 4)

CTAC_withctat_antisense_rep3 <- matrix( nrow = 10, ncol = 4)

TTACwithCTAC_antisense <- matrix( nrow = 10, ncol = 4)

TTACwithCTAC_antisense_rep2 <- matrix( nrow = 10, ncol = 4)

CTATwithCTAC_antisense <- matrix( nrow = 10, ncol = 4)

CTATwithCTAC_antisense_rep3 <- matrix( nrow = 10, ncol = 4)

CAACwithCTAC_antisense <- matrix( nrow = 10, ncol = 4)

CTAC_withctat_sense <- matrix( nrow = 10, ncol = 4)

CTAC_withctat_sense_rep1 <- matrix( nrow = 10, ncol = 4)

CTAT_withCTAC_sense <- matrix( nrow = 10, ncol = 4)

CTAT_withCTAC_sense_rep1 <- matrix( nrow = 10, ncol = 4)

#Counting occurrences of different intensities

CTAC_withcaac_antisense <- count_occurrences(data_source = data_CTAC_withcaac_antisense,data_well = CTAC_withcaac_antisense)

CTAC_withttac_antisense <- count_occurrences(data_source = data_CTAC_withttac_antisense,data_well = CTAC_withttac_antisense)

CTAC_withttac_antisense_rep2 <- count_occurrences(data_source = data_CTAC_withttac_antisense_rep2,data_well = CTAC_withttac_antisense_rep2)

CTAC_withctat_antisense <- count_occurrences(data_source = data_CTAC_withctat_antisense,data_well = CTAC_withctat_antisense)

CTAC_withctat_antisense_rep3 <- count_occurrences(data_source = data_CTAC_withctat_antisense_rep3,data_well = CTAC_withctat_antisense_rep3)

TTACwithCTAC_antisense <- count_occurrences(data_source = data_TTACwithCTAC_antisense,data_well = TTACwithCTAC_antisense )

TTACwithCTAC_antisense_rep2 <- count_occurrences(data_source = data_TTACwithCTAC_antisense_rep2,data_well = TTACwithCTAC_antisense_rep2 )

CTATwithCTAC_antisense <- count_occurrences(data_source = data_CTATwithCTAC_antisense,data_well = CTATwithCTAC_antisense)

CTATwithCTAC_antisense_rep3 <- count_occurrences(data_source = data_CTATwithCTAC_antisense_rep3,data_well = CTATwithCTAC_antisense_rep3)

CAACwithCTAC_antisense <- count_occurrences(data_source = data_CAACwithCTAC_antisense,data_well = CAACwithCTAC_antisense)

CTAC_withctat_sense <- count_occurrences(data_source = data_CTAC_withtat_sense,data_well = CTAC_withctat_sense)

CTAC_withctat_sense_rep1 <- count_occurrences(data_source = data_CTAC_withtat_sense_rep1,data_well = CTAC_withctat_sense_rep1)

CTAT_withCTAC_sense <- count_occurrences(data_source = data_CTAT_withCTAC_sense ,data_well = CTAT_withCTAC_sense)

CTAT_withCTAC_sense_rep1 <- count_occurrences(data_source = data_CTAT_withCTAC_sense_rep1 ,data_well = CTAT_withCTAC_sense_rep1)

matrix_names <- c("CTAC_withcaac_antisense", "CTAC_withttac_antisense", "CTAC_withttac_antisense_rep2",
                   "CTAC_withctat_antisense", "CTAC_withctat_antisense_rep3")

total_observations <- sum(sapply(mget(matrix_names), sum, na.rm = TRUE))

print(total_observations)

merged_matrix <- Reduce("+", mget(matrix_names))

count_per_row <- rowSums(merged_matrix, na.rm = TRUE)

stage_names <- paste("Stage", 1:(1 + length(count_per_row) - 1))
names(count_per_row) <- stage_names


print(count_per_row)


#rename columns

colnames(CTAC_withcaac_antisense) <- c("null", "low", "medium", "high")
colnames(CTAC_withttac_antisense) <- c("null", "low", "medium", "high")
colnames(CTAC_withttac_antisense_rep2) <- c("null", "low", "medium", "high")
colnames(CTAC_withctat_antisense) <- c("null", "low", "medium", "high")
colnames(CTAC_withctat_antisense_rep3) <- c("null", "low", "medium", "high")
colnames(TTACwithCTAC_antisense) <- c("null", "low", "medium", "high")
colnames(TTACwithCTAC_antisense_rep2) <- c("null", "low", "medium", "high")
colnames(CTATwithCTAC_antisense) <- c("null", "low", "medium", "high")
colnames(CTATwithCTAC_antisense_rep3) <- c("null", "low", "medium", "high")
colnames(CAACwithCTAC_antisense) <- c("null", "low", "medium", "high")
colnames(CTAC_withctat_sense) <- c("null", "low", "medium", "high")
colnames(CTAT_withCTAC_sense) <- c("null", "low", "medium", "high")
colnames(CTAC_withctat_sense_rep1) <- c("null", "low", "medium", "high")
colnames(CTAT_withCTAC_sense_rep1) <- c("null", "low", "medium", "high")

#Score calculation loop for each stage and each satellite

score_occurrence_data_combination <- data.frame(score_CTAC_withcaac_antisense = rep(0,10),
                                               
                                               score_CTAC_withttac_antisense = rep(0,10),
                                               
                                               score_CTAC_withctat_antisense = rep(0,10),
                                               
                                               score_TTACwithCTAC_antisense = rep(0,10),
                                               
                                               score_CTATwithCTAC_antisense = rep(0,10),
                                               
                                               score_CAACwithCTAC_antisense = rep(0,10), 
                                               
                                               score_CTAC_withctat_sense = rep(0,10), 
                                               
                                               score_CTAT_withCTAC_sense = rep(0,10),
                                               
                                               score_CTATwithCTAC_antisense_rep3 = rep(0,10), 
                                               
                                               score_CTAC_withctat_antisense_rep3 = rep(0,10),
                                               
                                               score_CTAC_withttat_antisense_rep2 = rep(0,10),
                                               
                                               score_TTACwithCTAC_antisense_rep2 = rep(0,10), 
                                               
                                               score_CTAC_withctat_sense_rep1 = rep(0,10), 
                                               
                                               score_CTAT_withCTAC_sense_rep1 = rep(0,10),
                                               
                                               stage = c(1:10),
                                               
                                               row.names = c("stage 1","stage 2","stage 3",
                                                             "stage 4","stage 5","stage 6",
                                                             "stage 7","stage 8","stage 9",
                                                             "stage 10"))


for ( i in c(1:10)){
  
  score_occurrence_data_combination[[i,1]] <- (CTAC_withcaac_antisense[i,2] + 2*CTAC_withcaac_antisense[i,3] + 3*CTAC_withcaac_antisense[i,4])
  
  score_occurrence_data_combination[[i,2]] <- (CTAC_withttac_antisense[i,2] + 2*CTAC_withttac_antisense[i,3] + 3*CTAC_withttac_antisense[i,4]) 
  
  score_occurrence_data_combination[[i,3]] <- (CTAC_withctat_antisense[i,2] + 2*CTAC_withctat_antisense[i,3] + 3*CTAC_withctat_antisense[i,4])
  
  score_occurrence_data_combination[[i,4]] <- (TTACwithCTAC_antisense[i,2] + 2*TTACwithCTAC_antisense[i,3] + 3*TTACwithCTAC_antisense[i,4])
  
  score_occurrence_data_combination[[i,5]] <- (CTATwithCTAC_antisense[i,2] + 2*CTATwithCTAC_antisense[i,3] + 3*CTATwithCTAC_antisense[i,4])
  
  score_occurrence_data_combination[[i,6]] <- (CAACwithCTAC_antisense[i,2] + 2*CAACwithCTAC_antisense[i,3] + 3*CAACwithCTAC_antisense[i,4])
  
  score_occurrence_data_combination[[i,7]] <- (CTAC_withctat_sense[i,2] + 2*CTAC_withctat_sense[i,3] + 3*CTAC_withctat_sense[i,4])
  
  score_occurrence_data_combination[[i,8]] <- (CTAT_withCTAC_sense[i,2] + 2*CTAT_withCTAC_sense[i,3] + 3*CTAT_withCTAC_sense[i,4])
  
  score_occurrence_data_combination[[i,9]] <- (CTATwithCTAC_antisense_rep3[i,2] + 2*CTATwithCTAC_antisense_rep3[i,3] + 3*CTATwithCTAC_antisense_rep3[i,4])
  
  score_occurrence_data_combination[[i,10]] <- (CTAC_withctat_antisense_rep3[i,2] + 2*CTAC_withctat_antisense_rep3[i,3] + 3*CTAC_withctat_antisense_rep3[i,4])
  
  score_occurrence_data_combination[[i,11]] <- (CTAC_withttac_antisense_rep2[i,2] + 2*CTAC_withttac_antisense_rep2[i,3] + 3*CTAC_withttac_antisense_rep2[i,4]) 
  
  score_occurrence_data_combination[[i,12]] <- (TTACwithCTAC_antisense_rep2[i,2] + 2*TTACwithCTAC_antisense_rep2[i,3] + 3*TTACwithCTAC_antisense_rep2[i,4])
  
  score_occurrence_data_combination[[i,13]] <- (CTAC_withctat_sense_rep1[i,2] + 2*CTAC_withctat_sense_rep1[i,3] + 3*CTAC_withctat_sense_rep1[i,4])
  
  score_occurrence_data_combination[[i,14]] <- (CTAT_withCTAC_sense_rep1[i,2] + 2*CTAT_withCTAC_sense_rep1[i,3] + 3*CTAT_withCTAC_sense_rep1[i,4])
  
  print(paste("Loop",i,"done"))
}

score_occurrence_data_normalize_combination <- score_occurrence_data_combination

for (i in c(1:10)){
  
  score_occurrence_data_normalize_combination[[i,1]] <- score_occurrence_data_normalize_combination[i,1]/(sum(CTAC_withcaac_antisense[i,1:4]) +1)
  
  score_occurrence_data_normalize_combination[[i,2]] <- score_occurrence_data_normalize_combination[i,2]/(sum(CTAC_withttac_antisense[i,1:4]) +1)
  
  score_occurrence_data_normalize_combination[[i,3]] <- score_occurrence_data_normalize_combination[i,3]/(sum(CTAC_withctat_antisense[i,1:4]) +1)
  
  score_occurrence_data_normalize_combination[[i,4]] <- score_occurrence_data_normalize_combination[i,4]/(sum(TTACwithCTAC_antisense[i,1:4]) +1)
  
  score_occurrence_data_normalize_combination[[i,5]] <- score_occurrence_data_normalize_combination[i,5]/(sum(CTATwithCTAC_antisense[i,1:4]) +1)
  
  score_occurrence_data_normalize_combination[[i,6]] <- score_occurrence_data_normalize_combination[i,6]/(sum(CAACwithCTAC_antisense[i,1:4]) +1)
  
  score_occurrence_data_normalize_combination[[i,7]] <- score_occurrence_data_normalize_combination[i,7]/(sum(CTAC_withctat_sense[i,1:4]) +1)
  
  score_occurrence_data_normalize_combination[[i,8]] <- score_occurrence_data_normalize_combination[i,8]/(sum(CTAT_withCTAC_sense[i,1:4]) +1)
  
  score_occurrence_data_normalize_combination[[i,9]] <- score_occurrence_data_normalize_combination[i,9]/(sum(CTATwithCTAC_antisense_rep3[i,1:4]) +1)
  
  score_occurrence_data_normalize_combination[[i,10]] <- score_occurrence_data_normalize_combination[i,10]/(sum(CTAC_withctat_antisense_rep3[i,1:4]) +1)
  
  score_occurrence_data_normalize_combination[[i,11]] <- score_occurrence_data_normalize_combination[i,11]/(sum(CTAC_withttac_antisense_rep2[i,1:4]) +1)
  
  score_occurrence_data_normalize_combination[[i,12]] <- score_occurrence_data_normalize_combination[i,12]/(sum(TTACwithCTAC_antisense_rep2[i,1:4]) +1)
  
  score_occurrence_data_normalize_combination[[i,13]] <- score_occurrence_data_normalize_combination[i,13]/(sum(CTAC_withctat_sense_rep1[i,1:4]) +1)
  
  score_occurrence_data_normalize_combination[[i,14]] <- score_occurrence_data_normalize_combination[i,14]/(sum(CTAT_withCTAC_sense_rep1[i,1:4]) +1)
  
}

#############################
#compilation of all CTAC
#CTAT sense and antisense
############################

CTAC_antisense_compile <- data.frame(Null = rep(0,10),
                                     
                                     Low = rep(0,10),
                                     
                                     Medium  = rep(0,10),
                                     
                                     High = rep(0,10), 
                                     
                                     row.names = c("stage 1","stage 2","stage 3",
                                                   "stage 4","stage 5","stage 6",
                                                   "stage 7","stage 8","stage 9",
                                                   "stage 10"))

CTAT_antisense_compile <- data.frame(Null = rep(0,10),
                                     
                                     Low = rep(0,10),
                                     
                                     Medium  = rep(0,10),
                                     
                                     High = rep(0,10), 
                                     
                                     row.names = c("stage 1","stage 2","stage 3",
                                                   "stage 4","stage 5","stage 6",
                                                   "stage 7","stage 8","stage 9",
                                                   "stage 10"))

CTAC_sense_compile <- data.frame(Null = rep(0,10),
                                 
                                 Low = rep(0,10),
                                 
                                 Medium  = rep(0,10),
                                 
                                 High = rep(0,10), 
                                 
                                 row.names = c("stage 1","stage 2","stage 3",
                                               "stage 4","stage 5","stage 6",
                                               "stage 7","stage 8","stage 9",
                                               "stage 10"))

CTAT_sense_compile <- data.frame(Null = rep(0,10),
                                 
                                 Low = rep(0,10),
                                 
                                 Medium  = rep(0,10),
                                 
                                 High = rep(0,10), 
                                 
                                 row.names = c("stage 1","stage 2","stage 3",
                                               "stage 4","stage 5","stage 6",
                                               "stage 7","stage 8","stage 9",
                                               "stage 10"))

score_occurrence_data_CTAC_antisense_compile <- data.frame(Score_CTAC_antisense_compile = rep(0,10), 
                                                          
                                                          stage = c(1:10),
                                                          
                                                          row.names = c("stage 1","stage 2","stage 3",
                                                                        "stage 4","stage 5","stage 6",
                                                                        "stage 7","stage 8","stage 9",
                                                                        "stage 10"))

score_occurrence_data_CTAT_antisense_compile <- data.frame(Score_CTAT_antisense_compile = rep(0,10), 
                                                          
                                                          stage = c(1:10),
                                                          
                                                          row.names = c("stage 1","stage 2","stage 3",
                                                                        "stage 4","stage 5","stage 6",
                                                                        "stage 7","stage 8","stage 9",
                                                                        "stage 10"))

score_occurrence_data_CTAT_sense_compile <- data.frame(Score_CTAT_sense_compile = rep(0,10), 
                                                      
                                                      stage = c(1:10),
                                                      
                                                      row.names = c("stage 1","stage 2","stage 3",
                                                                    "stage 4","stage 5","stage 6",
                                                                    "stage 7","stage 8","stage 9",
                                                                    "stage 10"))

score_occurrence_data_CTAC_sense_compile <- data.frame(Score_CTAC_sense_compile = rep(0,10), 
                                                      
                                                      stage = c(1:10),
                                                      
                                                      row.names = c("stage 1","stage 2","stage 3",
                                                                    "stage 4","stage 5","stage 6",
                                                                    "stage 7","stage 8","stage 9",
                                                                    "stage 10"))

for (i in c(1:4)) {
  for (j in c(1:10)) {
    
    CTAC_antisense_compile[[j,i]] <- mean(c(CTAC_withcaac_antisense[j,i], CTAC_withctat_antisense[j,i], 
                                         CTAC_withttac_antisense[j,i],CTAC_withctat_antisense_rep3[j,i], 
                                         CTAC_withttac_antisense_rep2[j,i]))
    
    CTAT_antisense_compile[[j,i]] <- mean(c(CTATwithCTAC_antisense[j,i], CTATwithCTAC_antisense_rep3[j,i]))
    
    CTAC_sense_compile[[j,i]] <- mean(c(CTAC_withctat_sense[j,i], CTAC_withctat_sense_rep1[j,i]))
    
    CTAT_sense_compile[[j,i]] <- mean(c(CTAT_withCTAC_sense[j,i], CTAT_withCTAC_sense_rep1[j,i]))
    
  }
}

for (i in c(1:10)){
  
  score_occurrence_data_CTAC_antisense_compile[[i,1]] <- (CTAC_antisense_compile[i,2] + 2*CTAC_antisense_compile[i,3] + 3*CTAC_antisense_compile[i,4])
  
  score_occurrence_data_CTAT_antisense_compile[[i,1]] <- (CTAT_antisense_compile[i,2] + 2*CTAT_antisense_compile[i,3] + 3*CTAT_antisense_compile[i,4])
  
  score_occurrence_data_CTAC_sense_compile[[i,1]] <- (CTAC_sense_compile[i,2] + 2*CTAC_sense_compile[i,3] + 3*CTAC_sense_compile[i,4])
  
  score_occurrence_data_CTAT_sense_compile[[i,1]] <- (CTAT_sense_compile[i,2] + 2*CTAT_sense_compile[i,3] + 3*CTAT_sense_compile[i,4])
  
}

score_occurrence_data_CTAC_antisense_compile_normalize <- score_occurrence_data_CTAC_antisense_compile

score_occurrence_data_CTAT_antisense_compile_normalize <- score_occurrence_data_CTAT_antisense_compile

score_occurrence_data_CTAC_sense_compile_normalize <- score_occurrence_data_CTAC_sense_compile

score_occurrence_data_CTAT_sense_compile_normalize <- score_occurrence_data_CTAT_sense_compile

for (i in c(1:10)) {
  
  score_occurrence_data_CTAC_antisense_compile_normalize[[i,1]] <- score_occurrence_data_CTAC_antisense_compile_normalize[i,1]/(sum(CTAC_antisense_compile[i,1:4]) +1)
  
  score_occurrence_data_CTAT_antisense_compile_normalize[[i,1]] <- score_occurrence_data_CTAT_antisense_compile_normalize[i,1]/(sum(CTAT_antisense_compile[i,1:4]) +1)
  
  score_occurrence_data_CTAC_sense_compile_normalize[[i,1]] <- score_occurrence_data_CTAC_sense_compile_normalize[i,1]/(sum(CTAC_sense_compile[i,1:4]) +1)
  
  score_occurrence_data_CTAT_sense_compile_normalize[[i,1]] <- score_occurrence_data_CTAT_sense_compile_normalize[i,1]/(sum(CTAT_sense_compile[i,1:4]) +1)
  
}


CTAC_sense_compile_mean <- data.frame(mean = rep(0,10),
                                      
                                      SD = rep(0,10),
                                 
                                      stage = c(1:10),
                                 
                                      row.names = c("stage 1","stage 2","stage 3",
                                               "stage 4","stage 5","stage 6",
                                               "stage 7","stage 8","stage 9",
                                               "stage 10"))

CTAT_sense_compile_mean <- data.frame(mean = rep(0,10),
                                      
                                      SD = rep(0,10),
                                      
                                      stage = c(1:10),
                                      
                                      row.names = c("stage 1","stage 2","stage 3",
                                                    "stage 4","stage 5","stage 6",
                                                    "stage 7","stage 8","stage 9",
                                                    "stage 10"))

CTAT_antisense_compile_mean <- data.frame(mean = rep(0,10),
                                          
                                          SD = rep(0,10),
                                      
                                          stage = c(1:10),
                                      
                                          row.names = c("stage 1","stage 2","stage 3",
                                                    "stage 4","stage 5","stage 6",
                                                    "stage 7","stage 8","stage 9",
                                                    "stage 10"))

CTAC_antisense_compile_mean <- data.frame(mean = rep(0,10),
                                          
                                          SD = rep(0,10),
                                          
                                          stage = c(1:10),
                                          
                                          row.names = c("stage 1","stage 2","stage 3",
                                                        "stage 4","stage 5","stage 6",
                                                        "stage 7","stage 8","stage 9",
                                                        "stage 10"))



for (j in c(1:10)) {
    
    CTAC_sense_compile_mean[[j,1]] <- mean(c(score_occurrence_data_normalize_combination$score_CTAC_withctat_sense[j], 
                                             score_occurrence_data_normalize_combination$score_CTAC_withctat_sense_rep1[j]))
    
    CTAC_sense_compile_mean[[j,2]] <- sd(c(score_occurrence_data_normalize_combination$score_CTAC_withctat_sense[j], 
                                             score_occurrence_data_normalize_combination$score_CTAC_withctat_sense_rep1[j]))
    
    CTAT_sense_compile_mean[[j,1]] <- mean(c(score_occurrence_data_normalize_combination$score_CTAT_withCTAC_sense[j], 
                                             score_occurrence_data_normalize_combination$score_CTAT_withCTAC_sense_rep1[j]))
    
    CTAT_sense_compile_mean[[j,2]] <- sd(c(score_occurrence_data_normalize_combination$score_CTAT_withCTAC_sense[j], 
                                             score_occurrence_data_normalize_combination$score_CTAT_withCTAC_sense_rep1[j]))
    
    CTAT_antisense_compile_mean[[j,1]] <- mean(c(score_occurrence_data_normalize_combination$score_CTATwithCTAC_antisense[j], 
                                             score_occurrence_data_normalize_combination$score_CTATwithCTAC_antisense_rep3[j]))
    
    CTAT_antisense_compile_mean[[j,2]] <- sd(c(score_occurrence_data_normalize_combination$score_CTATwithCTAC_antisense[j], 
                                                 score_occurrence_data_normalize_combination$score_CTATwithCTAC_antisense_rep3[j]))
    
    CTAC_antisense_compile_mean[[j,1]] <- mean(c(score_occurrence_data_normalize_combination$score_CTAC_withcaac_antisense[j],
                                                 score_occurrence_data_normalize_combination$score_CTAC_withttac_antisense[j],
                                                 score_occurrence_data_normalize_combination$score_CTAC_withctat_antisense[j],
                                                 score_occurrence_data_normalize_combination$score_CTAC_withctat_antisense_rep3[j],
                                                 score_occurrence_data_normalize_combination$score_CTAC_withttat_antisense_rep2[j]))
    
    CTAC_antisense_compile_mean[[j,2]] <- sd(c(score_occurrence_data_normalize_combination$score_CTAC_withcaac_antisense[j],
                                                 score_occurrence_data_normalize_combination$score_CTAC_withttac_antisense[j],
                                                 score_occurrence_data_normalize_combination$score_CTAC_withctat_antisense[j],
                                                 score_occurrence_data_normalize_combination$score_CTAC_withctat_antisense_rep3[j],
                                                 score_occurrence_data_normalize_combination$score_CTAC_withttat_antisense_rep2[j]))
    
    
}




########################
#Graphic with GGplot
######################


transparency = 0.2

line_size = 0.8


windowsFonts( A = windowsFont("Arial"))

graph <- ggplot(score_occurrence_data_normalize_combination,
                aes(x = stage)) +
  
  #CTAC lines
  
  geom_line(aes(y = score_CTAC_withcaac_antisense,
                color = "Forward AAACTAC"), size = line_size, alpha = transparency,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTAC_withttac_antisense,
                color = "Forward AAACTAC"), size = line_size,alpha = transparency,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTAC_withttat_antisense_rep2,
                color = "Forward AAACTAC"), size = line_size,alpha = transparency,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTAC_withctat_antisense,
                color = "Forward AAACTAC"), size = line_size,alpha = transparency,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTAC_withctat_antisense_rep3,
                color = "Forward AAACTAC"), size = line_size,alpha = transparency,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTAC_withctat_sense,
                color = "Reverse AAACTAC"), size = line_size, alpha = transparency,
            linetype = "dotdash") +
  
  geom_line(aes(y = score_CTAC_withctat_sense_rep1,
                color = "Reverse AAACTAC"), size = line_size, alpha = transparency,
            linetype = "dotdash") +
  
  #CTAT lines
  
  geom_line(aes(y = score_CTATwithCTAC_antisense,
                color = "Forward AAACTAT"), size = line_size, alpha = transparency,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTATwithCTAC_antisense_rep3,
                color = "Forward AAACTAT"), size = line_size, alpha = transparency,
            linetype = "solid") +
  
  geom_line(aes(y = score_CTAT_withCTAC_sense,
                color = "Reverse AAACTAT"), size = line_size, alpha = transparency,
            linetype = "dotdash") +
  
  geom_line(aes(y = score_CTAT_withCTAC_sense_rep1,
                color = "Reverse AAACTAT"), size = line_size, alpha = transparency,
            linetype = "dotdash") +
  
  #Compilation lines
  
  geom_line(data = CTAC_antisense_compile_mean,
            aes(y = mean,
                color = "Compilation foward AAACTAC"),
            linetype = "solid", size = 1.5) +
  
  geom_line(data = CTAT_antisense_compile_mean,
            aes(y = mean,
                color = "Compilation foward AAACTAT"),
            linetype = "solid", size = 1.5) +
  
  geom_line(data = CTAC_sense_compile_mean,
            aes(y = mean,
                color = "Compilation reverse AAACTAC"),
            linetype = "dotdash", size = 1.5) +
  
  geom_line(data = CTAT_sense_compile_mean,
            aes(y = mean,
                color = "Compilation reverse AAACTAT"),
            linetype = "dotdash", size = 1.5) +
  
  labs(
    x = "Egg chamber stages",
    y = "Normalized expression score",
    title = "The lncRNA expression levels of satellite DNA sequences \n across egg chamber development in Drosophila virllis",
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
    "Compilation foward AAACTAC" = "limegreen",
    "Compilation foward AAACTAT" = "magenta",
    "Compilation reverse AAACTAC" = "limegreen",
    "Compilation reverse AAACTAT" = "magenta"
  )) +
  
  theme_bw() +
  theme(
    legend.position = "right",
    legend.background = element_rect(fill = FALSE, color = FALSE), 
    legend.title = element_text(size = 12, face = "bold",family = "A"),
    legend.text = element_text(size = 9),
    text = element_text(size = 12, family = "A"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "A"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain", family = "A"),
    plot.margin = margin(20, 20, 20, 20)
  )

print(graph)


#############################################


#################################
graph_personalisable <- ggplot(score_occurrence_data_normalize_combination,
                               aes(x = stage)) +
  
  geom_line(aes(y = score_CTAC_withctat_sense,
                color = "Reverse AAACTAC"), size = 0.75,
            linetype = "dotdash") +
  
  geom_line(aes(y = score_CTAT_withCTAC_sense,
                color = "Reverse AAACTAT"), size = 0.75,
            linetype = "dotdash") +
  
  geom_line(data = score_occurrence_data_CTAC_antisense_compile_normalize,
            aes(y = Score_CTAC_antisense_compile,
                color = "Compilation AAACTAC"),
            linetype = "solid", size = 1.3) +
  
  geom_line(data = score_occurrence_data_CTAT_antisense_compile_normalize,
            aes(y = Score_CTAT_antisense_compile,
                color = "Compilation AAACTAT"),
            linetype = "solid", size = 1.3) +
  
  labs(
    x = "Egg chamber stages",
    y = "Normalized expression score",
    title = "The lncRNA expression levels of satellite DNA sequences \n across egg chamber development in Drosophila virllis",
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
    "Compilation AAACTAC" = "limegreen",
    "Compilation AAACTAT" = "magenta"
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

print(graph_personalisable)


############################################
#Data exploration
###########################################

par(mfrow=c(4,2), oma = c(0, 0, 3, 0))

for ( i in c(3:10)){
  
    plot(CTAC_withcaac_antisense[i,],
         xlab= " null = 1, low = 2, medium = 3, high = 4",
         ylab = "category occurrence",
         pch = 19,
         type = "o",
         main = paste0("Stage",i))
  
  tot <- sum(CTAC_withcaac_antisense[i,])
  
  text( x = 2.5, y = (0.75*max(CTAC_withcaac_antisense[i,])), label = paste0("N = ", tot))
  
}

mtext("CTAC_withcaac_antisense number of observations", outer = TRUE, cex = 1.5, font = 2)

par(mfrow=c(4,2), oma = c(0, 0, 3, 0))

for ( i in c(3:10)){
  
  plot(CTAC_withctat_antisense_rep3[i,],
       xlab= " null = 1, low = 2, medium = 3, high = 4",
       ylab = "category occurrence",
       pch = 19,
       type = "o",
       main = paste0("Stage",i))
  
  tot <- sum(CTAC_withctat_antisense_rep3[i,])
  
  text( x = 2.5, y = (0.75*max(CTAC_withctat_antisense_rep3[i,])), label = paste0("N = ", tot))
  
}

mtext("CTAC_withctat_antisense_rep3 number of observations", outer = TRUE, cex = 1.5, font = 2)

par(mfrow=c(4,2), oma = c(0, 0, 3, 0))

for ( i in c(3:10)){
  
  plot(CTAC_withctat_antisense[i,],
       xlab= " null = 1, low = 2, medium = 3, high = 4",
       ylab = "category occurrence",
       pch = 19,
       type = "o",
       main = paste0("Stage",i))
  
  tot <- sum(CTAC_withctat_antisense[i,])
  
  text( x = 2.5, y = (0.75*max(CTAC_withctat_antisense[i,])), label = paste0("N = ", tot))
  
}

mtext("CTAC_withctat_antisense number of observations", outer = TRUE, cex = 1.5, font = 2)

#GLM distribution Gamma test

#creation of data_frame mutate of sense and antisense CTAC and CTAT

score_CTAC_antisense <- score_occurrence_data_normalize_combination[,c(1,2,3,10,11,15)] %>%
  pivot_longer(
    cols = !stage,
    names_to = "strand",
    values_to = "score"
  )

score_CTAC_sense <- score_occurrence_data_normalize_combination[,c(7,13,15)] %>%
  pivot_longer(
    cols = !stage,
    names_to = "strand",
    values_to = "score"
  )

score_CTAT_antisense <- score_occurrence_data_normalize_combination[,c(5,9,15)] %>%
  pivot_longer(
    cols = !stage,
    names_to = "strand",
    values_to = "score"
  )

score_CTAT_sense <- score_occurrence_data_normalize_combination[,c(8,14,15)] %>%
  pivot_longer(
    cols = !stage,
    names_to = "strand",
    values_to = "score"
  )

hist(score_CTAC_antisense$score)
hist(score_CTAC_sense$score)
hist(score_CTAT_antisense$score)
hist(score_CTAT_sense$score)

boxplot(score_CTAC_antisense$score ~ score_CTAC_antisense$stage)
boxplot(score_CTAC_sense$score ~ score_CTAC_sense$stage)
boxplot(score_CTAT_antisense$score ~ score_CTAT_antisense$stage)
boxplot(score_CTAT_sense$score ~ score_CTAT_sense$stage)

concatenated <- qpcR:::cbind.na(score_CTAC_antisense[,3],score_CTAC_sense[,3],
                              score_CTAT_antisense[,3],score_CTAT_sense[,3])

pairs.panels(concatenated)



model_CTAC_antisense <- lm(data = score_CTAC_antisense, score ~ stage)
model_CTAC_sense <- lm(data = score_CTAC_sense, score ~ stage)
model_CTAT_antisense <- lm(data = score_CTAT_antisense, score ~ stage)
model_CTAT_sense <- lm(data = score_CTAT_sense, score ~ stage)


gvlma(model_CTAC_antisense)
check_model(model_CTAC_antisense)

gvlma(model_CTAC_sense)
check_model(model_CTAC_sense)


gvlma(model_CTAT_antisense)
check_model(model_CTAT_antisense)


gvlma(model_CTAT_sense)
check_model(model_CTAT_sense)


#verification of gam model as well as their application conditions
#####################################
#Curve comparison test

score_CTAC_antisense$strand <- "CTAC_antisense" 

score_CTAC_sense$strand <- "CTAC_sense" 

score_CTAT_antisense$strand <- "CTAT_antisense" 

score_CTAT_sense$strand <- "CTAT_sense" 

combined_data <- rbind(score_CTAC_antisense, score_CTAC_sense,score_CTAT_antisense,score_CTAT_sense)

#data export

combined_data$strand <- as.factor(combined_data$strand)

combined_model_tw <- gam(score ~ strand + s(stage, by = as.factor(strand), bs = "cr", k = 10), 
                      data = combined_data, 
                      method = "REML", family = tw(link = "log"),control = list(
                        maxit = 1000))

combined_model <- gam(list(score ~ strand + s(stage, by = strand, bs = "cr", k = 10), 
                           ~ s(stage)), 
                      data = combined_data,
                      optimizer = c("outer", "newton"),
                      method = "REML", 
                      family = gaulss(),
                      control = list(maxit = 1000)
)

summary(combined_model)
par(mfrow=c(2,2))
gam.check(combined_model)
k.check(combined_model)

concurvity(combined_model, full = FALSE)


model_performance(combined_model)
plot(combined_model)

anova <- anova.gam(combined_model)


stage_dense <- seq(min(combined_data$stage), 
                   max(combined_data$stage), 
                   length.out = 200)

df_predict_combine <- expand.grid(
  stage = stage_dense,
  strand = unique(combined_data$strand) 
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

# calculate R^2 by strand

combined_data$preds <- predict(combined_model, type = "response")[,1]

table_R2 <- combined_data %>%
  group_by(strand) %>%
  summarize(
    RSS = sum((score - preds)^2),                  
    TSS = sum((score - mean(score))^2),           
    R2  = round(1 - (RSS / TSS), 4)               
  )

print(table_R2)

legend_names <- c(
  "CTAC_antisense" = "Forward AAACTAC",
  "CTAC_sense"     = "Reverse AAACTAC",
  "CTAT_antisense" = "Foward AAACTAT",
  "CTAT_sense"     = "Reverse AAACTAT"
)

p_val_strand <- anova$pTerms.pv[1]

#####################################################

ploT_GAM <- ggplot() + 
  
  geom_point(data = combined_data, 
             aes(x = stage, y = score, colour = strand),
             position = position_jitter(width = 0.05)) + 
  
  geom_ribbon(data = df_predict_combine, 
              aes(x = stage, ymin = lower, ymax = upper, fill = strand), 
              alpha = 0.09) +
  
  geom_line(data = df_predict_combine, 
            aes(x = stage, y = fit, colour = strand, linetype = strand),
            linewidth = 1.2) +
  
  labs(
    x = "Ovariole developmental stages",
    y = "Normalize expression score",
    title = "Transcription evolution of forward and reverse LncRNA of AAACTAC and AAACTAT for D. virilis",
    subtitle = "Used probe were forward and reverse",
    colour = "LncRNA strand",
    linetype = "LncRNA strand" 
  ) +
  scale_color_manual(labels = legend_names,
                     values = c("CTAC_antisense" = "limegreen", 
                                "CTAC_sense"     = "limegreen", 
                                "CTAT_antisense" = "magenta", 
                                "CTAT_sense"     = "magenta")) +
  
  scale_fill_manual(labels = legend_names,
                    values = c("CTAC_antisense" = "limegreen", 
                               "CTAC_sense"     = "limegreen", 
                               "CTAT_antisense" = "magenta", 
                               "CTAT_sense"     = "magenta")) +
  
  scale_linetype_manual(labels = legend_names,
                        values = c("CTAC_antisense" = "solid", 
                                   "CTAC_sense"     = "dotted", 
                                   "CTAT_antisense" = "solid", 
                                   "CTAT_sense"     = "dotted")) +
  
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  ylim(c(-0.5, 5)) + 
  
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
  
  annotate("rect", xmin = 6.5, xmax = Inf, ymin = 4.2, ymax = 5,
           alpha = 1, fill = "white") +
  
  #annotate("text", x = Inf, y = Inf, label = paste0("Global ANOVA: p = ", format.pval(p_val_strand, digits = 2)), 
           #hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
  
  annotate("text", x = 8.5, y = 4.9, label = paste0("R² S.CTAC = ", round(table_R2[2,4], digits = 3), "***"), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic",
           color = "limegreen") + 
  
  annotate("text", x = 8.5, y = 4.6, label = paste0("R² A.CTAC = ", round(table_R2[1,4], digits = 3), "***" ), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic",
           color = "limegreen") +
  
  annotate("text", x = Inf, y = 4.9, label = paste0("R² S.CTAT = ", round(table_R2[4,4], digits = 3), "***"), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic",
           color = "magenta") + 
  
  annotate("text", x = Inf, y = 4.6, label = paste0("R² A.CTAT = ", round(table_R2[3,4], digits = 3), "***" ), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic",
           color = "magenta")
  

print(ploT_GAM)

combined_data_CTAC <- rbind(combined_data[combined_data$strand == "CTAC_antisense", ],
                              combined_data[combined_data$strand == "CTAC_sense", ])

df_predict_combine_CTAC <- rbind(df_predict_combine[df_predict_combine$strand == "CTAC_antisense", ],
                                  df_predict_combine[df_predict_combine$strand == "CTAC_sense", ])


plot_GAM_CTAC_only <- ggplot() + 
  
  geom_point(data = combined_data_CTAC, 
             aes(x = stage, y = score, colour = strand),
             position = position_jitter(width = 0.05)) + 
  
  geom_ribbon(data = df_predict_combine_CTAC, 
              aes(x = stage, ymin = lower, ymax = upper, fill = strand), 
              alpha = 0.09) +
  
  geom_line(data = df_predict_combine_CTAC, 
            aes(x = stage, y = fit, colour = strand, linetype = strand),
            linewidth = 1.2) +
  
  labs(
    x = "Ovariole developmental stages",
    y = "Normalize expression score",
    title = "Transcription evolution of forward and reverse LncRNA of AAACTAC and AAACTAT for D. virilis",
    subtitle = "Used probe were forward and reverse",
    colour = "LncRNA strand",
    linetype = "LncRNA strand" 
  ) +
  scale_color_manual(labels = legend_names,
                     values = c("CTAC_antisense" = "darkgreen", 
                                "CTAC_sense"     = "limegreen")) +
  
  scale_fill_manual(labels = legend_names,
                    values = c("CTAC_antisense" = "darkgreen", 
                               "CTAC_sense"     = "limegreen")) +
  
  scale_linetype_manual(labels = legend_names,
                        values = c("CTAC_antisense" = "solid", 
                                   "CTAC_sense"     = "dotted")) +
  
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  ylim(c(-0.5, 5)) + 
  
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.background = element_rect(fill = "transparent", color = NA), 
    legend.title = element_text(size = 12, face = "bold", family = "A"),
    legend.text = element_text(size = 12),
    text = element_text(size = 12, family = "A"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "A"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "plain", family = "A"),
    plot.margin = margin(20, 20, 20, 20),
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold")
  ) + 
  
  guides(fill = "none") +
  
  annotate("rect", xmin = 8.5, xmax = Inf, ymin = 4.2, ymax = 5,
           alpha = 1, fill = "white") +
  
  #annotate("text", x = Inf, y = Inf, label = paste0("Global ANOVA: p = ", format.pval(p_val_strand, digits = 2)), 
  #hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
  
  annotate("text", x = Inf, y = 4.9, label = paste0("R² F.CTAC = ", round(table_R2[2,4], digits = 3), "***"), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic",
           color = "darkgreen") + 
  
  annotate("text", x = Inf, y = 4.6, label = paste0("R² R.CTAC = ", round(table_R2[1,4], digits = 3), "***" ), 
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic",
           color = "limegreen") 

print(plot_GAM_CTAC_only)

############################################################


table_pvalues <- data.frame()

for (s in 3:10) {
  stage_data <- subset(combined_data, stage == s)
  
  if(length(unique(stage_data$strand)) > 1) {
    
    fit <- aov(score ~ strand, data = stage_data)
    
    tukey <- TukeyHSD(fit)$strand
    
    temp <- data.frame(
      Stage = s,
      Comparison = rownames(tukey),
      p_value = round(tukey[, "p adj"], 4)
    )
    table_pvalues <- rbind(table_pvalues, temp)
  }
}

print(table_pvalues)

signif_data <- table_pvalues %>%
  filter(p_value < 0.05) %>%
  mutate(
    label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )
