
library("ggplot2")

test <- read.csv("Values.csv",sep=",", dec=".", header=TRUE)
test2 <-read.csv("Values_2.csv",sep=",", dec=".", header=TRUE)
test3 <-read.csv("Values_3.csv",sep=",", dec=".", header=TRUE)
test4 <-read.csv("Values_4.csv",sep=",", dec=".", header=TRUE)

plot(test$C2.D.vir48_CTAC_CTAT_antisens_29052025.lif_._Image001 ~ test$C3.D.vir48_CTAC_CTAT_antisens_29052025.lif_._Image001,
     type = "p",
     pch = 19,
     col = "orange")

plot(test2$C2.D.vir48_CTAC_CTAT_antisens_29052025.lif_._Image004 ~ test2$C3.D.vir48_CTAC_CTAT_antisens_29052025.lif_._Image004,
     type = "p",
     pch = 19,
     col = "blue")


plot(test3$C2.D.vir48_CTAC_CTAT_antisens_29052025.lif_._Image005 ~ test3$C3.D.vir48_CTAC_CTAT_antisens_29052025.lif_._Image005,
     type = "p",
     pch = 19,
     col = "black")

plot(test4$C2.D.vir48_CTAC_CTAT_antisens_29052025.lif_._Image007 ~ test4$C3.D.vir48_CTAC_CTAT_antisens_29052025.lif_._Image007,
     type = "p",
     pch = 19,
     col = "black")

cor.test(test$C2.D.vir48_CTAC_CTAT_antisens_29052025.lif_._Image001, test$C3.D.vir48_CTAC_CTAT_antisens_29052025.lif_._Image001, method = "pearson")




cor.test(test4$C2.D.vir48_CTAC_CTAT_antisens_29052025.lif_._Image007, test4$C3.D.vir48_CTAC_CTAT_antisens_29052025.lif_._Image007, method = "pearson")

############################
#Graph d'essaie
############################

test4 <-read.csv("Values_4.csv",sep=",", dec=".", header=TRUE)
test_icq_blue <-read.csv("Values_icq_bluechanel.csv",sep=",", dec=".", header=TRUE)


test_icq_green <-read.csv("Values_icq_greenchanel.csv",sep=",", dec=".", header=TRUE)
test_van_steensel <-read.csv("Values_van_steensel_CCF.csv",sep=",", dec=".", header=TRUE)
test_coste_threshold <-read.csv("Values_coste_threshold_test.csv",sep=",", dec=".", header=TRUE)

plot(test_icq_green$C2.D.vir48_CTAC_CTAT_antisens_29052025.lif_._Image007 ~ test_icq_green$X.Ai.a..Bi.b.,
     type = "p",
     pch = 19,
     col = "black")

plot(test_icq_blue$C3.D.vir48_CTAC_CTAT_antisens_29052025.lif_._Image007 ~ test_icq_blue$X.Ai.a..Bi.b.,
     type = "p",
     pch = 19,
     col = "black")

plot(test_van_steensel$Y1 ~ test_van_steensel$dx,
     type = "l",
     lwd = 3,
     col = "black")

lines(test_van_steensel$Y3 ~test_van_steensel$dx,
         type = "l",
         lwd = 3,
         col = "blue",)

plot(test_coste_threshold$Pearson.s_coefficient_below ~ test_coste_threshold$ThrA,
     type = "p",
     lwd = 3,
     col = "black")














