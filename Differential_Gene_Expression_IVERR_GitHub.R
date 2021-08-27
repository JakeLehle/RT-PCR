setwd("C:/Users/Documents/UTSA/Lab/RT-PCR/PGC-LCs_Validation")
# Clear the global enviornment and check that you are using the most current version of R
rm(list = ls())
install.packages("installr")
library(installr)
updateR()
install.packages("reshape")
library(reshape)
install.packages("dplyr")
library(dplyr)
install.packages("ggplot2")
library(ggplot2)

library(dplyr)
Sertoli = c(0.029501909, 0.047830185, 0.039009393, 0.071400955, 0.064282306, 0.075950148, 0.022554098, 0.020139595, 0.033904228, 0.00645475, 0.00286432, 0.007288827, 28.78589252, 29.4473892, 27.54768639, 15.54947149, 16.64525723, 16.66174808, 0.08750307, 0.095257905, 0.081048079, 25.09662907, 24.96181882, 25.3741064, 0.003699351, 0.004785226, 0.004875832, 0.000113587, 2.15691E-05, 2.38527E-06)
Granulosa = c(0.011541335, 0.011928446, 0.011335075, 0.019347775, 0.029230452, 0.020836524, 1.128401248, 1.151289683, 1.108368412, 0.008196422, 0.007595136, 0.008084518, 0.071661493, 0.064669131, 0.069147486, 1.255868805, 1.306372848, 1.287600515, 0.01008212, 0.008428565, 0.007971158, 59.94157775, 61.3509779, 61.57704277, 1.28439E-05, 0.00049788, 0.000260008, 1.0164E-08, 6.32658E-06, 3.12229E-06)
ESC = c(0.004039193, 0.005385866, 0.002620471, 0.010002751, 0.005621464, 0.002438519, 0.02903798, 0.014772911, 0.022671169, 244.0631004, 243.5090617, 246.7824296, 0.025720165, 0.019658368, 0.02769145, 0.032944194, 0.025857567, 0.03120863, NA, NA, NA, 0.073008819, 0.076812625, 0.05660326, 757.0592424, 753.2494727, 417.2138525, 11.36168908, 11.36316123, 10.68206063)
iPSC = c(0.000164339, 0.001685181, 0.003749615, NA, 7.6432E-06, NA, NA, 0.009915155, 0.003140781, 103.7990192, 110.1199677, 109.0497534, 0.019122996, 0.020288182, 0.010471232, 0.150790604, 0.188399553, 0.165774859, NA, NA, NA, 0.07633671, 0.09000255, 0.084892901, 273.0146338, 274.9154617, 259.9116136, 3.739707647, 3.732797149, 3.948242791)
data = rbind(Sertoli,Granulosa,ESC,iPSC);
colnames(data) = c( "Ar", "Ar", "Ar", "Esr1", "Esr1", "Esr1", "Esr2", "Esr2", "Esr2", "Nanog", "Nanog", "Nanog", "Sox9", "Sox9", "Sox9", "Wt1", "Wt1", "Wt1", "Fshr", "Fshr", "Fshr", "Inha", "Inha", "Inha", "Pou5f1", "Pou5f1", "Pou5f1", "Sox2", "Sox2", "Sox2");
data
data_T <- t(data);
data_T
mdata <- melt(data_T, id=c("Sertoli", "Granulosa", "ESC", "iPSC"));
mdata
colnames(mdata)=c("Gene", "Cell", 'Fold_Expression');
mdata
#Lets sort the data by the gene name
mnewdata <- mdata[order(mdata$Gene),]
mnewdata
#You can see the row numbers are all jumbled now though so we need to reset them to the default
rownames(mnewdata) <- NULL
mnewdata

#Now you can just subset and pull out and plot whichever gene you want.
AR <- mnewdata[1:12,]
AR
AR_mean <- AR %>% group_by(Cell) %>% summarize(average = mean(Fold_Expression), sd = sd(Fold_Expression))
AR_mean

AR_p<- ggplot(data=AR_mean, aes(x=Cell, y=average, fill=Cell)) +
  theme(plot.title = element_text(hjust = 0.5, size= 30, face="bold"), text = element_text(size = 20), legend.position = "none")+
  labs(title="Ar", x="Cell Type", y=expression("Fold Expression")) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=average-sd, ymax=average+sd), width=.2,
                position=position_dodge(.9)) +
  geom_hline(yintercept = 0.007, linetype='dotted', col = 'red', size= 1.1 ) +
  scale_fill_brewer(palette="Set3")
AR_p
#///////////////

Esr1 <- mnewdata[13:24,]
Esr1
Esr1_mean <- Esr1 %>% group_by(Cell) %>% summarize(average = mean(Fold_Expression), sd = sd(Fold_Expression))
Esr1_mean[4,2] <- c(0.0000076432)


Esr1_p<- ggplot(data=Esr1_mean, aes(x=Cell, y=average, fill=Cell)) +
  theme(plot.title = element_text(hjust = 0.5, size= 30, face="bold"), text = element_text(size = 20), legend.position = "none")+
  labs(title="Esr1", x="Cell Type", y=expression("Fold Expression")) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=average-sd, ymax=average+sd), width=.2,
                position=position_dodge(.9)) +
  geom_hline(yintercept = 0.007, linetype='dotted', col = 'red', size= 1.1 ) +
  scale_fill_brewer(palette="Set3")
Esr1_p

#/////////////////////

Esr2 <- mnewdata[25:36,]
Esr2
Esr2_mean <- Esr2 %>% group_by(Cell) %>% summarize(average = mean(Fold_Expression), sd = sd(Fold_Expression))
Esr2_mean[4,2] <- c((0.009915155+0.003140781)/2) 


Esr2_p<- ggplot(data=Esr2_mean, aes(x=Cell, y=average, fill=Cell)) +
  theme(plot.title = element_text(hjust = 0.5, size= 30, face="bold"), text = element_text(size = 20), legend.position = "none")+
  labs(title="Esr2", x="Cell Type", y=expression("Fold Expression")) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=average-sd, ymax=average+sd), width=.2,
                position=position_dodge(.9)) +
  geom_hline(yintercept = 0.007, linetype='dotted', col = 'red', size= 1.1 ) +
  scale_fill_brewer(palette="Set3")
Esr2_p

install.packages("gg.gap")
library(gg.gap)

gg.gap(plot=Esr2_p,
       segments=c(0.04,0.9),
       tick_width = c(0.025, 0.1),
       margin = c(top = 0.1, right = 0.1, bottom = 0.1, left = 0.7),
       ylim=c(0,1.2))


#///////////////////////

NANOG <- mnewdata[37:48,]
NANOG
NANOG_mean <- NANOG %>% group_by(Cell) %>% summarize(average = mean(Fold_Expression), sd = sd(Fold_Expression))
NANOG_mean


NANOG_p<- ggplot(data=NANOG_mean, aes(x=Cell, y=average, fill=Cell)) +
  theme(plot.title = element_text(hjust = 0.5, size= 30, face="bold"), text = element_text(size = 20), legend.position = "none")+
  labs(title="Nanog", x="Cell Type", y=expression("Fold Expression")) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=average-sd, ymax=average+sd), width=.2,
                position=position_dodge(.9)) +
  geom_hline(yintercept = 0.007, linetype='dotted', col = 'red', size= 1.1 ) +
  scale_fill_brewer(palette="Set3")
NANOG_p

#///////////////////////

Sox9 <- mnewdata[49:60,]
Sox9
Sox9_mean <- Sox9 %>% group_by(Cell) %>% summarize(average = mean(Fold_Expression), sd = sd(Fold_Expression))
Sox9_mean


Sox9_p<- ggplot(data=Sox9_mean, aes(x=Cell, y=average, fill=Cell)) +
  theme(plot.title = element_text(hjust = 0.5, size= 30, face="bold"), text = element_text(size = 20), legend.position = "none")+
  labs(title="Sox9", x="Cell Type", y=expression("Fold Expression")) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=average-sd, ymax=average+sd), width=.2,
                position=position_dodge(.9)) +
  geom_hline(yintercept = 0.007, linetype='dotted', col = 'red', size= 1.1 ) +
  scale_fill_brewer(palette="Set3")
Sox9_p

#/////////////////////////

Wt1 <- mnewdata[61:72,]
Wt1
Wt1_mean <- Wt1 %>% group_by(Cell) %>% summarize(average = mean(Fold_Expression), sd = sd(Fold_Expression))
Wt1_mean


Wt1_p<- ggplot(data=Wt1_mean, aes(x=Cell, y=average, fill=Cell)) +
  theme(plot.title = element_text(hjust = 0.5, size= 30, face="bold"), text = element_text(size = 20), legend.position = "none")+
  labs(title="Wt1", x="Cell Type", y=expression("Fold Expression")) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=average-sd, ymax=average+sd), width=.2,
                position=position_dodge(.9)) +
  geom_hline(yintercept = 0.007, linetype='dotted', col = 'red', size= 1.1 ) +
  scale_fill_brewer(palette="Set3")
Wt1_p

#////////////////////////

Fshr <- mnewdata[73:84,]
Fshr
Fshr_mean <- Fshr %>% group_by(Cell) %>% summarize(average = mean(Fold_Expression), sd = sd(Fold_Expression))
Fshr_mean


Fshr_p<- ggplot(data=Fshr_mean, aes(x=Cell, y=average, fill=Cell)) +
  theme(plot.title = element_text(hjust = 0.5, size= 30, face="bold"), text = element_text(size = 20), legend.position = "none")+
  labs(title="Fshr", x="Cell Type", y=expression("Fold Expression")) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=average-sd, ymax=average+sd), width=.2,
                position=position_dodge(.9)) +
  geom_hline(yintercept = 0.007, linetype='dotted', col = 'red', size= 1.1 ) +
  scale_fill_brewer(palette="Set3")
Fshr_p

#/////////////////////////

Inha <- mnewdata[85:96,]
Inha
Inha_mean <- Inha %>% group_by(Cell) %>% summarize(average = mean(Fold_Expression), sd = sd(Fold_Expression))
Inha_mean


Inha_p<- ggplot(data=Inha_mean, aes(x=Cell, y=average, fill=Cell)) +
  theme(plot.title = element_text(hjust = 0.5, size= 30, face="bold"), text = element_text(size = 20), legend.position = "none")+
  labs(title="Inha", x="Cell Type", y=expression("Fold Expression")) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=average-sd, ymax=average+sd), width=.2,
                position=position_dodge(.9)) +
  geom_hline(yintercept = 0.007, linetype='dotted', col = 'red', size= 1.1 ) +
  scale_fill_brewer(palette="Set3")
Inha_p

#/////////////////////////

Pou5f1 <- mnewdata[97:108,]
Pou5f1
Pou5f1_mean <- Pou5f1 %>% group_by(Cell) %>% summarize(average = mean(Fold_Expression), sd = sd(Fold_Expression))
Pou5f1_mean


Pou5f1_p<- ggplot(data=Pou5f1_mean, aes(x=Cell, y=average, fill=Cell)) +
  theme(plot.title = element_text(hjust = 0.5, size= 30, face="bold"), text = element_text(size = 20), legend.position = "none")+
  labs(title="Pou5f1", x="Cell Type", y=expression("Fold Expression")) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=average-sd, ymax=average+sd), width=.2,
                position=position_dodge(.9)) +
  geom_hline(yintercept = 0.007, linetype='dotted', col = 'red', size= 1.1 ) +
  scale_fill_brewer(palette="Set3")
Pou5f1_p

#/////////////////////////

Sox2 <- mnewdata[109:120,]
Sox2
Sox2_mean <- Sox2 %>% group_by(Cell) %>% summarize(average = mean(Fold_Expression), sd = sd(Fold_Expression))
Sox2_mean


Sox2_p<- ggplot(data=Sox2_mean, aes(x=Cell, y=average, fill=Cell)) +
  theme(plot.title = element_text(hjust = 0.5, size= 30, face="bold"), text = element_text(size = 20), legend.position = "none")+
  labs(title="Sox2", x="Cell Type", y=expression("Fold Expression")) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=average-sd, ymax=average+sd), width=.2,
                position=position_dodge(.9)) +
  geom_hline(yintercept = 0.007, linetype='dotted', col = 'red', size= 1.1 ) +
  scale_fill_brewer(palette="Set3")
Sox2_p