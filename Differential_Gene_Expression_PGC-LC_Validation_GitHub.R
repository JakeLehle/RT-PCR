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

#////////////////////////////////////////////////////////////////////////////

iPSC = c(302.653075136779,	318.940557809098,	359.788737290081,	0.013122508871,	0.055199105607,	0.088972876418,	371.898785195213,	521.370910249840,	446.817087022375,	186.325413858015,	175.358074816729,	172.820846554677,	NA,	0.028507933983,	0.045632782034,	1.285810370790,	1.110299461968,	1.245347758373,	NA,	NA,	NA,	38.600591390428,	40.251751743720,	38.631172152179,	0.780147253325,	1.058635611129,	1.247707695185,	29.172062348138,	33.555092492165,	33.458168049654,	29.504193219549,	30.126104927702,	29.517887750882,	3.104868769397,	3.472388379535,	3.863224659555,	NA,	NA,	0.027446289807,	5.708132904482,	4.509886896732,	5.604823353928,	0.630999116243,	0.303820542786,	0.612582366367)
Epi_LC = c(215.285705232529,	203.905211273123,	197.104909388115,	0.139729396443,	0.181676531070,	0.129399664057,	12.084653625228,	12.916087902899,	12.817137697846,	16.794686233980,	17.363027560032,	15.626175222018,	3.321182686386,	3.118454033778,	3.009005626185,	4.508243560951,	4.008549717904,	3.927958475756,	0.036811770738,	0.039651393976,	NA,	4.705762161531,	4.391704193409,	4.262087949159,	0.090180258746,	0.091206620704,	0.249478655910,	0.648257775711,	0.772404275295,	0.498142640145,	0.036051842720,	0.036712816215,	NA,	0.199527794688,	0.091565452797,	0.211192976011,	NA,	NA,	NA,	3.776523972226,	3.868245153324,	4.162640396560,	0.356912730192,	0.379716585088,	0.282971390526)
PGC_LC = c(232.987855592034,	241.085110592765,	211.384292241974,	3.098272428167,	2.993008496858,	2.731197158510,	276.066027705101,	279.345718448808,	276.653173308305,	149.312273461751,	142.230420891442,	141.547790197161,	0.551090009533,	0.390305935385,	0.510632532736,	0.106488549802,	0.034326114271,	0.131484241546,	NA,	NA,	NA,	94.604009372577,	93.678272260391,	87.214205828088,	4.859665743779,	5.063736492536,	4.995685048874,	144.509660005648,	147.216306467003,	150.824874905132,	43.133318782645,	42.983886528673,	41.097183946650,	29.562461690663,	29.462264825630,	29.535701446472,	0.580093506876,	0.648054688027,	0.746176606095,	42.479453876410,	43.292942412150,	26.860879086843,	1.043201601960,	1.037023273759,	0.943322499084)
data = rbind(iPSC,Epi_LC,PGC_LC);

colnames(data) = c("Pou5f1",	"Pou5f1",	"Pou5f1",	"Wnt3",	"Wnt3",	"Wnt3",	"Nanog",	"Nanog",	"Nanog",	"Sox2",	"Sox2",	"Sox2",	"Fgf5",	"Fgf5",	"Fgf5",	"Dnmt3b",	"Dnmt3b",	"Dnmt3b",	"Sox17",	"Sox17",	"Sox17",	"Dppa3",	"Dppa3",	"Dppa3",	"Gata4",	"Gata4",	"Gata4",	"Tfap2c",	"Tfap2c",	"Tfap2c",	"Prdm14",	"Prdm14",	"Prdm14",	"Nanos3",	"Nanos3",	"Nanos3",	"Itgb3",	"Itgb3",	"Itgb3",	"Dnd1",	"Dnd1",	"Dnd1",	"Dazl",	"Dazl",	"Dazl")

data
data_T <- t(data)
data_T
mdata <- melt(data_T, id=c("iPSC", "Epi-LC", "PGC-LC"));
mdata
colnames(mdata)=c("Gene", "Cell", 'Fold_Expression');
mdata
#Lets sort the data by the gene name
mnewdata <- mdata[order(mdata$Gene),]
mnewdata
#You can see the row numbers are all jumbled now though so we need to reset them to the default
rownames(mnewdata) <- NULL
mnewdata

Pou5f1 <- mnewdata[1:9,]

df_mean_pou5f1 <- Pou5f1 %>% 
  group_by(Cell) %>% 
  summarize(average = mean(Fold_Expression))

df_mean_pou5f1
# Pou5f1 at stock concentration

Pou5f1_p <- ggplot(Pou5f1, aes(x=Cell, y=Fold_Expression, color=Cell)) + 
  theme(plot.title = element_text(hjust = 0.5, size= 20, face="bold"), text = element_text(size = 15), legend.position = "none") +
  labs(title="Pou5f1", x="Cell Type", y=expression("Fold Expression Change")) +
  ylim(0,400)+
  geom_boxplot(width=0.5) +
  geom_point(data = df_mean_pou5f1,
             mapping = aes(x = Cell, y = average),
             color="black") +
  geom_line(data = df_mean_pou5f1, 
            mapping = aes(x = Cell, y = average, group=1), 
            color= "black")

Pou5f1_p

#////////////////////////////////////////////////////////

Wnt3 <- mnewdata[10:18,]

df_mean_wnt3 <- Wnt3 %>% 
  group_by(Cell) %>% 
  summarize(average = mean(Fold_Expression))

df_mean_wnt3
# Wnt3 at stock concentration

Wnt3_p <- ggplot(Wnt3, aes(x=Cell, y=Fold_Expression, color=Cell)) + 
  theme(plot.title = element_text(hjust = 0.5, size= 20, face="bold"), text = element_text(size = 15), legend.position = "none") +
  labs(title="Wnt3", x="Cell Type", y=expression("Fold Expression Change")) +
  geom_boxplot(width=0.5) +
  geom_point(data = df_mean_wnt3,
             mapping = aes(x = Cell, y = average),
             color="black") +
  geom_line(data = df_mean_wnt3, 
            mapping = aes(x = Cell, y = average, group=1), 
            color= "black")

Wnt3_p

#////////////////////////////////////////////////////////

Nanog <- mnewdata[19:27,]

df_mean_nanog <- Nanog %>% 
  group_by(Cell) %>% 
  summarize(average = mean(Fold_Expression))

df_mean_nanog
# Nanog at stock concentration

Nanog_p <- ggplot(Nanog, aes(x=Cell, y=Fold_Expression, color=Cell)) + 
  theme(plot.title = element_text(hjust = 0.5, size= 20, face="bold"), text = element_text(size = 15), legend.position = "none") +
  labs(title="Nanog", x="Cell Type", y=expression("Fold Expression Change")) +
  geom_boxplot(width=0.5) +
  geom_point(data = df_mean_nanog,
             mapping = aes(x = Cell, y = average),
             color="black") +
  geom_line(data = df_mean_nanog, 
            mapping = aes(x = Cell, y = average, group=1), 
            color= "black")

Nanog_p

#///////////////////////////////////////////////////////

Sox2 <- mnewdata[28:36,]

df_mean_sox2 <- Sox2 %>% 
  group_by(Cell) %>% 
  summarize(average = mean(Fold_Expression))

df_mean_sox2
# Sox2 at stock concentration

Sox2_p <- ggplot(Sox2, aes(x=Cell, y=Fold_Expression, color=Cell)) + 
  theme(plot.title = element_text(hjust = 0.5, size= 20, face="bold"), text = element_text(size = 15), legend.position = "none") +
  labs(title="Sox2", x="Cell Type", y=expression("Fold Expression Change")) +
  geom_boxplot(width=0.5) +
  geom_point(data = df_mean_sox2,
             mapping = aes(x = Cell, y = average),
             color="black") +
  geom_line(data = df_mean_sox2, 
            mapping = aes(x = Cell, y = average, group=1), 
            color= "black")

Sox2_p

#/////////////////////////////////////////////////////////

Fgf5 <- mnewdata[37:45,]

df_mean_fgf5 <- Fgf5 %>% 
  group_by(Cell) %>% 
  summarize(average = mean(Fold_Expression))

df_mean_fgf5[1,2] <- c(0.037070358)
df_mean_fgf5
# Fgf5 at stock concentration

Fgf5_p <- ggplot(Fgf5, aes(x=Cell, y=Fold_Expression, color=Cell)) + 
  theme(plot.title = element_text(hjust = 0.5, size= 20, face="bold"), text = element_text(size = 15), legend.position = "none") +
  labs(title="Fgf5", x="Cell Type", y=expression("Fold Expression Change")) +
  geom_boxplot(width=0.5) +
  geom_point(data = df_mean_fgf5,
             mapping = aes(x = Cell, y = average),
             color="black") +
  geom_line(data = df_mean_fgf5, 
            mapping = aes(x = Cell, y = average, group=1), 
            color= "black")

Fgf5_p

#/////////////////////////////////////////////////////////

Dnmt3b <- mnewdata[46:54,]

df_mean_dnmt3b <- Dnmt3b %>% 
  group_by(Cell) %>% 
  summarize(average = mean(Fold_Expression))

df_mean_dnmt3b
# Dnmt3b at stock concentration

Dnmt3b_p <- ggplot(Dnmt3b, aes(x=Cell, y=Fold_Expression, color=Cell)) + 
  theme(plot.title = element_text(hjust = 0.5, size= 20, face="bold"), text = element_text(size = 15), legend.position = "none") +
  labs(title="Dnmt3b", x="Cell Type", y=expression("Fold Expression Change")) +
  geom_boxplot(width=0.5) +
  geom_point(data = df_mean_dnmt3b,
             mapping = aes(x = Cell, y = average),
             color="black") +
  geom_line(data = df_mean_dnmt3b, 
            mapping = aes(x = Cell, y = average, group=1), 
            color= "black")

Dnmt3b_p

#/////////////////////////////////////////////////////////

Sox17 <- mnewdata[55:63,]
Sox17
df_mean_sox17 <- Sox17 %>% 
  group_by(Cell) %>% 
  summarize(average = mean(Fold_Expression))
df_mean_sox17[2,2] <- c(0.03823158)
df_mean_sox17
# Sox17 at stock concentration

Sox17_p <- ggplot(Sox17, aes(x=Cell, y=Fold_Expression, color = "green")) +
  theme(plot.title = element_text(hjust = 0.5, size= 20, face="bold"), text = element_text(size = 15), legend.position = "none") +
  labs(title="Sox17", x="Cell Type", y=expression("Fold Expression Change")) +
  geom_boxplot(width=0.5)+
  ylim(0,0.041)+
  geom_point(data = df_mean_sox17,
             mapping = aes(x = Cell, y = average),
             color="black") +
  geom_line(data = df_mean_sox17, 
            mapping = aes(x = Cell, y = average, group=1), 
            color= "black")

Sox17_p 

#/////////////////////////////////////////////////////////

Dppa3 <- mnewdata[64:72,]

df_mean_dppa3 <- Dppa3 %>% 
  group_by(Cell) %>% 
  summarize(average = mean(Fold_Expression))

df_mean_dppa3
# Dppa3 at stock concentration

Dppa3_p <- ggplot(Dppa3, aes(x=Cell, y=Fold_Expression, color=Cell)) + 
  theme(plot.title = element_text(hjust = 0.5, size= 20, face="bold"), text = element_text(size = 15), legend.position = "none") +
  labs(title="Dppa3", x="Cell Type", y=expression("Fold Expression Change")) +
  geom_boxplot(width=0.5) +
  geom_point(data = df_mean_dppa3,
             mapping = aes(x = Cell, y = average),
             color="black") +
  geom_line(data = df_mean_dppa3, 
            mapping = aes(x = Cell, y = average, group=1), 
            color= "black")

Dppa3_p

#/////////////////////////////////////////////////////////

Gata4 <- mnewdata[73:81,]

df_mean_gata4 <- Gata4 %>% 
  group_by(Cell) %>% 
  summarize(average = mean(Fold_Expression))

df_mean_gata4
# Gata4 at stock concentration

Gata4_p <- ggplot(Gata4, aes(x=Cell, y=Fold_Expression, color=Cell)) + 
  theme(plot.title = element_text(hjust = 0.5, size= 20, face="bold"), text = element_text(size = 15), legend.position = "none") +
  labs(title="Gata4", x="Cell Type", y=expression("Fold Expression Change")) +
  geom_boxplot(width=0.5) +
  geom_point(data = df_mean_gata4,
             mapping = aes(x = Cell, y = average),
             color="black") +
  geom_line(data = df_mean_gata4, 
            mapping = aes(x = Cell, y = average, group=1), 
            color= "black")

Gata4_p

#//////////////////////////////////////////////////////////

Tfap2c <- mnewdata[82:90,]

df_mean_tfap2c <- Tfap2c %>% 
  group_by(Cell) %>% 
  summarize(average = mean(Fold_Expression))

df_mean_tfap2c
# Tfap2c at stock concentration

Tfap2c_p <- ggplot(Tfap2c, aes(x=Cell, y=Fold_Expression, color=Cell)) + 
  theme(plot.title = element_text(hjust = 0.5, size= 20, face="bold"), text = element_text(size = 15), legend.position = "none") +
  labs(title="Tfap2c", x="Cell Type", y=expression("Fold Expression Change")) +
  geom_boxplot(width=0.5) +
  geom_point(data = df_mean_tfap2c,
             mapping = aes(x = Cell, y = average),
             color="black") +
  geom_line(data = df_mean_tfap2c, 
            mapping = aes(x = Cell, y = average, group=1), 
            color= "black")

Tfap2c_p

#////////////////////////////////////////////////////////////

Prdm14 <- mnewdata[91:99,]

df_mean_prdm14 <- Prdm14 %>% 
  group_by(Cell) %>% 
  summarize(average = mean(Fold_Expression))

df_mean_prdm14[2,2] <- c(0.036382329467)
df_mean_prdm14
# Prdm14 at stock concentration

Prdm14_p <- ggplot(Prdm14, aes(x=Cell, y=Fold_Expression, color=Cell)) + 
  theme(plot.title = element_text(hjust = 0.5, size= 20, face="bold"), text = element_text(size = 15), legend.position = "none") +
  labs(title="Prdm14", x="Cell Type", y=expression("Fold Expression Change")) +
  geom_boxplot(width=0.5) +
  geom_point(data = df_mean_prdm14,
             mapping = aes(x = Cell, y = average),
             color="black") +
  geom_line(data = df_mean_prdm14, 
            mapping = aes(x = Cell, y = average, group=1), 
            color= "black")

Prdm14_p

#///////////////////////////////////////////////////////////

Nanos3 <- mnewdata[100:108,]

df_mean_nanos3 <- Nanos3 %>% 
  group_by(Cell) %>% 
  summarize(average = mean(Fold_Expression))

df_mean_nanos3
# Nanos3 at stock concentration

Nanos3_p <- ggplot(Nanos3, aes(x=Cell, y=Fold_Expression, color=Cell)) + 
  theme(plot.title = element_text(hjust = 0.5, size= 20, face="bold"), text = element_text(size = 15), legend.position = "none") +
  labs(title="Nanos3", x="Cell Type", y=expression("Fold Expression Change")) +
  geom_boxplot(width=0.5) +
  geom_point(data = df_mean_nanos3,
             mapping = aes(x = Cell, y = average),
             color="black") +
  geom_line(data = df_mean_nanos3, 
            mapping = aes(x = Cell, y = average, group=1), 
            color= "black")

Nanos3_p

#///////////////////////////////////////////////////////////

Itgb3 <- mnewdata[109:117,]
Itgb3
df_mean_itgb3 <- Itgb3 %>% 
  group_by(Cell) %>% 
  summarize(average = mean(Fold_Expression))

df_mean_itgb3[1,2] <- c(0.02744629)

# Itgb3 at stock concentration

Itgb3_p <- ggplot(Itgb3, aes(x=Cell, y=Fold_Expression, color=Cell)) + 
  theme(plot.title = element_text(hjust = 0.5, size= 20, face="bold"), text = element_text(size = 15), legend.position = "none") +
  labs(title="Itgb3", x="Cell Type", y=expression("Fold Expression Change")) +
  geom_boxplot(width=0.5) +
  ylim(0,0.8)+
  geom_point(data = df_mean_itgb3,
             mapping = aes(x = Cell, y = average),
             color="black") +
  geom_line(data = df_mean_itgb3, 
            mapping = aes(x = Cell, y = average, group=1), 
            color= "black")

Itgb3_p

#//////////////////////////////////////////////////////////

Dnd1 <- mnewdata[118:126,]

df_mean_dnd1 <- Dnd1 %>% 
  group_by(Cell) %>% 
  summarize(average = mean(Fold_Expression))

df_mean_dnd1
# Dnd1 at stock concentration

Dnd1_p <- ggplot(Dnd1, aes(x=Cell, y=Fold_Expression, color=Cell)) + 
  theme(plot.title = element_text(hjust = 0.5, size= 20, face="bold"), text = element_text(size = 15), legend.position = "none") +
  labs(title="Dnd1", x="Cell Type", y=expression("Fold Expression Change")) +
  geom_boxplot(width=0.5) +
  geom_point(data = df_mean_dnd1,
             mapping = aes(x = Cell, y = average),
             color="black") +
  geom_line(data = df_mean_dnd1, 
            mapping = aes(x = Cell, y = average, group=1), 
            color= "black")

Dnd1_p

#///////////////////////////////////////////////////////////

Dazl <- mnewdata[127:135,]

df_mean_dazl <- Dazl %>% 
  group_by(Cell) %>% 
  summarize(average = mean(Fold_Expression))

df_mean_dazl
# Dazl at stock concentration

Dazl_p <- ggplot(Dazl, aes(x=Cell, y=Fold_Expression, color=Cell)) + 
  theme(plot.title = element_text(hjust = 0.5, size= 20, face="bold"), text = element_text(size = 15), legend.position = "none") +
  labs(title="Dazl", x="Cell Type", y=expression("Fold Expression Change")) +
  geom_boxplot(width=0.5) +
  geom_point(data = df_mean_dazl,
             mapping = aes(x = Cell, y = average),
             color="black") +
  geom_line(data = df_mean_dazl, 
            mapping = aes(x = Cell, y = average, group=1), 
            color= "black")

Dazl_p
