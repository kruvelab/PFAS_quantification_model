library(readxl)
library(tidyr)
library(ggplot2)
library(plotly)
library(cowplot)
library(forcats)
library(gridExtra)
library(RColorBrewer)


setwd("C:/Users/MelanieLauria/Documents/GitHub/PFOA_semi_quant/results/Melanie_new_suspects")

#QCs plot (targeted-predicted in pg/ul)-----------------------------------------


QCs_data <- read_excel("310523_DataForFigures.xlsx",sheet = "QCs")
View(QCs)

QCs <- gather(QCs_data, Type, Value, -Compounds)

#order correctly (and not in alphabetical order) using factor
#QCs$type <- factor(QCs$type,
#                    levels=c("QC_A", "QC_B", "QC_C"))
#In the end I used fct_inorder() instead of using factors

plot_QC <- ggplot(QCs, aes(fill=Type, y=fct_inorder(Compounds), x=Value))+
  geom_bar(stat = "identity", position = "dodge")+
  labs(y="", x="pg/??l")+
  theme_classic()

plot_QC

#Fluorine mass balance plot-----------------------------------------------------


fmb_EOF <-read_excel("310523_DataForFigures.xlsx",sheet = "Sheet1")
fmb_PFAS <-read_excel("310523_DataForFigures.xlsx",sheet = "Sheet2")
#View(fmb_EOF)
#View(fmb_PFAS)
barwidth=0.3
my_labels <- c('SD2', 'SD3', 'SD1', 'PW3', 'PW1','PW2','PW4','PW5','GD4', 'GD3', 'GD2', 'GD1', 'GD5')

plot_fmb <- ggplot()+
  geom_bar(data = fmb_EOF,                                                       #EOF
           aes(fill=Analysis,
               y = Value, 
               x = Sample), 
           stat = "identity", 
           position = "stack",
           width = barwidth)+
  geom_errorbar(data = fmb_EOF, 
                mapping = aes(x=Sample, ymin=Value-Error, ymax=Value+Error),
                width=0.1, 
                linewidth=0.6)+
  
  geom_bar(data = fmb_PFAS,                                                      #stacked sumPFAS and sumSuspects
           aes(fill=fct_inorder(Analysis),
               y = Value, 
               x = Sample + barwidth), 
           stat = "identity", 
           position = "stack",
           width = barwidth)+
  geom_errorbar(data = fmb_PFAS, 
                mapping = aes(x = Sample + barwidth, ymin = Error_minus, ymax = Error_plus),
                width=0.1, 
                linewidth=0.6
                )+
  coord_flip()+                                                                  #Coordinates flipped
  scale_x_reverse(breaks= c(1:13), labels= my_labels)+
  xlab("")+
  ylab("ng F/g ww")+
  theme_classic()+
  facet_wrap(ncol=1, fct_inorder(Group)~., scales="free", strip.position="left")+
  scale_fill_brewer(palette="Blues")


plot_fmb




  
#only Swedish dolphins----------------------------------------------------------

plot_fmb_SD <- ggplot()+
  geom_bar(data = subset(fmb_EOF, Group== "Swedish Dolphins"), 
           aes(fill=Analysis,
               y = Value, 
               x = Sample), 
           stat = "identity", 
           position = "stack",
           width = barwidth)+
  geom_errorbar(data = subset(fmb_EOF, Group== "Swedish Dolphins"), 
                mapping = aes(x=Sample, ymin=Value-Error, ymax=Value+Error),
                width=0.1, 
                linewidth=0.6)+
  
  geom_bar(data = subset(fmb_PFAS, Group== "Swedish Dolphins"), 
           aes(fill=fct_inorder(Analysis),
               y = Value, 
               x = Sample + barwidth), 
           stat = "identity", 
           position = "stack",
           width = barwidth)+
  geom_errorbar(data = subset(fmb_PFAS, Group== "Swedish Dolphins"), 
                mapping = aes(x = Sample + barwidth, ymin = Error_minus, ymax = Error_plus),
                width=0.1, 
                linewidth=0.6
  )+
  coord_flip()

plot_fmb_SD

#only Pilot Whales--------------------------------------------------------------

plot_fmb_PW <- ggplot()+
  geom_bar(data = subset(fmb_EOF, Group== "Greenlandic Pilot Whales"), 
           aes(fill=Analysis,
               y = Value, 
               x = Sample), 
           stat = "identity", 
           position = "stack",
           width = barwidth)+
  geom_errorbar(data = subset(fmb_EOF, Group== "Greenlandic Pilot Whales"), 
                mapping = aes(x=Sample, ymin=Value-Error, ymax=Value+Error),
                width=0.1, 
                linewidth=0.6)+
  
  geom_bar(data = subset(fmb_PFAS, Group== "Greenlandic Pilot Whales"), 
           aes(fill=fct_inorder(Analysis),
               y = Value, 
               x = Sample + barwidth), 
           stat = "identity", 
           position = "stack",
           width = barwidth)+
  geom_errorbar(data = subset(fmb_PFAS, Group== "Greenlandic Pilot Whales"), 
                mapping = aes(x = Sample + barwidth, ymin = Error_minus, ymax = Error_plus),
                width=0.1, 
                linewidth=0.6
  )+
  coord_flip()

plot_fmb_PW

#Only Greenlandic Dolphins------------------------------------------------------

plot_fmb_GD <- ggplot()+
  geom_bar(data = subset(fmb_EOF, Group== "Greenlandic Dolphins"), 
           aes(fill=Analysis,
               y = Value, 
               x = Sample), 
           stat = "identity", 
           position = "stack",
           width = barwidth)+
  geom_errorbar(data = subset(fmb_EOF, Group== "Greenlandic Dolphins"), 
                mapping = aes(x=Sample, ymin=Value-Error, ymax=Value+Error),
                width=0.1, 
                linewidth=0.6)+
  
  geom_bar(data = subset(fmb_PFAS, Group== "Greenlandic Dolphins"), 
           aes(fill=fct_inorder(Analysis),
               y = Value, 
               x = Sample + barwidth), 
           stat = "identity", 
           position = "stack",
           width = barwidth)+
  geom_errorbar(data = subset(fmb_PFAS, Group== "Greenlandic Dolphins"), 
                mapping = aes(x = Sample + barwidth, ymin = Error_minus, ymax = Error_plus),
                width=0.1, 
                linewidth=0.6
  )+
  coord_flip()

plot_fmb_GD

#Combined plot of the 3---------------------------------------------------------
grid.arrange(plot_fmb_SD, plot_fmb_PW, plot_fmb_GD, nrow = 3)
