library(tidyverse)
library(plotly)
library(RColorBrewer)
setwd("~/GitHub/PFAS_semi_quant")
source("code/reading_excel.R")

data1 = read_delim("data/Sum PFAS ng_g results.csv",
                         delim = ",",
                         col_names = TRUE)

#col <- brewer.pal(5, "Set2")

data1 = data1 %>%
  mutate(order= case_when(Sample_type=="W. Greenland; White-beaked dolphin"~1,
                          Sample_type=="E. Greenland; Long-finned pilot whale"~2,
                          Sample_type=="Sweden; White-beaked dolphin"~3))

ggplot(data=data1, aes(x = reorder(Sample_type,order,fun=min()), y = Sum_PFAS))+
  geom_bar(aes(fill=Sample_type), position=position_dodge2(preserve = "single"), 
           stat="identity",
           col="black") +
  scale_fill_manual(values = c("#00C5CD","#054C70","#BBFFFF"))+
  geom_errorbar(aes(x=Sample_type, ymin = Sum_PFAS-SD, 
                    ymax = Sum_PFAS+SD),size=1.05,
                width=.9, stat = "identity",
                position = position_dodge2(preserve = "single"),)+
  geom_hline(yintercept = 0)+
  labs(x = "", y = "Sum PFAS ng/g",) +
  theme(text = element_text(size=25),
        panel.background = element_rect(fill = "white"),
        axis.ticks.x = element_line(size = 1),
        axis.title.x = element_text(size=20),
        axis.ticks.length.x = unit(.5, "cm"),
        axis.line.y = element_line(size = 1, color = "black"),
        legend.position = "none")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 750))
