# install.packages("extrafont")
# font_import(paths = "C:/Windows/Fonts")

# palette - https://coolors.co/palettes/trending

# fonts - https://fonts.google.com/ 
# (https://fonts.google.com/specimen/Khula)

# book - visualizations in R
# https://r4ds.had.co.nz/data-visualisation.html

# affinity designer for figures?


# I created this palette:
# https://coolors.co/274c77-6096ba-a3cef1-ee6c4d-8b8c89-293241

library(tidyverse)
library(scales)
library(cowplot)
library(grid)
library(gridExtra)

basecolor <- "#293241"
datapoints_color1 <- "#274c77"
datapoints_color2 <- "#6096ba"  
datapoints_color3 <- "#a3cef1"  
datapoints_color4 <- "#8b8c89"  
highlighter_color <- "#ee6c4d"
  
extrafont::loadfonts(device = "win")
#font <- extrafont::choose_font("Khula")
font <- extrafont::choose_font("Quicksand")
fontsize <- 12

        
#----theme----
my_theme <- theme(
  #remove the background of the plot
  plot.background = element_blank(),
  #and from the panel as well
  panel.background = element_blank(),
  #define the width and color of the axis on the plot
  axis.line = element_line(linewidth = 1,
                           color = basecolor),
  #if you use plot title you can specify parameters here
  #PS! use plot title only if you send or show the plot on its own 
  #for plots on the slide/thesis use slide title and figure caption 
  # plot.title = element_text(color = basecolor,
  #                           size = 14,
  #                           face = "bold"),
  #specify the size and style of the text on the plot, e.g. axis title
  text = element_text(family = font,
                      size = fontsize,
                      color = basecolor),
  legend.key = element_blank(),
  strip.background = element_blank(),
  #strip.text = element_blank(),
  strip.text = element_text(family = font,
                            size = fontsize,
                            color = basecolor),
  #to remove or adjust the position of the legend
  #"none" - is no legend; "top" "bottom", "right", "left";
  #or by coordinates. 
  #c(0, 0) corresponds to the "bottom left" 
  #and c(1, 1) corresponds to the "top right" position.
  legend.position = "none",
  #legend.position = c(0.9, 0.25),
  #if you have a legend title and text you can specify font size here
  #here it indicates no legend title
  legend.title = element_blank(), 
  legend.text = element_text(family = font,
                             size = fontsize,
                             color = basecolor),
  #specify axis marks text
  axis.text = element_text(family = font,
                           size = fontsize,
                           color = basecolor),
  #remove tick marks
  axis.ticks = element_blank(),
  #define the ratio of x and y axis
  #PS! for scatter plots it needs to be 1!
  #for predicted - measured plots also adjust the ranges!
  aspect.ratio = 1,
  #adjust the position of the axis title
  axis.title.x = element_text(hjust = c(1), vjust = c(0)),
  axis.title.y = element_text(hjust = c(1), vjust = c(1))
)        




#---- plot: correlation plot for predicted and measured ionization efficiencies ----

setwd("C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/PFOA_semi_quant")

# read in model
logIE_pred_model_PFAS_allData <- readRDS(file="models/230329_logIE_model_withPFAS_train_test.RData")

#correlation plot
IE_slope_cor = ggplot() +
  geom_point(data = logIE_pred_model_PFAS_allData$data$training_set,
             mapping = aes(logIE, logIE_predicted),
             color = "#274c77",
             alpha = 0.7,
             size = 2) +
  geom_point(data = logIE_pred_model_PFAS_allData$data$test_set,
             mapping = aes(logIE, logIE_predicted),
             color = "#ee6c4d",
             alpha = 0.7,
             size = 2) +
  geom_abline(intercept = -1, slope = 1) +
  geom_abline(intercept = 1, slope = 1) +
  geom_abline(intercept = 0, slope = 1) +
  ylab(substitute(paste("log", italic("IE"))["predicted"]))  +
  xlab(substitute(paste("log", italic("IE"))["measured"])) +
  theme(plot.title = element_text(size = fontsize),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.line.y = element_line(size = 1, color = basecolor),
        axis.line.x = element_line(size = 1, color = basecolor),
        axis.title.x = element_text(size=fontsize),
        axis.title.y = element_text(size=fontsize),
        aspect.ratio = 1,
        axis.text = element_text(family = font,
                                 size = fontsize,
                                 color = basecolor),
        legend.key = element_blank(),
        strip.background = element_blank(),
        text = element_text(family = font,
                            size = fontsize,
                            color = basecolor))+
  facet_wrap(~data_type) +
  annotation_logticks(colour = basecolor) +
  xlim(c(-4,5)) +
  ylim(c(-4,5)) +
  my_theme

IE_slope_cor
# ggsave(IE_slope_cor,  filename = "results/modelling_results/230329_model_PFAS_train_test_logIE.png", width=16, height=8, units = "cm", device = NULL)
# ggsave(IE_slope_cor, filename = "results/modelling_results/230329_model_PFAS_train_test_logIE.svg", width=16, height=8, units = "cm")

ggplotly(IE_slope_cor)

#---- plot: correlation plot: PFAS IE predictions with model without PFAS ----
setwd("C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/PFOA_semi_quant")

# read in model
logIE_pred_model_without_PFAS <- readRDS(file="models/230220_logIE_model_withoutPFAS_allData.RData")

#read in PFAS IE data
PFAS_data = read_delim("data_for_modelling/PFAS_logIE_anchored_PaDEL.csv")

#predict IE values for PFAS
PFAS_data = PFAS_data %>% 
  mutate(logIE_predicted = predict(logIE_pred_model_without_PFAS$model, newdata = PFAS_data)) %>% 
  select(logIE_predicted, everything())


#correlation plot between predicted and anchored IE values for PFAS
IE_slope_cor_model_without_PFAS = ggplot() +
  geom_point(data = PFAS_data,
             mapping = aes(logIE, logIE_predicted),
             color = "#ee6c4d",
             alpha = 0.7,
             size = 2) +
  geom_abline(intercept = -1, slope = 1) +
  geom_abline(intercept =  1, slope = 1) +
  geom_abline(intercept =  0, slope = 1) +
  ylab(substitute(paste("log", italic("IE"))["predicted"]))  +
  xlab(substitute(paste("log", italic("IE"))["measured"])) +
  theme(plot.title = element_text(size = fontsize),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.line.y = element_line(size = 1, color = basecolor),
        axis.line.x = element_line(size = 1, color = basecolor),
        axis.title.x = element_text(size=fontsize),
        axis.title.y = element_text(size=fontsize),
        aspect.ratio = 1,
        axis.text = element_text(family = font,
                                 size = fontsize,
                                 color = basecolor),
        legend.key = element_blank(),
        strip.background = element_blank(),
        text = element_text(family = font,
                            size = fontsize,
                            color = basecolor))+
  annotation_logticks(colour = basecolor) +
  xlim(c(-4,5)) +
  ylim(c(-4,5)) +
  my_theme

IE_slope_cor_model_without_PFAS
# ggsave(IE_slope_cor_model_without_PFAS,  filename = "results/modelling_results/230221_logIE_PFAS_with_model_withoutPFAS.png", width=8, height=8, units = "cm", device = NULL)
# ggsave(IE_slope_cor_model_without_PFAS, filename = "results/modelling_results/230221_logIE_PFAS_with_model_withoutPFAS.svg", width=8, height=8, units = "cm")

ggplotly(IE_slope_cor_model_without_PFAS)

rmse(PFAS_data$logIE, PFAS_data$logIE_predicted)



#---- plot: correlation plot: PFAS IE predictions with leave-one-out approach ----
setwd("C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/PFOA_semi_quant")

#read in PFAS pred IE data
PFAS_LOO_data = read_delim("results/modelling_results/PFAS_pred_logIEs_with_leave_one_out_approach.csv")

#correlation plot between predicted and anchored IE values for PFAS
IE_slope_cor_model_LOO_PFAS = ggplot() +
  geom_point(data = PFAS_LOO_data,
             mapping = aes(logIE, logIE_predicted),
             color = "#ee6c4d",
             alpha = 0.7,
             size = 2) +
  geom_abline(intercept = -1, slope = 1) +
  geom_abline(intercept =  1, slope = 1) +
  geom_abline(intercept =  0, slope = 1) +
  ylab(substitute(paste("log", italic("IE"))["predicted"]))  +
  xlab(substitute(paste("log", italic("IE"))["measured"])) +
  theme(plot.title = element_text(size = fontsize),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.line.y = element_line(size = 1, color = basecolor),
        axis.line.x = element_line(size = 1, color = basecolor),
        axis.title.x = element_text(size=fontsize),
        axis.title.y = element_text(size=fontsize),
        aspect.ratio = 1,
        axis.text = element_text(family = font,
                                 size = fontsize,
                                 color = basecolor),
        legend.key = element_blank(),
        strip.background = element_blank(),
        text = element_text(family = font,
                            size = fontsize,
                            color = basecolor))+
  annotation_logticks(colour = basecolor) +
  xlim(c(-4,5)) +
  ylim(c(-4,5)) +
  my_theme

IE_slope_cor_model_LOO_PFAS
# ggsave(IE_slope_cor_model_LOO_PFAS,  filename = "results/modelling_results/230221_logIE_PFAS_leaveOneOut.png", width=8, height=8, units = "cm", device = NULL)
# ggsave(IE_slope_cor_model_LOO_PFAS, filename = "results/modelling_results/230221_logIE_PFAS_leaveOneOut.svg", width=8, height=8, units = "cm")

ggplotly(IE_slope_cor_model_LOO_PFAS)

rmse(PFAS_LOO_data$logIE, PFAS_LOO_data$logIE_predicted)

#---- plot: cowplot - Liigand model vs leave-one-out approach predicted logIE for PFAS ---- 

joined_plot1 = plot_grid(IE_slope_cor_model_without_PFAS, IE_slope_cor_model_LOO_PFAS) #, labels = c('A', 'B')
# ggsave(joined_plot1,  filename = "results/modelling_results/230330_logIE_PFAS_oldModel_vs_LOO.png", width=16, height=7, units = "cm", device = NULL)
# ggsave(joined_plot1, filename = "results/modelling_results/230330_logIE_PFAS_oldModel_vs_LOO.svg", width=16, height=7, units = "cm")

#---- plot: homologue series quantification vs leave-one-out model quantification ---- 

summary_table_CF2 <- read.csv2("results/homologue_vs_IEmodel_results/homolgoue_series_conc_summaries/summary_table_CF2.csv")

quant_homologue_CF2 = ggplot() +
  geom_point(data = summary_table_CF2,
             mapping = aes(Theoretical_conc_uM, conc_homolog, color = pattern),
             alpha = 0.7,
             size = 2) +
  scale_color_manual(values=c("#274c77", "#a3cef1"))+
  geom_abline(intercept = -1, slope = 1) +
  geom_abline(intercept =  1, slope = 1) +
  geom_abline(intercept =  0, slope = 1) +
  ylab("Predicted concentration (uM)")  +
  xlab("Theoretical concentration (uM)") +
  theme(plot.title = element_text(size = fontsize),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.line.y = element_line(size = 1, color = basecolor),
        axis.line.x = element_line(size = 1, color = basecolor),
        axis.title.x = element_text(size=fontsize),
        axis.title.y = element_text(size=fontsize),
        aspect.ratio = 1,
        axis.text = element_text(family = font,
                                 size = fontsize,
                                 color = basecolor),
        legend.key = element_blank(),
        strip.background = element_blank(),
        text = element_text(family = font,
                            size = fontsize,
                            color = basecolor))+
  annotation_logticks(colour = basecolor) +
  scale_x_log10(limits  = c(10^-5, 10^-0), breaks = 10^(-5:0), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(limits  = c(10^-5, 10^-0), breaks = 10^(-5:0), labels = trans_format("log10", math_format(10^.x))) +
  my_theme

quant_model_LOO_CF2 = ggplot() +
  geom_point(data = summary_table_CF2,
             mapping = aes(Theoretical_conc_uM, conc_pred),
             color = "#274c77",
             alpha = 0.7,
             size = 2) +
  geom_abline(intercept = -1, slope = 1) +
  geom_abline(intercept =  1, slope = 1) +
  geom_abline(intercept =  0, slope = 1) +
  ylab("Predicted concentration (uM)")  +
  xlab("Theoretical concentration (uM)") +
  theme(plot.title = element_text(size = fontsize),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.line.y = element_line(size = 1, color = basecolor),
        axis.line.x = element_line(size = 1, color = basecolor),
        axis.title.x = element_text(size=fontsize),
        axis.title.y = element_text(size=fontsize),
        aspect.ratio = 1,
        axis.text = element_text(family = font,
                                 size = fontsize,
                                 color = basecolor),
        legend.key = element_blank(),
        strip.background = element_blank(),
        text = element_text(family = font,
                            size = fontsize,
                            color = basecolor))+
  annotation_logticks(colour = basecolor) +
  scale_x_log10(limits  = c(10^-5, 10^-0), breaks = 10^(-5:0), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(limits  = c(10^-5, 10^-0), breaks = 10^(-5:0), labels = trans_format("log10", math_format(10^.x))) +
  my_theme


#combine using cowplot
joined_plot2 <- plot_grid(quant_homologue_CF2, quant_model_LOO_CF2)

# #create common x and y labels
# y.grob <- textGrob("Predicted concentration (uM)", 
#                    gp=gpar(fontface=font, col=basecolor, fontsize=fontsize), rot=90)
# 
# x.grob <- textGrob("Theoretical concentration (uM)", 
#                    gp=gpar(fontface=font, col=basecolor, fontsize=fontsize))
# 
# #add to plot
# grid.arrange(arrangeGrob(joined_plot2, left = y.grob, bottom = x.grob))


# ggsave(joined_plot2,  filename = "results/homologue_vs_IEmodel_results/230221_conc_homolog_vs_model.png", width=16, height=7, units = "cm", device = NULL)
# ggsave(joined_plot2, filename = "results/homologue_vs_IEmodel_results/230221_conc_homolog_vs_model.svg", width=16, height=7, units = "cm")

summary_table_CF2CF2 <- read_delim("results/homologue_vs_IEmodel_results/homolgoue_series_conc_summaries/summary_table_CF2CF2_filtered.csv")

#---- plot: Fluorine mass balance (TF, EOF, PFAS (target + semi-quant)) ----

