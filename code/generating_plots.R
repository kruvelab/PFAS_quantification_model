

# palette - https://coolors.co/palettes/trending

# fonts - https://fonts.google.com/ 
# (https://fonts.google.com/specimen/Khula)

# book - visualizations in R
# https://r4ds.had.co.nz/data-visualisation.html

# affinity designer for figures?


# I created this palette:
# https://coolors.co/274c77-6096ba-a3cef1-ee6c4d-8b8c89-293241


# install.packages("extrafont")
# font_import(paths = "C:/Windows/Fonts")
library(tidyverse)
library(scales)
library(cowplot)
library(grid)
library(gridExtra)
library(patchwork)
library(lsa)
library(ggforce)

#specyfy your working directory here:
admin = "C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/PFOA_semi_quant"
setwd(admin)


basecolor <- "#293241"
datapoints_color1 <- "#274c77"
datapoints_color2 <- "#6096ba"  
datapoints_color3 <- "#a3cef1"  
datapoints_color4 <- "#8b8c89"  
highlighter_color <- "#ee6c4d"
  
extrafont::loadfonts(device = "win")
font <- extrafont::choose_font("Arial")
fontsize <- 9

        
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
  #legend.position = "none",
  legend.position = c(0.75, 0.25),
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
  #axis.ticks = element_blank(),
  #define the ratio of x and y axis
  #PS! for scatter plots it needs to be 1!
  #for predicted - measured plots also adjust the ranges!
  #aspect.ratio = 1,
  #adjust the position of the axis title
  axis.title.x = element_text(hjust = c(1), vjust = c(0)),
  axis.title.y = element_text(hjust = c(1), vjust = c(1))
)        




#---- plot: correlation plot for predicted and measured ionization efficiencies ----


# read in model
logIE_pred_model_PFAS_allData <- readRDS(file="models/230619_logIE_model_withPFAS_train_test.RData")

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
        #axis.ticks = element_line(color = basecolor),
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
  #annotation_logticks(colour = basecolor) +
  scale_x_continuous(limits = c(-4, 5), breaks=c(-4, -2, 0, 2, 4))+
  scale_y_continuous(limits = c(-4, 5), breaks=c(-4, -2, 0, 2, 4))+
  # xlim(c(-4,5)) +
  # ylim(c(-4,5)) +
  my_theme

IE_slope_cor
# ggsave(IE_slope_cor,  filename = "results/modelling_results/230703_model_PFAS_train_test_logIE.png", width=16, height=8, units = "cm", device = NULL)
# ggsave(IE_slope_cor, filename = "results/modelling_results/230703_model_PFAS_train_test_logIE.svg", width=16, height=8, units = "cm")
#ggsave(IE_slope_cor,  filename = "C:/Users/HelenSepman/OneDrive - Kruvelab/Helen_phd/presentations/figures/230411_presentation_nonpfas_test.png", width=8, height=8, units = "cm", device = NULL)


ggplotly(IE_slope_cor)

#---- plot: correlation plots separately: PFAS IE predictions with model without PFAS ----

IE_slope_cor_withoutPFAS = ggplot() +
  geom_point(data = logIE_pred_model_PFAS_allData$data$training_set %>% 
               filter(data_type == "non-PFAS"),
             mapping = aes(logIE, logIE_predicted),
             color = "#274c77",
             alpha = 0.7,
             size = 2) +
  geom_point(data = logIE_pred_model_PFAS_allData$data$test_set %>% 
               filter(data_type == "non-PFAS"),
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
        #axis.ticks = element_line(color = basecolor),
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
  #facet_wrap(~data_type) +
  #annotation_logticks(colour = basecolor) +
  scale_x_continuous(limits = c(-4, 5), breaks=c(-4, -2, 0, 2, 4))+
  scale_y_continuous(limits = c(-4, 5), breaks=c(-4, -2, 0, 2, 4))+
  # xlim(c(-4,5)) +
  # ylim(c(-4,5)) +
  my_theme


IE_slope_cor_withPFAS = ggplot() +
  geom_point(data = logIE_pred_model_PFAS_allData$data$training_set %>% 
               filter(data_type == "PFAS"),
             mapping = aes(logIE, logIE_predicted),
             color = "#274c77",
             alpha = 0.7,
             size = 2) +
  geom_point(data = logIE_pred_model_PFAS_allData$data$test_set %>% 
               filter(data_type == "PFAS"),
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
        #axis.ticks = element_line(color = basecolor),
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
  #facet_wrap(~data_type) +
  #annotation_logticks(colour = basecolor) +
  scale_x_continuous(limits = c(-4, 5), breaks=c(-4, -2, 0, 2, 4))+
  scale_y_continuous(limits = c(-4, 5), breaks=c(-4, -2, 0, 2, 4))+
  # xlim(c(-4,5)) +
  # ylim(c(-4,5)) +
  my_theme


#---- plot: correlation plot: PFAS IE predictions with model without PFAS ----
setwd("C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/PFOA_semi_quant")

# read in model
logIE_pred_model_without_PFAS <- readRDS(file="models/230703_logIE_model_withoutPFAS_allData.RData")

#read in PFAS IE data
PFAS_data = read_delim("data_for_modelling/230703_PFAS_IE_anchored_PaDEL.csv")

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
  #annotation_logticks(colour = basecolor) +
  scale_x_continuous(limits = c(-4, 5), breaks=c(-4, -2, 0, 2, 4))+
  scale_y_continuous(limits = c(-4, 5), breaks=c(-4, -2, 0, 2, 4))+
  # xlim(c(-4,5)) +
  # ylim(c(-4,5)) +
  my_theme

IE_slope_cor_model_without_PFAS
# ggsave(IE_slope_cor_model_without_PFAS,  filename = "results/modelling_results/230411_logIE_PFAS_with_model_withoutPFAS.png", width=8, height=8, units = "cm", device = NULL)
# ggsave(IE_slope_cor_model_without_PFAS, filename = "results/modelling_results/230411_logIE_PFAS_with_model_withoutPFAS.svg", width=8, height=8, units = "cm")

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
  scale_x_continuous(limits = c(-4, 5), breaks=c(-4, -2, 0, 2, 4))+
  scale_y_continuous(limits = c(-4, 5), breaks=c(-4, -2, 0, 2, 4))+
  #annotation_logticks(colour = basecolor) +
  # xlim(c(-4,5)) +
  # ylim(c(-4,5)) +
  my_theme

IE_slope_cor_model_LOO_PFAS
# ggsave(IE_slope_cor_model_LOO_PFAS,  filename = "results/modelling_results/230411_logIE_PFAS_leaveOneOut.png", width=8, height=8, units = "cm", device = NULL)
# ggsave(IE_slope_cor_model_LOO_PFAS, filename = "results/modelling_results/230411_logIE_PFAS_leaveOneOut.svg", width=8, height=8, units = "cm")

ggplotly(IE_slope_cor_model_LOO_PFAS)

rmse(PFAS_LOO_data$logIE, PFAS_LOO_data$logIE_predicted)

#---- plot: joined plot - Liigand model vs leave-one-out approach predicted logIE for PFAS ---- 
#put plots together
joined_plot1 = IE_slope_cor_model_without_PFAS + IE_slope_cor_model_LOO_PFAS 

#take off some axis
joined_plot1[[2]] = joined_plot1[[2]] + theme(axis.text.y = element_blank(),
                                              axis.ticks.y = element_blank(),
                                              axis.title.y = element_blank() )

joined_plot1[[1]] = joined_plot1[[1]] + theme(axis.ticks.x = element_blank(),
                                              axis.title.x = element_blank() )
#annotate plots
joined_plot1 = joined_plot1 + plot_annotation(tag_levels = "A")


# joined_plot1 = plot_grid(IE_slope_cor_model_without_PFAS, IE_slope_cor_model_LOO_PFAS) #, labels = c('A', 'B')
# ggsave(joined_plot1,  filename = "results/modelling_results/230331_logIE_PFAS_oldModel_vs_LOO.png", width=16, height=8, units = "cm", device = NULL)
# ggsave(joined_plot1, filename = "results/modelling_results/230331_logIE_PFAS_oldModel_vs_LOO.svg", width=16, height=8, units = "cm")

#---- plot: joined plot - homologue series quantification vs leave-one-out model quantification ---- 

summary_table_CF2 <- read_delim("results/homologue_vs_IEmodel_results/230703_summary_table_CF2_concentrations.csv")

quant_homologue_CF2 = ggplot() +
  geom_point(data = summary_table_CF2,
             mapping = aes(Theoretical_conc_uM/10^6, conc_homolog_uM/10^6, color = pattern),
             alpha = 0.7,
             size = 2) +
  scale_color_manual(values=c("#274c77", "#a3cef1"))+
  geom_abline(intercept = -1, slope = 1) +
  geom_abline(intercept =  1, slope = 1) +
  geom_abline(intercept =  0, slope = 1) +
  ylab("Predicted concentration (M)")  +
  xlab("Theoretical concentration (M)") +
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
        axis.text.y = element_text(hjust = 0), # fix y axis number to align them
        legend.key = element_blank(),
        strip.background = element_blank(),
        text = element_text(family = font,
                            size = fontsize,
                            color = basecolor))+
  annotation_logticks(colour = basecolor) +
  scale_x_log10(limits  = c(10^-11, 10^-6), breaks = 10^(-10:6), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(limits  = c(10^-11, 10^-6), breaks = 10^(-10:6), labels = trans_format("log10", math_format(10^.x))) +
  my_theme

quant_model_LOO_CF2 = ggplot() +
  geom_point(data = summary_table_CF2,
             mapping = aes(Theoretical_conc_uM/10^6, conc_pred_uM/10^6),
             color = "#274c77",
             alpha = 0.7,
             size = 2) +
  geom_abline(intercept = -1, slope = 1) +
  geom_abline(intercept =  1, slope = 1) +
  geom_abline(intercept =  0, slope = 1) +
  ylab("Predicted concentration (M)")  +
  xlab("Theoretical concentration (M)") +
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
        axis.text.y = element_text(hjust = 0), # fix y axis number to align them
        legend.key = element_blank(),
        strip.background = element_blank(),
        text = element_text(family = font,
                            size = fontsize,
                            color = basecolor))+
  annotation_logticks(colour = basecolor) +
  scale_x_log10(limits  = c(10^-11, 10^-6), breaks = 10^(-10:6), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(limits  = c(10^-11, 10^-6), breaks = 10^(-10:6), labels = trans_format("log10", math_format(10^.x))) +
  my_theme


#put plots together
joined_plot2 = quant_homologue_CF2 + quant_model_LOO_CF2 

#take off some axis
joined_plot2[[2]] = joined_plot2[[2]] + theme(axis.text.y = element_blank(),
                                              axis.ticks.y = element_blank(),
                                              axis.title.y = element_blank() )

joined_plot2[[1]] = joined_plot2[[1]] + theme(axis.ticks.x = element_blank(),
                                              axis.title.x = element_blank() )
#annotate plots
joined_plot2 = joined_plot2 + plot_annotation(tag_levels = "A")



# ggsave(joined_plot2,  filename = "results/homologue_vs_IEmodel_results/230331_conc_homolog_vs_model.png", width=16, height=8, units = "cm", device = NULL)
# ggsave(joined_plot2, filename = "results/homologue_vs_IEmodel_results/230331_conc_homolog_vs_model.svg", width=16, height=8, units = "cm")

#---- plot: joined plot - CF2CF2 homologue series quantification vs leave-one-out model quantification ---- 

summary_table_CF2CF2 <- read_delim("results/homologue_vs_IEmodel_results/230703_summary_table_CF2CF2_concentrations.csv")

quant_homologue_CF2CF2 = ggplot() +
  geom_point(data = summary_table_CF2CF2,
             mapping = aes(Theoretical_conc_uM/10^6, conc_homolog_uM/10^6, color = pattern),
             alpha = 0.7,
             size = 2) +
  scale_color_manual(values=c("#274c77", "#a3cef1"))+
  geom_abline(intercept = -1, slope = 1) +
  geom_abline(intercept =  1, slope = 1) +
  geom_abline(intercept =  0, slope = 1) +
  ylab("Predicted concentration (M)")  +
  xlab("Theoretical concentration (M)") +
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
        axis.text.y = element_text(hjust = 0), # fix y axis number to align them
        legend.key = element_blank(),
        strip.background = element_blank(),
        text = element_text(family = font,
                            size = fontsize,
                            color = basecolor))+
  annotation_logticks(colour = basecolor) +
  scale_x_log10(limits  = c(10^-11, 10^-6), breaks = 10^(-10:6), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(limits  = c(10^-11, 10^-6), breaks = 10^(-10:6), labels = trans_format("log10", math_format(10^.x))) +
  my_theme

quant_model_LOO_CF2CF2 = ggplot() +
  geom_point(data = summary_table_CF2CF2,
             mapping = aes(Theoretical_conc_uM/10^6, conc_pred_uM/10^6),
             color = "#274c77",
             alpha = 0.7,
             size = 2) +
  geom_abline(intercept = -1, slope = 1) +
  geom_abline(intercept =  1, slope = 1) +
  geom_abline(intercept =  0, slope = 1) +
  ylab("Predicted concentration (M)")  +
  xlab("Theoretical concentration (M)") +
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
        axis.text.y = element_text(hjust = 0), # fix y axis number to align them
        legend.key = element_blank(),
        strip.background = element_blank(),
        text = element_text(family = font,
                            size = fontsize,
                            color = basecolor))+
  annotation_logticks(colour = basecolor) +
  scale_x_log10(limits  = c(10^-11, 10^-6), breaks = 10^(-10:6), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(limits  = c(10^-11, 10^-6), breaks = 10^(-10:6), labels = trans_format("log10", math_format(10^.x))) +
  my_theme


#put plots together
joined_plot4 = quant_homologue_CF2CF2 + quant_model_LOO_CF2CF2

#take off some axis
joined_plot4[[2]] = joined_plot4[[2]] + theme(axis.text.y = element_blank(),
                                              axis.ticks.y = element_blank(),
                                              axis.title.y = element_blank() )

joined_plot4[[1]] = joined_plot4[[1]] + theme(axis.ticks.x = element_blank(),
                                              axis.title.x = element_blank() )
#annotate plots
joined_plot4 = joined_plot4 + plot_annotation(tag_levels = "A")



# ggsave(joined_plot4,  filename = "results/homologue_vs_IEmodel_results/230331_conc_homolog_vs_model_C2F4.png", width=16, height=8, units = "cm", device = NULL)
# ggsave(joined_plot4, filename = "results/homologue_vs_IEmodel_results/230331_conc_homolog_vs_model_C2F4.svg", width=16, height=8, units = "cm")


#---- plot: all previous combined ----

joined_plot_all1 = (IE_slope_cor_withoutPFAS + IE_slope_cor_withPFAS)
joined_plot_all2 = (IE_slope_cor_model_without_PFAS + IE_slope_cor_model_LOO_PFAS)
joined_plot_all3 = (quant_homologue_CF2 + quant_model_LOO_CF2)

# G1 + plot_spacer() + G2 + plot_layout(widths = c(4, -1.1 ,4.5),guides = "collect")& theme(legend.position = "top")
#take off some axis
joined_plot_all1[[2]] = joined_plot_all1[[2]] + theme(axis.text.y = element_blank(),
                                              axis.ticks.y = element_blank(),
                                              axis.title.y = element_blank() ) 
joined_plot_all1[[1]] = joined_plot_all1[[1]] + theme(axis.ticks.x = element_blank(),
                                              axis.title.x = element_blank() )

joined_plot_all2[[2]] = joined_plot_all2[[2]] + theme(axis.text.y = element_blank(),
                                                      axis.ticks.y = element_blank(),
                                                      axis.title.y = element_blank() ) 
joined_plot_all2[[1]] = joined_plot_all2[[1]] + theme(axis.ticks.x = element_blank(),
                                                      axis.title.x = element_blank() )

joined_plot_all3[[2]] = joined_plot_all3[[2]] + theme(axis.text.y = element_blank(),
                                                      axis.ticks.y = element_blank(),
                                                      axis.title.y = element_blank() ) 
joined_plot_all3[[1]] = joined_plot_all3[[1]] + theme(axis.ticks.x = element_blank(),
                                                      axis.title.x = element_blank() )


joined_plot_all = (joined_plot_all1 + plot_layout(tag_level = "new"))/
  (joined_plot_all2 + plot_layout(tag_level = "new"))/
  (joined_plot_all3 + plot_layout(tag_level = "new"))  +  plot_annotation(tag_levels = c('A', '1')) +
  theme(plot.tag = element_text(size = fontsize))

joined_plot_all = (joined_plot_all1)/
  (joined_plot_all2)/
  (joined_plot_all3)  +  plot_annotation(tag_levels = c('A')) +
  theme(plot.tag = element_text(size = fontsize))

# ggsave(joined_plot_all, filename = "results/modelling_results/230703_all_modelling_comparison.svg", width=16, height=24, units = "cm")
# ggsave(joined_plot_all, filename = "results/modelling_results/230703_all_modelling_comparison.png", width=16, height=24, units = "cm")

#---- TEST plot: all previous combined ----

joined_plot_all11 = (IE_slope_cor_withoutPFAS)
joined_plot_all12 = (IE_slope_cor_withPFAS)
joined_plot_all21 = (IE_slope_cor_model_without_PFAS)
joined_plot_all22 = (IE_slope_cor_model_LOO_PFAS)
joined_plot_all31 = (quant_homologue_CF2)
joined_plot_all32 = (quant_model_LOO_CF2)

#take off some axis
joined_plot_all12 = joined_plot_all12 + theme(axis.ticks.y = element_blank(),
                                                      axis.title.y = element_blank() ) 
joined_plot_all11 = joined_plot_all11 + theme(axis.ticks.x = element_blank(),
                                                      axis.title.x = element_blank(),
                                                      axis.text.x = element_blank())

joined_plot_all22 = joined_plot_all22 + theme(axis.ticks.y = element_blank(),
                                                      axis.title.y = element_blank() ) 
joined_plot_all21 = joined_plot_all21 + theme(axis.ticks.x = element_blank(),
                                                      axis.title.x = element_blank(),
                                                      axis.text.x = element_blank())

joined_plot_all32 = joined_plot_all32 + theme(axis.ticks.y = element_blank(),
                                                      axis.title.y = element_blank() ) 
joined_plot_all31 = joined_plot_all31 + theme(axis.ticks.x = element_blank(),
                                                      axis.title.x = element_blank(),
                                                      axis.text.x = element_blank())



design <- "AACCEE
           AACCEE
           BBDDFF
           BBDDFF"


joined_plot_all = joined_plot_all11 + joined_plot_all12 +
  joined_plot_all21 + joined_plot_all22 +
  joined_plot_all31 + joined_plot_all32 +
  plot_layout(design = design) + plot_annotation(tag_levels = c('A'))
  
  (joined_plot_all1 + plot_layout(tag_level = "new"))/
  (joined_plot_all2 + plot_layout(tag_level = "new", byrow = F))/
  (joined_plot_all3 + plot_layout(tag_level = "new"))  +  plot_annotation(tag_levels = c('A', '1')) +
  theme(plot.tag = element_text(size = fontsize))

joined_plot_all = (joined_plot_all1) +
  (joined_plot_all2) +
  (joined_plot_all3)  +  plot_annotation(tag_levels = c('A')) +
  theme(plot.tag = element_text(size = fontsize)) 
# 
# ggsave(joined_plot_all, filename = "results/modelling_results/230811_all_modelling_comparison.svg", width=18, height=10, units = "cm")
# ggsave(joined_plot_all, filename = "results/modelling_results/230811_all_modelling_comparison.png", width=18, height=10, units = "cm")


#---- plot: Fluorine mass balance (EOF vs quantified PFAS (target + suspect model quant)) ----
fluorine_massbalance_data <- read_excel("results/Melanie_new_suspects/310523_DataForFigures.xlsx",sheet = "FluorineMassBalance")

fluorine_massbalance_data = fluorine_massbalance_data %>% 
  mutate(bars_separation = case_when(Analysis == "EOF" ~ "EOF",
                                     TRUE ~ "quantified")) %>% 
  group_by(Sample) %>% 
  mutate(sample_type = str_split(Sample, "-")[[1]][1]) %>% 
  ungroup() %>% 
  mutate(sample_type = case_when(sample_type == "SL" ~ "Swedish Dolphins",
                                 sample_type == "PL"~ "Greenlandic Pilot Whales",
                                 TRUE ~ "Greenlandic Dolphins")) %>% 
  mutate(sample_type = factor(sample_type, levels = c( "Swedish Dolphins", "Greenlandic Pilot Whales","Greenlandic Dolphins")),
         Analysis = factor(Analysis, levels = c( "EOF","Suspects", "∑PFAS")))

F_balance_barplot =  ggplot(fluorine_massbalance_data, aes(fill=Analysis, y=Sample, x=Value))+
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_manual(values=c("#ee6c4d", "#274c77", "#a3cef1" )) +
  labs(y="", x="ng F-/g w w")+
  my_theme



ggplot() + 
  geom_bar(data = fluorine_massbalance_data %>%  filter(bars_separation != "EOF"), #, sample_type  == "SL"),
           mapping = aes(y=Value , x=Sample, fill=Analysis),
           stat ="identity",
           position = "stack",
           width = 0.3,
           just = -0.5) +
  geom_bar(data = fluorine_massbalance_data %>%  filter(bars_separation == "EOF"), #, sample_type  == "SL"),
           mapping = aes(y=Value , x=Sample, fill=Analysis),
           stat ="identity",
           position = "stack",
           width = 0.3,
           just = 0.5) +
  scale_fill_manual(breaks=c("EOF", "∑PFAS", "Suspects"), values=c("#274c77", "#a3cef1", "#ee6c4d" )) +
  facet_grid(~sample_type, scales = "free",space ="free") +
  my_theme +
  theme(
    strip.text.x = element_blank()
  )


# ---- model vs target - QC plots ----
#modified the names in the file:
QC_all = read_delim("results/modelling_results/QC_target_quant_LOO_quant_for_plot.csv", delim = ",") %>%  
  select(Compound, Filename, conc_pg_uL, quant_type)

plot_QC <- ggplot(QC_all %>% filter(grepl("AL", Filename)), aes(fill=quant_type, y=fct_inorder(Compound), x=(conc_pg_uL)))+
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_manual(values=c("#ee6c4d", "#274c77")) +
  labs(y="", x="pg/µl")+
  #scale_x_log10() +
  my_theme


# ggsave(plot_QC,  filename = "results/homologue_vs_IEmodel_results/smaller_230811_QC_model_target_sampleA.png", width=12, height=14, units = "cm", device = NULL)
# ggsave(plot_QC, filename = "results/homologue_vs_IEmodel_results/smaller_230811_QC_model_target_sampleA.svg", width=12, height=14, units = "cm")

QC_all= QC_all %>% 
  pivot_wider(names_from = quant_type, values_from = conc_pg_uL) %>% 
  mutate(quant_difference = target-model)

plot_QC2 <- ggplot(QC_all %>% filter(grepl("AL", Filename)), aes(fill=Filename, y=fct_inorder(Compound), x=quant_difference))+
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_manual(values=c("#274c77")) +
  labs(y="", x="Δpg/µl")+
  my_theme

# ggsave(plot_QC2,  filename = "results/homologue_vs_IEmodel_results/230704_QC_target_minus_model_sampleA.png", width=12, height=14, units = "cm", device = NULL)
# ggsave(plot_QC2, filename = "results/homologue_vs_IEmodel_results/230704_QC_target_minus_model_sampleA.svg", width=12, height=14, units = "cm")

# ---- model vs target - QC plots to SI (QC-B and QC-C) ----
#modified the names in the file:
QC_all = read_delim("results/modelling_results/QC_target_quant_LOO_quant_for_plot.csv", delim = ",") %>%  
  select(Compound, Filename, conc_pg_uL, quant_type)

plot_QC0 <- ggplot(QC_all %>% filter(grepl("AL", Filename)), aes(fill=quant_type, y=fct_inorder(Compound), x=(conc_pg_uL)))+
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_manual(values=c("#ee6c4d", "#274c77")) +
  labs(y="", x="pg/µl")+
  #scale_x_log10() +
  my_theme

plot_QC1 <- ggplot(QC_all %>% filter(grepl("BL", Filename)), aes(fill=quant_type, y=fct_inorder(Compound), x=(conc_pg_uL)))+
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_manual(values=c("#ee6c4d", "#274c77")) +
  labs(y="", x="pg/µl")+
  #scale_x_log10() +
  my_theme

plot_QC2 <- ggplot(QC_all %>% filter(grepl("CL", Filename)), aes(fill=quant_type, y=fct_inorder(Compound), x=(conc_pg_uL)))+
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_manual(values=c("#ee6c4d", "#274c77")) +
  labs(y="", x="pg/µl")+
  #scale_x_log10() +
  my_theme

#put plots together
joined_plot5 = plot_QC0 + plot_QC1 + plot_QC2 

#take off some axis
joined_plot5[[2]] = joined_plot5[[2]] + theme(axis.text.y = element_blank(),
                                              axis.ticks.y = element_blank(),
                                              axis.title.y = element_blank(),
                                              axis.title.x = element_blank(),
                                              legend.position = "none")

joined_plot5[[1]] = joined_plot5[[1]] + theme(axis.title.x = element_blank(),
                                              legend.position = "none")

joined_plot5[[3]] = joined_plot5[[3]] + theme(axis.text.y = element_blank(),
                                              axis.ticks.y = element_blank(),
                                              axis.title.y = element_blank())


#annotate plots
joined_plot5 = joined_plot5 + plot_annotation(tag_levels = "A")


# ggsave(joined_plot5,  filename = "results/homologue_vs_IEmodel_results/230811_QC-ABC.png", width=20, height=14, units = "cm", device = NULL)
# ggsave(joined_plot5, filename = "results/homologue_vs_IEmodel_results/230811_QC-ABC.svg", width=20, height=14, units = "cm")


