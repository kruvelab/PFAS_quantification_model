

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
library(readxl)


#specify your working directory here:
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
# ggsave(IE_slope_cor,  filename = "results/modelling_results/230331_model_PFAS_train_test_logIE.png", width=16, height=8, units = "cm", device = NULL)
# ggsave(IE_slope_cor, filename = "results/modelling_results/230331_model_PFAS_train_test_logIE.svg", width=16, height=8, units = "cm")
#ggsave(IE_slope_cor,  filename = "C:/Users/HelenSepman/OneDrive - Kruvelab/Helen_phd/presentations/figures/230411_presentation_nonpfas_test.png", width=8, height=8, units = "cm", device = NULL)


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

summary_table_CF2 <- read_delim("results/homologue_vs_IEmodel_results/summary_table_CF2_concentrations.csv")

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

summary_table_CF2CF2 <- read_delim("results/homologue_vs_IEmodel_results/homolgoue_series_conc_summaries/summary_table_CF2CF2_filtered.csv")

#---- plot: QCs ----

QCs_target <- read_excel("results/Melanie_new_suspects/310523_DataForFigures.xlsx",sheet = "QCs", skip = 1)
QCs_target = QCs_target %>% 
  gather(Filename, target_conc, -Compound)

QCs_model = read_delim("results/modelling_results/targets_qc_model_concentrations_with_LOO_changing_names.csv")

QCs_all = QCs_target %>% 
  left_join(QCs_model) %>% 
  mutate(difference =  conc_pred_pg_uL*2 - target_conc) %>% 
  drop_na(difference)


plot_QC <- ggplot(QCs_all, aes(fill=Filename, y=fct_inorder(Compound), x=difference))+
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_manual(values=c("#274c77", "#6096ba", "#a3cef1"))+
  labs(y="", x="pg/µl")+
  my_theme

plot_QC

ggsave(plot_QC, filename = "results/Melanie_new_suspects/QCs_target_vs_model_quant.svg", width=14, height=16, units = "cm")


# Need to multiply model resutls with 2 due to dilution!



#---- plot: Fluorine mass balance (TF, EOF, PFAS (target + semi-quant)) ----


#---- PCA of Targets and suspects based on PaDEL decriptors that ended up in the model after cleaning the descriptors ----

# get the descriptor names form the model
logIE_pred_model_train_test = readRDS(file="models/230329_logIE_model_withPFAS_train_test.RData")
descriptors = colnames(data_clean) 
descriptors = descriptors[! descriptors%in% c("logIE", "pH.aq.", "polarity_index", "viscosity", "NH4", "name", "data_type", "SMILES")]

#suspect SMILES
suspects_SMILES = read_delim("results/Melanie_new_suspects/suspects_smiles_melanie_updated_semicolon2.csv")

suspects_SMILES = suspects_SMILES %>% 
  select(ID, SMILES) %>% 
  rename(Compound = ID) %>% 
  mutate(type = "suspect")

#target SMILES
target_SMILES = read_SMILES(filename = "data_for_modelling/Smiles_for_Target_PFAS_semicolon.csv",
                          compounds_to_be_removed_as_list = c("HFPO-DA", "MeFOSE", "EtFOSE", "10:2 mono PAP", "4:2 mono PAP", "6:2 mono PAP", "8:2 mono PAP"))

target_SMILES = target_SMILES %>% 
  select(Compound, SMILES) %>% 
  mutate(type = "target")

SMILES_all = target_SMILES %>% 
  bind_rows(suspects_SMILES)

#calculate PaDEL descriptors
SMILES_all = PaDEL_original(SMILES_all)

#select only the PaDEL descr that are in the model (also removed eluent descriptors here)
SMILES_all = SMILES_all %>% 
  select(Compound, SMILES, type, descriptors)


#Remove zero variance columns as otherwise cannot scale the data and do PCA
SMILES_all= SMILES_all %>% 
  select(Compound, SMILES, type) %>% 
  bind_cols(SMILES_all[ , which(apply(SMILES_all, 2, var) != 0)])

#create the PCA vectors
chemicals_pca = prcomp(SMILES_all %>%
                         select(-c(Compound, SMILES, type)), 
                       center = TRUE, 
                       scale = TRUE) 
summary(chemicals_pca)

loadings_pca = as_tibble(chemicals_pca$rotation)
scores_pca = as_tibble(predict(chemicals_pca))

scores_pca = scores_pca %>%
  bind_cols(SMILES_all %>%
              select(Compound, SMILES, type)) #add the column with clusters and fruit names


plot_pc1_pc2 = ggplot() + 
  geom_point(mapping = aes(x = PC1, 
                           y = PC2, 
                           color = type),
             data = scores_pca,
             size = 2,
             alpha = 0.75) + 
  scale_color_manual(values=c("#274c77", "#ee6c4d"))+
  labs(x = "PC1 (17.0%)", y = "PC2 (15.8%)") +
  my_theme +
  theme(legend.position = "right")


plot_pc1_pc2


plot_pc1_pc3 = ggplot() + 
  geom_point(mapping = aes(x = PC1, 
                           y = PC3, 
                           color = type),
             data = scores_pca,
             size = 2,
             alpha = 0.75) + 
  scale_color_manual(values=c("#274c77", "#ee6c4d"))+
labs(x = "PC1 (17.0%)", y = "PC3 (11.8%)") +
  my_theme +
  theme(legend.position = "right")

# PC3 (11.8%)
plot_pc1_pc3

joined_plot3 = plot_pc1_pc2 + plot_pc1_pc3
joined_plot3[[1]] = joined_plot3[[1]] + theme(legend.position = "none")

# ggsave(joined_plot3,  filename = "results/Melanie_new_suspects/230411_PCA_targets_suspects.png", width=18, height=10, units = "cm", device = NULL)
# ggsave(joined_plot3, filename = "results/Melanie_new_suspects/230411_PCA_targets_suspects.svg", width=18, height=10, units = "cm")



# 3D plot
# car::scatter3d(data = scores_pca,
#                x = scores_pca$PC1,
#                y = scores_pca$PC2,
#                z = scores_pca$PC3,
#                surface.col = 1:14,
#                main="3D PCA",
#                xlab = "PC1",
#                ylab = "PC2",
#                zlab = "PC3",
#                surface = FALSE,
#                groups = as.factor(scores_pca$type),
#                axis.col = c("black","black","black"),
#                axis.size = 2,
#                sphere.size = 3)


# # cosine similarities
# # first - names of compounds to col names
# SMILES_all_similarity = SMILES_all %>% 
#   select(-c(SMILES, type)) %>% 
#   unique()
# col_names = SMILES_all_similarity$Compound
# 
# SMILES_all_similarity = t(SMILES_all_similarity %>% select(-Compound))
# colnames(SMILES_all_similarity) = col_names
# 
# 
# similarity_scores = cosine(as.matrix(SMILES_all_similarity))
# 


#---- t-SNE of Targets and suspects based on PaDEL decriptors that ended up in the model after cleaning the descriptors ----

# load your omic data here as mydata
install.packages("tsne")
trn <- data.matrix(SMILES_all[-c(1,2,3)])

library(tsne)

cols <- rainbow(10)

# this is the epoch callback function used by tsne. 
# x is an NxK table where N is the number of data rows passed to tsne, and K is the dimension of the map. 
# Here, K is 2, since we use tsne to map the rows to a 2D representation (map).
ecb = function(x, y){ plot(x, t='n'); text(x, labels=trn[,65], col=cols[trn[,65] +1]); }

tsne_res = tsne(trn[,1:64], epoch_callback = ecb, perplexity=50, epoch=50)

