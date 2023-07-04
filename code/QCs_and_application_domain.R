setwd("C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/PFOA_semi_quant")
source("code/functions.R")
library(tidyverse)
library(readxl)
library(Rtsne)


###  QC sample concentrations with LOO models

## ---- Reading in LC-MS data of calibration solutions ----
Orbitrap_dataset_raw = read_excel_allsheets(filename = "data_for_modelling/Batch 1 Semi Quant w frag.xlsx")

Orbitrap_dataset_raw = Orbitrap_dataset_raw %>% 
  group_by(Compound) %>%
  mutate(Theoretical_amt = case_when(
    Filename == "2020071205-cal21" ~ mean(Theoretical_amt[Filename=="2020071205-cal22"]),
    TRUE ~ Theoretical_amt))%>%
  ungroup()

SMILES_data = read_SMILES(filename = "data_for_modelling/Smiles_for_Target_PFAS_semicolon.csv",
                          compounds_to_be_removed_as_list = c("HFPO-DA", "MeFOSE", "EtFOSE", "10:2 mono PAP", "4:2 mono PAP", "6:2 mono PAP", "8:2 mono PAP", "PFOcDA"))

data_real_conc = Orbitrap_dataset_raw %>%
  mutate(Theoretical_amt = replace(Theoretical_amt , grepl("NaN", Theoretical_amt, fixed = TRUE), NA)) %>%
  group_by(Compound) %>%
  mutate(Theoretical_amt = case_when(
    Filename == "2020071205-cal21" ~ mean(Theoretical_amt[Filename=="2020071205-cal22"]),
    TRUE ~ Theoretical_amt))%>%
  ungroup() %>%
  left_join(SMILES_data) %>%
  drop_na(SMILES) %>%
  mutate(RT = as.numeric(RT),
         area_IC = Area*IC,
         Theoretical_conc_uM = Theoretical_amt/Molecular_weight) 

#---- calculate concentrations ----

table_with_PFAS_LOO_pred_logIEs = read_delim("results/modelling_results/PFAS_pred_logIEs_with_leave_one_out_approach.csv")
directory_with_LOO_models = "C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/PFOA_semi_quant/models/leave_one_out_approach"
data_detected_PFAS = data_real_conc

LOO_model_pred_QC_conc = concentration_forPFAS_pretrained_models(SMILES_names_with_homologues = data_real_conc %>% select(Compound, SMILES) %>%  unique(), # use all SMILES
                                                                    table_with_PFAS_LOO_pred_logIEs = table_with_PFAS_LOO_pred_logIEs,
                                                                    directory_with_LOO_models = directory_with_LOO_models,
                                                                    data_detected_PFAS = data_detected_PFAS)


QC_conc = LOO_model_pred_QC_conc$predicted_conc %>% 
  filter(grepl("QC", Filename))


#write_delim(QC_conc, "results/modelling_results/targets_qc_model_concentrations_with_LOO.csv", delim = ",")

# ---- combine target and model quantification for QCs so that they could be made into a plot ----

QC_model_conc = QC_conc %>% 
  select(Compound, Filename, conc_pred_pg_uL) %>% 
  mutate(quant_type = "model") %>% 
  rename(conc_pg_uL = conc_pred_pg_uL)

QCs_target_data <- read_excel("results/Melanie_new_suspects/310523_DataForFigures.xlsx", sheet = "QCs", skip = 1) %>% 
  gather(key = "Filename", value = "conc_pg_uL", - c(Compound)) %>% 
  mutate(quant_type = "target")

QC_all = QC_model_conc %>% 
  bind_rows(QCs_target_data)

#write_delim(QC_all, "results/modelling_results/QC_target_quant_LOO_quant_for_plot.csv", delim = ",")

# ---- QC potential plots ----
#modified the names in the file:
QC_all = read_delim("results/modelling_results/QC_target_quant_LOO_quant_for_plot.csv", delim = ",") %>%  
  select(Compound, Filename, conc_pg_uL, quant_type)

plot_QC <- ggplot(QC_all, aes(fill=quant_type, y=fct_inorder(Compound), x=conc_pg_uL))+
  geom_bar(stat = "identity", position = "dodge")+
  labs(y="", x="pg/µl")+
  theme_classic()

# fold-differences?

QC_all= QC_all %>% 
  pivot_wider(names_from = quant_type, values_from = conc_pg_uL) %>% 
  mutate(quant_fold_difference = case_when(model > target ~ model/target,
                                           TRUE ~ target/model)) 
  

QC_all %>% 
  group_by(Filename) %>% 
  summarize(mean_fold_difference = mean(quant_fold_difference)) %>% 
  ungroup()


#---- PCA of Targets and suspects based on PaDEL decriptors that ended up in the model after cleaning the descriptors ----

# get the descriptor names form the model
logIE_pred_model_train_test = readRDS(file="models/230619_logIE_model_withPFAS_train_test.RData")
descriptors = colnames(data_clean) 
descriptors = descriptors[! descriptors%in% c("logIE", "pH.aq.", "polarity_index", "viscosity", "NH4", "name", "data_type", "SMILES")]

#target SMILES
target_SMILES <- data_all_binded %>% 
  filter(data_type == "PFAS") %>% 
  rename(Compound = name) %>% 
  select(Compound, SMILES) %>% 
  mutate(type = "Targets") %>% 
  unique()

#suspect SMILES
suspects_SMILES = read_delim("results/Melanie_new_suspects/suspects_smiles_melanie_updated_semicolon2.csv")

suspects_SMILES = suspects_SMILES %>% 
  select(ID, SMILES) %>% 
  left_join(target_SMILES) %>% 
  filter(is.na(Compound)) %>% 
  select(-Compound) %>% 
  rename(Compound = ID) %>% 
  mutate(type = "Suspects")


#all old data smiles that were used for modelling
old_data <- data_all_binded %>% 
  filter(data_type == "non-PFAS") %>% 
  rename(Compound = name) %>%
  select(Compound, SMILES) %>% 
  mutate(type = "Liigand et al.") %>% 
  unique()


SMILES_all = target_SMILES %>% 
  bind_rows(suspects_SMILES) %>% 
  bind_rows(old_data)

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
             alpha = 0.5) + 
  scale_color_manual(values=c("#274c77", "#ee6c4d","#8b8c89"))+
  labs(x = "PC1 (21.6%)", y = "PC2 (8.1%)") +
  my_theme +
  theme(legend.position = "bottom")


plot_pc1_pc2


plot_pc1_pc3 = ggplot() + 
  geom_point(mapping = aes(x = PC1, 
                           y = PC3, 
                           color = type),
             data = scores_pca,
             size = 2,
             alpha = 0.75) + 
  scale_color_manual(values=c("#274c77", "#ee6c4d","#8b8c89"))+
  labs(x = "PC1 (21.6%)", y = "PC3 (6.9%)") +
  my_theme +
  theme(legend.position = "right")

# PC3 (11.8%)
plot_pc1_pc3

joined_plot3 = plot_pc1_pc2 + plot_pc1_pc3
joined_plot3[[1]] = joined_plot3[[1]] + theme(legend.position = "none") 

joined_plot3 = joined_plot3 + plot_annotation(tag_levels = "A")

# ggsave(joined_plot3,  filename = "results/modelling_results/230704_PCA_targets_suspects.png", width=18, height=10, units = "cm", device = NULL)
# ggsave(joined_plot3, filename = "results/modelling_results/230704_PCA_targets_suspects.svg", width=18, height=10, units = "cm")



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

##---- t-SNE ----

data_for_tsne = SMILES_all %>%
  unique()

data_matrix  = as.matrix(data_for_tsne %>% select(-c(SMILES, Compound, type)))
tsne_out <- Rtsne(data_matrix)

# Conversion of matrix to dataframe
tsne_plot <- data.frame(x = tsne_out$Y[,1],
                        y = tsne_out$Y[,2]) %>% 
  bind_cols(type = data_for_tsne$type)

# Plotting the plot using ggplot() function
tsne_plot = ggplot2::ggplot(tsne_plot)+ geom_point(aes(x=x,y=y, color = type), alpha = 0.6, size = 2.5) + 
  scale_color_manual(values=c("#a3cef1","#274c77", "#ee6c4d" )) + 
  my_theme + theme(legend.position = "bottom")

# ggsave(tsne_plot,  filename = "results/modelling_results/230704_tsne.png", width=12, height=13, units = "cm", device = NULL)
# ggsave(tsne_plot, filename = "results/modelling_results/230704_tsne.svg", width=12, height=13, units = "cm")






