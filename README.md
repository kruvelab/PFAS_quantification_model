# PFAS quantification model
This GitHub repository includes all supplementary data, codes, and models for the article:

**"Closing the organofluorine mass balance in marine mammals using suspect screening and machine learning-based quantification"**
by Mélanie Z. Lauria, Helen Sepman, Thomas Ledbetter, Merle Plassmann, Anna Roos, Malene Simon, Jonathan P. Benskin, and Anneli Kruve.

An ionization efficiency prediction model was trained on a joined dataset of 33 target PFAS measured in this work and 100 chemicals measured in negative mode by Liigand et al. (01_data_for_modelling/IE_training_data). For training the models, Pharmaceutical Data Exploration Laboratory (PaDEL)​​ molecular descriptors and eluent descriptors (aqueous pH, polarity index, viscosity, surface tension and presence of NH4) were used. The RMSEs over all data points in training and test sets were 0.43 and 0.79 logIE units, respectively. This corresponds to fold errors of 2.7× and 6.2×. 

### How can (tentatively) identified chemicals be quantified with the ionization efficiency prediction approach?
The predicted logIE indicates how the response factor of a respective chemical relates to other chemicals and is instrument independent. However, response factors are instrument and lab-specific, meaning that the magnitude of RF depends on the source design, instrument vendor and to some extent even on the software used for data integration. Therefore, predicted ionization efficiency needs to be converted to a measurement-specific response factor before quantification. For this, calibrants need to be used (33 target PFAS were used here) and they are measured together with the suspects and their response factors (slopes of the calibration graphs) are calculated. For these calibrants, a calibration graph between measurement-specific response factors and predicted ionization efficiencies is constructed. Using this calibration graph, the logIE is predicted for suspect chemicals, converted into a predicted response factor and used to calculate the concentration (Figure 1).  

![modelling_workflow](https://github.com/kruvelab/PFAS_quantification_model/assets/48623628/7ec7fba3-979d-4813-a077-798b7644b376)
Figure 1. Workflow to predict concentration using predicted log*IE* approach.

### How to use the PFAS quantification model developed here?
All developed models can be accessed on this GitHub page. The work is ongoing to implement the final model in an R-package *MS2Quant* for easier accessibility. 

For now, to use the quantification model developed, the file *"05_suspect_quantification/suspect_quantification.R"* can be used as a template. For predicting the concentrations the following inputs are needed:
* cal_filename_data - a TraceFinder file with integrated signal areas and theoretical amounts specified in the file for creating calibration curves.
* cal_filename_smiles - SMILES for calibrants; the Excel sheet name of the calibrant should be provided as *ID* for merging the information
* sus_filename_data - a TraceFinder file with integrated signal areas for the suspects
* sus_filename_smiles - SMILES for (tentatively) identified suspects; the Excel sheet name of the suspect should be provided as *ID* for merging the information
* logIE_pred_model - a final model that can be read in from "03_models/230619_logIE_model_withPFAS_allData.Rdata"
* filename_eluent = a .csv file with the gradient program for calculating eluent descriptors; see "05_suspect_quantification/eluent.csv"
* organic_modifier - organic modifier used (either "MeCN" or "MeOH"; default "MeCN")
* pH.aq. - pH of the aqueous phase; default "7.0"
* NH4 - binary; presence (1) or absence (0) of NH4 if additive containing NH4 was used; default "1"
* compounds_to_be_removed_as_list - optional; if some chemicals must be removed from the analysis.
