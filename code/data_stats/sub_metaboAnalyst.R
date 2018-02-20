#Script for metaboAnalyst

source("../function_definitions.R")

#data to character for input
#human_sepsis_data_char <- capture.output(write.csv(x = human_sepsis_data[,c(-1,-2,-5)], quote = TRUE))
#human_sepsis_data_char <- paste(human_sepsis_data_char, sep = "\n", collapse = "\n")
human_sepsis_data <- get_human_sepsis_data()
human_sepsis_data <- human_sepsis_data[human_sepsis_data$`CAP / FP` != "-",]
write.csv(x = subset(human_sepsis_data, Day == 0, c(-1:-3,-5)), file = "metaboAnalyst_sepsis_temp.csv")

human_sepsis_legend <- get_human_sepsis_legend()

#Init
mSet<-InitDataObjects("conc", "pathinteg", FALSE)
mSet<-Read.TextData(mSet, "metaboAnalyst_sepsis_temp.csv", "rowu", "disc") #might fail bc. of path requirement

#Set up mSetObj with the list of compounds
mSet<-Setup.MapData(mSet, colnames(human_sepsis_data)[-1:-5]);
mSet<-Setup.MapData(mSet, human_sepsis_legend[-1:-5,2])

#Sanity checks
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet)
mSet<-IsSmallSmplSize(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "MeanCenter", ref=NULL, ratio=FALSE, ratioNum=20)

###Overrepresentation erichment test

# Cross reference list of compounds against libraries (hmdb, pubchem, chebi, kegg, metlin)
mSet<-CrossReferencing(mSet, "name");

# Create the mapping results table
mSet<-CreateMappingResultTable(mSet)

# Set the metabolite filter
mSet<-SetMetabolomeFilter(mSet, T);

# Select metabolite set library, refer to 
mSet<-SetCurrentMsetLib(mSet, "pathway", 0);

# Calculate hypergeometric score, results table generated in your working directory
mSet<-CalculateHyperScore(mSet)

# Plot the ORA, bar-graph
mSet<-PlotORA(mSet, "ora_0_", "bar", "png", 72, width=NA)


###Biomarker identification

mSet<-SetAnalysisMode(mSet, "univ")

# Prepare data for biomarker analysis
mSet<-PrepareROCData(mSet)

# Perform univariate ROC curve analysis, 
mSet<-Perform.UnivROC(mSet, feat.nm = "C4", imgName = "C4", "png", dpi=300, isAUC=F, isOpt=T, optMethod="closest.topleft", isPartial=F, measure="sp", cutoff=0.2)

# Perform univariate ROC curve analysis, resulting in a partial AUC with a 95% CI band 
mSet<-Perform.UnivROC(mSet, feat.nm = "Tyr", imgName = "Tyr", "png", dpi=300, isAUC=T, isOpt=T, optMethod="closest.topleft", isPartial=T, measure="se", cutoff=0.2)

# Perform calculation of feature importance (AUC, p value, fold change)
mSet<-CalculateFeatureRanking(mSet, 4)


##multivariate biomarker identification

mSet<-SetAnalysisMode(mSet, "explore")

# Prepare data for biomarker analysis
mSet<-PrepareROCData(mSet)

# Perform multivariate ROC curve analysis, using SVM classification and ranking
mSet<-PerformCV.explore(mSet, cls.method = "rf", rank.method = "rf", lvNum = 2)

# Comparison plot of ROC curves of all models
mSet<-PlotROCMulti(mSet, imgName = "ROC_all_models", format = "png", dpi = 300, mdl.inx= 0, avg.method = "threshold", show.conf = 0, show.holdout = 0, focus="fpr", cutoff=0.5)

# Plot the ROC curve of a single selected model, in this case model 1 and display the confidence interval
mSet<-PlotROCMulti(mSet, imgName = "ROC_model1", format = "png", dpi = 300, mdl.inx = 1, avg.method = "threshold", show.conf = 1, 0, "fpr", 0.2)

# Plot predicted class probabilities for each sample for a selected model, not showing labels of wrongly classified samples
mSet<-PlotProbView(mSet, imgName = "multi_roc_prob", format = "png", dpi = 300, mdl.inx = -1, show = 0, showPred = 0)

# Plot the predictive accuracy of models with increasing number of features
mSet<-PlotAccuracy(mSet, imgName = "multi_roc_accuracy", format = "png", dpi = 300)

# Plot the most important features of a selected model ranked from most to least important
mSet<-PlotImpVars(mSet, imgName = "multi_roc_impvar", format="png", dpi=300, mdl.inx = -1, measure="freq", feat.num=15)


####Model selection

# Set the biomarker analysis mode to perform ROC Curve Based Model Creation and Evaluation ("test")
mSet<-SetAnalysisMode(mSet, "test")

# Prepare data for biomarker analysis
mSet<-PrepareROCData(mSet)

# Perform calculation of feature importance (AUC, p value, fold change)
mSet<-CalculateFeatureRanking(mSet)

# Manually select a subset of features for ROC analysis to build a classifier 
selected.cmpds <- c("C4", "Tyr")

# Manually select a subset of samples for ROC analysis hold-out data for validation purposes
selected.smpls <- c(rownames(mSet$dataSet$orig)[sample(nrow(mSet$dataSet$orig), 5)])

# Prepare the custom data for model creation and sample hold-out 
mSet<-SetCustomData(mSet, selected.cmpds, selected.smpls)

# Perform ROC curve analysis, using SVM classification 
mSet<-PerformCV.test(mSet, method = "rf", lvNum = 2)

# Plot the ROC curve for the created model
mSet<-PlotROCTest(mSet, imgName = "cls_roc_0_", format="png",  dpi=300, mdl.inx = 0, avg.method = "threshold", 0, 0, "fpr", 0.5)

# Plot the predicted class probabilities for each sample using the user-created classifier, not showing labels of wrongly classified samples
mSet<-PlotProbViewTest(mSet, imgName = "cls_prob_0_", format="png",  dpi=300, mdl.inx =-1, show=0, showPred= 0)

# Plot the predictive accuracy of the model with increasing number of features
mSet<-PlotAccuracy(mSet, imgName = "cls_accu_0_", format="png",  dpi=300)
