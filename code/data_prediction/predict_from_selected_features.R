#Load libraries
library(data.table)
library(reshape2)
library(ggplot2)
library(matrixStats)
library(ranger)
library(missRanger)
library(kernlab)
library(ropls)
library(pROC)
library(caTools)
library(AUC)

#Import central functions
source("../function_definitions.R")

#Set path
out_dir <- "../../results"
out_dir_pred_metab <- "../../results/data_pred_metab"
out_dir_pred <- "../../results/data_pred_featsel"

#Make sure paths exist
if (!dir.exists(out_dir))
  dir.create(out_dir)
if (!dir.exists(out_dir_pred))
  dir.create(out_dir_pred)

#Import data
##Import clinical data
human_sepsis_data <- get_human_sepsis_data()

##Import experiment data
rat_sepsis_data <- get_rat_sepsis_data()

##Import viterbi path data
viterbi_path <- fread(input = "../HMM/tram_viterbi_paths_all_septic_samples_2LJHMM.csv", header = TRUE, data.table = FALSE)

##Filter out patients without sepsis
human_nonsepsis_data <- human_sepsis_data[human_sepsis_data$`CAP / FP` == "-", ]
human_sepsis_data <- human_sepsis_data[human_sepsis_data$`CAP / FP` != "-", ]

#Get significantly different classes
human_sig_diff_res <- human_sig_diffs_along_days(human_sepsis_data, corr_fdr = FALSE)
day_sig_u_diff <- human_sig_diff_res$day_sig_u_diff
day_sig_t_diff <- human_sig_diff_res$day_sig_t_diff

#Get sig vars
sig_u_class <- na.omit(colnames(day_sig_u_diff[,-1])[colAnys(day_sig_u_diff[, -1] <= 0.05)])
sig_t_class <- na.omit(colnames(day_sig_t_diff[,-1])[colAnys(day_sig_t_diff[, -1] <= 0.05)])

#Try prediction
##On human
###Build data set for "vanilla" prediction
human_sepsis_data_ml <- subset(x = human_sepsis_data, select = c(2,4,6:ncol(human_sepsis_data)))
#strip_start <- strip_start <- which(colnames(human_sepsis_data) == "Urea")
#human_sepsis_data_ml <- subset(human_sepsis_data_ml, select = -which(colnames(human_sepsis_data_ml) %in% colnames(human_sepsis_data)[strip_start:ncol(human_sepsis_data)]))
human_sepsis_data_ml$Survival <- as.numeric(as.factor(human_sepsis_data_ml$Survival)) - 1 #Dependent variable transformation
###Add CAP and FP as independent variables
#human_sepsis_data_ml <- cbind(human_sepsis_data_ml[,1:2], data.frame(CAP = human_sepsis_data$`CAP / FP` == "CAP"), data.frame(FP = human_sepsis_data$`CAP / FP` == "FP"), human_sepsis_data_ml[,-1:-2])
###Strip phenomenological variables
#strip_start <- which(colnames(human_sepsis_data_ml) == "Urea")
#human_sepsis_data_ml <- human_sepsis_data_ml[,-strip_start:-ncol(human_sepsis_data_ml)]
human_sepsis_data_ml <- subset(human_sepsis_data_ml, select = c(colnames(human_sepsis_data_ml)[1:2], sig_t_class))
###No, strip metabolic variables
#human_sepsis_data_ml <- human_sepsis_data_ml[,c(1:2,strip_start:ncol(human_sepsis_data_ml))]
data_ml_subset <- apply(human_sepsis_data_ml, 1, function(x){ sum(!is.na(x)) > length(x) - 9})
human_sepsis_data_ml <- human_sepsis_data_ml[data_ml_subset,]
#human_sepsis_data_ml$`CAP / FP` <- as.factor(human_sepsis_data$`CAP / FP`)
colnames(human_sepsis_data_ml) <- make.names(colnames(human_sepsis_data_ml))
###Impute missing values
human_sepsis_data_ml_full <- missRanger(human_sepsis_data_ml, pmm.k = 3, num.trees = 100)
###Scale values
human_sepsis_data_ml_full[, -1:-2] <- scale(human_sepsis_data_ml_full[, -1:-2])
##Set important parameters
fml <- Survival ~ .
num_folds <- 5
rg.num.trees <- 1000
##Set up iteration
num_repeats <- 20
day_tab <- table(human_sepsis_data$Day[data_ml_subset])
day_set <- rep(as.numeric(names(day_tab)), each = num_repeats)
var_range <- 2:30
var_set <- rep(var_range, each = num_repeats)
var_set_name_list <- list()
rg.npr.repeat.df <- data.frame(var = var_set, tpr = 0, tnr = 0, fpr = 0, fnr = 0)
ks.npr.repeat.df <- data.frame(var = var_set, tpr = 0, tnr = 0, fpr = 0, fnr = 0)
lm.npr.repeat.df <- data.frame(var = var_set, tpr = 0, tnr = 0, fpr = 0, fnr = 0)
ol.npr.repeat.df <- data.frame(var = var_set, tpr = 0, tnr = 0, fpr = 0, fnr = 0)
rg.red.npr.repeat.df <- data.frame(var = var_set, tpr = 0, tnr = 0, fpr = 0, fnr = 0)
ks.red.npr.repeat.df <- data.frame(var = var_set, tpr = 0, tnr = 0, fpr = 0, fnr = 0)
lm.red.npr.repeat.df <- data.frame(var = var_set, tpr = 0, tnr = 0, fpr = 0, fnr = 0)
ol.red.npr.repeat.df <- data.frame(var = var_set, tpr = 0, tnr = 0, fpr = 0, fnr = 0)
rg.red.auc.repeat.df <- data.frame(var = var_set, auc = 0)
ks.red.auc.repeat.df <- data.frame(var = var_set, auc = 0)
lm.red.auc.repeat.df <- data.frame(var = var_set, auc = 0)
ol.red.auc.repeat.df <- data.frame(var = var_set, auc = 0)
##Do the big CV loop
d <- 1
{
  ###Select data with Day <= x, remove Patient ID
  human_sepsis_data_ml <- subset(human_sepsis_data_ml_full, subset = human_sepsis_data$Day[data_ml_subset] <= day_set[d], select = -1)
  ###Get Patient ID seperately
  patient_number <- subset(human_sepsis_data_ml_full, subset = human_sepsis_data$Day[data_ml_subset] <= day_set[d], select = Patient)[[1]]
  ###Test performance with cross-validation
  v.rg.npr <- list()
  v.ks.npr <- list()
  v.lm.npr <- list()
  v.ol.npr <- list()
  v.rg.auc <- list()
  v.ks.auc <- list()
  v.lm.auc <- list()
  yPredRG <- rep(0, nrow(human_sepsis_data_ml))
  yPredKS <- rep(0, nrow(human_sepsis_data_ml))
  yPredLM <- rep(0, nrow(human_sepsis_data_ml))
  #fold_set <- ml.split.folds.quasistrat(num_folds = num_folds, class = human_sepsis_data_ml$Survival, non_strat_class = patient_number)
  fold_set <- ml.split.folds.strat(num_folds = num_folds, class = human_sepsis_data_ml$Survival)
  for (fold in 1:num_folds){
    fold_learn_set <- human_sepsis_data_ml[-fold_set[[fold]],]
    fold_test_set <- human_sepsis_data_ml[fold_set[[fold]],]
    r_fold_learn_set <- fold_learn_set
    r_fold_test_set <- fold_test_set
    r_fold_learn_set$Survival <- as.factor(r_fold_learn_set$Survival)
    r_fold_test_set$Survival <- as.factor(r_fold_test_set$Survival)
    rgCV <- ranger(data = r_fold_learn_set, dependent.variable.name = "Survival", num.trees = rg.num.trees, write.forest = T, save.memory = F, probability = TRUE)
    ksCV <- ksvm(fml, fold_learn_set, type = "C-svc", kernel = "vanilladot", scaled = FALSE)
    lmCV <- lm(formula = fml, data = fold_learn_set)
    olCV <- opls(x = fold_learn_set[, -1], y = fold_learn_set[, 1], predI = 1, orthoI = 1, printL = F, plotL = F)
    v.rg.npr[[fold]] <- ml.npr(predict(rgCV, r_fold_test_set)$predictions[, 2] > 0.5, fold_test_set$Survival)
    v.ks.npr[[fold]] <- ml.npr(predict(ksCV, fold_test_set), fold_test_set$Survival)
    v.lm.npr[[fold]] <- ml.npr(predict(lmCV, fold_test_set[,-1]), fold_test_set$Survival)
    v.ol.npr[[fold]] <- ml.npr(predict(olCV, fold_test_set[,-1]), fold_test_set$Survival)
    v.rg.auc[[fold]] <- ml.auc(ref = fold_test_set$Survival, conf = predict(rgCV, r_fold_test_set)$predictions[, 2])
    v.rg.auc[[fold]] <- ml.auc(ref = fold_test_set$Survival, conf = predict(ksCV, fold_test_set, type = "decision"))
    v.lm.auc[[fold]] <- ml.auc(ref = fold_test_set$Survival, conf = predict(lmCV, fold_test_set[,-1])) #Can not compute ROC for OPLS-DA because prediction confidence is not available
    yPredRG[fold_set[[fold]]] <- predict(rgCV, fold_test_set)$predictions
    yPredKS[fold_set[[fold]]] <- predict(ksCV, fold_test_set)
    yPredLM[fold_set[[fold]]] <- predict(lmCV, fold_test_set)
  }
  rg.npr.repeat.df[d,-1] <- colMeans(rbindlist(v.rg.npr))
  ks.npr.repeat.df[d,-1] <- colMeans(rbindlist(v.ks.npr))
  lm.npr.repeat.df[d,-1] <- colMeans(rbindlist(v.lm.npr))
  ol.npr.repeat.df[d,-1] <- colMeans(rbindlist(v.ol.npr))

  for (v in seq_along(var_set)){
    ###Get variable importance and validate on special data subsets
    rg <- ranger(data = human_sepsis_data_ml, dependent.variable.name = "Survival", num.trees = rg.num.trees * 3, write.forest = T, save.memory = F, classification = TRUE, importance = "permutation")
    var_importance <- rg$variable.importance
    ###Reduce data set to important variables
    human_sepsis_data_ml_red <- human_sepsis_data_ml
    human_sepsis_data_ml_red[, 1 + which(var_importance < sort(var_importance, decreasing = TRUE)[var_set[v]])] <- NULL
    var_set_name_list[[v]] <- names(var_importance[var_importance < sort(var_importance, decreasing = TRUE)[var_set[v]]])
    fold_set <- ml.split.folds.strat(num_folds = num_folds, class = human_sepsis_data_ml_red$Survival)
    v.rg.npr <- list()
    v.ks.npr <- list()
    v.lm.npr <- list()
    v.ol.npr <- list()
    v.rg.auc <- list()
    v.ks.auc <- list()
    v.lm.auc <- list()
    yPredRG <- rep(0, nrow(human_sepsis_data_ml_red))
    yPredKS <- rep(0, nrow(human_sepsis_data_ml_red))
    yPredLM <- rep(0, nrow(human_sepsis_data_ml_red))
    for (fold in 1:num_folds){
      fold_learn_set <- human_sepsis_data_ml_red[-fold_set[[fold]],]
      fold_test_set <- human_sepsis_data_ml_red[fold_set[[fold]],]
      r_fold_learn_set <- fold_learn_set
      r_fold_test_set <- fold_test_set
      r_fold_learn_set$Survival <- as.factor(r_fold_learn_set$Survival)
      r_fold_test_set$Survival <- as.factor(r_fold_test_set$Survival)
      rgCV <- ranger(data = r_fold_learn_set, dependent.variable.name = "Survival", num.trees = rg.num.trees, write.forest = T, save.memory = F, classification = TRUE, probability = TRUE)
      ksCV <- ksvm(fml, fold_learn_set, type = "C-svc", kernel = "vanilladot", scaled = FALSE)
      lmCV <- lm(formula = fml, data = fold_learn_set)
      olCV <- opls(x = fold_learn_set[, -1], y = fold_learn_set[, 1], predI = 1, orthoI = 1, printL = F, plotL = F)
      v.rg.npr[[fold]] <- ml.npr(predict(rgCV, r_fold_test_set)$predictions[, 2] > 0.5, fold_test_set$Survival)
      v.ks.npr[[fold]] <- ml.npr(predict(ksCV, fold_test_set), fold_test_set$Survival)
      v.lm.npr[[fold]] <- ml.npr(predict(lmCV, fold_test_set[,-1]), fold_test_set$Survival)
      v.ol.npr[[fold]] <- ml.npr(predict(olCV, fold_test_set[,-1]), fold_test_set$Survival)
      v.rg.auc[[fold]] <- ml.auc(ref = fold_test_set$Survival, conf = predict(rgCV, r_fold_test_set)$predictions[, 2])
      v.ks.auc[[fold]] <- ml.auc(ref = fold_test_set$Survival, conf = predict(ksCV, fold_test_set, type = "decision"))
      v.lm.auc[[fold]] <- ml.auc(ref = fold_test_set$Survival, conf = predict(lmCV, fold_test_set[,-1])) #Can not compute ROC for OPLS-DA because prediction confidence is not available
      yPredRG[fold_set[[fold]]] <- predict(rgCV, human_sepsis_data_ml_red[fold_set[[fold]],])$predictions
      yPredKS[fold_set[[fold]]] <- predict(ksCV, human_sepsis_data_ml_red[fold_set[[fold]],])
      yPredLM[fold_set[[fold]]] <- predict(lmCV, human_sepsis_data_ml_red[fold_set[[fold]],])
    }
    rg.red.npr.repeat.df[v,-1] <- colMeans(rbindlist(v.rg.npr))
    ks.red.npr.repeat.df[v,-1] <- colMeans(rbindlist(v.ks.npr))
    lm.red.npr.repeat.df[v,-1] <- colMeans(rbindlist(v.lm.npr))
    ol.red.npr.repeat.df[v,-1] <- colMeans(rbindlist(v.ol.npr))
    
    rg.red.auc.repeat.df[v,-1] <- mean(unlist(v.rg.auc))
    ks.red.auc.repeat.df[v,-1] <- mean(unlist(v.ks.auc))
    lm.red.auc.repeat.df[v,-1] <- mean(unlist(v.lm.auc))
  }
}

rg.npr.repeat.df$Method <- "RF"
ks.npr.repeat.df$Method <- "SVM"
lm.npr.repeat.df$Method <- "LM"
ol.npr.repeat.df$Method <- "OPLSDA"
rg.red.npr.repeat.df$Method <- "RF"
ks.red.npr.repeat.df$Method <- "SVM"
lm.red.npr.repeat.df$Method <- "LM"
ol.red.npr.repeat.df$Method <- "OPLSDA"

rg.red.auc.repeat.df$Method <- "RF"
ks.red.auc.repeat.df$Method <- "SVM"
lm.red.auc.repeat.df$Method <- "LM"

#Calculate accuracy
rg.npr.repeat.df$accuracy <- (rg.npr.repeat.df$tpr + rg.npr.repeat.df$tnr) / rowSums(rg.npr.repeat.df[, c(-1,-6)])
ks.npr.repeat.df$accuracy <- (ks.npr.repeat.df$tpr + ks.npr.repeat.df$tnr) / rowSums(ks.npr.repeat.df[, c(-1,-6)])
lm.npr.repeat.df$accuracy <- (lm.npr.repeat.df$tpr + lm.npr.repeat.df$tnr) / rowSums(lm.npr.repeat.df[, c(-1,-6)])
ol.npr.repeat.df$accuracy <- (ol.npr.repeat.df$tpr + ol.npr.repeat.df$tnr) / rowSums(ol.npr.repeat.df[, c(-1,-6)])
rg.red.npr.repeat.df$accuracy <- (rg.red.npr.repeat.df$tpr + rg.red.npr.repeat.df$tnr) / rowSums(rg.red.npr.repeat.df[, c(-1,-6)])
ks.red.npr.repeat.df$accuracy <- (ks.red.npr.repeat.df$tpr + ks.red.npr.repeat.df$tnr) / rowSums(ks.red.npr.repeat.df[, c(-1,-6)])
lm.red.npr.repeat.df$accuracy <- (lm.red.npr.repeat.df$tpr + lm.red.npr.repeat.df$tnr) / rowSums(lm.red.npr.repeat.df[, c(-1,-6)])
ol.red.npr.repeat.df$accuracy <- (ol.red.npr.repeat.df$tpr + ol.red.npr.repeat.df$tnr) / rowSums(ol.red.npr.repeat.df[, c(-1,-6)])

#Melt to long table
npr_by_day_data <- melt(rbindlist(list(rg.npr.repeat.df, ks.npr.repeat.df, lm.npr.repeat.df)), id.vars = c("var", "Method"))
npr_red_by_day_data <- melt(rbindlist(list(rg.red.npr.repeat.df, ks.red.npr.repeat.df, lm.red.npr.repeat.df)), id.vars = c("var", "Method"))
auc_red_data <- melt(rbindlist(list(rg.red.auc.repeat.df, ks.red.auc.repeat.df, lm.red.auc.repeat.df)), id.vars = c("var", "Method"))

##Remove redundancy
npr_by_day_data <- subset(npr_by_day_data, subset = variable %in% c("fpr", "fnr", "accuracy"))
npr_red_by_day_data <- subset(npr_red_by_day_data, subset = variable %in% c("fpr", "fnr", "accuracy"))

#Combine for nice plot
npr_by_day_data$dimensions <- "All variables"
npr_red_by_day_data$dimensions <- "Important variables"
npr_data <- rbind(npr_by_day_data, npr_red_by_day_data)

#Plot FPR and FNR
npr_by_feat_num_plot <- ggplot(subset(npr_data, variable %in% c("fnr", "fpr")), aes(x = var, y = value, color = Method)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", size = 0.5) +
  facet_grid(variable ~ dimensions) +
  ylim(0, 1) +
  xlab("Features included") +
  ggtitle("") +
  theme_bw()
ggsave(filename = "npr_by_feat_num.png", plot = npr_by_feat_num_plot, path = out_dir_pred, width = 5, height = 5, units = "in")

#Plot accuracy
acc_by_feat_num_plot <- ggplot(subset(npr_data, variable == "accuracy"), aes(x = var, y = value, color = Method)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean") +
  facet_grid( ~ dimensions) +
  ylim(0, 1) +
  xlab("Features included") +
  ylab("Accuracy") +
  #ggtitle("Normal stratified Cross Validation may shows overfitting") +
  #ggtitle("Quasi-stratified Cross Validation prevents overfitting") +
  theme_bw()
ggsave(filename = "acc_by_feat_num.png", plot = acc_by_feat_num_plot, path = out_dir_pred, width = 6, height = 4, units = "in")

#Plot AUC
auc_by_feat_num_plot <- ggplot(subset(auc_red_data, variable == "auc"), aes(x = var, y = value, color = Method)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean") +
  #facet_grid( ~ dimensions) +
  ylim(0, 1) +
  xlab("Features included") +
  ylab("AUC") +
  #ggtitle("Normal stratified Cross Validation may shows overfitting") +
  #ggtitle("Quasi-stratified Cross Validation prevents overfitting") +
  theme_bw()
ggsave(filename = "auc_by_feat_num.png", plot = auc_by_feat_num_plot, path = out_dir_pred, width = 6, height = 3.5, units = "in")
