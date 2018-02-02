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

##Import clinical validation data
human_validation_data <- get_Ferrario_validation_data()

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
human_sepsis_data_ml <- subset(human_sepsis_data_ml, select = c(colnames(human_sepsis_data_ml)[1:2], intersect(sig_t_class, colnames(human_validatio_data)))
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
num_folds <- 4
rg.num.trees <- 1000
##Set up iteration
num_repeats <- 20
day_tab <- table(human_sepsis_data$Day[data_ml_subset])
day_set <- rep(as.numeric(names(day_tab)), each = num_repeats)
var_range <- 2:22
var_set <- rep(var_range, each = num_repeats)
var_set_name_list <- list()
tot_n_var <- ncol(human_sepsis_data_ml_full)-2
rg.npr.repeat.df <- data.frame(var = rep(tot_n_var, num_repeats), tpr = 0, tnr = 0, fpr = 0, fnr = 0, ppv = 0, npv = 0)
ks.npr.repeat.df <- data.frame(var = rep(tot_n_var, num_repeats), tpr = 0, tnr = 0, fpr = 0, fnr = 0, ppv = 0, npv = 0)
lm.npr.repeat.df <- data.frame(var = rep(tot_n_var, num_repeats), tpr = 0, tnr = 0, fpr = 0, fnr = 0, ppv = 0, npv = 0)
ol.npr.repeat.df <- data.frame(var = rep(tot_n_var, num_repeats), tpr = 0, tnr = 0, fpr = 0, fnr = 0, ppv = 0, npv = 0)
rg.acc.repeat.df <- data.frame(var = rep(tot_n_var, num_repeats), acc = 0)
ks.acc.repeat.df <- data.frame(var = rep(tot_n_var, num_repeats), acc = 0)
lm.acc.repeat.df <- data.frame(var = rep(tot_n_var, num_repeats), acc = 0)
ol.acc.repeat.df <- data.frame(var = rep(tot_n_var, num_repeats), acc = 0)
rg.auc.repeat.df <- data.frame(var = rep(tot_n_var, num_repeats), auc = 0)
ks.auc.repeat.df <- data.frame(var = rep(tot_n_var, num_repeats), auc = 0)
lm.auc.repeat.df <- data.frame(var = rep(tot_n_var, num_repeats), auc = 0)
ol.auc.repeat.df <- data.frame(var = rep(tot_n_var, num_repeats), auc = 0)
rg.red.npr.repeat.df <- data.frame(var = var_set, tpr = 0, tnr = 0, fpr = 0, fnr = 0, ppv = 0, npv = 0)
ks.red.npr.repeat.df <- data.frame(var = var_set, tpr = 0, tnr = 0, fpr = 0, fnr = 0, ppv = 0, npv = 0)
lm.red.npr.repeat.df <- data.frame(var = var_set, tpr = 0, tnr = 0, fpr = 0, fnr = 0, ppv = 0, npv = 0)
ol.red.npr.repeat.df <- data.frame(var = var_set, tpr = 0, tnr = 0, fpr = 0, fnr = 0, ppv = 0, npv = 0)
rg.red.acc.repeat.df <- data.frame(var = var_set, acc = 0)
ks.red.acc.repeat.df <- data.frame(var = var_set, acc = 0)
lm.red.acc.repeat.df <- data.frame(var = var_set, acc = 0)
ol.red.acc.repeat.df <- data.frame(var = var_set, acc = 0)
rg.red.auc.repeat.df <- data.frame(var = var_set, auc = 0)
ks.red.auc.repeat.df <- data.frame(var = var_set, auc = 0)
lm.red.auc.repeat.df <- data.frame(var = var_set, auc = 0)
ol.red.auc.repeat.df <- data.frame(var = var_set, auc = 0)
##Do the big CV loop
d <- 1
r_count <- 1
{
  ###Select data with Day <= x, remove Patient ID ("select = -1")
  human_sepsis_data_ml <- subset(human_sepsis_data_ml_full, subset = human_sepsis_data$Day[data_ml_subset] <= day_set[d], select = -1)
  ###Get Patient ID seperately
  patient_number <- subset(human_sepsis_data_ml_full, subset = human_sepsis_data$Day[data_ml_subset] <= day_set[d], select = Patient)[[1]]
  ###Repeat num_repeat times
  for (v in 1:num_repeats){
    ###Test performance with cross-validation
    v.rg.npr <- list()
    v.ks.npr <- list()
    v.lm.npr <- list()
    v.ol.npr <- list()
    v.rg.acc <- list()
    v.ks.acc <- list()
    v.lm.acc <- list()
    v.ol.acc <- list()
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
      lmCV <- glm(formula = fml, data = fold_learn_set, family = binomial(link = "probit"))
      #olCV <- opls(x = fold_learn_set[, -1], y = fold_learn_set[, 1], predI = 1, orthoI = 1, printL = F, plotL = F)
      v.rg.npr[[fold]] <- ml.npr(predict(rgCV, r_fold_test_set)$predictions[, 2] > 0.5, fold_test_set$Survival)
      v.ks.npr[[fold]] <- ml.npr(predict(ksCV, fold_test_set), fold_test_set$Survival)
      v.lm.npr[[fold]] <- ml.npr(predict(lmCV, fold_test_set[,-1]), fold_test_set$Survival)
      #v.ol.npr[[fold]] <- ml.npr(predict(olCV, fold_test_set[,-1]), fold_test_set$Survival)
      v.rg.acc[[fold]] <- ml.acc(predict(rgCV, r_fold_test_set)$predictions[, 2] > 0.5, fold_test_set$Survival)
      v.ks.acc[[fold]] <- ml.acc(predict(ksCV, fold_test_set), fold_test_set$Survival)
      v.lm.acc[[fold]] <- ml.acc(predict(lmCV, fold_test_set[,-1]), fold_test_set$Survival)
      #v.ol.acc[[fold]] <- ml.acc(predict(olCV, fold_test_set[,-1]), fold_test_set$Survival)
      v.rg.auc[[fold]] <- ml.auc(ref = fold_test_set$Survival, conf = predict(rgCV, r_fold_test_set)$predictions[, 2])
      v.rg.auc[[fold]] <- ml.auc(ref = fold_test_set$Survival, conf = predict(ksCV, fold_test_set, type = "decision"))
      v.lm.auc[[fold]] <- ml.auc(ref = fold_test_set$Survival, conf = predict(lmCV, fold_test_set[,-1])) #Can not compute ROC for OPLS-DA because prediction confidence is not available
      yPredRG[fold_set[[fold]]] <- predict(rgCV, fold_test_set)$predictions
      yPredKS[fold_set[[fold]]] <- predict(ksCV, fold_test_set)
      yPredLM[fold_set[[fold]]] <- predict(lmCV, fold_test_set)
    }
    rg.npr.repeat.df[v,-1] <- colMeans(rbindlist(v.rg.npr))
    ks.npr.repeat.df[v,-1] <- colMeans(rbindlist(v.ks.npr))
    lm.npr.repeat.df[v,-1] <- colMeans(rbindlist(v.lm.npr))
    #ol.npr.repeat.df[v,-1] <- colMeans(rbindlist(v.ol.npr))
    
    rg.acc.repeat.df[v,-1] <- mean(unlist(v.rg.acc))
    ks.acc.repeat.df[v,-1] <- mean(unlist(v.ks.acc))
    lm.acc.repeat.df[v,-1] <- mean(unlist(v.lm.acc))
    #ol.acc.repeat.df[v,-1] <- mean(unlist(v.ol.acc))
    
    rg.auc.repeat.df[v,-1] <- mean(unlist(v.rg.auc))
    ks.auc.repeat.df[v,-1] <- mean(unlist(v.ks.auc))
    lm.auc.repeat.df[v,-1] <- mean(unlist(v.lm.auc))
  }

  ###Repeat CV with selected features
  for (v in seq_along(var_set)){
    human_sepsis_data_ml_red <- human_sepsis_data_ml #kept as dummy to avoid variable renaming
    rg <- ranger(data = human_sepsis_data_ml, dependent.variable.name = "Survival", num.trees = rg.num.trees * 3, write.forest = T, save.memory = F, classification = TRUE, importance = "permutation")
    var_importance <- rg$variable.importance
    ###Reduce data set to important variables
    human_sepsis_data_ml_red[, 1 + which(var_importance < sort(var_importance, decreasing = TRUE)[var_set[v]])] <- NULL
    ###Save variable set for later analysis
    var_set_name_list[[v]] <- colnames(human_sepsis_data_ml_red)[-1]
    fold_set <- ml.split.folds.strat(num_folds = num_folds, class = human_sepsis_data_ml$Survival)
    v.rg.npr <- list()
    v.ks.npr <- list()
    v.lm.npr <- list()
    v.ol.npr <- list()
    v.rg.acc <- list()
    v.ks.acc <- list()
    v.lm.acc <- list()
    v.ol.acc <- list()
    v.rg.auc <- list()
    v.ks.auc <- list()
    v.lm.auc <- list()
    yPredRG <- rep(0, nrow(human_sepsis_data_ml))
    yPredKS <- rep(0, nrow(human_sepsis_data_ml))
    yPredLM <- rep(0, nrow(human_sepsis_data_ml))
    for (fold in 1:num_folds){
      ###Build learn and test set
      fold_learn_set <- human_sepsis_data_ml_red[-fold_set[[fold]],]
      fold_test_set <- human_sepsis_data_ml_red[fold_set[[fold]],]
      ###Get variable importance and validate on special data subsets
      r_count <- r_count + 1
      ####Continue with the other stuff
      r_fold_learn_set <- fold_learn_set
      r_fold_test_set <- fold_test_set
      r_fold_learn_set$Survival <- as.factor(r_fold_learn_set$Survival)
      r_fold_test_set$Survival <- as.factor(r_fold_test_set$Survival)
      rgCV <- ranger(data = r_fold_learn_set, dependent.variable.name = "Survival", num.trees = rg.num.trees, write.forest = T, save.memory = F, classification = TRUE, probability = TRUE)
      ksCV <- ksvm(fml, fold_learn_set, type = "C-svc", kernel = "vanilladot", scaled = FALSE)
      lmCV <- glm(formula = fml, data = fold_learn_set, family = binomial(link = "probit"))
      #olCV <- opls(x = fold_learn_set[, -1], y = fold_learn_set[, 1], predI = 1, orthoI = 1, printL = F, plotL = F)
      v.rg.npr[[fold]] <- ml.npr(predict(rgCV, r_fold_test_set)$predictions[, 2] > 0.5, fold_test_set$Survival)
      v.ks.npr[[fold]] <- ml.npr(predict(ksCV, fold_test_set), fold_test_set$Survival)
      v.lm.npr[[fold]] <- ml.npr(predict(lmCV, fold_test_set[,-1]), fold_test_set$Survival)
      #v.ol.npr[[fold]] <- ml.npr(predict(olCV, fold_test_set[,-1]), fold_test_set$Survival)
      v.rg.acc[[fold]] <- ml.acc(predict(rgCV, r_fold_test_set)$predictions[, 2] > 0.5, fold_test_set$Survival)
      v.ks.acc[[fold]] <- ml.acc(predict(ksCV, fold_test_set), fold_test_set$Survival)
      v.lm.acc[[fold]] <- ml.acc(predict(lmCV, fold_test_set[,-1]), fold_test_set$Survival)
      #v.ol.acc[[fold]] <- ml.acc(predict(olCV, fold_test_set[,-1]), fold_test_set$Survival)
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
    #ol.red.npr.repeat.df[v,-1] <- colMeans(rbindlist(v.ol.npr))
    
    rg.red.acc.repeat.df[v,-1] <- mean(unlist(v.rg.acc))
    ks.red.acc.repeat.df[v,-1] <- mean(unlist(v.ks.acc))
    lm.red.acc.repeat.df[v,-1] <- mean(unlist(v.lm.acc))
    #ol.red.acc.repeat.df[v,-1] <- mean(unlist(v.ol.acc))
    
    rg.red.auc.repeat.df[v,-1] <- mean(unlist(v.rg.auc))
    ks.red.auc.repeat.df[v,-1] <- mean(unlist(v.ks.auc))
    lm.red.auc.repeat.df[v,-1] <- mean(unlist(v.lm.auc))
  }
}

rg.npr.repeat.df$Method <- "RF"
ks.npr.repeat.df$Method <- "SVM"
lm.npr.repeat.df$Method <- "Probit"
ol.npr.repeat.df$Method <- "OPLSDA"
rg.red.npr.repeat.df$Method <- "RF"
ks.red.npr.repeat.df$Method <- "SVM"
lm.red.npr.repeat.df$Method <- "Probit"
ol.red.npr.repeat.df$Method <- "OPLSDA"

rg.auc.repeat.df$Method <- "RF"
ks.auc.repeat.df$Method <- "SVM"
lm.auc.repeat.df$Method <- "Probit"
rg.red.auc.repeat.df$Method <- "RF"
ks.red.auc.repeat.df$Method <- "SVM"
lm.red.auc.repeat.df$Method <- "Probit"

#Calculate accuracy
rg.npr.repeat.df$accuracy <- rg.acc.repeat.df$acc
ks.npr.repeat.df$accuracy <- ks.acc.repeat.df$acc
lm.npr.repeat.df$accuracy <- lm.acc.repeat.df$acc
ol.npr.repeat.df$accuracy <- ol.acc.repeat.df$acc
rg.red.npr.repeat.df$accuracy <- rg.red.acc.repeat.df$acc
ks.red.npr.repeat.df$accuracy <- ks.red.acc.repeat.df$acc
lm.red.npr.repeat.df$accuracy <- lm.red.acc.repeat.df$acc
ol.red.npr.repeat.df$accuracy <- ol.red.acc.repeat.df$acc

#Calculate Negative/Positive Predictive Value

#Melt to long table
npr_by_day_data <- melt(rbindlist(list(rg.npr.repeat.df, ks.npr.repeat.df, lm.npr.repeat.df)), id.vars = c("var", "Method"))
npr_red_by_day_data <- melt(rbindlist(list(rg.red.npr.repeat.df, ks.red.npr.repeat.df, lm.red.npr.repeat.df)), id.vars = c("var", "Method"))
auc_data <- melt(rbindlist(list(rg.auc.repeat.df, ks.auc.repeat.df, lm.auc.repeat.df)), id.vars = c("var", "Method"))
auc_red_data <- melt(rbindlist(list(rg.red.auc.repeat.df, ks.red.auc.repeat.df, lm.red.auc.repeat.df)), id.vars = c("var", "Method"))

##Remove redundancy
npr_by_day_data <- subset(npr_by_day_data, subset = variable %in% c("ppv", "npv", "accuracy"))
npr_red_by_day_data <- subset(npr_red_by_day_data, subset = variable %in% c("ppv", "npv", "accuracy"))

#Combine for nice plot
npr_by_day_data$dimensions <- "All variables"
npr_red_by_day_data$dimensions <- "Important variables"
npr_data <- rbind(npr_by_day_data, npr_red_by_day_data)
auc_data$dimensions <- "All variables"
auc_red_data$dimensions <- "Important variables"
auc_data_full <- rbind(auc_data, auc_red_data)

#Plot FPR and FNR
npr_by_feat_num_plot <- ggplot(subset(npr_data, variable %in% c("ppv", "npv")), aes(x = var, y = value, color = Method)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", size = 0.5) +
  facet_grid(variable ~ dimensions, scales = "free_x") +
  ylim(0, 1) +
  scale_x_discrete(limits = c(var_range[seq(1, length(var_range), 4)], tot_n_var)) + 
  xlab("Features included") +
  ggtitle("") +
  theme_bw()
ggsave(filename = "npr_by_feat_num.png", plot = npr_by_feat_num_plot, path = out_dir_pred, width = 8, height = 4, units = "in")

#Plot accuracy
acc_by_feat_num_plot <- ggplot(subset(npr_data, variable == "accuracy"), aes(x = var, y = value, color = Method)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean") +
  facet_grid( ~ dimensions, scales = "free_x") +
  ylim(0, 1) +
  scale_x_discrete(limits = c(var_range[seq(1, length(var_range), 4)], tot_n_var)) +
  xlab("Features included") +
  ylab("Accuracy") +
  #ggtitle("Normal stratified Cross Validation may shows overfitting") +
  #ggtitle("Quasi-stratified Cross Validation prevents overfitting") +
  theme_bw()
ggsave(filename = "acc_by_feat_num.png", plot = acc_by_feat_num_plot, path = out_dir_pred, width = 8, height = 2.5, units = "in")

#Plot AUC
auc_by_feat_num_plot <- ggplot(subset(auc_data_full, variable == "auc"), aes(x = var, y = value, color = Method)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean") +
  facet_grid( ~ dimensions, scales = "free_x") +
  scale_x_discrete(limits = c(var_range[seq(1, length(var_range), 4)], tot_n_var)) + 
  ylim(0, 1) +
  xlab("Features included") +
  ylab("AUC") +
  #ggtitle("Normal stratified Cross Validation may shows overfitting") +
  #ggtitle("Quasi-stratified Cross Validation prevents overfitting") +
  theme_bw()
ggsave(filename = "auc_by_feat_num.png", plot = auc_by_feat_num_plot, path = out_dir_pred, width = 8, height = 2.5, units = "in")

#Plot metabs against metabs for six most important variables
# selection <- unique(na.omit(unlist(var_set_name_list[1:10])))
# plot(as.formula(paste0("~ ", paste0(selection, collapse = "+"))), data = human_sepsis_data_ml[,selection], col = human_sepsis_data_ml$Survival + 1)
# 
# hsiclasso_set <- c('Leptin', 'PC aa C30:2', 'TSH', 'PC aa C30:0', 'FT4', 'C5-M-DC', 'C4', 'C16', 'GH', 't4-OH-Pro', 'PC aa C38:3', 'Prolactin', 'DHEA', 'Leu', 'SM C24:1')
# selection <- na.omit(match(hsiclasso_set, colnames(human_sepsis_data)))
# hsd <- subset(human_sepsis_data, Day == 0)
# colnames(hsd) <- make.names(colnames(hsd))
# plot(as.formula(paste0("~ ", paste0(colnames(hsd)[selection], collapse = "+"))), data = hsd[,selection], col = as.numeric(as.factor(hsd$Survival)) + 1)
