#Load libraries
library(data.table)
library(reshape2)
library(ggplot2)
library(matrixStats)
library(ranger)
library(caret)
library(missRanger)
library(kernlab)
library(C50)
library(e1071)
#library(mixOmics)
library(corpcor)
library(tictoc)
#library(ropls)
#library(pROC)
#library(caTools)
#library(AUC)

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

##Import clinicla legend
human_sepsis_legend <- get_human_sepsis_legend()

##Filter out patients without sepsis
human_nonsepsis_data <- human_sepsis_data[human_sepsis_data$`CAP / FP` == "-", ]
human_sepsis_data <- human_sepsis_data[human_sepsis_data$`CAP / FP` != "-", ]

##only keep non-lipid metabolites
#human_sepsis_data <- subset(human_sepsis_data, select = c(1:5, which(human_sepsis_legend %in% c("acylcarnitine", "amino acid", "biogenic amine", "sugar"))))

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
human_sepsis_data_ml <- human_sepsis_data[, 1:which(colnames(human_sepsis_data) == "H1")]
human_sepsis_data_ml <- subset(human_sepsis_data_ml, Day == 0)
human_sepsis_data_ml <- subset(human_sepsis_data_ml, select = !(colnames(human_sepsis_data_ml) %in% c("Histamin", "PC aa C36:0", "DOPA"))) #remove metabs with questionable range
human_sepsis_data_ml$Survival <- as.numeric(factor(human_sepsis_data_ml$Survival, levels = sort(unique(human_sepsis_data_ml$Survival)))) - 1 #Dependent variable transformation
human_sepsis_data_ml <- human_sepsis_data_ml[, c(1:5, 5 + which(!colAnys(human_sepsis_data_ml == 0)[-1:-5]))]
###Strip phenomenological variables
human_sepsis_data_ml <- subset(human_sepsis_data_ml, select = c(colnames(human_sepsis_data_ml)[1:5], intersect(sig_t_class, colnames(human_sepsis_data_ml))))
colnames(human_sepsis_data_ml) <- make.names(colnames(human_sepsis_data_ml))
###Impute missing values
human_sepsis_data_ml[, -1:-5] <- scale(missRanger(human_sepsis_data_ml[, -1:-5], pmm.k = 3, num.trees = 100))

#Parralel nested TPLOCV-RFE
##Ranger
tic()
rg_tlpocv_res <- tlpocv_rfe_parallel(data_x = human_sepsis_data_ml[, -1:-5], 
                                     data_y = human_sepsis_data_ml["Survival"], 
                                     mc.cores = 7)
toc()
##Logit
lm_tr_fun <- function(tr_x, tr_y){
  data <- data.frame(tr_y = tr_y[[1]], tr_x)
  glm(formula = tr_y ~ ., data = data, family = binomial(link = "gaussian"))
}
lm_prob_fun <- function(classifier, te_x){
  predict(classifier, te_x, type = "response")
}
lm_varimp_fun <- function(classifier){
  abs(classifier$coefficients)
}
tic()
lm_tlpocv_res <- tlpocv_rfe_parallel(data_x = human_sepsis_data_ml[-1:-5], 
                                     data_y = human_sepsis_data_ml["Survival"], 
                                     mc.cores = 7, 
                                     train_fun = lm_tr_fun, 
                                     prob_fun = lm_prob_fun, 
                                     varimp_fun = lm_varimp_fun)
toc()
##SVM
sv_tr_fun <- function(tr_x, tr_y){
  data <- data.frame(tr_y = tr_y[[1]], tr_x)
  svm(formula = tr_y ~ ., data = data, scale = FALSE, type = "C-classification", kernel = "linear")
}
sv_prob_fun <- function(classifier, te_x){
  p <- predict(classifier, te_x, decision.values = TRUE)
  attr(p, "decision.value")
}
sv_varimp_fun <- function(classifier){
  d <- data.frame(diag(1, nrow = length(classifier$scaled)))
  colnames(d) <- attr(classifier$terms, "term.labels")
  abs(sapply(d, function(r) attr(predict(classifier, t(r), decision.values = TRUE), "decision.values")[1]))
}
tic()
sv_tlpocv_res <- tlpocv_rfe_parallel(data_x = human_sepsis_data_ml[, -1:-5], 
                                     data_y = human_sepsis_data_ml["Survival"], 
                                     mc.cores = 40, 
                                     train_fun = sv_tr_fun, 
                                     prob_fun = sv_prob_fun, 
                                     varimp_fun = sv_varimp_fun)
toc()
save.image()

#-----------------------------


rg.npr.repeat.df$Method <- "RF"
ks.npr.repeat.df$Method <- "SVM"
lm.npr.repeat.df$Method <- "Probit"
ol.npr.repeat.df$Method <- "OPLSDA"
rg.red.npr.repeat.df$Method <- "RF"
ks.red.npr.repeat.df$Method <- "SVM"
lm.red.npr.repeat.df$Method <- "Probit"
ol.red.npr.repeat.df$Method <- "OPLSDA"
rg.red.npr.FVal.df$Method <- "RF"
ks.red.npr.FVal.df$Method <- "SVM"
lm.red.npr.FVal.df$Method <- "Probit"

rg.auc.repeat.df$Method <- "RF"
ks.auc.repeat.df$Method <- "SVM"
lm.auc.repeat.df$Method <- "Probit"
rg.red.auc.repeat.df$Method <- "RF"
ks.red.auc.repeat.df$Method <- "SVM"
lm.red.auc.repeat.df$Method <- "Probit"
rg.red.auc.FVal.df$Method <- "RF"
ks.red.auc.FVal.df$Method <- "SVM"
lm.red.auc.FVal.df$Method <- "Probit"

#Calculate accuracy
rg.npr.repeat.df$accuracy <- rg.acc.repeat.df$acc
ks.npr.repeat.df$accuracy <- ks.acc.repeat.df$acc
lm.npr.repeat.df$accuracy <- lm.acc.repeat.df$acc
ol.npr.repeat.df$accuracy <- ol.acc.repeat.df$acc
rg.red.npr.repeat.df$accuracy <- rg.red.acc.repeat.df$acc
ks.red.npr.repeat.df$accuracy <- ks.red.acc.repeat.df$acc
lm.red.npr.repeat.df$accuracy <- lm.red.acc.repeat.df$acc
ol.red.npr.repeat.df$accuracy <- ol.red.acc.repeat.df$acc
rg.red.npr.FVal.df$accuracy <- rg.red.acc.FVal.df$acc
ks.red.npr.FVal.df$accuracy <- ks.red.acc.FVal.df$acc
lm.red.npr.FVal.df$accuracy <- lm.red.acc.FVal.df$acc
ol.red.npr.FVal.df$accuracy <- ol.red.acc.FVal.df$acc

#Calculate Negative/Positive Predictive Value

#Melt to long table
npr_by_day_data <- melt(rbindlist(list(rg.npr.repeat.df, ks.npr.repeat.df, lm.npr.repeat.df)), id.vars = c("var", "Method"))
npr_red_by_day_data <- melt(rbindlist(list(rg.red.npr.repeat.df, ks.red.npr.repeat.df, lm.red.npr.repeat.df)), id.vars = c("var", "Method"))
npr_FVal_data <- melt(rbindlist(list(rg.red.npr.FVal.df, ks.red.npr.FVal.df, lm.red.npr.FVal.df)), id.vars = c("var", "Method"))
auc_data <- melt(rbindlist(list(rg.auc.repeat.df, ks.auc.repeat.df, lm.auc.repeat.df)), id.vars = c("var", "Method"))
auc_red_data <- melt(rbindlist(list(rg.red.auc.repeat.df, ks.red.auc.repeat.df, lm.red.auc.repeat.df)), id.vars = c("var", "Method"))
auc_FVal_data <- melt(rbindlist(list(rg.red.auc.FVal.df, ks.red.auc.FVal.df, lm.red.auc.FVal.df)), id.vars = c("var", "Method"))

##Remove redundancy
npr_by_day_data <- subset(npr_by_day_data, subset = variable %in% c("ppv", "npv", "accuracy"))
npr_red_by_day_data <- subset(npr_red_by_day_data, subset = variable %in% c("ppv", "npv", "accuracy"))
npr_FVal_data <- subset(npr_FVal_data, subset = variable %in% c("ppv", "npv", "accuracy"))

#Combine for nice plot
npr_by_day_data$dimensions <- "All variables"
npr_red_by_day_data$dimensions <- "Important variables"
npr_data <- rbind(npr_by_day_data, npr_red_by_day_data)
npr_red_by_day_data$dimensions <- "UKJ data"
npr_FVal_data$dimensions <- "Ferrario data (validation)"
npr_red_FVal_data <- rbind(npr_red_by_day_data, npr_FVal_data)
auc_data$dimensions <- "All variables"
auc_red_data$dimensions <- "Important variables"
auc_data_full <- rbind(auc_data, auc_red_data)
auc_red_data$dimensions <- "UK data"
auc_FVal_data$dimensions <- "Ferrario data (validation)"
auc_red_FVal_data <- rbind(auc_red_data, auc_FVal_data)

#Plot FPR and FNR
npr_by_feat_num_plot <- ggplot(subset(npr_data, variable %in% c("ppv", "npv")), aes(x = var, y = value, color = Method)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", size = 0.5, position = position_dodge(width = 0.4)) +
  facet_grid(variable ~ dimensions, scales = "free_x") +
  ylim(0, 1) +
  scale_x_discrete(limits = c(var_range[seq(1, length(var_range), 4)], tot_n_var)) + 
  xlab("Features included") +
  ggtitle("") +
  theme_bw()
ggsave(filename = "npr_by_feat_num_isect_Ferrario.png", plot = npr_by_feat_num_plot, path = out_dir_pred, width = 8, height = 4, units = "in")

#Plot accuracy
acc_by_feat_num_plot <- ggplot(subset(npr_data, variable == "accuracy"), aes(x = var, y = value, color = Method)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", position = position_dodge(width = 0.4)) +
  facet_grid( ~ dimensions, scales = "free_x") +
  ylim(0, 1) +
  scale_x_discrete(limits = c(var_range[seq(1, length(var_range), 4)], tot_n_var)) +
  xlab("Features included") +
  ylab("Accuracy") +
  #ggtitle("Normal stratified Cross Validation may shows overfitting") +
  #ggtitle("Quasi-stratified Cross Validation prevents overfitting") +
  theme_bw()
ggsave(filename = "acc_by_feat_num_isect_Ferrario.png", plot = acc_by_feat_num_plot, path = out_dir_pred, width = 8, height = 2.5, units = "in")

#Plot AUC
auc_by_feat_num_plot <- ggplot(subset(auc_data_full, variable == "auc"), aes(x = var, y = value, color = Method)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", position = position_dodge(width = 0.4)) +
  facet_grid( ~ dimensions, scales = "free_x") +
  scale_x_discrete(limits = c(var_range[seq(1, length(var_range), 4)], tot_n_var)) + 
  ylim(0, 1) +
  xlab("Features included") +
  ylab("AUC") +
  #ggtitle("Normal stratified Cross Validation may shows overfitting") +
  #ggtitle("Quasi-stratified Cross Validation prevents overfitting") +
  theme_bw()
ggsave(filename = "auc_by_feat_num_isect_Ferrario.png", plot = auc_by_feat_num_plot, path = out_dir_pred, width = 8, height = 2.5, units = "in")

#Plot npr on validation data from Ferrario et al.
npr_by_feat_num_FVal_plot <- ggplot(subset(npr_FVal_data, variable %in% c("ppv", "npv")), aes(x = var, y = value, color = Method)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", size = 0.5, position = position_dodge(width = 0.4)) +
  facet_grid( ~ variable, scales = "free_x") +
  ylim(0, 1) +
  scale_x_discrete(limits = c(var_range[seq(1, length(var_range), 4)], tot_n_var)) + 
  xlab("Features included") +
  ggtitle("") +
  theme_bw()
ggsave(filename = "npr_by_feat_num_Ferrario_validation.png", plot = npr_by_feat_num_FVal_plot, path = out_dir_pred, width = 4, height = 4, units = "in")

#Plot accuracy on validation data from Ferrario et al.
acc_by_feat_num_FVal_plot <- ggplot(subset(npr_FVal_data, variable == "accuracy"), aes(x = var, y = value, color = Method)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", position = position_dodge(width = 0.4)) +
  ylim(0, 1) +
  scale_x_discrete(limits = c(var_range[seq(1, length(var_range), 4)], tot_n_var)) +
  xlab("Features included") +
  ylab("Accuracy") +
  theme_bw()
ggsave(filename = "acc_by_feat_num_Ferrario_validation.png", plot = acc_by_feat_num_FVal_plot, path = out_dir_pred, width = 4, height = 2.5, units = "in")

#Plot auc on validation data from Ferrario et al.
auc_by_feat_num_FVal_plot <- ggplot(subset(auc_FVal_data, variable == "auc"), aes(x = var, y = value, color = Method)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", position = position_dodge(width = 0.4)) +
  scale_x_discrete(limits = c(var_range[seq(1, length(var_range), 4)], tot_n_var)) + 
  ylim(0, 1) +
  xlab("Features included") +
  ylab("AUC") +
  theme_bw()
ggsave(filename = "auc_by_feat_num_Ferrario_validation.png", plot = auc_by_feat_num_FVal_plot, path = out_dir_pred, width = 4, height = 2.5, units = "in")

#Plot npr and acc on validation and test data
npr_red_FVal_plot <- ggplot(subset(npr_red_FVal_data, variable %in% c("accuracy", "ppv", "npv")), aes(x = var, y = value, color = Method)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", size = 0.5, position = position_dodge(width = 0.4)) +
  facet_grid(variable ~ dimensions) +
  ylim(0, 1) +
  xlab("Features included") +
  ggtitle("") +
  theme_bw()
ggsave(filename = "acc_npr_test_and_Ferrario_validation.png", plot = npr_red_FVal_plot, path = out_dir_pred, width = 6, height = 6, units = "in")

#Plot metabs against metabs for six most important variables
# selection <- unique(na.omit(unlist(var_set_name_list[1:10])))
# plot(as.formula(paste0("~ ", paste0(selection, collapse = "+"))), data = human_sepsis_data_ml[,selection], col = human_sepsis_data_ml$Survival + 1)
# 
# hsiclasso_set <- c('Leptin', 'PC aa C30:2', 'TSH', 'PC aa C30:0', 'FT4', 'C5-M-DC', 'C4', 'C16', 'GH', 't4-OH-Pro', 'PC aa C38:3', 'Prolactin', 'DHEA', 'Leu', 'SM C24:1')
# selection <- na.omit(match(hsiclasso_set, colnames(human_sepsis_data)))
# hsd <- subset(human_sepsis_data, Day == 0)
# colnames(hsd) <- make.names(colnames(hsd))
# plot(as.formula(paste0("~ ", paste0(colnames(hsd)[selection], collapse = "+"))), data = hsd[,selection], col = as.numeric(as.factor(hsd$Survival)) + 1)
