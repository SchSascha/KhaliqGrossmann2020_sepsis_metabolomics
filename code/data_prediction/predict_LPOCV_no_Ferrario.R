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
# ##Logit
# lm_tr_fun <- function(tr_x, tr_y){
#   data <- data.frame(tr_y = tr_y[[1]], tr_x)
#   glm(formula = tr_y ~ ., data = data, family = binomial(link = "gaussian"))
# }
# lm_prob_fun <- function(classifier, te_x){
#   predict(classifier, te_x, type = "response")
# }
# lm_varimp_fun <- function(classifier){
#   abs(classifier$coefficients)
# }
# tic()
# lm_tlpocv_res <- tlpocv_rfe_parallel(data_x = human_sepsis_data_ml[-1:-5], 
#                                      data_y = human_sepsis_data_ml["Survival"], 
#                                      mc.cores = 7, 
#                                      train_fun = lm_tr_fun, 
#                                      prob_fun = lm_prob_fun, 
#                                      varimp_fun = lm_varimp_fun)
# toc()
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
sv_tlpocv_res <- tlpocv_rfe_parallel(data_x = human_sepsis_data_ml[-1:-5], 
                                     data_y = human_sepsis_data_ml["Survival"], 
                                     mc.cores = 6, 
                                     train_fun = sv_tr_fun, 
                                     prob_fun = sv_prob_fun, 
                                     varimp_fun = sv_varimp_fun)
toc()
save.image()

#Plot AUCs and ROCs

rg_AUC_data <- data.frame(Reduce("rbind", rg_tlpocv_res$inner_AUC))
rg_AUC_data <- rbind(rg_AUC_data, data.frame(t(rg_tlpocv_res$outer_AUC)))
colnames(rg_AUC_data) <- (length(rg_tlpocv_res$inner_AUC[[1]]) + 1):2
rg_AUC_data$stage <- c(rep("Test data", length(rg_tlpocv_res$inner_AUC)), "Validation data")
rg_AUC_data_long <- melt(rg_AUC_data, id.vars = c("stage"))
rg_AUC_data_long$variable <- as.numeric(rg_AUC_data_long$variable)
rg_AUC_data_long$value <- 1 - rg_AUC_data_long$value
p <- ggplot(data = rg_AUC_data_long, mapping = aes(x = variable, y = value, color = stage, fill = stage)) + 
  stat_summary(fun.y = median, fun.ymax = function(x) quantile(x, p = 0.75), fun.ymin = function(x) quantile(x, p = 0.25), geom = "ribbon", alpha = 0.5, colour = NA) + 
  stat_summary(fun.y = median, geom = "line") + 
  ylim(0, 1) + 
  ylab("Median AUC") +
  xlab("Nubmer of features") +
  ggtitle("RF: Test and Validation AUC of nested TLPOCV-RFE") +
  theme_bw()
ggsave(filename = "RF_TLPOCV_RFE_AUC.png", path = out_dir_pred, plot = p, width = 10, height = 5, units = "in")

sv_AUC_data <- data.frame(Reduce("rbind", sv_tlpocv_res$inner_AUC))
sv_AUC_data <- rbind(sv_AUC_data, data.frame(t(sv_tlpocv_res$outer_AUC)))
colnames(sv_AUC_data) <- (length(sv_tlpocv_res$inner_AUC[[1]]) + 1):2
sv_AUC_data$stage <- c(rep("Test data", length(sv_tlpocv_res$inner_AUC)), "Validation data")
sv_AUC_data_long <- melt(sv_AUC_data, id.vars = c("stage"))
sv_AUC_data_long$variable <- as.numeric(sv_AUC_data_long$variable)
p <- ggplot(data = sv_AUC_data_long, mapping = aes(x = variable, y = value, color = stage, fill = stage)) + 
  stat_summary(fun.y = median, fun.ymax = function(x) quantile(x, p = 0.75), fun.ymin = function(x) quantile(x, p = 0.25), geom = "ribbon", alpha = 0.5, colour = NA) + 
  stat_summary(fun.y = median, geom = "line") + 
  ylim(0, 1) + 
  ylab("Median AUC") +
  xlab("Nubmer of features") +
  ggtitle("SVM: Test and Validation AUC of nested TLPOCV-RFE") +
  theme_bw()
ggsave(filename = "SVM_TLPOCV_RFE_AUC.png", path = out_dir_pred, plot = p, width = 10, height = 5, units = "in")

sv_int_ROC_data <- lapply(sv_tlpocv_res$int_sample_ranking[[length(sv_tlpocv_res$int_sample_ranking)]], function(e) ml.roc(ref = e[, 1], conf = e[, 2]))
sv_int_ROC_data <- lapply(sv_int_ROC_data, function(e) aggregate(x = e[, 2], by = list(FPR = e[, 1]), FUN = min))
sv_int_ROC_data <- data.frame(Reduce("rbind", sv_int_ROC_data))
if (min(sv_int_ROC_data[, 1]) != 0)
  sv_int_ROC_data <- rbind(sv_int_ROC_data, c(0, 0))
if (max(sv_int_ROC_data[, 2]) != 1)
  sv_int_ROC_data <- rbind(sv_int_ROC_data, c(1, 1))
colnames(sv_int_ROC_data) <- c("FPR", "TPR")
sv_int_ROC_data$stage <- "Test data"
ranks_2feat <- sv_tlpocv_res$ext_sample_ranking[[length(sv_tlpocv_res$ext_sample_ranking)]]
sv_ext_ROC_data <- data.frame(ml.roc(ref = ranks_2feat[, 1], conf = ranks_2feat[, 2]))
sv_ext_ROC_data <- aggregate(x = sv_ext_ROC_data[, 2], by = list(FPR = sv_ext_ROC_data[, 1]), FUN = min)
if (min(sv_ext_ROC_data[, 1]) != 0)
  sv_ext_ROC_data <- rbind(sv_ext_ROC_data, c(0, 0))
if (max(sv_ext_ROC_data[, 2]) != 1)
  sv_ext_ROC_data <- rbind(sv_ext_ROC_data, c(1, 1))
colnames(sv_ext_ROC_data) <- c("FPR", "TPR")
sv_ext_ROC_data$stage <- "Validation data"
sv_ROC_data <- rbind(sv_int_ROC_data, sv_ext_ROC_data)
p <- ggplot(data = sv_ROC_data, mapping = aes(x = FPR, y = TPR, color = stage, fill = stage)) + 
  stat_summary(fun.y = median, fun.ymax = function(x) quantile(x, p = 0.75), fun.ymin = function(x) quantile(x, p = 0.25), geom = "ribbon", alpha = 0.5, colour = NA) + 
  stat_summary(fun.y = median, geom = "line") + 
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  ggtitle("SVM: Test and Validation ROC of nested TLPOCV RFE\nbest 2 feature set") +
  theme_bw() + 
  theme(panel.grid = element_blank())
ggsave(filename = "SVM_TLPOCV_RFE_ROC_2feat.png", path = out_dir_pred, plot = p, width = 6, height = 5, units = "in")

rg_int_ROC_data <- lapply(rg_tlpocv_res$int_sample_ranking[[length(rg_tlpocv_res$int_sample_ranking)]], function(e) ml.roc(ref = e[, 1], conf = e[, 2]))
rg_int_ROC_data <- lapply(rg_int_ROC_data, function(e) aggregate(x = e[, 2], by = list(FPR = e[, 1]), FUN = min))
rg_int_ROC_data <- data.frame(Reduce("rbind", rg_int_ROC_data))
if (min(rg_int_ROC_data[, 1]) != 0)
  rg_int_ROC_data <- rbind(rg_int_ROC_data, c(0, 0))
if (max(rg_int_ROC_data[, 2]) != 1)
  rg_int_ROC_data <- rbind(rg_int_ROC_data, c(1, 1))
colnames(rg_int_ROC_data) <- c("FPR", "TPR")
rg_int_ROC_data$stage <- "Test data"
ranks_2feat <- rg_tlpocv_res$ext_sample_ranking[[length(rg_tlpocv_res$ext_sample_ranking)]]
rg_ext_ROC_data <- data.frame(ml.roc(ref = ranks_2feat[, 1], conf = ranks_2feat[, 2]))
rg_ext_ROC_data <- aggregate(x = rg_ext_ROC_data[, 2], by = list(FPR = rg_ext_ROC_data[, 1]), FUN = min)
if (min(rg_ext_ROC_data[, 1]) != 0)
  rg_ext_ROC_data <- rbind(rg_ext_ROC_data, c(0, 0))
if (max(rg_ext_ROC_data[, 2]) != 1)
  rg_ext_ROC_data <- rbind(rg_ext_ROC_data, c(1, 1))
colnames(rg_ext_ROC_data) <- c("FPR", "TPR")
rg_ext_ROC_data$stage <- "Validation data"
rg_ROC_data <- rbind(rg_int_ROC_data, rg_ext_ROC_data)
p <- ggplot(data = rg_ROC_data, mapping = aes(x = FPR, y = TPR, color = stage, fill = stage)) + 
  stat_summary(fun.y = median, fun.ymax = function(x) quantile(x, p = 0.75), fun.ymin = function(x) quantile(x, p = 0.25), geom = "ribbon", alpha = 0.5, colour = NA) + 
  stat_summary(fun.y = median, geom = "line") + 
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  ggtitle("RF: Test and Validation ROC of nested TLPOCV RFE\nbest 2 feature set") +
  theme_bw() + 
  theme(panel.grid = element_blank())
ggsave(filename = "RF_TLPOCV_RFE_ROC_2feat.png", path = out_dir_pred, plot = p, width = 6, height = 5, units = "in")
