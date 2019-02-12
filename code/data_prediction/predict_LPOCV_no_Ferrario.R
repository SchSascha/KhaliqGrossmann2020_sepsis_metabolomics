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
human_sepsis_data_ml <- human_sepsis_data_ml[, c(1:6, 6 + which(!colAnys(human_sepsis_data_ml == 0)[-1:-6]))]
###Strip phenomenological variables
#human_sepsis_data_ml <- subset(human_sepsis_data_ml, select = c(colnames(human_sepsis_data_ml)[1:5], intersect(sig_t_class, colnames(human_sepsis_data_ml))))
colnames(human_sepsis_data_ml) <- make.names(colnames(human_sepsis_data_ml))
###Impute missing values
human_sepsis_data_ml[, -1:-6] <- scale(missRanger(human_sepsis_data_ml[, -1:-6], pmm.k = 3, num.trees = 100))

#Parralel nested TPLOCV-RFE
##Ranger
tic()
rg_tlpocv_res <- tlpocv_rfe_parallel(data_x = human_sepsis_data_ml[-1:-6],
                                     data_y = human_sepsis_data_ml["Survival"],
                                     mc.cores = 30)
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
sv_tlpocv_res <- tlpocv_rfe_parallel(data_x = human_sepsis_data_ml[-1:-6],
                                     data_y = human_sepsis_data_ml["Survival"],
                                     mc.cores = 30,
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
rg_AUC_data_long$variable <- as.numeric(as.character(rg_AUC_data_long$variable))
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

rg_int_ROC_data <- lapply(lapply(rg_tlpocv_res$int_sample_ranking, `[[`, length(rg_tlpocv_res$int_sample_ranking[[1]])), function(e) ml.roc(ref = 1 - e[, 1], conf = e[, 2]))
rg_int_ROC_data <- lapply(rg_int_ROC_data, function(e) stepfun(e[, 1], c(0, e[, 2])))
rg_int_ROC_data <- lapply(rg_int_ROC_data, function(f) data.frame(FPR = seq(0, 1, by = 0.01), TPR = f(seq(0, 1, by = 0.01))))
rg_int_ROC_data <- data.frame(Reduce("rbind", rg_int_ROC_data))
if (min(rg_int_ROC_data[, 1]) != 0)
  rg_int_ROC_data <- rbind(c(0, 0), rg_int_ROC_data)
if (max(rg_int_ROC_data[, 1]) != 1)
  rg_int_ROC_data <- rbind(rg_int_ROC_data, c(1, 1))
colnames(rg_int_ROC_data) <- c("FPR", "TPR")
rg_int_ROC_data$stage <- "Test data"
ranks_2feat <- rg_tlpocv_res$ext_sample_ranking[[length(rg_tlpocv_res$ext_sample_ranking)]]
rg_ext_ROC_data <- data.frame(ml.roc(ref = 1 - ranks_2feat[, 1], conf = ranks_2feat[, 2]))
if (min(rg_ext_ROC_data[, 2]) != 0)
  rg_ext_ROC_data <- rbind(c(0, 0), rg_ext_ROC_data)
if (max(rg_ext_ROC_data[, 1]) != 1)
  rg_ext_ROC_data <- rbind(rg_ext_ROC_data, c(1, 1))
colnames(rg_ext_ROC_data) <- c("FPR", "TPR")
rg_ext_ROC_data$stage <- "Validation data"
rg_ROC_data <- rbind(rg_int_ROC_data, rg_ext_ROC_data)
p <- ggplot(data = rg_ROC_data, mapping = aes(x = FPR, y = TPR, color = stage, fill = stage)) + 
  stat_summary(data = rg_int_ROC_data, 
               fun.y = mean, fun.ymax = function(x) quantile(x, p = 0.75), fun.ymin = function(x) quantile(x, p = 0.25), 
               geom = "ribbon", alpha = 0.5, colour = NA) + 
  stat_summary(data = rg_int_ROC_data, 
               fun.y = median, geom = "line") +
  geom_line(data = rg_ext_ROC_data) +
  geom_segment(x = 0, y = 0, xend = 1, yend = 1, size = 0.25, inherit.aes = FALSE) +
  geom_text(x = 0.95, y = 0.15, hjust = "right", color = scales::hue_pal()(2)[1],
            label = paste0("Training AUC = ", format(mean(subset(rg_AUC_data_long, variable == 2 & stage == "Test data", "value")[[1]]), digits = 3))) +
  geom_text(x = 0.95, y = 0.10, hjust = "right", color = scales::hue_pal()(2)[2],
            label = paste0("Validation AUC = ", format(subset(rg_AUC_data_long, variable == 2 & stage == "Validation data", "value"), digits = 3))) +
  scale_x_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
  ggtitle("RF: Test and Validation ROC of nested TLPOCV RFE\nbest 2 feature set") +
  theme_bw() + 
  theme(panel.grid = element_blank())
ggsave(filename = "RF_TLPOCV_RFE_ROC_2feat.png", path = out_dir_pred, plot = p, width = 6, height = 5, units = "in")

p <- ggplot(data = subset(human_sepsis_data_ml, Day == 0, c("Survival", rg_tlpocv_res$best_features[[75]])), mapping = aes_string(x = rg_tlpocv_res$best_features[[75]][1], y = rg_tlpocv_res$best_features[[75]][2], color = quote(factor(Survival)))) + 
  geom_point() + 
  theme_bw()
ggsave(plot = p, filename = "RF_class_separation_2feat.png", width = 5, height = 4, units = "in", path = out_dir_pred)

sv_AUC_data <- data.frame(Reduce("rbind", sv_tlpocv_res$inner_AUC))
sv_AUC_data <- rbind(sv_AUC_data, data.frame(t(sv_tlpocv_res$outer_AUC)))
colnames(sv_AUC_data) <- (length(sv_tlpocv_res$inner_AUC[[1]]) + 1):2
sv_AUC_data$stage <- c(rep("Test data", length(sv_tlpocv_res$inner_AUC)), "Validation data")
sv_AUC_data_long <- melt(sv_AUC_data, id.vars = c("stage"))
sv_AUC_data_long$variable <- as.numeric(as.character(sv_AUC_data_long$variable))
sv_AUC_data_long$value <- 1 - sv_AUC_data_long$value
p <- ggplot(data = sv_AUC_data_long, mapping = aes(x = variable, y = value, color = stage, fill = stage)) +
  stat_summary(fun.y = median, fun.ymax = function(x) quantile(x, p = 0.75), fun.ymin = function(x) quantile(x, p = 0.25), geom = "ribbon", alpha = 0.5, colour = NA) +
  stat_summary(fun.y = median, geom = "line") +
  ylim(0, 1) +
  ylab("Median AUC") +
  xlab("Nubmer of features") +
  ggtitle("SVM: Test and Validation AUC of nested TLPOCV-RFE") +
  theme_bw()
ggsave(filename = "SVM_TLPOCV_RFE_AUC.png", path = out_dir_pred, plot = p, width = 10, height = 5, units = "in")

sv_int_ROC_data <- lapply(lapply(sv_tlpocv_res$int_sample_ranking, `[[`, length(sv_tlpocv_res$int_sample_ranking[[1]])), function(e) ml.roc(ref = 1 - e[, 1], conf = e[, 2]))
sv_int_ROC_data <- lapply(sv_int_ROC_data, function(e) stepfun(e[, 1], c(0, e[, 2])))
sv_int_ROC_data <- lapply(sv_int_ROC_data, function(f) data.frame(FPR = seq(0, 1, by = 0.01), TPR = f(seq(0, 1, by = 0.01))))
sv_int_ROC_data <- data.frame(Reduce("rbind", sv_int_ROC_data))
if (min(sv_int_ROC_data[, 1]) != 0)
  sv_int_ROC_data <- rbind(c(0, 0), sv_int_ROC_data)
if (max(sv_int_ROC_data[, 1]) != 1)
  sv_int_ROC_data <- rbind(sv_int_ROC_data, c(1, 1))
colnames(sv_int_ROC_data) <- c("FPR", "TPR")
sv_int_ROC_data$stage <- "Test data"
ranks_2feat <- sv_tlpocv_res$ext_sample_ranking[[length(sv_tlpocv_res$ext_sample_ranking)]]
sv_ext_ROC_data <- data.frame(ml.roc(ref = 1 - ranks_2feat[, 1], conf = ranks_2feat[, 2]))
if (min(sv_ext_ROC_data[, 2]) != 0)
  sv_ext_ROC_data <- rbind(c(0, 0), sv_ext_ROC_data)
if (max(sv_ext_ROC_data[, 1]) != 1)
  sv_ext_ROC_data <- rbind(sv_ext_ROC_data, c(1, 1))
colnames(sv_ext_ROC_data) <- c("FPR", "TPR")
sv_ext_ROC_data$stage <- "Validation data"
sv_ROC_data <- rbind(sv_int_ROC_data, sv_ext_ROC_data)
p <- ggplot(data = sv_ROC_data, mapping = aes(x = FPR, y = TPR, color = stage, fill = stage)) +
  stat_summary(data = sv_int_ROC_data,
               fun.y = mean, fun.ymax = function(x) quantile(x, p = 0.75), fun.ymin = function(x) quantile(x, p = 0.25),
               geom = "ribbon", alpha = 0.5, colour = NA) +
  stat_summary(data = sv_int_ROC_data,
               fun.y = median, geom = "line") +
  geom_line(data = sv_ext_ROC_data) +
  geom_segment(x = 0, y = 0, xend = 1, yend = 1, size = 0.25, inherit.aes = FALSE) +
  geom_text(x = 0.95, y = 0.15, hjust = "right", color = scales::hue_pal()(2)[1],
            label = paste0("Training AUC = ", format(mean(subset(sv_AUC_data_long, variable == 2 & stage == "Test data", "value")[[1]]), digits = 3))) +
  geom_text(x = 0.95, y = 0.10, hjust = "right", color = scales::hue_pal()(2)[2],
            label = paste0("Validation AUC = ", format(subset(sv_AUC_data_long, variable == 2 & stage == "Validation data", "value"), digits = 3))) +
  scale_x_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
  ggtitle("SVM: Test and Validation ROC of nested TLPOCV RFE\nbest 2 feature set") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(filename = "SVM_TLPOCV_RFE_ROC_2feat.png", path = out_dir_pred, plot = p, width = 6, height = 5, units = "in")

p <- ggplot(data = subset(human_sepsis_data_ml, Day == 0, c("Survival", sv_tlpocv_res$best_features[[75]])), mapping = aes_string(x = sv_tlpocv_res$best_features[[75]][1], y = sv_tlpocv_res$best_features[[75]][2], color = quote(factor(Survival)))) + 
  geom_point() + 
  theme_bw()
ggsave(plot = p, filename = "SVM_class_separation_2feat.png", width = 5, height = 4, units = "in", path = out_dir_pred)
