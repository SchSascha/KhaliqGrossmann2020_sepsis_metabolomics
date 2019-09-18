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
library(pheatmap)
library(cowplot)
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
human_data <- get_human_sepsis_data()

##Import clinicla legend
human_sepsis_legend <- get_human_sepsis_legend()

##Filter out patients without sepsis
human_nonsepsis_data <- human_data[human_data$`CAP / FP` == "-", ]
human_sepsis_data <- human_data[human_data$`CAP / FP` != "-", ]

##only keep non-lipid metabolites
#human_sepsis_data <- subset(human_sepsis_data, select = c(1:5, which(human_sepsis_legend %in% c("acylcarnitine", "amino acid", "biogenic amine", "sugar"))))

#Try prediction
##On human
###Build data set for "vanilla" prediction
human_sepsis_data_ml <- human_sepsis_data[, 1:which(colnames(human_sepsis_data) == "H1")]
human_sepsis_data_ml <- subset(human_sepsis_data_ml, Day == 0)
human_sepsis_data_ml <- subset(human_sepsis_data_ml, select = !(colnames(human_sepsis_data_ml) %in% c("Histamin", "PC aa C36:0", "DOPA"))) #remove metabs with questionable range
human_sepsis_data_ml$Survival <- as.numeric(factor(human_sepsis_data_ml$Survival, levels = sort(unique(human_sepsis_data_ml$Survival)))) - 1 #Dependent variable transformation
col_all_nonzeros <- 6 + which(!colAnys(human_sepsis_data_ml == 0)[-1:-6])
human_sepsis_data_ml <- human_sepsis_data_ml[, c(1:6, col_all_nonzeros)]
human_nonsepsis_data_ml <- human_nonsepsis_data[, 1:which(colnames(human_sepsis_data) == "H1")]
human_nonsepsis_data_ml <- subset(human_nonsepsis_data_ml, Day == 0)
human_nonsepsis_data_ml <- subset(human_nonsepsis_data_ml, select = !(colnames(human_nonsepsis_data_ml) %in% c("Histamin", "PC aa C36:0", "DOPA"))) #remove metabs with questionable range
human_nonsepsis_data_ml$Survival <- as.numeric(factor(human_nonsepsis_data_ml$Survival, levels = sort(unique(human_nonsepsis_data_ml$Survival)))) - 1 #Dependent variable transformation
human_nonsepsis_data_ml <- human_nonsepsis_data_ml[, c(1:6, col_all_nonzeros)]
###Strip phenomenological variables
#human_sepsis_data_ml <- subset(human_sepsis_data_ml, select = c(colnames(human_sepsis_data_ml)[1:5], intersect(sig_t_class, colnames(human_sepsis_data_ml))))
colnames(human_sepsis_data_ml) <- make.names(colnames(human_sepsis_data_ml))
colnames(human_nonsepsis_data_ml) <- make.names(colnames(human_nonsepsis_data_ml))
###Impute missing values
scaled <- scale(missRanger(human_sepsis_data_ml[, -1:-6], pmm.k = 3, num.trees = 100))
sepsis_data_scales <- lapply(c("scaled:center", "scaled:scale"), function(w) attr(scaled, which = w))
human_sepsis_data_ml[, -1:-6] <- scaled
human_nonsepsis_data_ml[, -1:-6] <- scale(missRanger(human_nonsepsis_data_ml[, -1:-6], pmm.k = 3, num.trees = 100), 
                                          center = sepsis_data_scales[[1]], 
                                          scale = sepsis_data_scales[[2]])

##Set number of features to try in RFE
set_sizes <- c(1:20, 25, 30, 40, 50, 60, 80, 100, 120, 150, 180)
set_sizes <- set_sizes[set_sizes < (ncol(human_sepsis_data_ml) - 6)]

#Parralel nested TPLOCV-RFE
##Ranger
tic()
rg_tlpocv_res <- tlpocv_rfe_parallel(data_x = human_sepsis_data_ml[-1:-6],
                                     data_y = human_sepsis_data_ml["Survival"],
                                     set_sizes = set_sizes,
                                     mc.cores = 6)
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
  rownames(d) <- colnames(d)
  empty <- data.frame(d[1, ])
  empty[1:length(empty)] <- 0
  colnames(empty) <- colnames(d)
  intercept <- attr(predict(classifier, empty, decision.values = TRUE), "decision.values")[1]
  t(abs(attr(predict(classifier, d, decision.values = TRUE), "decision.values") + intercept))
}
tic()
sv_tlpocv_res <- tlpocv_rfe_parallel(data_x = human_sepsis_data_ml[-1:-6],
                                     data_y = human_sepsis_data_ml["Survival"],
                                     set_sizes = set_sizes,
                                     mc.cores = 6,
                                     train_fun = sv_tr_fun,
                                     prob_fun = sv_prob_fun,
                                     varimp_fun = sv_varimp_fun)
toc()
save.image()

#Test performance on nonseptic patients
num_ns_reps <- 60
num_best_feats <- 6
##Random Forest
l <- length(rg_tlpocv_res$best_features)
rg_ns_feats <- rg_tlpocv_res$best_features[(l-1):(l - num_best_feats + 1)]
rg_cf_pred <- list()
for (n in seq_along(rg_ns_feats)){
  rg_cf_pred[[n]] <- list()
  for (num in 1:num_ns_reps){
    rg_cf <- ranger(data = human_sepsis_data_ml[, c("Survival", rg_ns_feats[[n]])], write.forest = TRUE, probability = TRUE, save.memory = FALSE, dependent.variable.name = "Survival", num.trees = 500)
    rg_cf_pred[[n]][[num]] <- predict(object = rg_cf, data = human_nonsepsis_data_ml, type = "response")
  }
}
rg_cf_pred <- lapply(rg_cf_pred, lapply, `[[`, "predictions")
rg_cf_pred <- lapply(rg_cf_pred, lapply, function(e) e[, 1])
rg_cf_auc <- lapply(rg_cf_pred, sapply, function(p) ml.auc(ref = 1 - human_nonsepsis_data_ml$Survival, conf = p))
rg_cf_auc <- sapply(rg_cf_auc, mean)
rg_cf_auc
rg_cf_roc <- lapply(rg_cf_pred, lapply, function(p) ml.roc(ref = 1 - human_nonsepsis_data_ml$Survival, conf = p))
##SVM
l <- length(sv_tlpocv_res$best_features)
sv_ns_feats <- sv_tlpocv_res$best_features[(l-1):(l - num_best_feats + 1)]
sv_cf_pred <- list()
for (n in seq_along(sv_ns_feats)){
  sv_cf_pred[[n]] <- list()
  for (num in 1:num_ns_reps){
    sv_cf <- sv_tr_fun(tr_x = human_sepsis_data_ml[, sv_ns_feats[[n]]], tr_y = human_sepsis_data_ml["Survival"])
    sv_cf_pred[[n]][[num]] <- sv_prob_fun(classifier = sv_cf, te_x = human_nonsepsis_data_ml[sv_ns_feats[[n]]])
  }
}
sv_cf_auc <- lapply(sv_cf_pred, sapply, function(p) ml.auc(ref = 1 - human_nonsepsis_data_ml$Survival, conf = p))
sv_cf_auc <- sapply(sv_cf_auc, mean)
sv_cf_auc
sv_cf_roc <- lapply(sv_cf_pred, lapply, function(p) ml.roc(ref = 1 - human_nonsepsis_data_ml$Survival, conf = p))

#Plot AUCs and ROCs

rg_AUC_data <- data.frame(Reduce("rbind", rg_tlpocv_res$inner_AUC))
rg_AUC_data <- rbind(rg_AUC_data, data.frame(t(rg_tlpocv_res$outer_AUC)))
colnames(rg_AUC_data) <- sapply(rg_tlpocv_res$best_features, length)
rg_AUC_data$Stage <- c(rep("Test data", length(rg_tlpocv_res$inner_AUC)), "Validation data")
rg_AUC_data_long <- melt(rg_AUC_data, id.vars = c("Stage"))
rg_AUC_data_long$variable <- as.numeric(as.character(rg_AUC_data_long$variable))
rg_AUC_data_long$value <- 1 - rg_AUC_data_long$value
p <- ggplot(data = rg_AUC_data_long, mapping = aes(x = variable, y = value, color = Stage, fill = Stage)) + 
  stat_summary(fun.y = median, fun.ymax = function(x) quantile(x, p = 0.75), fun.ymin = function(x) quantile(x, p = 0.25), geom = "ribbon", alpha = 0.5, colour = NA) + 
  stat_summary(fun.y = median, geom = "line") + 
  ylim(0, 1) + 
  ylab("Median AUC") +
  xlab("Number of features") +
  ggtitle("RF: Test and Validation AUC of nested TLPOCV-RFE") +
  theme_bw()
ggsave(filename = "RF_TLPOCV_RFE_AUC.png", path = out_dir_pred, plot = p, width = 10, height = 5, units = "in")

f2_pos <- which(sapply(rg_tlpocv_res$best_features, length) == 2)
rg_int_ROC_data <- lapply(lapply(rg_tlpocv_res$int_sample_ranking, `[[`, f2_pos), function(e) ml.roc(ref = 1 - e[, 1], conf = e[, 2]))
rg_int_ROC_data <- lapply(rg_int_ROC_data, function(e) stepfun(e[, 1], c(0, e[, 2])))
rg_int_ROC_data <- lapply(rg_int_ROC_data, function(f) data.frame(FPR = seq(0, 1, by = 0.01), TPR = f(seq(0, 1, by = 0.01))))
rg_int_ROC_data <- data.frame(Reduce("rbind", rg_int_ROC_data))
if (min(rg_int_ROC_data[, 1]) != 0)
  rg_int_ROC_data <- rbind(c(0, 0), rg_int_ROC_data)
if (max(rg_int_ROC_data[, 1]) != 1)
  rg_int_ROC_data <- rbind(rg_int_ROC_data, c(1, 1))
colnames(rg_int_ROC_data) <- c("FPR", "TPR")
rg_int_ROC_data$Stage <- "Test data"
ranks_2feat <- rg_tlpocv_res$ext_sample_ranking[[f2_pos]]
rg_ext_ROC_data <- data.frame(ml.roc(ref = 1 - ranks_2feat[, 1], conf = ranks_2feat[, 2]))
if (min(rg_ext_ROC_data[, 2]) != 0)
  rg_ext_ROC_data <- rbind(c(0, 0), rg_ext_ROC_data)
if (max(rg_ext_ROC_data[, 1]) != 1)
  rg_ext_ROC_data <- rbind(rg_ext_ROC_data, c(1, 1))
colnames(rg_ext_ROC_data) <- c("FPR", "TPR")
rg_ext_ROC_data$Stage <- "Validation data"
rg_ROC_data <- rbind(rg_int_ROC_data, rg_ext_ROC_data)
rsize <- rel(3.5)
p_roc_rg_2feat <- ggplot(data = rg_ROC_data, mapping = aes(x = FPR, y = TPR, color = Stage, fill = Stage)) + 
  stat_summary(data = rg_int_ROC_data, 
               fun.y = mean, fun.ymax = function(x) quantile(x, p = 0.75), fun.ymin = function(x) quantile(x, p = 0.25), 
               geom = "ribbon", alpha = 0.5, colour = NA) + 
  stat_summary(data = rg_int_ROC_data, 
               fun.y = median, geom = "line", size = rel(1.5)) +
  geom_line(data = rg_ext_ROC_data, size = rel(1.5)) +
  geom_segment(x = 0, y = 0, xend = 1, yend = 1, size = 0.5, inherit.aes = FALSE) +
  geom_text(x = 0.95, y = 0.15, hjust = "right", color = scales::hue_pal()(2)[1], size = rel(4.5),
            label = paste0("Test AUC = ", format(mean(subset(rg_AUC_data_long, variable == 2 & Stage == "Test data", "value")[[1]]), digits = 3, nsmall = 3))) +
  geom_text(x = 0.95, y = 0.10, hjust = "right", color = scales::hue_pal()(2)[2], size = rel(4.5),
            label = paste0("Validation AUC = ", format(subset(rg_AUC_data_long, variable == 2 & Stage == "Validation data", "value"), digits = 3, nsmall = 3))) +
  scale_x_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
  #ggtitle("RF: Test and Validation ROC of nested TLPOCV RFE\nbest 2 feature set") +
  geom_text(label = "ROC curve\nRandom Forest\n2 features", x = 0.5, y = 0.85, size = rel(4.5), inherit.aes = FALSE) +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = "none", text = element_text(size = rsize), plot.margin = margin(r = 10, l = 10))
ggsave(filename = "RF_TLPOCV_RFE_ROC_2feat.png", path = out_dir_pred, plot = p_roc_rg_2feat, width = 6, height = 5, units = "in")

f2 <- rg_tlpocv_res$best_features[[which(sapply(rg_tlpocv_res$best_features, length) == 2)]]
p <- ggplot(data = subset(human_sepsis_data_ml, Day == 0, c("Survival", f2)), mapping = aes_string(x = f2[1], y = f2[2], color = quote(factor(Survival)))) + 
  geom_point() + 
  theme_bw()
ggsave(plot = p, filename = "RF_class_separation_2feat.png", width = 5, height = 4, units = "in", path = out_dir_pred)

sv_AUC_data <- data.frame(Reduce("rbind", sv_tlpocv_res$inner_AUC))
sv_AUC_data <- rbind(sv_AUC_data, data.frame(t(sv_tlpocv_res$outer_AUC)))
colnames(sv_AUC_data) <- sapply(sv_tlpocv_res$best_features, length)
sv_AUC_data$Stage <- c(rep("Test data", length(sv_tlpocv_res$inner_AUC)), "Validation data")
sv_AUC_data_long <- melt(sv_AUC_data, id.vars = c("Stage"))
sv_AUC_data_long$variable <- as.numeric(as.character(sv_AUC_data_long$variable))
sv_AUC_data_long$value <- 1 - sv_AUC_data_long$value
p <- ggplot(data = sv_AUC_data_long, mapping = aes(x = variable, y = value, color = Stage, fill = Stage)) +
  stat_summary(fun.y = median, fun.ymax = function(x) quantile(x, p = 0.75), fun.ymin = function(x) quantile(x, p = 0.25), geom = "ribbon", alpha = 0.5, colour = NA) +
  stat_summary(fun.y = median, geom = "line") +
  ylim(0, 1) +
  ylab("Median AUC") +
  xlab("Nubmer of features") +
  ggtitle("SVM: Test and Validation AUC of nested TLPOCV-RFE") +
  theme_bw()
ggsave(filename = "SVM_TLPOCV_RFE_AUC.png", path = out_dir_pred, plot = p, width = 10, height = 5, units = "in")

f2_pos <- which(sapply(sv_tlpocv_res$best_features, length) == 2)
sv_int_ROC_data <- lapply(lapply(sv_tlpocv_res$int_sample_ranking, `[[`, f2_pos), function(e) ml.roc(ref = 1 - e[, 1], conf = e[, 2]))
sv_int_ROC_data <- lapply(sv_int_ROC_data, function(e) stepfun(e[, 1], c(0, e[, 2])))
sv_int_ROC_data <- lapply(sv_int_ROC_data, function(f) data.frame(FPR = seq(0, 1, by = 0.01), TPR = f(seq(0, 1, by = 0.01))))
sv_int_ROC_data <- data.frame(Reduce("rbind", sv_int_ROC_data))
if (min(sv_int_ROC_data[, 1]) != 0)
  sv_int_ROC_data <- rbind(c(0, 0), sv_int_ROC_data)
if (max(sv_int_ROC_data[, 1]) != 1)
  sv_int_ROC_data <- rbind(sv_int_ROC_data, c(1, 1))
colnames(sv_int_ROC_data) <- c("FPR", "TPR")
sv_int_ROC_data$Stage <- "Test data"
ranks_2feat <- sv_tlpocv_res$ext_sample_ranking[[f2_pos]]
sv_ext_ROC_data <- data.frame(ml.roc(ref = 1 - ranks_2feat[, 1], conf = ranks_2feat[, 2]))
if (min(sv_ext_ROC_data[, 2]) != 0)
  sv_ext_ROC_data <- rbind(c(0, 0), sv_ext_ROC_data)
if (max(sv_ext_ROC_data[, 1]) != 1)
  sv_ext_ROC_data <- rbind(sv_ext_ROC_data, c(1, 1))
colnames(sv_ext_ROC_data) <- c("FPR", "TPR")
sv_ext_ROC_data$Stage <- "Validation data"
sv_ROC_data <- rbind(sv_int_ROC_data, sv_ext_ROC_data)
rsize <- rel(3.5)
p_roc_svm_2feat <- ggplot(data = sv_ROC_data, mapping = aes(x = FPR, y = TPR, color = Stage, fill = Stage)) +
  stat_summary(data = sv_int_ROC_data,
               fun.y = mean, fun.ymax = function(x) quantile(x, p = 0.75), fun.ymin = function(x) quantile(x, p = 0.25),
               geom = "ribbon", alpha = 0.5, colour = NA) +
  stat_summary(data = sv_int_ROC_data,
               fun.y = median, geom = "line", size = rel(1.5)) +
  geom_line(data = sv_ext_ROC_data, size = rel(1.5)) +
  geom_segment(x = 0, y = 0, xend = 1, yend = 1, size = 0.5, inherit.aes = FALSE) +
  geom_text(x = 0.95, y = 0.15, hjust = "right", color = scales::hue_pal()(2)[1], size = rel(4.5),
            label = paste0("Test AUC = ", format(mean(subset(sv_AUC_data_long, variable == 2 & Stage == "Test data", "value")[[1]]), digits = 3, nsmall = 3))) +
  geom_text(x = 0.95, y = 0.10, hjust = "right", color = scales::hue_pal()(2)[2], size = rel(4.5),
            label = paste0("Validation AUC = ", format(subset(sv_AUC_data_long, variable == 2 & Stage == "Validation data", "value"), digits = 3, nsmall = 3))) +
  scale_x_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
  #ggtitle("SVM: Test and Validation ROC of nested TLPOCV RFE\nbest 2 feature set") +
  geom_text(label = "ROC curve\nSupport Vector Machine\n2 features", x = 0.5, y = 0.85, size = rel(4.5), inherit.aes = FALSE) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom", legend.direction = "horizontal", legend.text = element_text(size = rsize), text = element_text(size = rsize), plot.margin = margin(r = 20, l = 10))
ggsave(filename = "SVM_TLPOCV_RFE_ROC_2feat.png", path = out_dir_pred, plot = p_roc_svm_2feat, width = 6, height = 5, units = "in")

f2 <- sv_tlpocv_res$best_features[[which(sapply(sv_tlpocv_res$best_features, length) == 2)]]
p <- ggplot(data = subset(human_sepsis_data_ml, Day == 0, c("Survival", f2)), mapping = aes_string(x = f2[1], y = f2[2], color = quote(factor(Survival)))) + 
  geom_point() + 
  theme_bw()
ggsave(plot = p, filename = "SVM_class_separation_2feat.png", width = 5, height = 4, units = "in", path = out_dir_pred)

rg_nsep_roc_dfs <- lapply(rg_cf_roc, lapply, as.data.frame)
rg_nsep_roc_dfs <- lapply(rg_nsep_roc_dfs, lapply, function(e) stepfun(e[, 1], c(0, e[, 2])))
rg_nsep_roc_dfs <- lapply(rg_nsep_roc_dfs, lapply, function(f) data.frame(FPR = seq(0, 1, by = 0.01), TPR = f(seq(0, 1, by = 0.01))))
rg_nsep_roc_dfs <- lapply(rg_nsep_roc_dfs, lapply, function(df) return(rbind(data.frame(FPR = 0, TPR = 0), df)))
rg_nsep_roc_dfs <- lapply(lapply(rg_cf_roc, lapply, as.data.frame), lapply, function(df) return(rbind(data.frame(TPR = 0, FPR = 0), df)))
rg_nonsep_melt_roc <- melt(rg_nsep_roc_dfs, id.vars = 1:2)
rg_nonsep_melt_roc$Num_features <- factor(sapply(rg_ns_feats, length)[rg_nonsep_melt_roc$L1])
rg_nonsep_melt_roc$Group <- interaction(rg_nonsep_melt_roc$L2, rg_nonsep_melt_roc$Num_features)
rg_nonsep_roc <- ggplot(data = rg_nonsep_melt_roc, mapping = aes(x = FPR, y = TPR, color = Num_features, fill = Num_features)) +
  stat_summary(data = rg_nonsep_melt_roc,
               fun.y = mean, fun.ymax = function(x) quantile(x, p = 0.75), fun.ymin = function(x) quantile(x, p = 0.25),
               geom = "ribbon", alpha = 0.1, colour = NA) +
  stat_summary(data = rg_nonsep_melt_roc,
               fun.y = median, geom = "line") +
  geom_segment(x = 0, y = 0, xend = 1, yend = 1, size = 0.25, inherit.aes = FALSE) +
  scale_x_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
  ggtitle("RF: performance on nonseptic patients\ntrained on septic patients") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(plot = rg_nonsep_roc, filename = "RF_train_sep_test_nonsep_roc.png", width = 5, height = 4, units = "in", path = out_dir_pred)

sv_nsep_roc_dfs <- lapply(sv_cf_roc, lapply, as.data.frame)
sv_nsep_roc_dfs <- lapply(sv_nsep_roc_dfs, lapply, function(e) stepfun(e[, 1], c(0, e[, 2])))
sv_nsep_roc_dfs <- lapply(sv_nsep_roc_dfs, lapply, function(f) data.frame(FPR = seq(0, 1, by = 0.01), TPR = f(seq(0, 1, by = 0.01))))
sv_nsep_roc_dfs <- lapply(sv_nsep_roc_dfs, lapply, function(df) return(rbind(data.frame(FPR = 0, TPR = 0), df)))
sv_nonsep_melt_roc <- melt(sv_nsep_roc_dfs, id.vars = 1:2)
sv_nonsep_melt_roc$Num_features <- factor(sapply(sv_ns_feats, length)[sv_nonsep_melt_roc$L1])
sv_nonsep_melt_roc$Group <- interaction(sv_nonsep_melt_roc$L2, sv_nonsep_melt_roc$Num_features)
sv_nonsep_roc <- ggplot(data = sv_nonsep_melt_roc, mapping = aes(x = FPR, y = TPR, color = Num_features, fill = Num_features)) +
  stat_summary(data = sv_nonsep_melt_roc,
               fun.y = mean, fun.ymax = function(x) quantile(x, p = 0.75), fun.ymin = function(x) quantile(x, p = 0.25),
               geom = "ribbon", alpha = 0.1, colour = NA) +
  stat_summary(data = sv_nonsep_melt_roc,
               fun.y = median, geom = "line") +
  geom_segment(x = 0, y = 0, xend = 1, yend = 1, size = 0.25, inherit.aes = FALSE) +
  scale_x_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
  ggtitle("SVM: performance on nonseptic patients\ntrained on septic patients") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(plot = sv_nonsep_roc, filename = "SVM_train_sep_test_nonsep_roc.png", width = 5, height = 4, units = "in", path = out_dir_pred)

#make figure panel including heatmap of significant differences in metabolites & biochemical parameters (Fig 2A)
phm <- readRDS(file = "../data_stats/sig_feat_all_pats_heatmap.RData")
phm$gtable$grobs[[7]]$label <- "Metab.\ngroup"
phm$gtable$grobs[[7]]$gp$lineheight <- 0.7
panel2 <- plot_grid(phm$gtable, 
                    plot_grid(p_roc_rg_2feat, p_roc_svm_2feat, ncol = 2, nrow = 1, align = "h", axis = "tb", labels = c("", "C")), 
                    ncol = 1, nrow = 2, rel_heights = c(1, 0.4), labels = c("A", "B"))
ggsave(file = "all_features_sig_heatmap_ROC.png", path = out_dir_pred, plot = panel2, width = 9, height = 10 + 3.5, units = "in")
ggsave(file = "all_features_sig_heatmap_ROC.svg", path = out_dir_pred, plot = panel2, width = 9, height = 10 + 3.5, units = "in")
  