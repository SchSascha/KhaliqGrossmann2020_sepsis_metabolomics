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
tic()
tlpocv_res <- tlpocv_rfe_parallel(data_x = human_sepsis_data_ml[, -1:-5], data_y = human_sepsis_data_ml["Survival"], mc.cores = 7)
toc()

##Build sample pair list
tlpo_s_df <- t(combn(x = 1:nrow(human_sepsis_data_ml), m = 2))
##Run Feature Selection
feature_list <- list()
feature_list[[1]] <- list()
feature_list[[1]][[1]] <- colnames(human_sepsis_data_ml)[-1:-5]
inner_AUC <- list()
inner_AUC[[1]] <- rep(-100, length(feature_list[[1]]) - 1)
inner_AUC[2:nrow(tlpo_s_df)] <- inner_AUC[1]
outer_AUC <- rep(0, length(feature_list[[1]]) - 1)
pair_dir <- list()
pair_ext_dir <- list()
tic()
###External LPOCV
for (p_ext in 1:nrow(tlpo_s_df)){ #for each validation pair
  tlpo_int_s_df <- tlpo_s_df[(!tlpo_s_df[, 1] %in% tlpo_s_df[p_ext, ]) & !(tlpo_s_df[, 2] %in% tlpo_s_df[p_ext, ]), ] #exclude samples in validation pair from training
  ####Internal LPOCV/RFE
  feature_list[[p_ext]] <- list()
  pair_dir[[p_ext]] <- list()
  pair_ext_dir[[p_ext]] <- list()
  feature_list[[p_ext]][[1]] <- colnames(human_sepsis_data_ml)[-1:-5]
  feat_start <- feature_list[[p_ext]][[1]] #start with all features
  for (feat_count in seq_along(feat_start[-length(feat_start)])){ #the hard way: remove features one by one
    feat_sel <- feature_list[[p_ext]][[feat_count]] #update feature selection
    pair_dir[[p_ext]][[feat_count]] <- rep(-100, nrow(tlpo_int_s_df)) #edge direction in dominance graph
    pair_ext_dir[[p_ext]][[feat_count]] <- rep(-100, nrow(tlpo_int_s_df)) #edge direction in dominance graph for validation pair
    feat_imp <- matrix(0, nrow = nrow(tlpo_int_s_df), ncol = length(feat_sel)) #feature importance matrix for this run
    for (p_int in 1:nrow(tlpo_int_s_df)){ #for each test pair
      tr_x <- human_sepsis_data_ml[-c(tlpo_int_s_df[p_int, ], tlpo_s_df[p_ext, ]), feat_sel]
      tr_y <- human_sepsis_data_ml[-c(tlpo_int_s_df[p_int, ], tlpo_s_df[p_ext, ]), "Survival"]
      te_x <- human_sepsis_data_ml[tlpo_int_s_df[p_int, ], feat_sel]
      te_y <- human_sepsis_data_ml[tlpo_int_s_df[p_int, ], "Survival"]
      rg_te <- ranger(data = human_sepsis_data_ml[-c(tlpo_int_s_df[p_int, ], tlpo_s_df[p_ext, ]), c("Survival", feat_sel)], 
                      dependent.variable.name = "Survival", 
                      write.forest = T, 
                      probability = T, 
                      save.memory = F, 
                      importance = "impurity")
      pred_te <- predict(rg_te, te_x)
      pair_dir[[p_ext]][[feat_count]][p_int] <- heaviside(pred_te$predictions[1, 1], pred_te$predictions[2, 1]) #determine edge direction
      feat_imp[p_int, ] <- rg_te$variable.importance #collect variable importance
    }
    ####Calculate AUC as in Perez et al., 2018
    #####Get out degree for each vertex (sample)
    out_degree <- rep(0, max(tlpo_int_s_df[, 1:2]))
    d1 <- aggregate(x = pair_dir[[p_ext]][[feat_count]], by = list(v = tlpo_int_s_df[, 1]), FUN = sum) 
    d2 <- aggregate(x = 1 - pair_dir[[p_ext]][[feat_count]], by = list(v = tlpo_int_s_df[, 2]), FUN = sum)
    out_degree[d1$v] <- d1$x
    out_degree[d2$v] <- out_degree[d2$v] + d2$x
    tlpo_od <- apply(tlpo_int_s_df[, 1:2], 1:2, function(x) out_degree[x])
    tlpo_h <- heaviside(tlpo_od[, 1], tlpo_od[, 2])
    #####use only the positive-negative pairs for AUC (positive first sample, negative second sample!)
    pos1_idx <- human_sepsis_data_ml$Survival[tlpo_int_s_df[, 1]] == 1
    pos2_idx <- human_sepsis_data_ml$Survival[tlpo_int_s_df[, 2]] == 1
    neg1_idx <- human_sepsis_data_ml$Survival[tlpo_int_s_df[, 1]] == 0
    neg2_idx <- human_sepsis_data_ml$Survival[tlpo_int_s_df[, 2]] == 0
    sum_h <- sum(tlpo_h[pos1_idx & neg2_idx]) + sum(1 - tlpo_h[pos2_idx & neg1_idx])
    inner_AUC[[p_ext]][feat_count] <- sum_h / (sum(pos1_idx & neg2_idx) + sum(pos2_idx & neg1_idx))
    feat_imp <- colMeans(feat_imp) #rank features by importance
    feature_list[[p_ext]][[feat_count + 1]] <- feat_sel[feat_imp > sort(feat_imp)[1]] #down-select feature set
  }
}
toc()
###Find best feature set/most common feature combination for each feature set size
best_feat_set <- list()
for (feat_count in seq_along(feature_list[[1]])[-length(feature_list[[1]])]){
  fls <- unlist(lapply(feature_list, `[[`, feat_count)) #gets all feature sets for a fixed set size from all validation pairs
  fl_hist <- table(fls)
  fl_hist <- sort(fl_hist, decreasing = TRUE) #beware! feature names are not in original order
  best_feat_set[[feat_count]] <- names(fl_hist)[1:(length(feature_list[[1]][[feat_count]]))]
}
###Evaluate feature sets on external validation pairs
for (p_ext in 1:nrow(tlpo_s_df)){ #for each validation pair
  pair_ext_dir[[p_ext]] <- rep(0, length(best_feat_set))
  tlpo_int_s_df <- tlpo_s_df[(!tlpo_s_df[, 1] %in% tlpo_s_df[p_ext, ]) & !(tlpo_s_df[, 2] %in% tlpo_s_df[p_ext, ]), ] #exclude samples in validation pair from training
  for (feat_count in seq_along(best_feat_set)){
    feat_sel <- best_feat_set[[feat_count]]
    tr_x <- human_sepsis_data_ml[-tlpo_s_df[p_ext, ], feat_sel]
    tr_y <- human_sepsis_data_ml[-tlpo_s_df[p_ext, ], "Survival"]
    va_x <- human_sepsis_data_ml[tlpo_s_df[p_ext, ], feat_sel]
    va_y <- human_sepsis_data_ml[tlpo_s_df[p_ext, ], "Survival"]
    rg_va <- ranger(data = human_sepsis_data_ml[-tlpo_s_df[p_ext, ], c("Survival", feat_sel)], 
                    dependent.variable.name = "Survival", 
                    write.forest = T, 
                    probability = T, 
                    save.memory = F, 
                    importance = "impurity")
    pred_va <- predict(rg_va, va_x)
    pair_ext_dir[[p_ext]][feat_count] <- heaviside(pred_va$predictions[1, 1], pred_va$predictions[2, 1]) #determine edge direction
  }
}
####Calculate AUC as in Perez et al., 2018
#####Get out degree for each vertex (sample)
for (feat_count in seq_along(pair_ext_dir[[1]])){
  out_degree <- rep(0, max(tlpo_s_df[, 1:2]))
  d1 <- aggregate(x = sapply(pair_ext_dir, `[[`, feat_count), by = list(v = tlpo_s_df[, 1]), FUN = sum) 
  d2 <- aggregate(x = 1 - sapply(pair_ext_dir, `[[`, feat_count), by = list(v = tlpo_s_df[, 2]), FUN = sum)
  out_degree[d1$v] <- d1$x
  out_degree[d2$v] <- out_degree[d2$v] + d2$x
  tlpo_od <- apply(tlpo_int_s_df[, 1:2], 1:2, function(x) out_degree[x])
  tlpo_h <- heaviside(tlpo_od[, 1], tlpo_od[, 2])
  #####use only the positive-negative pairs for AUC (positive first sample, negative second sample!)
  pos1_idx <- human_sepsis_data_ml$Survival[tlpo_int_s_df[, 1]] == 1
  pos2_idx <- human_sepsis_data_ml$Survival[tlpo_int_s_df[, 2]] == 1
  neg1_idx <- human_sepsis_data_ml$Survival[tlpo_int_s_df[, 1]] == 0
  neg2_idx <- human_sepsis_data_ml$Survival[tlpo_int_s_df[, 2]] == 0
  sum_h <- sum(tlpo_h[pos1_idx & neg2_idx]) + sum(1 - tlpo_h[pos2_idx & neg1_idx])
  outer_AUC[feat_count] <- sum_h / (sum(pos1_idx & neg2_idx) + sum(pos2_idx & neg1_idx))
}


#-----------------------------

#Validate variable set only on Ferrario et al. data, use recent leave pair out cross validation for AUC and ROC generation (Perez et al., 2018)
##Build sample pair list
#tlpo_s_df <- as.data.frame(t(combn(x = 1:nrow(human_sepsis_data_ml), m = 2)))
tlpo_s_df <- expand.grid(s1 = which(human_validation_data$Survival == 1), s2 = which(human_validation_data$Survival == 0))
##Run classification
u_var_set_name_list <- unique(var_set_name_list)
for (n in seq_along(u_var_set_name_list)){ #only unique sets, bc. mostly deterministic results
  dir_key <- paste0("dir_pset", n)
  tlpo_s_df[[dir_key]] <- 0
  for (p in 1:nrow(tlpo_s_df)){
    hvd <- human_sepsis_data_ml[c("Survival", u_var_set_name_list[[n]])]
    hvd_tr <- hvd[-as.numeric(tlpo_s_df[p, 1:2]), ]
    hvd_te <- hvd[as.numeric(tlpo_s_df[p, 1:2]), ]
    rg_tlpo <- ranger(data = hvd_tr, dependent.variable.name = "Survival", num.trees = 500, write.forest = T, save.memory = F, probability = TRUE)
    rg_pred <- predict(rg_tlpo, hvd_te)
    rg_tlpo_s_df[[dir_key]][p] <- rg_pred$predictions[1, 1] > rg_pred$predictions[2, 1]
  }
}
##Calculate AUC as in Airola et al., 2011
tournament_score <- list()
for (n in 3:ncol(tlpo_s_df)){
  out_degree <- rep(0, max(tlpo_s_df[, 1:2]))
  pos_out <- tapply(X = as.numeric(tlpo_s_df[, n] == 1), FUN = sum, INDEX = tlpo_s_df[, 1])
  neg_out <- tapply(X = as.numeric(tlpo_s_df[, n] == 0), FUN = sum, INDEX = tlpo_s_df[, 2])
  out_degree[as.numeric(names(pos_out))] <- pos_out
  out_degree[as.numeric(names(neg_out))] <- neg_out
  tlpo_od <- apply(tlpo_s_df[, 1:2], 1:2, function(x) out_degree[x])
  tlpo_h <- heaviside(tlpo_od[, 1], tlpo_od[, 2])
  tournament_score[[n-2]] <- sum(tlpo_h) / nrow(tlpo_s_df)
}

##Calculate tournament scores
tournament_score <- list()
for (n in 3:ncol(tlpo_s_df)){
  out_degree <- c(table(tlpo_s_df[, c(1, n)])[, 2], 0) + c(0, table(tlpo_s_df[, c(2, n)])[, 1]) #shift with c() is necessary bc. first counts 1..n-1 and second 2..n
  tlpo_od <- apply(tlpo_s_df[, 1:2], 1:2, function(x) out_degree[x])
  tlpo_h <- heaviside(tlpo_od[, 1], tlpo_od[, 2])
  tournament_score[[n-2]] <- 1 / (length(tlpo_h) + 1) * sum(tlpo_h)
}




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
