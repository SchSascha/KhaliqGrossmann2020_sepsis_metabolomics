#Load libraries
library(data.table)
library(reshape2)
library(ggplot2)
library(matrixStats)
library(ranger)
library(missRanger)
library(kernlab)
library(heatmaply)

#Import central functions
source("../function_definitions.R")

#Set path
out_dir <- "../../results"
out_dir_stats <- "../../results/data_stats"
out_dir_pred <- "../../results/data_pred_pheno"

#Make sure paths exist
if (!dir.exists(out_dir))
  dir.create(out_dir)
if (!dir.exists(out_dir_stats))
  dir.create(out_dir_stats)
if (!dir.exists(out_dir_pred))
  dir.create(out_dir_pred)

#Import data
##Import clinical data
human_sepsis_data <- get_human_sepsis_data()

##Import experiment data
rat_sepsis_data <- get_rat_sepsis_data()


#Get significantly different classes
human_sig_diff_res <- human_sig_diffs_along_days(human_sepsis_data, corr_fdr = TRUE)
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
strip_start <- which(colnames(human_sepsis_data_ml) == "Urea")
#human_sepsis_data_ml <- human_sepsis_data_ml[,-strip_start:-ncol(human_sepsis_data_ml)]
###No, strip metabolic variables
human_sepsis_data_ml <- human_sepsis_data_ml[,c(1:2,strip_start:ncol(human_sepsis_data_ml))]
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
rg.num.trees <- 500

##Try Survival forest
###Build specific data set from first measurement
human_surv_subset_ml <- cbind(data.frame(Day = human_sepsis_data$Day[data_ml_subset]), human_sepsis_data_ml_full)
surv_pats_max_days <- tapply(X = human_surv_subset_ml$Day, INDEX = human_surv_subset_ml$Patient, FUN = max)
surv_pats_min_days <- tapply(X = human_surv_subset_ml$Day, INDEX = human_surv_subset_ml$Patient, FUN = max)
human_surv_subset_ml <- human_surv_subset_ml[match(names(surv_pats_min_days), human_surv_subset_ml$Patient), ]
#human_sep_surv_ml$Day <- surv_pats_max_days[match(human_sep_surv_ml$Patient, names(surv_pats_max_days))] #For "time to fail" coding
human_surv_subset_ml$MaxDay <- surv_pats_max_days[match(human_surv_subset_ml$Patient, names(surv_pats_max_days))] #For interval coding: Day is interval start, MaxDay is end, 1 - Survival is event
human_surv_full_dat <- human_surv_subset_ml
human_survival_col_dat <- human_surv_subset_ml$Survival
human_surv_subset_ml <- subset(human_surv_subset_ml, Survival == 0, select = -1:-3)
human_surv_subset_ml$MaxDay <- human_surv_subset_ml$MaxDay + 1
human_surv_subset_ml$MaxDay[human_surv_subset_ml$Survival == 1] <- 100
human_surv_obj <- Surv(human_surv_subset_ml$MaxDay)
###Build RF
surv.rg <- ranger(human_surv_obj ~ ., data = human_surv_subset_ml, write.forest = TRUE, num.trees = rg.num.trees, importance = "permutation")
sp <- predict(surv.rg, cbind(data.frame(MaxDay = 0), subset(human_sepsis_data_ml_full, select = -1:-2)))
heatmaply(x = sp$survival, row_side_colors = data.frame(Survival = human_sepsis_data_ml_full$Survival))

##Set up iteration
num_repeats <- 6
day_tab <- table(human_sepsis_data$Day[data_ml_subset])
day_set <- rep(as.numeric(names(day_tab)), each = num_repeats)
rg.npr.repeat.df <- data.frame(day = day_set, tpr = 0, tnr = 0, fpr = 0, fnr = 0)
ks.npr.repeat.df <- data.frame(day = day_set, tpr = 0, tnr = 0, fpr = 0, fnr = 0)
lm.npr.repeat.df <- data.frame(day = day_set, tpr = 0, tnr = 0, fpr = 0, fnr = 0)
rg.red.npr.repeat.df <- data.frame(day = day_set, tpr = 0, tnr = 0, fpr = 0, fnr = 0)
ks.red.npr.repeat.df <- data.frame(day = day_set, tpr = 0, tnr = 0, fpr = 0, fnr = 0)
lm.red.npr.repeat.df <- data.frame(day = day_set, tpr = 0, tnr = 0, fpr = 0, fnr = 0)
##Try prediction on last day for each patient
# pat_num <- unique(human_sepsis_data$Patient)
# full_index <- match(pat_num, human_sepsis_data$Patient[order(human_sepsis_data$`Sample ID`, decreasing = TRUE)])
# full_index <- nrow(human_sepsis_data) - full_index + 1
# fold_set <- ml.split.folds.strat(num_folds = num_folds, class = human_sepsis_data_ml_full$Survival[full_index])
# last_day_ml_data <- human_sepsis_data_ml_full[full_index,]
# v.rg.npr <- list()
# for (fold in 1:num_folds){
#   fold_learn_set <- last_day_ml_data[-fold_set[[fold]],]
#   rgCV <- ranger(data = fold_learn_set, dependent.variable.name = "Survival", num.trees = rg.num.trees, write.forest = T, save.memory = F, classification = TRUE)
#   v.rg.npr[[fold]] <- ml.npr(predict(rgCV, last_day_ml_data[fold_set[[fold]],])$predictions, last_day_ml_data$Survival[fold_set[[fold]]])
# }
# colMeans(rbindlist(v.rg.npr)) #Only a little better than first day.
##Do the big CV loop
for (d in seq_along(day_set)){
  ###Select data with Day <= x, remove Patient ID
  human_sepsis_data_ml <- subset(human_sepsis_data_ml_full, subset = human_sepsis_data$Day[data_ml_subset] <= day_set[d], select = -1)
  ###Get Patient ID seperately
  patient_number <- subset(human_sepsis_data_ml_full, subset = human_sepsis_data$Day[data_ml_subset] <= day_set[d], select = Patient)[[1]]
  {
  ###Build gram matrix; only on learn set! build on test set seperately x = learn_set, y = test_set
  #num_data <- as.matrix(human_sepsis_data_ml[,-1])
  #num_data <- num_data[, !colAnyNAs(num_data)]
  #g_num_data <- kernelMatrix(kernel = vanilladot(), x = num_data, y = num_data)
  ###Supplant original data with gram matrix
  #human_sepsis_data_ml <- cbind(data.frame(Survival = human_sepsis_data_ml$Survival), g_num_data)
  }
  ###Test performance with cross-validation
  v.rg.npr <- list()
  v.ks.npr <- list()
  v.lm.npr <- list()
  yPredRG <- rep(0, nrow(human_sepsis_data_ml))
  yPredKS <- rep(0, nrow(human_sepsis_data_ml))
  yPredLM <- rep(0, nrow(human_sepsis_data_ml))
  fold_set <- ml.split.folds.quasistrat(num_folds = num_folds, class = human_sepsis_data_ml$Survival, non_strat_class = patient_number)
  #fold_set <- ml.split.folds.strat(num_folds = num_folds, class = human_sepsis_data_ml$Survival)
  for (fold in 1:num_folds){
    fold_learn_set <- human_sepsis_data_ml[-fold_set[[fold]],]
    rgCV <- ranger(data = fold_learn_set, dependent.variable.name = "Survival", num.trees = rg.num.trees, write.forest = T, save.memory = F, classification = TRUE)
    #ksCV <- ksvm(fml, fold_learn_set, type = "C-svc", kernel = "vanilladot")
    ksCV <- ksvm(fml, fold_learn_set, type = "C-svc", kernel = "rbfdot", kpar = "automatic", C = 10)
    lmCV <- lm(formula = fml, data = fold_learn_set)
    v.rg.npr[[fold]] <- ml.npr(predict(rgCV, human_sepsis_data_ml[fold_set[[fold]],])$predictions, human_sepsis_data_ml$Survival[fold_set[[fold]]])
    v.ks.npr[[fold]] <- ml.npr(predict(ksCV, human_sepsis_data_ml[fold_set[[fold]],]), human_sepsis_data_ml$Survival[fold_set[[fold]]])
    v.lm.npr[[fold]] <- ml.npr(predict(lmCV, human_sepsis_data_ml[fold_set[[fold]],]), human_sepsis_data_ml$Survival[fold_set[[fold]]])
    yPredRG[fold_set[[fold]]] <- predict(rgCV, human_sepsis_data_ml[fold_set[[fold]],])$predictions
    yPredKS[fold_set[[fold]]] <- predict(ksCV, human_sepsis_data_ml[fold_set[[fold]],])
    yPredLM[fold_set[[fold]]] <- predict(lmCV, human_sepsis_data_ml[fold_set[[fold]],])
  }
  # print("Classifier performance:")
  # print("Random Forest")
  # print(colMeans(rbindlist(v.rg.npr)))
  # print("C-SVM")
  # print(colMeans(rbindlist(v.ks.npr)))
  # print("Linear Model")
  # print(colMeans(rbindlist(v.lm.npr)))
  rg.npr.repeat.df[d,-1] <- colMeans(rbindlist(v.rg.npr))
  ks.npr.repeat.df[d,-1] <- colMeans(rbindlist(v.ks.npr))
  lm.npr.repeat.df[d,-1] <- colMeans(rbindlist(v.lm.npr))

  ###Get variable importance and validate on special data subsets
  #print("Most important variables according to Random Forest internal ranking:")
  rg <- ranger(data = human_sepsis_data_ml, dependent.variable.name = "Survival", num.trees = rg.num.trees * 3, write.forest = T, save.memory = F, classification = TRUE, importance = "permutation")
  #print(sort(rg$variable.importance, decreasing = TRUE)[1:10])
  var_importance <- rg$variable.importance
  #print(names(var_importance)[which(var_importance > quantile(x = var_importance, p = c(0.95))) + 1])
  # print("Random Forest validation on set of first day measurements when learnt on non-first day measurements:")
  # rg <- ranger(data = human_sepsis_data_ml[human_sepsis_data$Day > 0, ], dependent.variable.name = "Survival", num.trees = 1500, write.forest = T, save.memory = F, classification = TRUE, importance = "permutation")
  # print(t(ml.npr(predict(rg, human_sepsis_data_ml[human_sepsis_data$Day == 0, ])$predictions, human_sepsis_data_ml$Survival[human_sepsis_data$Day == 0])))
  # print("Random Forest validation on set of non-first day measurements when learnt on first day measurements:")
  # rg <- ranger(data = human_sepsis_data_ml[human_sepsis_data$Day == 0, ], dependent.variable.name = "Survival", num.trees = 1500, write.forest = T, save.memory = F, classification = TRUE, importance = "permutation")
  # print(t(ml.npr(predict(rg, human_sepsis_data_ml[human_sepsis_data$Day > 0, ])$predictions, human_sepsis_data_ml$Survival[human_sepsis_data$Day > 0])))
  
  ###Reduce data set to important variables
  human_sepsis_data_ml_red <- human_sepsis_data_ml
  human_sepsis_data_ml_red[, which(var_importance < quantile(x = var_importance, p = c(0.5))) + 1] <- NULL
  v.rg.npr <- list()
  v.ks.npr <- list()
  v.lm.npr <- list()
  yPredRG <- rep(0, nrow(human_sepsis_data_ml_red))
  yPredKS <- rep(0, nrow(human_sepsis_data_ml_red))
  yPredLM <- rep(0, nrow(human_sepsis_data_ml_red))
  for (fold in 1:num_folds){
    fold_learn_set <- human_sepsis_data_ml_red[-fold_set[[fold]],]
    rgCV <- ranger(data = fold_learn_set, dependent.variable.name = "Survival", num.trees = rg.num.trees, write.forest = T, save.memory = F, classification = TRUE)
    #ksCV <- ksvm(fml, fold_learn_set, type = "C-svc", kernel = "vanilladot")
    ksCV <- ksvm(fml, fold_learn_set, type = "C-svc", kernel = "rbfdot", kpar = "automatic", C = 10)
    lmCV <- lm(formula = fml, data = fold_learn_set)
    v.rg.npr[[fold]] <- ml.npr(predict(rgCV, human_sepsis_data_ml_red[fold_set[[fold]],])$predictions, human_sepsis_data_ml_red$Survival[fold_set[[fold]]])
    v.ks.npr[[fold]] <- ml.npr(predict(ksCV, human_sepsis_data_ml_red[fold_set[[fold]],]), human_sepsis_data_ml_red$Survival[fold_set[[fold]]])
    v.lm.npr[[fold]] <- ml.npr(predict(lmCV, human_sepsis_data_ml_red[fold_set[[fold]],]), human_sepsis_data_ml_red$Survival[fold_set[[fold]]])
    yPredRG[fold_set[[fold]]] <- predict(rgCV, human_sepsis_data_ml_red[fold_set[[fold]],])$predictions
    yPredKS[fold_set[[fold]]] <- predict(ksCV, human_sepsis_data_ml_red[fold_set[[fold]],])
    yPredLM[fold_set[[fold]]] <- predict(lmCV, human_sepsis_data_ml_red[fold_set[[fold]],])
  }
  # print("Classifier performance on data set of important variables:")
  # print("Random Forest")
  # print(colMeans(rbindlist(v.rg.npr)))
  # print("C-SVM")
  # print(colMeans(rbindlist(v.ks.npr)))
  # print("Linear Model")
  # print(colMeans(rbindlist(v.lm.npr)))
  rg.red.npr.repeat.df[d,-1] <- colMeans(rbindlist(v.rg.npr))
  ks.red.npr.repeat.df[d,-1] <- colMeans(rbindlist(v.ks.npr))
  lm.red.npr.repeat.df[d,-1] <- colMeans(rbindlist(v.lm.npr))
}

rg.npr.repeat.df$Method <- "RF"
ks.npr.repeat.df$Method <- "SVM"
lm.npr.repeat.df$Method <- "LM"
rg.red.npr.repeat.df$Method <- "RF"
ks.red.npr.repeat.df$Method <- "SVM"
lm.red.npr.repeat.df$Method <- "LM"

#Calculate accuracy
rg.npr.repeat.df$accuracy <- (rg.npr.repeat.df$tpr + rg.npr.repeat.df$tnr) / rowSums(rg.npr.repeat.df[, c(-1,-6)])
ks.npr.repeat.df$accuracy <- (ks.npr.repeat.df$tpr + ks.npr.repeat.df$tnr) / rowSums(ks.npr.repeat.df[, c(-1,-6)])
lm.npr.repeat.df$accuracy <- (lm.npr.repeat.df$tpr + lm.npr.repeat.df$tnr) / rowSums(lm.npr.repeat.df[, c(-1,-6)])
rg.red.npr.repeat.df$accuracy <- (rg.red.npr.repeat.df$tpr + rg.red.npr.repeat.df$tnr) / rowSums(rg.red.npr.repeat.df[, c(-1,-6)])
ks.red.npr.repeat.df$accuracy <- (ks.red.npr.repeat.df$tpr + ks.red.npr.repeat.df$tnr) / rowSums(ks.red.npr.repeat.df[, c(-1,-6)])
lm.red.npr.repeat.df$accuracy <- (lm.red.npr.repeat.df$tpr + lm.red.npr.repeat.df$tnr) / rowSums(lm.red.npr.repeat.df[, c(-1,-6)])

#Melt to long table
perf_by_day_data <- melt(rbindlist(list(rg.npr.repeat.df, ks.npr.repeat.df, lm.npr.repeat.df)), id.vars = c("day", "Method"))
perf_red_by_day_data <- melt(rbindlist(list(rg.red.npr.repeat.df, ks.red.npr.repeat.df, lm.red.npr.repeat.df)), id.vars = c("day", "Method"))

##Remove redundancy
perf_by_day_data <- subset(perf_by_day_data, subset = variable %in% c("fpr", "fnr", "accuracy"))
perf_red_by_day_data <- subset(perf_red_by_day_data, subset = variable %in% c("fpr", "fnr", "accuracy"))

#Combine for nice plot
perf_by_day_data$dimensions <- "All variables"
perf_red_by_day_data$dimensions <- "Important variables"
perf_data <- rbind(perf_by_day_data, perf_red_by_day_data)

perf_by_day_plot <- ggplot(subset(perf_data, variable %in% c("fnr", "fpr")), aes(x = day, y = value, color = Method)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean") +
  facet_grid(variable ~ dimensions) +
  ylim(0, 1) +
  xlab("Day included") +
  #ggtitle("Normal stratified Cross Validation may shows overfitting") +
  ggtitle("Quasi-stratified Cross Validation prevents overfitting") +
  theme_bw()
ggsave(filename = "prediction_performance_by_day_inclusion.png", plot = perf_by_day_plot, path = out_dir_pred, width = 5, height = 5, units = "in")

perf_by_day_acc_plot <- ggplot(subset(perf_data, variable == "accuracy"), aes(x = day, y = value, color = Method)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean") +
  facet_grid( ~ dimensions) +
  ylim(0, 1) +
  xlab("Day included") +
  #ggtitle("Normal stratified Cross Validation may shows overfitting") +
  #ggtitle("Quasi-stratified Cross Validation prevents overfitting") +
  theme_bw()
ggsave(filename = "prediction_performance_acc_by_day_inclusion.png", plot = perf_by_day_acc_plot, path = out_dir_pred, width = 5, height = 3, units = "in")
