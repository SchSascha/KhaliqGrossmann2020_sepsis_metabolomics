#Load libraries
library(data.table)
library(reshape2)
library(ggplot2)
library(matrixStats)
library(ranger)
library(missRanger)
library(kernlab)
#library(mixOmics)
library(corpcor)
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

##Import clinical validation data
human_validation_data <- get_Ferrario_validation_data()

##Filter out patients without sepsis
human_nonsepsis_data <- human_sepsis_data[human_sepsis_data$`CAP / FP` == "-", ]
human_sepsis_data <- human_sepsis_data[human_sepsis_data$`CAP / FP` != "-", ]

##Filter out unimportant rows and columns from validation data
human_validation_data <- subset(human_validation_data, Day == 0, c(-1, -2, -4))
colnames(human_validation_data)[grep(x = colnames(human_validation_data), pattern = "Survival")] <- "Survival"

##Filter outlier in Ferrario data
human_validation_data <- subset(human_validation_data, subset = C4 < 400)

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
human_validation_data$Survival <- as.numeric(factor(human_validation_data$Survival, levels = sort(unique(human_validation_data$Survival)))) - 1
human_sepsis_data_ml <- human_sepsis_data_ml[rowSums(is.na(human_sepsis_data_ml)) < 9, ]
human_sepsis_data_ml <- human_sepsis_data_ml[, c(1:5, 5 + which(!colAnys(human_sepsis_data_ml == 0)[-1:-5]))]
###Strip phenomenological variables
human_sepsis_data_ml <- subset(human_sepsis_data_ml, select = c(colnames(human_sepsis_data_ml)[1:5], intersect(sig_t_class, colnames(human_sepsis_data_ml))))
colnames(human_sepsis_data_ml) <- make.names(colnames(human_sepsis_data_ml))
colnames(human_validation_data) <- make.names(colnames(human_validation_data))
###Impute missing values
human_sepsis_data_ml[, -1:-5] <- missRanger(human_sepsis_data_ml[, -1:-5], pmm.k = 3, num.trees = 100)
#human_validation_data[, -1] <- missRanger(human_validation_data[, -1], pmm.k = 3, num.trees = 100)
###Keep only common metabolites
shared_metabs <- intersect(colnames(human_sepsis_data_ml)[-1:-2], colnames(human_validation_data)[-1])
human_sepsis_data_ml_shared <- human_sepsis_data_ml[, c(colnames(human_sepsis_data_ml)[1:5], shared_metabs)]
human_validation_data <- human_validation_data[, c("Survival", shared_metabs)]
#human_sepsis_data_ml_shared[, -1:-5] <- human_sepsis_data_ml_shared[, 5 + order(colnames(human_sepsis_data_ml_shared[, -1:-5]))]
#human_validation_data[, -1] <- human_validation_data[, 1 + order(colnames(human_validation_data[, -1]))]
###Do CORAL
####Get mean and cov without outliers
tr_box_stats <- lapply(human_sepsis_data_ml_shared[, -1:-5], boxplot.stats)
vl_box_stats <- lapply(human_validation_data[, -1], boxplot.stats)
tr_in_idx <- lapply(1:(ncol(human_sepsis_data_ml_shared) - 5),
                     function(x){
                       d <- human_sepsis_data_ml_shared[, x + 5] 
                       which(!(d %in% tr_box_stats[[x]]$out))
                     })
vl_in_idx <- lapply(1:(ncol(human_validation_data) - 1),
                     function(x){
                       d <- human_validation_data[, x + 1] 
                       which(!(d %in% vl_box_stats[[x]]$out))
                     })
tr_means <- sapply(1:(ncol(human_sepsis_data_ml_shared) - 5), 
                   function(x){ 
                     d <- human_sepsis_data_ml_shared[, x + 5] 
                     d <- d[tr_in_idx[[x]]]
                     mean(d) 
                   })
vl_means <- sapply(1:(ncol(human_validation_data) - 1), 
                   function(x){ 
                     d <- human_validation_data[, x + 1]
                     d <- d[vl_in_idx[[x]]]
                     mean(d) 
                   })
tr_sds <- sapply(1:(ncol(human_sepsis_data_ml_shared) - 5), 
                   function(x){ 
                     d <- human_sepsis_data_ml_shared[, x + 5] 
                     d <- d[tr_in_idx[[x]]]
                     sd(d) 
                   })
vl_sds <- sapply(1:(ncol(human_validation_data) - 1), 
                   function(x){ 
                     d <- human_validation_data[, x + 1]
                     d <- d[vl_in_idx[[x]]]
                     sd(d) 
                   })
human_validation_data[, -1] <- scale(human_validation_data[, -1], center = vl_means, scale = vl_sds)
human_sepsis_data_ml_shared[, -1:-5] <- scale(human_sepsis_data_ml_shared[, -1:-5], center = tr_means, scale = tr_sds)
# human_validation_data[, -1] <- scale(human_validation_data[, -1], center = TRUE, scale = TRUE)
# human_sepsis_data_ml_shared[, -1:-5] <- scale(human_sepsis_data_ml_shared[, -1:-5], center = TRUE, scale = TRUE)
tr_cov <- lapply(seq_along(tr_in_idx),
                  function(x){
                    sapply(seq_along(tr_in_idx), 
                           function(y){
                             in_idx <- intersect(tr_in_idx[[x]], tr_in_idx[[y]])
                             cov(x = human_sepsis_data_ml_shared[in_idx, x + 5], y = human_sepsis_data_ml_shared[in_idx, y + 5])
                           })
                  })
vl_cov <- lapply(seq_along(vl_in_idx),
                  function(x){
                    sapply(seq_along(vl_in_idx), 
                           function(y){
                             in_idx <- intersect(vl_in_idx[[x]], vl_in_idx[[y]])
                             cov(x = human_validation_data[in_idx, x + 1], y = human_validation_data[in_idx, y + 1])
                           })
                  })

hsdms <- human_sepsis_data_ml_shared
# C_source <- Reduce("rbind", tr_cov)
# # C_source <- as.matrix(cov.shrink(hsdms[, -1:-5]))
# svd_C_source <- svd(C_source)
# isqrt_C_source <- with(svd_C_source, u %*% diag(1/sqrt(d)) %*% t(v))
# C_target <- Reduce("rbind", vl_cov)
# # C_target <- cov.shrink(human_validation_data[, -1])
# svd_C_target <- svd(C_target)
# sqrt_C_target <- with(svd_C_target, sqrt_C_target <- u %*% diag(sqrt(d)) %*% t(v))
# human_sepsis_data_ml_shared[, -1:-5] <- as.matrix(hsdms[, -1:-5]) %*% isqrt_C_source %*% sqrt_C_target
##Set important parameters
fml <- Survival ~ .
num_folds <- 5
rg.num.trees <- 500
##Set up iteration
num_repeats <- 20
day_tab <- table(human_sepsis_data_ml_shared$Day)
day_set <- rep(as.numeric(names(day_tab)), each = num_repeats)
var_range <- c(2:min(6, length(shared_metabs)))
var_set <- rep(var_range, each = num_repeats)
var_set_name_list <- list()
tot_n_var <- ncol(human_sepsis_data_ml_shared)-5
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
rg.red.npr.FVal.df <- data.frame(var = var_set, tpr = 0, tnr = 0, fpr = 0, fnr = 0, ppv = 0, npv = 0)
ks.red.npr.FVal.df <- data.frame(var = var_set, tpr = 0, tnr = 0, fpr = 0, fnr = 0, ppv = 0, npv = 0)
lm.red.npr.FVal.df <- data.frame(var = var_set, tpr = 0, tnr = 0, fpr = 0, fnr = 0, ppv = 0, npv = 0)
ol.red.npr.FVal.df <- data.frame(var = var_set, tpr = 0, tnr = 0, fpr = 0, fnr = 0, ppv = 0, npv = 0)
rg.red.acc.FVal.df <- data.frame(var = var_set, acc = 0)
ks.red.acc.FVal.df <- data.frame(var = var_set, acc = 0)
lm.red.acc.FVal.df <- data.frame(var = var_set, acc = 0)
ol.red.acc.FVal.df <- data.frame(var = var_set, acc = 0)
rg.red.auc.FVal.df <- data.frame(var = var_set, auc = 0)
ks.red.auc.FVal.df <- data.frame(var = var_set, auc = 0)
lm.red.auc.FVal.df <- data.frame(var = var_set, auc = 0)
ol.red.auc.FVal.df <- data.frame(var = var_set, auc = 0)
##Do the big CV loop
d <- 1
r_count <- 1
{
  ###Select data with Day <= x, remove Patient ID ("select = -1")
  human_sepsis_data_ml <- subset(human_sepsis_data_ml_shared, subset = Day <= day_set[d])
  hsdm_prop_cols <- which(colnames(human_sepsis_data_ml)[1:5] != "Survival")
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
      fold_learn_set <- human_sepsis_data_ml[-fold_set[[fold]], -hsdm_prop_cols]
      fold_test_set <- human_sepsis_data_ml[fold_set[[fold]], -hsdm_prop_cols]
      r_fold_learn_set <- fold_learn_set
      r_fold_test_set <- fold_test_set
      r_fold_learn_set$Survival <- factor(r_fold_learn_set$Survival, levels = sort(unique(r_fold_learn_set$Survival)))
      r_fold_test_set$Survival <- factor(r_fold_test_set$Survival, levels = sort(unique(r_fold_learn_set$Survival)))
      class_count <- table(fold_learn_set$Survival)
      case_weights <- rep(class_count[1]/class_count[2], nrow(fold_learn_set))
      case_weights[fold_learn_set$Survival == 0] <- class_count[2]/class_count[1]
      case_weights_int <- class_count[match(fold_learn_set$Survival, names(class_count))]
      rgCV <- ranger(data = r_fold_learn_set, dependent.variable.name = "Survival", num.trees = rg.num.trees, write.forest = T, save.memory = F, probability = TRUE, case.weights = case_weights)
      ksCV <- ksvm(fml, fold_learn_set, type = "C-svc", kernel = "vanilladot", scaled = FALSE)
      lmCV <- glm(formula = fml, data = fold_learn_set, family = binomial(link = "probit"), weights = case_weights_int)
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
    rg <- ranger(data = human_sepsis_data_ml[, -hsdm_prop_cols], dependent.variable.name = "Survival", num.trees = rg.num.trees * 3, write.forest = T, save.memory = F, classification = TRUE, importance = "permutation")
    var_importance <- rg$variable.importance
    ###Reduce data set to important variables
    human_sepsis_data_ml_red[, 5 + which(var_importance < sort(var_importance, decreasing = TRUE)[var_set[v]])] <- NULL
    ###Save variable set for later analysis
    var_set_name_list[[v]] <- colnames(human_sepsis_data_ml_red)[-1:-5]
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
      fold_learn_set <- human_sepsis_data_ml_red[-fold_set[[fold]], -hsdm_prop_cols]
      fold_test_set <- human_sepsis_data_ml_red[fold_set[[fold]], -hsdm_prop_cols]
      ###Get variable importance and validate on special data subsets
      r_count <- r_count + 1
      ####Continue with the other stuff
      r_fold_learn_set <- fold_learn_set
      r_fold_test_set <- fold_test_set
      r_fold_learn_set$Survival <- factor(r_fold_learn_set$Survival, levels = sort(unique(r_fold_learn_set$Survival)))
      r_fold_test_set$Survival <- factor(r_fold_test_set$Survival, levels = sort(unique(r_fold_learn_set$Survival)))
      class_count <- table(fold_learn_set$Survival)
      case_weights <- rep(class_count[1]/class_count[2], nrow(fold_learn_set))
      case_weights[fold_learn_set$Survival == 0] <- class_count[2]/class_count[1]
      case_weights_int <- class_count[match(fold_learn_set$Survival, names(class_count))]
      rgCV <- ranger(data = r_fold_learn_set, dependent.variable.name = "Survival", num.trees = rg.num.trees, write.forest = T, save.memory = F, probability = TRUE, case.weights = case_weights)
      ksCV <- ksvm(fml, fold_learn_set, type = "C-svc", kernel = "vanilladot", scaled = FALSE)
      lmCV <- glm(formula = fml, data = fold_learn_set, family = binomial(link = "probit"), weights = case_weights_int)
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
    
    ###Validate result on data from Ferrario et al.
    class_count <- table(human_sepsis_data_ml_red$Survival)
    case_weights <- rep(class_count[1]/class_count[2], nrow(human_sepsis_data_ml_red))
    case_weights[human_sepsis_data_ml_red$Survival == 0] <- class_count[2]/class_count[1]
    case_weights_int <- class_count[match(human_sepsis_data_ml_red$Survival, names(class_count))]
    rgFV <- ranger(data = human_sepsis_data_ml_red[, -hsdm_prop_cols], dependent.variable.name = "Survival", num.trees = rg.num.trees, write.forest = T, save.memory = F, probability = TRUE)
    ksFV <- ksvm(fml, human_sepsis_data_ml_red[, -hsdm_prop_cols], type = "C-svc", kernel = "vanilladot", scaled = FALSE)
    lmFV <- glm(formula = fml, data = human_sepsis_data_ml_red[, -hsdm_prop_cols], family = binomial(link = "probit"), control = glm.control(maxit = 100))
    
    rg.red.npr.FVal.df[v,-1] <- ml.npr(predict(rgFV, human_validation_data)$predictions[, 2] > 0.5, human_validation_data$Survival)
    ks.red.npr.FVal.df[v,-1] <- ml.npr(predict(ksFV, human_validation_data), human_validation_data$Survival)
    lm.red.npr.FVal.df[v,-1] <- ml.npr(predict(lmFV, human_validation_data[,-1]), human_validation_data$Survival)

    rg.red.acc.FVal.df[v,-1] <- ml.acc(predict(rgFV, human_validation_data)$predictions[, 2] > 0.5, human_validation_data$Survival)
    ks.red.acc.FVal.df[v,-1] <- ml.acc(predict(ksFV, human_validation_data), human_validation_data$Survival)
    lm.red.acc.FVal.df[v,-1] <- ml.acc(predict(lmFV, human_validation_data[,-1]), human_validation_data$Survival)
    
    rg.red.auc.FVal.df[v,-1] <- ml.auc(ref = human_validation_data$Survival, conf = predict(rgFV, human_validation_data)$predictions[, 2])
    ks.red.auc.FVal.df[v,-1] <- ml.auc(ref = human_validation_data$Survival, conf = predict(ksCV, human_validation_data, type = "decision"))
    lm.red.auc.FVal.df[v,-1] <- ml.auc(ref = human_validation_data$Survival, conf = predict(lmCV, human_validation_data[,-1]))
  }
}

#Validate variable set only on Ferrario et al. data, use recent leave pair out cross validation for AUC and ROC generation (Perez et al., 2018)
##Build sample pair list
tlpo_s_df <- expand.grid(s1 = 1:nrow(human_validation_data), s2 = 1:nrow(human_validation_data))
tlpo_s_df <- tlpo_s_df[tlpo_s_df$s1 != tlpo_s_df$s2, ]
##Run classification
tlpo_res <- list()
u_var_set_name_list <- unique(var_set_name_list)
for (n in seq_along(u_var_set_name_list)){ #only unique sets, bc. mostly deterministic results
  dir_key <- paste0("dir_pset", n)
  tlpo_s_df[[dir_key]] <- 0
  for (p in 1:nrow(tlpo_s_df)){
    hvd <- human_validation_data[c("Survival", u_var_set_name_list[[n]])]
    hvd_tr <- hvd[-as.numeric(tlpo_s_df[p, 1:2]), ]
    hvd_te <- hvd[as.numeric(tlpo_s_df[p, 1:2]), ]
    rg_tlpo <- ranger(data = hvd_tr, dependent.variable.name = "Survival", num.trees = rg.num.trees, write.forest = T, save.memory = F, probability = TRUE)
    rg_pred <- predict(rg_tlpo, hvd_te)
    tlpo_s_df[[dir_key]][p] <- rg_pred$predictions[1, 1] > rg_pred$predictions[2, 1]
  }
}
##Calculate tournament scores
#TODO: find error, AUC is > 1
tournament_score <- list()
for (n in 3:ncol(tlpo_s_df)){
  out_degree <- table(tlpo_s_df[, c(1, n)])[, 2]
  tlpo_od <- apply(tlpo_s_df[, 1:2], 1:2, function(x) out_degree[x])
  tlpo_h <- heaviside(tlpo_od[, 1], tlpo_od[, 2])
  tournament_score[[n-2]] <- 1 / (sum(human_validation_data$Survival == 0) * sum(human_validation_data$Survival == 1)) * sum(tlpo_h)
}

# hsd <- subset(human_sepsis_data, Day == 0)
# hvd <- human_validation_data
# cns1 <- make.names(colnames(human_sepsis_data))
# cns2 <- make.names(colnames(human_validation_data))
# colnames(hsd) <- cns1
# colnames(hvd) <- cns2
# hsd[, -1:-5] <- missRanger(hsd[, -1:-5])
# hvd[, -1] <- missRanger(hvd[, -1])
# smet <- intersect(cns1, cns2)
# mxOmX <- rbind(hsd[smet][, -1], hvd[smet][, -1])
# mxOmY <- c(hsd$Survival, c("NS", "S")[1 + hvd$Survival])
# study <- rep(1:2, times = c(nrow(hsd), nrow(hvd)))
# mxo_res <- mixOmics(X = as.matrix(mxOmX), Y = factor(mxOmY), study = study, ncomp = 3, scale = TRUE, keepX = c(60, 30, 5))
# mxo_perf <- perf(object = mxo_res, validation = "MFold", folds = 10, nrepeat = 5, auc = TRUE)
# lapply(lapply(mxo_perf$auc, `[[`, 1), `[[`, 1)

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
