#Load libraries
library(data.table)
library(reshape2)
library(ggplot2)
library(matrixStats)
library(ranger)
library(missRanger)
library(kernlab)

#Import central functions
source("../function_definitions.R")

#Set path
out_dir <- "../../results/data_stats"

#Make sure paths exist
if (!dir.exists(out_dir))
  dir.create(out_dir)

#Import data
##Import clinical data
human_sepsis_data <- get_human_sepsis_data()

##Import experiment data
rat_sepsis_data <- get_rat_sepsis_data()

#Try prediction
##On human
###Build data set
human_sepsis_data_ml <- subset(x = human_sepsis_data, select = c(4,6:ncol(human_sepsis_data)))
human_sepsis_data_ml$Survival <- as.numeric(as.factor(human_sepsis_data_ml$Survival)) - 1
colnames(human_sepsis_data_ml) <- make.names(colnames(human_sepsis_data_ml))
###Impute missing values
human_sepsis_data_ml <- missRanger(human_sepsis_data_ml, pmm.k = 3, num.trees = 100)
###Test performance with cross-validation
fml <- Survival ~ .
num_folds <- 5
v.rg.npr <- list()
v.ks.npr <- list()
v.lm.npr <- list()
yPredRG <- rep(0, nrow(human_sepsis_data_ml))
yPredKS <- rep(0, nrow(human_sepsis_data_ml))
yPredLM <- rep(0, nrow(human_sepsis_data_ml))
fold_set <- ml.split.folds.strat(num_folds = num_folds, class = human_sepsis_data_ml$Survival)
for (fold in 1:num_folds){
  fold_learn_set <- human_sepsis_data_ml[-fold_set[[fold]],]
  rgCV <- ranger(data = fold_learn_set, dependent.variable.name = "Survival", num.trees = 500, write.forest = T, save.memory = F, classification = TRUE)
  ksCV <- ksvm(fml, fold_learn_set, type = "C-svc", kernel = "rbfdot", kpar = "automatic", C = 10)
  lmCV <- lm(formula = fml, data = fold_learn_set)
  v.rg.npr[[fold]] <- ml.npr(predict(rgCV, human_sepsis_data_ml[fold_set[[fold]],])$predictions, human_sepsis_data_ml$Survival[fold_set[[fold]]])
  v.ks.npr[[fold]] <- ml.npr(predict(ksCV, human_sepsis_data_ml[fold_set[[fold]],]), human_sepsis_data_ml$Survival[fold_set[[fold]]])
  v.lm.npr[[fold]] <- ml.npr(predict(lmCV, human_sepsis_data_ml[fold_set[[fold]],]), human_sepsis_data_ml$Survival[fold_set[[fold]]])
  yPredRG[fold_set[[fold]]] <- predict(rgCV, human_sepsis_data_ml[fold_set[[fold]],])$predictions
  yPredKS[fold_set[[fold]]] <- predict(ksCV, human_sepsis_data_ml[fold_set[[fold]],])
  yPredLM[fold_set[[fold]]] <- predict(lmCV, human_sepsis_data_ml[fold_set[[fold]],])
}
colMeans(rbindlist(v.rg.npr))
colMeans(rbindlist(v.ks.npr))
colMeans(rbindlist(v.lm.npr))

rg <- ranger(data = human_sepsis_data_ml, formula = fml, num.trees = 1500, write.forest = T, save.memory = F, classification = TRUE, importance = "permutation")
#rg$variable.importance[rg$variable.importance > 0.004]
