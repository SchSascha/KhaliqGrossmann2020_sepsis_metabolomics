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
out_dir <- "../../results"
out_dir_stats <- "../../results/data_stats"

#Make sure paths exist
if (!dir.exists(out_dir))
  dir.create(out_dir)
if (!dir.exists(out_dir_stats))
  dir.create(out_dir_stats)

#Import data
##Import clinical data
human_sepsis_data <- get_human_sepsis_data()

##Import experiment data
rat_sepsis_data <- get_rat_sepsis_data()

#Try prediction
##On human
###Build data set
human_sepsis_data_ml <- subset(x = human_sepsis_data, select = c(4,6:ncol(human_sepsis_data)))
human_sepsis_data_ml$Survival <- as.numeric(as.factor(human_sepsis_data_ml$Survival)) - 1 #Dependent variable transformation
#human_sepsis_data_ml$`CAP / FP` <- as.factor(human_sepsis_data$`CAP / FP`)
colnames(human_sepsis_data_ml) <- make.names(colnames(human_sepsis_data_ml))
###Impute missing values
human_sepsis_data_ml <- missRanger(human_sepsis_data_ml, pmm.k = 3, num.trees = 100)
{#The following block greatly increases SVM performance and greatly decreases RF performance
###Scale values
#human_sepsis_data_ml[, -1] <- scale(human_sepsis_data_ml[, -1])
###Build gram matrix
#num_data <- as.matrix(human_sepsis_data_ml[,-1])
#num_data <- num_data[, !colAnyNAs(num_data)]
#g_num_data <- kernelMatrix(kernel = vanilladot(), x = num_data, y = num_data)
###Supplant original data with gram matrix
#human_sepsis_data_ml <- cbind(data.frame(Survival = human_sepsis_data_ml$Survival), g_num_data)
}
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
print("Classifier performance:")
print("Random Forest")
print(colMeans(rbindlist(v.rg.npr)))
print("C-SVM")
print(colMeans(rbindlist(v.ks.npr)))
print("Linear Model")
print(colMeans(rbindlist(v.lm.npr)))

###Get variable importance and validate on special data subsets
print("Most important variables according to Random Forest internal ranking:")
rg <- ranger(data = human_sepsis_data_ml, dependent.variable.name = "Survival", num.trees = 1500, write.forest = T, save.memory = F, classification = TRUE, importance = "permutation")
print(sort(rg$variable.importance, decreasing = TRUE)[1:10])
var_importance <- rg$variable.importance
print("Random Forest validation on set of first day measurements when learnt on non-first day measurements:")
rg <- ranger(data = human_sepsis_data_ml[human_sepsis_data$Day > 0, ], dependent.variable.name = "Survival", num.trees = 1500, write.forest = T, save.memory = F, classification = TRUE, importance = "permutation")
print(t(ml.npr(predict(rg, human_sepsis_data_ml[human_sepsis_data$Day == 0, ])$predictions, human_sepsis_data_ml$Survival[human_sepsis_data$Day == 0])))
print("Random Forest validation on set of non-first day measurements when learnt on first day measurements:")
rg <- ranger(data = human_sepsis_data_ml[human_sepsis_data$Day == 0, ], dependent.variable.name = "Survival", num.trees = 1500, write.forest = T, save.memory = F, classification = TRUE, importance = "permutation")
print(t(ml.npr(predict(rg, human_sepsis_data_ml[human_sepsis_data$Day > 0, ])$predictions, human_sepsis_data_ml$Survival[human_sepsis_data$Day > 0])))

###Reduce data set to important variables
human_sepsis_data_ml_red <- human_sepsis_data_ml
human_sepsis_data_ml[, which(var_importance < quantile(x = var_importance, p = c(0.25))) + 1] <- NULL
v.rg.npr <- list()
v.ks.npr <- list()
v.lm.npr <- list()
yPredRG <- rep(0, nrow(human_sepsis_data_ml_red))
yPredKS <- rep(0, nrow(human_sepsis_data_ml_red))
yPredLM <- rep(0, nrow(human_sepsis_data_ml_red))
for (fold in 1:num_folds){
  fold_learn_set <- human_sepsis_data_ml_red[-fold_set[[fold]],]
  rgCV <- ranger(data = fold_learn_set, dependent.variable.name = "Survival", num.trees = 500, write.forest = T, save.memory = F, classification = TRUE)
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
print("Classifier performance on data set of important variables:")
print("Random Forest")
print(colMeans(rbindlist(v.rg.npr)))
print("C-SVM")
print(colMeans(rbindlist(v.ks.npr)))
print("Linear Model")
print(colMeans(rbindlist(v.lm.npr)))