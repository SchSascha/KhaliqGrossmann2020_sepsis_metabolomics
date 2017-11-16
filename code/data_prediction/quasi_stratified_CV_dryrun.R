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

human_sepsis_data_ml <- subset(x = human_sepsis_data, select = c(2,4,6:ncol(human_sepsis_data)))
human_sepsis_data_ml$Survival <- as.numeric(as.factor(human_sepsis_data_ml$Survival)) - 1 #Dependent variable transformation
###Strip phenomenological variables
strip_start <- which(colnames(human_sepsis_data_ml) == "Urea")
human_sepsis_data_ml <- human_sepsis_data_ml[,-strip_start:-ncol(human_sepsis_data_ml)]
###No, strip metabolic variables
#human_sepsis_data_ml <- human_sepsis_data_ml[,c(1,strip_start:ncol(human_sepsis_data_ml))]
data_ml_subset <- apply(human_sepsis_data_ml, 1, function(x){ sum(!is.na(x)) > length(x) - 9})
human_sepsis_data_ml <- human_sepsis_data_ml[data_ml_subset,]
#human_sepsis_data_ml$`CAP / FP` <- as.factor(human_sepsis_data$`CAP / FP`)
colnames(human_sepsis_data_ml) <- make.names(colnames(human_sepsis_data_ml))
###Impute missing values
human_sepsis_data_ml_full <- missRanger(human_sepsis_data_ml, pmm.k = 3, num.trees = 100)
###Scale values
human_sepsis_data_ml_full[, -1:-2] <- scale(human_sepsis_data_ml_full[, -1:-2])
##Set important variables
num_folds <- 5
num_repeats <- 6
day_tab <- table(human_sepsis_data$Day[data_ml_subset])
day_set <- rep(as.numeric(names(day_tab)), each = num_repeats)

for (d in seq_along(day_set)){
  ###Select data with Day <= x, remove Patient ID
  human_sepsis_data_ml <- subset(human_sepsis_data_ml_full, subset = human_sepsis_data$Day[data_ml_subset] <= day_set[d], select = -1)
  ###Get Patient ID seperately
  patient_number <- subset(human_sepsis_data_ml_full, subset = human_sepsis_data$Day[data_ml_subset] <= day_set[d], select = Patient)[[1]]
  
  ###Test quasi-stratified cross-validation
  fold_set <- ml.split.folds.quasistrat(num_folds = num_folds, class = human_sepsis_data_ml$Survival, non_strat_class = patient_number, num_sampling_steps = 10000)
}