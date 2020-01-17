#' Read the clinical data set on sepsis cases
#'
#' @return the data.frame with data, column names are identical to file content (non-standard for R)
#' @export
#'
#' @examples
get_human_sepsis_data <- function(){
  require(matrixStats)
  data <- read.csv(file = "../../data/measurements/Summary human sample data_2019_10_24.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, dec = ",", check.names = FALSE)
  data <- data[, -(which(colSds(as.matrix(data[, -1:-5])) == 0) + 5)] #remove columns with fixed values
  group <- data[["CAP / FP"]] #add group column with combinations of Survival and Sepsis/Nonsepsis
  group[group == "-"] <- "non-Septic"
  group[group %in% c("CAP", "FP")] <- "Septic"
  group <- paste0(group, "-", data$Survival)
  data <- cbind(data[, 1:5], data.frame(Group = group, stringsAsFactors = FALSE), data[, -1:-5])
  pat_exclude <- c(2, 5, 10, 19, 27, 29, 37, 40) 
  data <- subset(data, !Patient %in% pat_exclude) #remove patients we later found to be conceptually too close to sepsis
  data <- remove_bile_acids_from_measurements(data) #remove bile acids because wrong mesaurement technique for them
  return(data)
}

#' Read the clinical data set on sepsis shock cases from Ferrario et al., Scientific Reports, 2016
#'
#' @return data.frame with measurements for day1 and day7, formatted the same way as data from Mervyn
#' @export
#'
#' @examples
get_Ferrario_validation_data <- function(){
  #Read data
  data_day1 <- read.csv(file = "../../data/measurements/Ferrario_validation_day1.csv", header = FALSE, skip = 3, sep = "\t", stringsAsFactors = FALSE, dec = ",", check.names = FALSE)
  data_day7 <- read.csv(file = "../../data/measurements/Ferrario_validation_day7.csv", header = FALSE, skip = 3, sep = "\t", stringsAsFactors = FALSE, dec = ",", check.names = FALSE)
  header_day1 <- read.csv(file = "../../data/measurements/Ferrario_validation_day1.csv", header = FALSE, sep = "\t", stringsAsFactors = FALSE, dec = ",", check.names = FALSE, nrows = 3)
  header_day7 <- read.csv(file = "../../data/measurements/Ferrario_validation_day7.csv", header = FALSE, sep = "\t", stringsAsFactors = FALSE, dec = ",", check.names = FALSE, nrows = 3)
  #Build Patient characteristic
  data <- data.frame(Patient = unlist(c(header_day1[1, -1], header_day7[1, -1])), Survival28 = unlist(c(header_day1[2, -1], header_day7[2, -1])), Survival90 = unlist(c(header_day1[3, -1], header_day7[3, -1])), Day = rep(c(0,6), each = 20), stringsAsFactors = FALSE)
  #Transpose concentration data
  rownames(data_day1) <- data_day1$V1
  rownames(data_day7) <- data_day7$V1
  data_day1$V1 <- NULL
  data_day7$V1 <- NULL
  data_day1 <- t(data_day1)
  data_day7 <- t(data_day7)
  #Combine and edit metabolite names
  mets <- union(colnames(data_day1), colnames(data_day7))
  m1 <- match(colnames(data_day1), mets)
  m2 <- match(colnames(data_day7), mets)
  mets <- trimws(mets)
  mets <- sub(pattern = "LPC", replacement = "lysoPC", x = mets)
  mets <- sub(pattern = "SM OH", replacement = "SM (OH)", x = mets)
  mets[mets == "Sugars"] <- "H1"
  mets[mets == "Met:SO"] <- "Met-SO"
  mets[mets == "Ac Orn"] <- "Ac-Orn"
  mets[mets == "Alpha AAA"] <- "alpha-AAA"
  mets[mets == "T4 OH Pro"] <- "t4-OH-Pro"
  mets[mets == "Total DMA"] <- "total DMA"
  #Copy concentrations into data.frame
  data <- cbind(data, matrix(NA, ncol = length(mets), nrow = nrow(data)))
  colnames(data)[-1:-4] <- mets
  data[1:nrow(data_day1), m1 + 4] <- data_day1
  data[(nrow(data_day1) + 1):nrow(data), m2 + 4] <- data_day7
  data <- data[order(data$Patient),]
  return(data)
}

#' Read the legend to the clinical data set on sepsis cases. The legend has two group assignments of different coarseness.
#'
#' @return data.frame, first column stores the column names in the measurement data table, second the fine grained grouping, third the coarse grained grouping
#' @export
#'
#' @examples
get_human_sepsis_legend <- function(){
  data <- read.csv(file = "../../data/measurements/Legend human sample data pheno class.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, dec = ",", check.names = FALSE)
  data <- rbind(data[1:5, ], c("Group", "Combined Sepsis status and Survival", rep("", ncol(data)-2)), data[-1:-5, ])
  data$name[data$name == "Carnitine "] <- "Carnitine" #beware the random space
  data <- remove_bile_acids_from_legend(data)
  return(data)
}

#' Read the grouping of phenomenological variables.
#'
#' @return data.frame, table of associated groups
#' @export
#'
#' @examples
get_human_pheno_var_groups <- function(){
  data <- read.csv(file = "../../data/measurements/Human phenomenological groups.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, dec = ",", check.names = FALSE)
  return(data)
}

#' Read patient meta data table
#'
#' @return data.frame with patient meta data
#' @export
#'
#' @examples
get_human_meta_data <- function(){
  library(stringi)
  data <- read.csv(file = "../../data/measurements/Raw data Germany Human samples_2019_10_24.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, dec = ",", check.names = FALSE)
  data$`Patient ID` <- stri_trim_both(data$`Patient ID`)
  data$Diagnosis <- stri_trim_both(data$Diagnosis)
  return(data)
}

#' Read patient illenss data - why they where in the ICU
#'
#' @return
#' @export
#'
#' @examples
get_human_illness_data <- function(){
  library(stringi)
  data <- read.csv(file = "../../data/measurements/Vic_patients_downselected.csv", skip = 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  data$Diagnoses <- stri_trim_both(data$Diagnoses)
  data$Diagnoses[data$Diagnoses == "fracture femur"] <- "Fracture femur"
  data$Diagnoses[data$Diagnoses == "fp"] <- "Intra-abdominal sepsis"
  data$Diagnoses[data$Diagnoses == "pn"] <- "Community-acquired sepsis"
  data$Number <- stri_trim_both(data$Number)
  return(data)
}

#' Read the data set on provoked sepsis and controls in rat
#'
#' @return the data.frame with data, column names are identical to file content (non-standard for R)
#' @export
#'
#' @examples
get_rat_sepsis_data <- function(){
  data <- read.csv(file = "../../data/measurements/Summary rat sample data 61 edit.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, dec = ",", check.names = FALSE, blank.lines.skip = TRUE, strip.white = TRUE)
  data <- data[apply(data, 1, function(x){ sum(is.na(x)) }) < ncol(data) - 20, ] #Strip rows without any non-NA value
  data <- data[, apply(data, 2, function(x) { sum(is.na(x)) }) < nrow(data)] #Strip columns without any non-NA value
  data <- data[, apply(data, 2, function(x){ length(unique(x))}) > 1] #Strip columns with constant values
  data$BE <- data$BE - (min(data$BE, na.rm = TRUE) - 0.1) # without -0.1 later normalization will return -Inf for the lowest value
  data <- remove_bile_acids_from_measurements(data)
  colns <- colnames(data)
  colns[colns == "Tot Cholesterol"] <- "Total Cholesterol"
  colns[colns == "LDL/VLDL"] <- "LDL"
  colnames(data) <- colns
  grp <- data$group
  grp[grp == "control"] <- "Sham"
  grp[grp == "septic survivor"] <- "Septic-S"
  grp[grp == "septic non-survivor"] <- "Septic-NS"
  data$group <- grp
  return(data)
}

#' Read the plasma concentration range data on healthy french volunteers. Lipid concentrations are from FIA-MS/MS, hence not comparable to our data
#'
#' @return data.frame with mean, median
#' @export
#'
#' @examples
get_french_normal_data <- function(){
  data <- fread(input = "../../data/measurements/normal_plasma_concentrations/normal_all_metabs_combined.csv", sep = "\t", data.table = FALSE)
  data <- data[data$`Mean ± SD (µmol/L)` != "ND", ]
  iqr <- data$`Inter-quartile Range`
  iqr <- stri_sub(str = iqr, from = 2, to = stri_length(iqr) - 1)
  iqrv <- strsplit(x = iqr, split = ";", fixed = TRUE)
  iqrv <- sapply(iqrv, as.numeric)
  data$LQ <- iqrv[1, ]
  data$UQ <- iqrv[2, ]
  msd <- data$`Mean ± SD (µmol/L)`
  msd <- strsplit(x = msd, split = substr(msd[1], 6, 6))
  msd <- sapply(msd, as.numeric)
  data$Mean <- msd[1, ]
  data$SD <- msd[2, ]
  wol_met <- substr(data$Metabolite, 1, 2) == "L-"
  wol_met[data$Metabolite == "N-Acetylornithine"] <- TRUE
  data$Metabolite[wol_met] <- substring(data$Metabolite[wol_met], first = 3)
  data$Metabolite[data$Metabolite == "Glutamic acid"] <- "Glutamate"
  data$Metabolite[data$Metabolite == "Hexanoylcarnitine  (Fumarylcarnitine)"] <- "Hexanoylcarnitine (Fumarylcarnitine)"
  data$Median <- as.numeric(data$Median)
  data <- na.omit(data)
  return(data)
}

#' Read and return Copasi kinetic model parameter estimation results
#'
#' @return list of data.frames with parameter esimation results per day
#' @export
#'
#' @examples
get_human_kinetic_model_fit_results <- function(){
  #Read manually exported result files
  kin_fit_res_files <- list.files(path = "../../results/kinetic_model/", pattern = "stepwise_day._ES_..\\.txt", full.names = TRUE)
  kin_fit_res_file_con <- lapply(kin_fit_res_files, file, open = "r", blocking = TRUE) 
  raw_res_data <- lapply(kin_fit_res_file_con, readLines)
  dummy <- lapply(kin_fit_res_file_con, close)
  res_empty_lines <- lapply(lapply(raw_res_data, sapply, stri_cmp_eq, e2 = ""), which)
  res_data <- lapply(1:length(raw_res_data), function(rn) raw_res_data[[rn]][res_empty_lines[[rn]][1]:res_empty_lines[[rn]][2]])
  res_tables <- lapply(lapply(res_data, paste, collapse = "\n"), fread, data.table = FALSE)
  names(res_tables) <- stri_extract_last_regex(kin_fit_res_files, "day[0-9]")
  #Read Copasi's automatic result report files
  kin_fit_report_files <- list.files(path = "../../results/kinetic_model/", pattern = "stepwise_day._ES_.._report\\.txt", full.names = TRUE)
  kin_fit_report_file_con <- lapply(kin_fit_report_files, file, open = "r", blocking = TRUE) 
  raw_report_data <- lapply(kin_fit_report_file_con, readLines)
  dummy <- lapply(kin_fit_report_file_con, close)
  report_empty_lines <- lapply(lapply(raw_report_data, sapply, stri_cmp_eq, e2 = ""), which)
  report_data <- lapply(1:length(raw_report_data), function(rn) raw_report_data[[rn]][report_empty_lines[[rn]][11]:report_empty_lines[[rn]][12]])
  report_tables <- lapply(lapply(report_data, paste, collapse = "\n"), fread, data.table = FALSE)
  names(report_tables) <- stri_extract_last_regex(kin_fit_report_files, "day[0-9]")
  #Combine and return
  comb_results <- res_tables
  comb_results <- c(comb_results, report_tables[setdiff(names(report_tables), names(res_tables))])
  comb_results <- comb_results[order(names(comb_results))]
  comb_results <- lapply(comb_results, `[`, -1)
  return(comb_results)
}

#' Remove all bile acids from the metabolite legend table
#'
#' @param df the metabolite legend
#'
#' @return the legend without bile acids
#' @export
#'
#' @examples
remove_bile_acids_from_legend <- function(df){
  df[df[, 3] != "bile acid", ] #remove because of unsuitable measurement method
}

#' Remove all bile acids from the measurements table, select by colname
#'
#' @param df the measurement table
#'
#' @return the measurement table without bile acids
#' @export
#'
#' @examples
remove_bile_acids_from_measurements <- function(df){
  df[, -grep(pattern = ".*CA$", x = colnames(df))] #remove because of unsuitable measurement method
}

#' Read the data set on provoked sepsis and controls in rat
#'
#' @return the data.frame with data, column names are identical to file content (non-standard for R)
#' @export
#'
#' @examples
get_rat_sepsis_data_unedit <- function(){
  data <- read.csv(file = "../../data/measurements/Summary rat sample data.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, dec = ",", check.names = FALSE, blank.lines.skip = TRUE, strip.white = TRUE)
  data <- data[apply(data, 1, function(x){ sum(is.na(x)) }) < ncol(data) - 20, ] #Strip rows without any non-NA value
  data <- data[, apply(data, 2, function(x) { sum(is.na(x)) }) < nrow(data)] #Strip columns without any non-NA value
  data <- data[, apply(data, 2, function(x){ length(unique(x))}) > 1] #Strip columns with constant values
  data$BE <- data$BE - min(data$BE, na.rm = TRUE)
  data <- remove_bile_acids_from_measurements(data)
  colns <- colnames(data)
  colns[colns == "Tot Cholesterol"] <- "Total Cholesterol"
  colns[colns == "LDL/VLDL"] <- "LDL"
  colnames(data) <- colns
  return(data)
}

#' Read the legend to the clinical data set on sepsis and controls in rat. The legend has two group assignments of different coarseness.
#'
#' @return data.frame, first column stores the column names in the measurement data table, second the fine grained grouping, third the coarse grained grouping
#' @export
#'
#' @examples
get_rat_sepsis_legend <- function(){
  data <- read.csv(file = "../../data/measurements/Legend rat sample data.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, dec = ",", check.names = FALSE)
  data[[1]][data[[1]] == "Tot Cholesterol"] <- "Total Cholesterol"
  data[[1]][data[[1]] == "LDL/VLDL"] <- "LDL"
  data <- remove_bile_acids_from_legend(data)
  data <- data[!(data[, 1] %in% c("c4-OH-Pro", "Dopamine", "Nitro-Tyr", "PEA")), ] #excluded bc. of all NA or all 0
  return(data)
}

#' Get time points of significant differences
#'
#' @param diff_data data.frame as provided by rat_sig_diffs_along_time
#' @param alpha significance niveau, e.g. 0.05 (default)
#' @param time_var name of the time point column
#'
#' @return long format data.frame with time points at which a significant difference was found
#' @export
#'
#' @examples
get_sig_var_pos <- function(diff_data, alpha = 0.05, time_var = "Time"){
  #Get sig var pos
  sig_diff_pos <- lapply(diff_data[,-1], function(x, t_var, a){ diff_data[[t_var]][which(x <= a)] }, t_var = time_var, a = alpha)
  ##Clean
  sig_diff_pos <- sig_diff_pos[lapply(sig_diff_pos, length) > 0]
  ##Transform to long table
  for (n in seq_along(sig_diff_pos)){
    sig_diff_pos[[n]] <- data.frame(Time = sig_diff_pos[[n]], variable = names(sig_diff_pos)[n])
    colnames(sig_diff_pos[[n]])[1] = time_var
  }
  
  sig_diff_pos_long <- rbindlist(sig_diff_pos)
  return(sig_diff_pos_long)
}

#' Simulate measurements and count deviations for group of interest. Deviation means a value outside of the range defined by another group of samples.
#'
#' @param number dummy variable for lapply
#' @param n number of samples
#' @param d number of features per sample
#' @param dev_idx index to samples to count deviations for
#' @param sample_groups which samples belong together, e.g. to the same patient
#'
#' @return numeric vector of deviation counts per sample of interest
#' @export
#'
#' @examples
sim_dev <- function(number = 1, n, d, dev_idx, sample_groups){
  mat <- matrix(rlnorm(n = n * d), nrow = n, ncol = d)
  dev_max <- colMaxs(mat[-dev_idx, ])
  dev_min <- colMins(mat[-dev_idx, ])
  udev <- mat > matrix(dev_max, ncol = d, nrow = n, byrow = TRUE)
  ldev <- mat < matrix(dev_min, ncol = d, nrow = n, byrow = TRUE)
  sdev <- aggregate(udev | ldev, by = list(Sample = sample_groups), FUN = max) #count same metabolite at different time points as one deviation
  dev_score <- data.frame(Sample = sdev$Sample, score = rowSums(sdev[, -1]))
  res <- dev_score$score[match(unique(sample_groups[dev_idx]), dev_score$Sample)]
  names(res) <- unique(sample_groups[dev_idx])
  return(res)
}

#' Simulate measurements and count deviations in a group of interest for each column. Deviation means a value outside of the range defined by another group of samples.
#'
#' @param number dummy variable for lapply
#' @param n number of samples
#' @param d number of features per sample
#' @param dev_idx index to samples to count deviations for
#' @param sample_groups which samples belong together, e.g. to the same patient
#'
#' @return numeric vector of deviating sample group counts per column. Has length = d
#' @export
#'
#' @examples
sim_dev_met <- function(number = 1, n, d, dev_idx, sample_groups){
  mat <- matrix(rlnorm(n = n * d), nrow = n, ncol = d)
  dev_max <- colMaxs(mat[-dev_idx, ])
  dev_min <- colMins(mat[-dev_idx, ])
  udev <- mat > matrix(dev_max, ncol = d, nrow = n, byrow = TRUE)
  ldev <- mat < matrix(dev_min, ncol = d, nrow = n, byrow = TRUE)
  met_dev <- aggregate(udev | ldev, by = list(Sample = sample_groups), FUN = max) #count same metabolite at different time points as one deviation
  samples_per_met_dev <- sapply(met_dev[, -1], function(met) c(sum(met > 0), sum(met == 0)))
  return(samples_per_met_dev[1, ])
}

#' Simulate measurements and count the number of triples whose deviations cover all samples
#'
#' @param number dummy variable for lapply
#' @param n number of samples
#' @param d number of features per sample
#' @param dev_idx index to samples to count deviations for
#' @param sample_groups which samples belong together, e.g. to the same patient
#'
#' @return numeric vector of the indizes where you have triples
#' @export
#'
#' @examples
sim_dev_tuple <- function(number = 1, n, d, dev_idx, sample_groups){
  mat <- matrix(rlnorm(n = n * d), nrow = n, ncol = d)
  dev_max <- colMaxs(mat[-dev_idx, ])
  dev_min <- colMins(mat[-dev_idx, ])
  udev <- mat > matrix(dev_max, ncol = d, nrow = n, byrow = TRUE)
  ldev <- mat < matrix(dev_min, ncol = d, nrow = n, byrow = TRUE)
  met_dev <- aggregate(udev | ldev, by = list(Sample = sample_groups), FUN = max) #count same metabolite at different time points as one deviation
  met_dev <- met_dev[, c(1, 1 + which(colAnys(as.matrix(met_dev[, -1]))))]
  dev_idx_per_met <- lapply(lapply(met_dev[, -1], as.logical), which)
  dev_idx_met_crossover <- lapply(dev_idx_per_met, function(e) lapply(dev_idx_per_met, union, e))
  full_idx <- Reduce("rbind", lapply(lapply(dev_idx_met_crossover, lapply, length), sapply, `==`, length(unique(sample_groups[dev_idx]))))
  return((sum(full_idx) - sum(diag(full_idx))) / 2)
}

#' Scale values in a matrix or data.frame by dividing through their column-wise maximum. No min-max scaling, unless the minimum value in each column is 0!
#'
#' @param x matrix or data.frame
#' @param subset boolean or numerical index. For example 5:100, -1:-4 or c(F, F, F, F, T, ... , T). Default is all columns.
#'
#' @return scaled input matrix or data.frame such that the max value in each column is 1
#' @export
#'
#' @examples
max_norm <- function(x, subset = 1:ncol(x)){
  require(matrixStats)
  x_max_norm <- x
  x_max_norm[, subset] <- t(t(x_max_norm[, subset]) / colMaxs(as.matrix(x_max_norm[, subset]), na.rm = TRUE))
  return(x_max_norm)
}

#' Compute the [true/false]-[positive/negative]-rate given a prediction and reference from a binary classification task
#'
#' @param pred predicted class, coded in [-1, 1], but really any [<=0, >0] counts
#' @param ref true class, see pred for info
#'
#' @return a list with elements "tpr", "tnr", "fpr" and "fnr"
#' @export
#'
#' @examples
ml.npr <- function(pred, ref){
  res <- list()
  tp <- sum(pred > 0 & ref > 0)
  tn <- sum(pred <= 0 & ref <= 0)
  fp <- sum(pred > 0 & ref <= 0)
  fn <- sum(pred <= 0 & ref > 0)
  res[["tpr"]] <- sum(pred > 0 & ref > 0)/sum(ref > 0)
  res[["tnr"]] <- sum(pred <= 0 & ref <= 0)/sum(ref <= 0)
  res[["fpr"]] <- sum(pred > 0 & ref <= 0)/sum(ref <= 0)
  res[["fnr"]] <- sum(pred <= 0 & ref > 0)/sum(ref > 0)
  res[["ppv"]] <- tp / sum(pred > 0)
  res[["npv"]] <- tn / sum(pred <= 0)
  return(res)
}

#' Compute the accuracy measure for classifier predictions. Works for only for binary classes.
#'
#' @param pred predicted class
#' @param ref true class
#'
#' @return accuracy, single value
#' @export
#'
#' @examples
ml.acc <- function(pred, ref){
  return(sum(pred > 0 & ref > 0 | pred <= 0 & ref <= 0) / length(ref))
}

#' Calculate the points of the Receiver Operatoring Curve
#'
#' @param ref predicted class, coded in [-1, 1], but really any [<=0, >0] counts
#' @param conf classification confidence
#'
#' @return a matrix of column vectors for x and y coordinates
#' @export
#'
#' @examples
ml.roc <- function(ref, conf){
  conf_order <- order(conf, decreasing = TRUE)
  ref <- ref[conf_order]
  roc <- matrix(0, ncol = 2, nrow = length(conf))
  colnames(roc) <- c("FPR", "TPR")
  tp <- ref > 0
  fp <- ref <= 0
  tpr <- cumsum(tp) / sum(tp)
  fpr <- cumsum(fp) / sum(fp)
  roc[,1] <- fpr
  roc[,2] <- tpr
  return(roc)
}

#' Calculate the Area Under the receiver-operator-Curve
#'
#' @param ref predicted class, coded in [-1, 1], but really any [<=0, >0] counts
#' @param conf classification confidence
#'
#' @return the AUC value
#' @export
#'
#' @examples
ml.auc <- function(ref, conf){
  roc <- ml.roc(ref, conf)
  fpr_diff <- diff(c(roc[,1], 1))
  tpr_diff <- diff(c(roc[,2], 1))
  auc <- sum(fpr_diff * roc[,2] + 0.5 * fpr_diff * tpr_diff)
  return(auc)
}

#' Compute the Heaviside function H(x1, x2). Works on scalars and vectors
#'
#' @param x1 number or numeric vector
#' @param x2 number or numeric vector of same length as x1
#'
#' @return numeric vector with length of x1 (or x2 if x1 is a number) of Heaviside function values
#' @export
#'
#' @examples
heaviside <- function(x1, x2){
  return(as.numeric((x1 < x2) + 0.5 * (x1 == x2)))
}

#' Perform (nested) Tournament Leave Pair Out Cross Validation (TLPOCV) with Recursive Feature Elimination (RFE) in parralel via parallel::mclapply. Requires at least 5 samples from each class. Default model is from Random Forest package 'ranger'.
#' TLPOCV is implemented after Perez et al., 2018:
#' @Article{Perez2018,
#'   author    = {Ileana Montoya Perez and Antti Airola and Peter J Boström and Ivan Jambor and Tapio Pahikkala},
#'   title     = {Tournament leave-pair-out cross-validation for receiver operating characteristic analysis},
#'   journal   = {Statistical Methods in Medical Research},
#'   year      = {2018},
#'   pages     = {096228021879519},
#'   doi       = {10.1177/0962280218795190},
#'   publisher = {{SAGE} Publications},
#' }
#'
#' @param data_x data.frame of predictors, samples in rows and features in columns
#' @param data_y data.frame, univariate, classes coded as [0, 1]
#' @param train_fun function with arguments tr_x (predictors) and tr_y (classes), returns a model trained on tr_x and tr_y
#' @param prob_fun function with arguments classifier (model returned by train_fun) and te_x (test set predictors), returns the class probabilities or class scores of samples for the first class as a numeric vector
#' @param varimp_fun function with argument classifier (model returned by train_fun), returns the variable importance in classifier as a numeric vector
#' @param mc.cores number of cores to use in parallel processing of pairs, cf. parallel::mclapply
#' @param set_sizes numerical vector of feature set sizes to try, default is all possible set sizes
#'
#' @return list of lists: AUC values of the inner LPOCV iterations, AUC values of the outer TLPOCV iterations, selected best features for each feature set size with the order corresponding to that in the outer AUC list, sample ranks combined with classes for each feature set size corresponding to validation, sample ranks combined with classes for each feature set size corresponding to testing
#' @export
#'
#' @examples
tlpocv_rfe_parallel <- function(data_x, data_y, train_fun, prob_fun, varimp_fun, set_sizes = ncol(data_x):2, mc.cores = 7){
  if (missing(train_fun) | missing(prob_fun) | missing(varimp_fun)){
    warning("train_fun, prob_fun or varimp_fun is missing, using ranger")
    library(ranger)
    train_fun <- function(tr_x, tr_y){
      ranger(data = cbind(tr_y, tr_x), 
             dependent.variable.name = colnames(tr_y), 
             write.forest = T, 
             probability = T, 
             save.memory = F, 
             importance = "impurity", 
             num.threads = 1)
    }
    prob_fun <- function(classifier, te_x){
      predict(classifier, te_x)$predictions[, 1]
    }
    varimp_fun <- function(classifier){
      classifier$variable.importance
    }
  }
  if (max(set_sizes) >= ncol(data_x)){
    set_sizes <- set_sizes[set_sizes < ncol(data_x)]
  }
  set_sizes <- c(1, set_sizes) #include 1 as "buffer" in RFE loop even if already in set_sizes
  set_sizes <- sort(set_sizes, decreasing = TRUE)
  library(parallel)
  ##Build sample pair list
  tlpo_s_df <- t(combn(x = 1:nrow(data_x), m = 2))
  ##Run Feature Selection
  outer_AUC <- rep(0, length(ncol(data_x)))
  ###Function for internal TLPOCV with RFE
  int_tlpocv_rfe <- function(p_ext, tlpo_s_df, data_x, data_y, set_sizes){
    tlpo_int_s_df <- tlpo_s_df[(!tlpo_s_df[, 1] %in% tlpo_s_df[p_ext, ]) & !(tlpo_s_df[, 2] %in% tlpo_s_df[p_ext, ]), ] #exclude samples in validation pair from training
    feature_list <- list()
    feature_list[[1]] <- colnames(data_x) #start with all features
    pair_dir <- list()
    inner_AUC <- rep(-100, length(set_sizes) - 1)
    int_sample_ranking <- list()
    for (feat_count in seq_along(set_sizes)){
      feat_sel <- feature_list[[feat_count]] #update feature selection
      pair_dir[[feat_count]] <- rep(-100, nrow(tlpo_int_s_df)) #edge direction in dominance graph
      feat_imp <- matrix(0, nrow = nrow(tlpo_int_s_df), ncol = length(feat_sel)) #feature importance matrix for this run
      for (p_int in 1:nrow(tlpo_int_s_df)){ #for each test pair
        tr_x <- data.frame(data_x[-c(tlpo_int_s_df[p_int, ], tlpo_s_df[p_ext, ]), feat_sel])
        tr_y <- data.frame(data_y[-c(tlpo_int_s_df[p_int, ], tlpo_s_df[p_ext, ]), 1])
        te_x <- data.frame(data_x[tlpo_int_s_df[p_int, ], feat_sel])
        te_y <- data.frame(data_y[tlpo_int_s_df[p_int, ], 1])
        colnames(tr_x) <- feat_sel
        colnames(te_x) <- feat_sel
        colnames(tr_y) <- colnames(data_y)
        colnames(te_y) <- colnames(data_y)
        rg_te <- train_fun(tr_x = tr_x, tr_y = tr_y)
        pred_te <- prob_fun(classifier = rg_te, te_x = te_x)
        pair_dir[[feat_count]][p_int] <- heaviside(pred_te[2], pred_te[1]) #determine edge direction
        feat_imp[p_int, ] <- varimp_fun(classifier = rg_te) #collect variable importance
      }
      ####Calculate AUC as in Perez et al., 2018
      #####Get out degree for each vertex (sample)
      out_degree <- rep(0, nrow(data_y))
      d1 <- aggregate(x = pair_dir[[feat_count]], by = list(v = tlpo_int_s_df[, 1]), FUN = sum) 
      d2 <- aggregate(x = 1 - pair_dir[[feat_count]], by = list(v = tlpo_int_s_df[, 2]), FUN = sum)
      out_degree[d1$v] <- d1$x
      out_degree[d2$v] <- out_degree[d2$v] + d2$x
      tlpo_od <- apply(tlpo_int_s_df[, 1:2], 1:2, function(x) out_degree[x])
      tlpo_h <- heaviside(tlpo_od[, 2], tlpo_od[, 1])
      int_sample_ranking[[feat_count]] <- cbind(data_y, out_degree)[-tlpo_s_df[p_ext, ],]
      #####use only the positive-negative pairs for AUC (positive first sample, negative second sample!)
      pos1_idx <- data_y[tlpo_int_s_df[, 1], 1] == 1
      pos2_idx <- data_y[tlpo_int_s_df[, 2], 1] == 1
      neg1_idx <- data_y[tlpo_int_s_df[, 1], 1] == 0
      neg2_idx <- data_y[tlpo_int_s_df[, 2], 1] == 0
      sum_h <- sum(tlpo_h[pos1_idx & neg2_idx]) + sum(1 - tlpo_h[pos2_idx & neg1_idx])
      inner_AUC[feat_count] <- sum_h / (sum(pos1_idx & neg2_idx) + sum(pos2_idx & neg1_idx))
      feat_imp <- colMeans(feat_imp) #rank features by importance
      feature_list[[feat_count + 1]] <- feat_sel[order(feat_imp, decreasing = TRUE)[1:set_sizes[feat_count]]] #down-select feature set
    }
    feature_list[[length(feature_list)]] <- NULL
    return(list(inner_AUC = inner_AUC, feature_list = feature_list, int_sample_ranking = int_sample_ranking))
  }
  ###Internal TLPOCV
  int_tlpocv_res <- mclapply(1:nrow(tlpo_s_df), int_tlpocv_rfe, tlpo_s_df = tlpo_s_df, data_x = data_x, data_y = data_y, set_sizes = set_sizes, mc.cores = mc.cores)
  ###Extract fields
  inner_AUC <- lapply(int_tlpocv_res, `[[`, "inner_AUC")
  feature_list <- lapply(int_tlpocv_res, `[[`, "feature_list")
  int_sample_ranking <- lapply(int_tlpocv_res, `[[`, "int_sample_ranking")
  ###Find best feature set/most common feature combination for each feature set size
  best_feat_set <- list()
  for (feat_count in seq_along(feature_list[[1]])){
    fls <- unlist(lapply(feature_list, `[[`, feat_count)) #gets all feature sets for a fixed set size from all validation pairs
    fl_hist <- table(fls)
    fl_hist <- sort(fl_hist, decreasing = TRUE) #beware! feature names are not in original order
    best_feat_set[[feat_count]] <- names(fl_hist)[1:(length(feature_list[[1]][[feat_count]]))]
  }
  ###Function for external TLPOCV
  ext_tlpocv_eval <- function(p_ext, tlpo_s_df, data_x, data_y, best_feat_set){
    pair_ext_dir <- rep(0, length(best_feat_set))
    tlpo_int_s_df <- tlpo_s_df[(!tlpo_s_df[, 1] %in% tlpo_s_df[p_ext, ]) & !(tlpo_s_df[, 2] %in% tlpo_s_df[p_ext, ]), ] #exclude samples in validation pair from training
    for (feat_count in seq_along(best_feat_set)){
      feat_sel <- best_feat_set[[feat_count]]
      tr_x <- data.frame(data_x[-tlpo_s_df[p_ext, ], feat_sel])
      tr_y <- data.frame(data_y[-tlpo_s_df[p_ext, ], 1])
      va_x <- data.frame(data_x[tlpo_s_df[p_ext, ], feat_sel])
      va_y <- data.frame(data_y[tlpo_s_df[p_ext, ], 1])
      colnames(tr_x) <- feat_sel
      colnames(va_x) <- feat_sel
      colnames(tr_y) <- colnames(data_y)
      colnames(va_y) <- colnames(data_y)
      rg_va <- train_fun(tr_x = tr_x, tr_y = tr_y)
      pred_va <- prob_fun(classifier = rg_va, te_x = va_x)
      pair_ext_dir[feat_count] <- heaviside(pred_va[2], pred_va[1]) #determine edge direction
    }
    return(list(pair_ext_dir = pair_ext_dir))
  }
  ###Evaluate feature sets on external validation pairs
  ext_tlpocv_res <- mclapply(1:nrow(tlpo_s_df), ext_tlpocv_eval, tlpo_s_df = tlpo_s_df, data_x = data_x, data_y = data_y, best_feat_set = best_feat_set, mc.cores = mc.cores)
  pair_ext_dir <- lapply(ext_tlpocv_res, `[[`, "pair_ext_dir")
  ####Calculate AUC as in Perez et al., 2018
  #####Get out degree for each vertex (sample)
  ext_sample_ranking <- list()
  for (feat_count in seq_along(pair_ext_dir[[1]])){
    out_degree <- rep(0, nrow(data_y))
    d1 <- aggregate(x = sapply(pair_ext_dir, `[[`, feat_count), by = list(v = tlpo_s_df[, 1]), FUN = sum) 
    d2 <- aggregate(x = 1 - sapply(pair_ext_dir, `[[`, feat_count), by = list(v = tlpo_s_df[, 2]), FUN = sum)
    out_degree[d1$v] <- d1$x
    out_degree[d2$v] <- out_degree[d2$v] + d2$x
    tlpo_od <- apply(tlpo_s_df[, 1:2], 1:2, function(x) out_degree[x])
    tlpo_h <- heaviside(tlpo_od[, 2], tlpo_od[, 1])
    ext_sample_ranking[[feat_count]] <- out_degree
    #####use only the positive-negative pairs for AUC (positive first sample, negative second sample!)
    pos1_idx <- data_y[tlpo_s_df[, 1], 1] == 1
    pos2_idx <- data_y[tlpo_s_df[, 2], 1] == 1
    neg1_idx <- data_y[tlpo_s_df[, 1], 1] == 0
    neg2_idx <- data_y[tlpo_s_df[, 2], 1] == 0
    sum_h <- sum(tlpo_h[pos1_idx & neg2_idx]) + sum(1 - tlpo_h[pos2_idx & neg1_idx])
    outer_AUC[feat_count] <- sum_h / (sum(pos1_idx & neg2_idx) + sum(pos2_idx & neg1_idx))
  }
  ext_sample_ranking <- lapply(ext_sample_ranking, function(e) return(cbind(data_y, e)))
  return(list(inner_AUC = inner_AUC, outer_AUC = outer_AUC, best_features = best_feat_set, ext_sample_ranking = ext_sample_ranking, int_sample_ranking = int_sample_ranking))
}

#' Plot the training and validation AUC from nested TLPOCV-RFE along feature set sizes. If a feature set is of size 2 plot training and validation ROC curves and a dot graph of both features.
#'
#' @param tlpocv_rfe_res list output from tlpocv_rfe_parallel()
#' @param file_prefix prefix to append to plot files
#' @param plot_dir directory to plot into
#'
#' @return built ggplot plots ready for saving to file
#' @export
#'
#' @examples
plot_tlpocv_rfe_res <- function(tlpocv_rfe_res, file_prefix, plot_dir){
  library(ggplot2)
  res <- list()
  AUC_data <- data.frame(Reduce("rbind", tlpocv_rfe_res$inner_AUC))
  AUC_data <- rbind(AUC_data, data.frame(t(tlpocv_rfe_res$outer_AUC)))
  colnames(AUC_data) <- sapply(tlpocv_rfe_res$best_features, length)
  AUC_data$stage <- c(rep("Test data", length(tlpocv_rfe_res$inner_AUC)), "Validation data")
  AUC_data_long <- melt(AUC_data, id.vars = c("stage"))
  AUC_data_long$variable <- as.numeric(as.character(AUC_data_long$variable))
  AUC_data_long$value <- 1 - AUC_data_long$value
  p_AUC_RFE <- ggplot(data = AUC_data_long, mapping = aes(x = variable, y = value, color = stage, fill = stage)) + 
    stat_summary(fun.y = median, fun.ymax = function(x) quantile(x, p = 0.75), fun.ymin = function(x) quantile(x, p = 0.25), geom = "ribbon", alpha = 0.5, colour = NA) + 
    stat_summary(fun.y = median, geom = "line") + 
    ylim(0, 1) + 
    ylab("Median AUC") +
    xlab("Nubmer of features") +
    ggtitle("Test and Validation AUC of nested TLPOCV-RFE") +
    theme_bw()
  ggsave(filename = paste0(file_prefix, "TLPOCV_RFE_AUC.png"), path = plot_dir, plot = p_AUC_RFE, width = 10, height = 5, units = "in")
  
  res[[1]] <- ggplot_build(p_AUC_RFE)
  
  set_sizes <- sapply(tlpocv_rfe_res$best_features, length)
  if (2 %in% set_sizes){
    f2_pos <- which(sapply(tlpocv_rfe_res$best_features, length) == 2)
    int_ROC_data <- lapply(lapply(tlpocv_rfe_res$int_sample_ranking, `[[`, f2_pos), function(e) ml.roc(ref = 1 - e[, 1], conf = e[, 2]))
    int_ROC_data <- lapply(int_ROC_data, function(e) stepfun(e[, 1], c(0, e[, 2])))
    int_ROC_data <- lapply(int_ROC_data, function(f) data.frame(FPR = seq(0, 1, by = 0.01), TPR = f(seq(0, 1, by = 0.01))))
    int_ROC_data <- data.frame(Reduce("rbind", int_ROC_data))
    if (min(int_ROC_data[, 1]) != 0)
      int_ROC_data <- rbind(c(0, 0), int_ROC_data)
    if (max(int_ROC_data[, 1]) != 1)
      int_ROC_data <- rbind(int_ROC_data, c(1, 1))
    colnames(int_ROC_data) <- c("FPR", "TPR")
    int_ROC_data$stage <- "Test data"
    ranks_2feat <- tlpocv_rfe_res$ext_sample_ranking[[f2_pos]]
    ext_ROC_data <- data.frame(ml.roc(ref = 1 - ranks_2feat[, 1], conf = ranks_2feat[, 2]))
    if (min(ext_ROC_data[, 2]) != 0)
      ext_ROC_data <- rbind(c(0, 0), ext_ROC_data)
    if (max(ext_ROC_data[, 1]) != 1)
      ext_ROC_data <- rbind(ext_ROC_data, c(1, 1))
    colnames(ext_ROC_data) <- c("FPR", "TPR")
    ext_ROC_data$stage <- "Validation data"
    ROC_data <- rbind(int_ROC_data, ext_ROC_data)
    p <- ggplot(data = ROC_data, mapping = aes(x = FPR, y = TPR, color = stage, fill = stage)) + 
      stat_summary(data = int_ROC_data, 
                   fun.y = mean, fun.ymax = function(x) quantile(x, p = 0.75), fun.ymin = function(x) quantile(x, p = 0.25), 
                   geom = "ribbon", alpha = 0.5, colour = NA) + 
      stat_summary(data = int_ROC_data, 
                   fun.y = median, geom = "line") +
      geom_line(data = ext_ROC_data) +
      geom_segment(x = 0, y = 0, xend = 1, yend = 1, size = 0.25, inherit.aes = FALSE) +
      geom_text(x = 0.95, y = 0.15, hjust = "right", color = scales::hue_pal()(2)[1],
                label = paste0("Training AUC = ", format(mean(subset(AUC_data_long, variable == 2 & stage == "Test data", "value")[[1]]), digits = 3))) +
      geom_text(x = 0.95, y = 0.10, hjust = "right", color = scales::hue_pal()(2)[2],
                label = paste0("Validation AUC = ", format(subset(AUC_data_long, variable == 2 & stage == "Validation data", "value"), digits = 3))) +
      scale_x_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
      scale_y_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
      ggtitle("Test and Validation ROC of nested TLPOCV RFE\nbest 2 feature set") +
      theme_bw() + 
      theme(panel.grid = element_blank())
    ggsave(filename = paste0(file_prefix, "TLPOCV_RFE_ROC_2feat.png"), path = plot_dir, plot = p, width = 6, height = 5, units = "in")
    
    res[[2]] <- ggplot_build(p)
    
    f2 <- tlpocv_rfe_res$best_features[[which(sapply(tlpocv_rfe_res$best_features, length) == 2)]]
    p <- ggplot(data = subset(human_sepsis_data_ml, Day == 0, c("Survival", f2)), mapping = aes_string(x = f2[1], y = f2[2], color = quote(factor(Survival)))) + 
      geom_point() + 
      theme_bw()
    ggsave(plot = p, filename = paste0(file_prefix, "class_separation_2feat.png"), width = 5, height = 4, units = "in", path = plot_dir)
    
    res[[3]] <- ggplot_build(p)
  }
  return(res)
}

#' Generate folds for stratified cross-validation in a binary classification scenario
#'
#' @param num_folds number of folds
#' @param class class membership, numerical coding in [<=0, >0]
#'
#' @return vectors of indices combined in a list
#' @export
#'
#' @examples
ml.split.folds.strat <- function(num_folds, class){
  pos_idx <- class > 0
  neg_idx <- class <= 0
  pos_sample <- sample(x = which(pos_idx), size = sum(pos_idx), replace = FALSE)
  neg_sample <- sample(x = which(neg_idx), size = sum(neg_idx), replace = FALSE)
  pos_sample_from <- ceiling(length(pos_sample)*(1:num_folds)/num_folds)
  pos_sample_to <- ceiling(length(pos_sample)*(0:(num_folds-1))/num_folds)+1
  neg_sample_from <- ceiling(length(neg_sample)*(1:num_folds)/num_folds)
  neg_sample_to <- ceiling(length(neg_sample)*(0:(num_folds-1))/num_folds)+1
  res <- list()
  for (fold in 1:num_folds){
    res[[fold]] <- c(pos_sample[pos_sample_from[fold]:pos_sample_to[fold]],
                     neg_sample[neg_sample_from[fold]:neg_sample_to[fold]])
  }
  return(res)
}

#' Generate folds for stratified cross-validation in a binary classification scenario, but fully seperate the second class set into learn and test set.
#' Randomly samples partitions into folds to find one that has approximately equal number of cases and approx. equal class proportions (for stratified cross-validation).
#' Currently only supports binary stratified class and the optimization following the random sampling is not implemented.
#'
#' @param num_folds number of folds
#' @param class class membership, numerical coding in [<=0, >0]
#' @param non_strat_class class membership, for full seperation of classes
#' @param num_opt_steps number of optimisation steps following random sampling of partitions
#' @param num_sampling_steps number of random samples to generate
#'
#' @return
#' @export
#'
#' @examples
ml.split.folds.quasistrat <- function(num_folds, class, non_strat_class, num_opt_steps = 500, num_sampling_steps = 5000){
  #Get fields of stratified class
  class_set <- unique(class)
  num_non_strat_class1 <- length(table(non_strat_class[class == class_set[1]]))
  num_non_strat_class2 <- length(table(non_strat_class[class == class_set[2]]))
  #Get data stats
  class1_tab <- table(non_strat_class[class == class_set[1]])
  class2_tab <- table(non_strat_class[class == class_set[2]])
  #Find ideal partition parameters
  class1_part_count <- sum(class1_tab) / num_folds
  class2_part_count <- sum(class2_tab) / num_folds
  #Set random partition
  class1_assig <- sample.int(n = num_folds, size = num_non_strat_class1, replace = TRUE)
  class2_assig <- sample.int(n = num_folds, size = num_non_strat_class2, replace = TRUE)
  #Get goodness impression
  class1_num_cases <- sapply(1:num_folds, function(x){ sum(class1_tab[class1_assig == x]) })
  class2_num_cases <- sapply(1:num_folds, function(x){ sum(class2_tab[class2_assig == x]) })
  best_class1_assig_score <- sum((class1_part_count - class1_num_cases)^2)
  best_class2_assig_score <- sum((class2_part_count - class2_num_cases)^2)
  #Prepare sampling
  best_class1_assig <- class1_assig
  best_class2_assig <- class2_assig
  #Run random sampling to improve
  for (n in 1:num_sampling_steps){
    #Set new random partition
    class1_assig <- sample.int(n = num_folds, size = num_non_strat_class1, replace = TRUE)
    class2_assig <- sample.int(n = num_folds, size = num_non_strat_class2, replace = TRUE)
    #Get goodness impression
    class1_num_cases <- sapply(1:num_folds, function(x){ sum(class1_tab[class1_assig == x]) })
    class2_num_cases <- sapply(1:num_folds, function(x){ sum(class2_tab[class2_assig == x]) })
    class1_assig_score <- sum((class1_part_count - class1_num_cases)^2)
    class2_assig_score <- sum((class2_part_count - class2_num_cases)^2)
    if (class1_assig_score < best_class1_assig_score){
      best_class1_assig_score <- class1_assig_score
      best_class1_assig <- class1_assig
    }
    if (class2_assig_score < best_class2_assig_score){
      best_class2_assig_score <- class2_assig_score
      best_class2_assig <- class2_assig
    }
  }
  #Combine class1 and class2 partitions into folds
  ##Get members as index list
  class1_indizes <- lapply(1:num_folds, function(x){ which(non_strat_class %in% names(class1_tab)[which(best_class1_assig == x)]) })
  class2_indizes <- lapply(1:num_folds, function(x){ which(non_strat_class %in% names(class2_tab)[which(best_class2_assig == x)]) })
  ##Sort first class folds ascending
  class1_num_cases <- sapply(1:num_folds, function(x){ sum(class1_tab[best_class1_assig == x]) })
  class2_num_cases <- sapply(1:num_folds, function(x){ sum(class2_tab[best_class2_assig == x]) })
  class1_indizes <- class1_indizes[order(class1_num_cases, decreasing = TRUE)]
  class2_indizes <- class2_indizes[order(class2_num_cases, decreasing = FALSE)]
  ##Combine
  res <- lapply(1:num_folds, function(x){ c(class1_indizes[[x]], class2_indizes[[x]]) })
  ##Feedback about partition goodness
  #print(c(best_class1_assig_score, best_class2_assig_score, class1_num_cases, class2_num_cases))
  ##Return
  return(res)
}

#' Report p values for differences between survivors and non-survivors for each metabolite at each day seperately.
#' Only days with at least 3 samples are analysed. For testing we employ the U-test (wilcoxon rank sum test) and the two-sample t-test.
#'
#' @param data The human metabolome time series data set
#' @param corr_fdr Correct the p-values for FDR (default: TRUE)
#' @param time_var name of the time point column
#' @param status_var name of the status column
#' @param case_var name of the case column
#' @param descr_till_col column number up to which the case is described (time, status, case, ...)
#'
#' @return a list with members day_sig_u_diff and day_sig_t_diff. These are data.frames of p-values from U-test and t-test respectively.
#' @export
#'
#' @examples
human_sig_diffs_along_days <- function(data, corr_fdr = TRUE, time_var = "Day", status_var = "Survival", case_var = "Patient", descr_till_col = 6){
  day_survival_tab <- table(data[c(time_var, status_var)])
  ##Keep days with large enough sample count
  day_survival_tab <- day_survival_tab[rowMins(day_survival_tab) >= 3,]
  ##Get significant differences for each day
  day_sig_u_diff <- data.frame(Day = as.numeric(rownames(day_survival_tab)), matrix(0, nrow = nrow(day_survival_tab), ncol = ncol(data)-descr_till_col))
  colnames(day_sig_u_diff) <- c(time_var, colnames(data)[-1:-descr_till_col])
  day_sig_t_diff <- day_sig_u_diff
  non_short_stay_pats <- unique(subset(data, data[[time_var]] >= 0, case_var))[[1]]
  for (d in seq_along(day_sig_u_diff[[time_var]])){
    #Reduce to significant differences
    day <- day_sig_u_diff[[time_var]][d]
    for (n in seq_along(day_sig_u_diff[1,-1])){
      g1 <- subset(data, subset = data[[status_var]] == "S" & data[[time_var]] == day & data[[case_var]] %in% non_short_stay_pats, select = n + descr_till_col)
      g2 <- subset(data, subset = data[[status_var]] == "NS" & data[[time_var]] == day & data[[case_var]] %in% non_short_stay_pats, select = n + descr_till_col)
      g1 <- na.omit(g1)[[1]]
      g2 <- na.omit(g2)[[1]]
      if (length(g1) > 1 && length(g2) > 1 && length(table(g1)) > 1 && length(table(g2)) > 1){
        u_res <- wilcox.test(x = g1, y = g2)
        day_sig_u_diff[d, n + 1] <- u_res$p.value
        t_res <- t.test(x = g1, y = g2, var.equal = FALSE)
        day_sig_t_diff[d, n + 1] <- t_res$p.value
      }
      else{
        day_sig_u_diff[d, n + 1] <- NA
        day_sig_t_diff[d, n + 1] <- NA
      }
    }
  }
  if (corr_fdr){
    #correct the p-vals for FDR
    #apply() may be faster for this, but it rearranges the result; a reshape with matrix() would be necessary, which is badly readable
    for (d in 1:nrow(day_sig_u_diff)){
      day_sig_u_diff[d, -1] <- p.adjust(p = day_sig_u_diff[d,-1], method = "fdr")
      day_sig_t_diff[d, -1] <- p.adjust(p = day_sig_t_diff[d,-1], method = "fdr")
    }
  }
  res <- list(day_sig_u_diff = day_sig_u_diff, day_sig_t_diff = day_sig_t_diff)
}


rat_sig_diffs_along_time <- function(data, corr_fdr = TRUE){
  time_survival_tab <- table(data[c("time point", "group")])
  ##Keep times with large enough sample count
  time_survival_tab <- time_survival_tab[rowMins(time_survival_tab) >= 2,]
  ##Get significant differences for each time
  time_sig_u_diff <- data.frame(Time = rownames(time_survival_tab), matrix(1, nrow = nrow(time_survival_tab), ncol = ncol(data)-4))
  colnames(time_sig_u_diff)[-1] <- colnames(data)[-1:-4]
  time_sig_t_diff <- time_sig_u_diff
  time_mag_g1_t_diff <- time_sig_u_diff
  time_mag_g2_t_diff <- time_sig_u_diff
  time_fold_change <- time_sig_u_diff
  groups <- unique(data$group)
  for (d in seq_along(time_sig_u_diff$Time)){
    #Reduce to significant differences
    time <- time_sig_u_diff$Time[d]
    for (n in seq_along(time_sig_u_diff[1,-1])){
      g1 <- subset(data, subset = group == groups[1] & data$`time point` == time, select = n + 4)
      g2 <- subset(data, subset = group == groups[2] & data$`time point` == time, select = n + 4)
      g1 <- na.omit(g1)[[1]]
      g2 <- na.omit(g2)[[1]]
      if (length(g1) > 1 && length(g2) > 1 && length(table(g1)) > 1 && length(table(g2)) > 1){
        u_res <- wilcox.test(x = g1, y = g2)
        time_sig_u_diff[d, n + 1] <- u_res$p.value
        t_res <- t.test(x = g1, y = g2, var.equal = F)
        time_sig_t_diff[d, n + 1] <- t_res$p.value
      }
      else{
        time_sig_u_diff[d, n + 1] <- NA
        time_sig_t_diff[d, n + 1] <- NA
      }
      time_mag_g1_t_diff[d, n + 1] <- mean(g1)
      time_mag_g2_t_diff[d, n + 1] <- mean(g2)
      time_fold_change[d, n + 1] <- log(mean(g1) / mean(g2))
    }
  }
  if (corr_fdr){
    #correct the p-vals for FDR
    #apply() may be faster for this, but it rearranges the result; a reshape with matrix() would be necessary, which is badly readable
    for (d in 1:nrow(time_sig_u_diff)){
      time_sig_u_diff[d, -1] <- p.adjust(p = time_sig_u_diff[d,-1], method = "fdr")
      time_sig_t_diff[d, -1] <- p.adjust(p = time_sig_t_diff[d,-1], method = "fdr")
    }
  }
  res <- list(time_sig_u_diff = time_sig_u_diff, time_sig_t_diff = time_sig_t_diff, time_mag_g1_t_diff = time_mag_g1_t_diff, time_mag_g2_t_diff = time_mag_g2_t_diff, time_fold_change = time_fold_change)
}

#' Import Anna's significantly differentially concentrated metabolite names (acronyms like "C10" or "SM OH C16:0")
#'
#' @return chr vector of metabolite names
#' @export
#'
#' @examples
get_annas_rat_sig_diffs <- function(){
  mets = fread(input = "../../data/measurements/annas_sig_metabs.csv", sep = "\t", header = TRUE, data.table = FALSE)
  mets = subset(mets, isPresent == 1, "Variable")
  return(mets[[1]])
}

#' Fit a linear model with nlme::lme() and apply car::Anova().
#'
#' @param data data, cases in rows
#' @param formula formula connecting the variablees
#' @param dvar name of dependent variable in data
#' @param id.vars names of independent and grouping/misc. variables in data
#' @param random random effect specification, e.g. ~1|subject
#' @param use.corAR use corAR1()?
#'
#' @return function call to fit linear model with lme() from nlme package
#' @export
#'
#' @examples
fit_lin_mod_lme <- function(dvar, data, formula, id.vars, random, use.corAR = FALSE, control = lmeControl()){
  test_data <- data[, c(id.vars, dvar)]
  colnames(test_data)[ncol(test_data)] <- as.character(formula[[2]])
  for (s in id.vars){
    test_data[[s]] <- factor(test_data[[s]]) #to re-level
  }
  res <- NULL
  #fit linear model
  if (use.corAR){
    model.pre <- lme(fixed = formula, random = random, data = test_data)
    ac_val <- ACF(model.pre)[2,2] #autocorrelation value for a 1 day-lag
    ac_val <- max(min(ac_val, 0.99), -0.99) #value is in some cases >= 1
    try(res <- bquote(lme(fixed = .(formula), random = .(random), correlation = corAR1(form = .(random), value = .(ac_val)), data = .(test_data), method = "REML", keep.data = TRUE, control = .(control))))
    #try(res <- lmer(formula = formula, data = test_data, REML = TRUE))
  }
  else{
    try(res <- bquote(lme(fixed = .(formula), random = .(random), data = .(test_data), method = "REML", keep.data = TRUE, control = .(control))))
  }
  return(res)
}

#' Return model matrix for linear models fitted with lme() from the nlme package when the original data is not in the current environment
#'
#' @param object linear model build with lme()
#' @param data data used to build the linear model. Requires keep.data = TRUE in call to lme() if not supplied here.
#'
#' @return model matrix for lme objects
#' @export
#'
#' @examples
model.matrix.lme <- function(object, data){
  if (missing(data))
    data <- object$data
  model.matrix(object = object$terms, data = data, contrast.arg = object$contrasts)
}

#' Perform type III ANOVA, optionally with repeated measures, on each column of a data.frame indicated by col.set
#'
#' @param data data.frame
#' @param random random effect specification
#' @param formula data formula object
#' @param use.corAR use autocorrelated residual structure
#' @param col.set character vector of column names to iterate over
#' @param id.vars id variables in data.frame to use in random effect and formula
#'
#' @return list with fitted linear models, shapiro test p value and ANOVA p values
#' @export
#'
#' @examples
t3ANOVA <- function(data, random, formula, use.corAR, col.set, id.vars, control = lmeControl()){
  library(parallel)
  lin.models <- mclapply(X = col.set, FUN = fit_lin_mod_lme, data = data, formula = formula, id.vars = id.vars, random = random, use.corAR = use.corAR, control = control)
  lin.models <- mclapply(lin.models, function(m) try(eval(m)))
  lin.models <- lin.models[sapply(lin.models, function(m) class(m) == "lme")]
  model.normality.p <- unlist(lapply(lin.models, function(e) shapiro.test(residuals(e))$p.value))
  anova.terms <- rownames(Anova(lin.models[[1]], type = 3))
  anova.ps <- sapply(lin.models, function(x){ Anova(x, type = 3, white.adjust = TRUE)[[3]] }, USE.NAMES = TRUE)
  rownames(anova.ps) <- anova.terms
  res <- list(models = lin.models, normality.p = model.normality.p, ps = anova.ps)
  return(res)
}


#' Generate contrast matrix specific for Sepsis vs Control and Sepsis S vs Sepsis NS
#'
#' @param num_seq measurement day vector with as many entries as factor levels for Day effect in linear model
#'
#' @return
#' @export
#'
#' @examples
get_S_NS_C_contrmat <- function(num_seq){
  ld <- length(num_seq)
  sd <- seq_along(num_seq)
  contr.m <- matrix(0, ncol = ld * 2, nrow = ld * 3)
  colnames(contr.m) <- c(paste0("Seps vs Comp, D", num_seq), paste0("S vs NS, D", num_seq))
  #####Sepsis vs Nonsepsis
  contr.m[sd, sd] <- diag(-2, nrow = ld, ncol = ld)
  contr.m[sd + ld, sd] <- diag(1, nrow = ld, ncol = ld)
  contr.m[sd + 2 * ld, sd] <- diag(1, nrow = ld, ncol = ld)
  #####Septic survivors vs septic nonsurvivors
  contr.m[sd + ld, sd + ld] <- diag(-1, nrow = ld, ncol = ld)
  contr.m[sd + 2 * ld, sd + ld] <- diag(1, nrow = ld, ncol = ld)
  
  return(contr.m)
}

#' Return column and row coordinates of TRUE values in matrix
#'
#' @param matrix logical or 0/1 matrix
#'
#' @return data.frame with x-location vector and y-location vector
#' @export
#'
#' @examples
which.xy <- function(matrix){
  dims <- dim(matrix)
  locs <- which(matrix)
  loc.x <- (locs - 1) %% dims[2] + 1
  loc.y <- (locs - 1) %/% dims[1] + 1
  return(data.frame(x = loc.x, y = loc.y))
}

#' Compute the correlation and p value column vs column, remove NAs as in cov()'s use="pairwise.complete.obs"
#'
#' @param data a matrix or data.frame with numeric data
#' @param cmethod method of correlation test, either "Pearson" (default), "Spearman" or "Kendall"
#' @param adjust method of p value correction, see p.adjust; default is "fdr"
#'
#' @return a list with a 
#' @export
#'
#' @examples
my.corr.test <- function(x, method = "pearson", adjust = "fdr"){
  cmat <- matrix(0, nrow = ncol(x), ncol = ncol(x))
  pmat <- matrix(0, nrow = ncol(x), ncol = ncol(x))
  for (n in 1:ncol(x)){
    for (m in n:ncol(x)){
      rowset <- !is.na(x[, n]) & !is.na(x[, m])
      ctres <- cor.test(x = x[rowset, n], y = x[rowset, m], method = method)
      cmat[n,m] <- ctres$estimate
      cmat[m,n] <- ctres$estimate
      pmat[n,m] <- ctres$p.value
      pmat[m,n] <- ctres$p.value
    }
  }
  upmat <- upper.tri(pmat, diag = FALSE)
  pmat[upmat] <- p.adjust(p = pmat[upmat], method = adjust)
  #pmat <- matrix(p.adjust(p = pmat[upmat], method = adjust), ncol = ncol(data)) #checked for correct order of values in matrix
  colnames(pmat) <- colnames(x)
  rownames(pmat) <- colnames(x)
  colnames(cmat) <- colnames(x)
  rownames(cmat) <- colnames(x)
  return(list(r = cmat, p = pmat))
}


#' Compute the cubic root of x
#'
#' @param x real number x
#'
#' @return
#' @export
#'
#' @examples
cbrt_fun <- function(x){
  x ^ 1/3
}

#' Compute the cubic power of x, hence the inverse of the cubic root
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
cbrt_inv_fun <- function(x){
  x ^ 3
}

cbrt_trans <- scales::trans_new(name = "cbrt", transform = cbrt_fun, inverse = cbrt_inv_fun)

#' Compute one bootstrap sample of concentration correlations in Survivors. This is a dirty hack to reduce overhead in parralelization.
#'
#' @param n_S number of Survivor patients to sample
#' @param corr_dat list, probably with fields NS_sig_pairs and NS_grouped_sig_pairs and the corresponding coefficients
#'
#' @return the list corr_dat extended by fields for Survivors from one random sample of survivor patients
#' @export
#'
#' @examples
bootstrap_S_corr_fun <- function(r, n_S, corr_dat){
  set.seed(r)
  rand_subset <- sample.int(n = n_S, size = n_NS)
  S_corr <- my.corr.test(x = subset(subset(human_sepsis_data_normal, Survival == "S" & Day == d, select = -1:-6), 1:n_S %in% rand_subset), adjust = "fdr")
  S_grouped_corr <- my.corr.test(x = subset(human_sepsis_data_normal_grouped, Survival == "S" & Day == d, select = -1:-6), adjust = "fdr")
  ###Get significant metabolite pairs in survivors
  xy <- which.xy(S_grouped_corr$p <= 0.05)
  xy <- subset(xy, x < y)
  if (nrow(xy) > 0) {
    coeffs <- apply(xy, c(1), function(cc){ S_grouped_corr$r[cc[1], cc[2]] })
    cdir <- ifelse(coeffs > 0, "+", "-")
    corr_dat[["S_grouped_sig_pairs"]] <- paste0(cols_grouped_metab[xy$x], " ~(", cdir, ") ", cols_grouped_metab[xy$y])
    corr_dat[["S_grouped_sig_coeff"]] <- apply(xy, c(1), function(cc){ S_grouped_corr$coeff[cc[1], cc[2]] })
  }
  else{
    corr_dat[["S_grouped_sig_pairs"]] <- c()
    corr_dat[["S_grouped_sig_coeff"]] <- c()
  }
  xy <- which.xy(S_corr$p <= 0.05)
  xy <- subset(xy, x < y)
  if (nrow(xy) > 0){
    coeffs <- apply(xy, c(1), function(cc){ S_corr$r[cc[1], cc[2]] })
    cdir <- ifelse(coeffs > 0, "+", "-")
    corr_dat[["S_sig_pairs"]] <- paste0(cols_metab[xy$x], " ~(", cdir, ") ", cols_metab[xy$y])
    corr_dat[["S_sig_coeff"]] <- coeffs
  }
  else{
    corr_dat[["S_sig_pairs"]] <- c()
    corr_dat[["S_sig_coeff"]] <- c()
  }
  return(corr_dat)
}

col.na.omit <- function(data){
  return(t(na.omit(t(data))))
}

#' Create a manual color scale to be used in a ggplot for a given set of levels where one is assigned a black color. Default levels are S, NS, Control and Healthy_Fr.
#' Colors are generated as in ggplot with n-1 levels + black.
#'
#' @param name the title of the legend
#' @param levels the levels for which to generate the colors
#' @param black_pos the position in the levels vector that will be assigned black color
#'
#' @return an instance of scale_colour_manual()
#' @export
#'
#' @examples
human_col_scale <- function(name = "Group", levels = c("Septic-NS", "non-Septic-S", "Septic-S", "non-Septic-NS", "Healthy_Fr"), black_pos = 5, black_color = "grey50", ...){
  library(scales)
  color_set <- hue_pal()(length(levels) - 1)
  color_set <- c(color_set, black_color) #last is black (or grey)
  names(color_set) <- c(levels[seq_along(levels)[-black_pos]], levels[black_pos]) #assign correct level to black (last)
  return(scale_colour_manual(name = name, values = color_set, ...))
}

#' Build arrows (feature influence on biplot points). 
#' Arrows will be sorted by length. 
#' Arrow and label positions will be rescaled to fit into the plot limits and further rescaled by scale_by to deal with labels that overlap a plot border. 
#' Arrows are shortened by shorten_arr_by; the labels indicate the exact magnitude of the influence.
#'
#' @param x tranformed features
#' @param w original features, all >= 0
#' @param num number of arrows to show, beginning from the longest
#' @param xmax max x coordinate of plot
#' @param xmin min x coordinate of plot
#' @param ymax max y coordinate of plot
#' @param ymin min y coordinate of plot
#' @param scale_by (inverse) scaling factor, increase if labels overlap plot border
#' @param shorten_arw_by the absolute length to shorten an arrow by before scaling to plot limits
#'
#' @return a list of labels and arrow coordinates to feed into ggplot2::geom_path
#' @export
#'
#' @examples
make_ordination_arrows <- function(x, w, num = 5, xmax, xmin, ymax, ymin, scale_by = 1.0, shorten_arw_by = 1){
  library(vegan)
  ph_ars <- wascores(x = x, w = w, expand = TRUE)
  ph_ars_len <- apply(ph_ars, 1, function(row) sqrt(sum(row^2)))
  ph_ars_len_full <- ph_ars_len
  ph_ars <- as.data.frame(ph_ars[order(ph_ars_len, decreasing = TRUE)[1:num], ])
  ph_ars_len <- ph_ars_len[order(ph_ars_len, decreasing = TRUE)[1:num]]
  ph_ars_names <- ph_ars
  ph_ars$Group <- 1:nrow(ph_ars)
  ph_ars <- rbind(data.frame(PC1 = rep(0, nrow(ph_ars)), PC2 = rep(0, nrow(ph_ars)), Group = ph_ars$Group), ph_ars)
  ph_ars_s <- ph_ars
  ph_ars_names_s <- ph_ars_names
  xyscale <- max(ph_ars$PC1/xmax, ph_ars$PC1/xmin, ph_ars$PC2/ymax, ph_ars$PC2/ymin) * scale_by
  ph_ars_names_s[, 1:2] <- ph_ars_names[, 1:2] / xyscale
  ph_ars_s[, 1:2] <- ph_ars[, 1:2] * ((ph_ars_len - shorten_arw_by) / ph_ars_len) / xyscale
  ph_ars_names_s$label <- rownames(ph_ars_names_s)
  return(list(ph_ars_names_s = ph_ars_names_s, ph_ars_s = ph_ars_s, ph_ars_len_full = ph_ars_len_full))
}

get_lipid_isobars <- function(){
  library(data.table)
  library(stringi)
  #read isobars table
  lipid_table <- fread(input = "../../data/id_maps/tabula-List-of-Isobaric-and-Isomeric-Lipid-Species_v1_2018.csv", data.table = FALSE)
  #clean
  lipid_table <- lipid_table[-grep(pattern = "^\\[", x = lipid_table[, 2]), ]
  lipid_table <- lipid_table[!(rowSums(lipid_table == "") == ncol(lipid_table)), ]
  lipid_table <- lipid_table[!duplicated(lipid_table, MARGIN = 1:ncol(lipid_table)), ]
  #fill in blanks
  biocrates_anno <- rle(lipid_table[, 1])
  biocrates_anno_blanks <- which(biocrates_anno$values == "")
  biocrates_anno$values[biocrates_anno_blanks] <- biocrates_anno$values[biocrates_anno_blanks - 1]
  lipid_table[, 1] <- inverse.rle(biocrates_anno)
  isobar_anno <- rle(lipid_table[, 2])
  isobar_anno_blanks <- which(isobar_anno$values == "")
  isobar_anno$values[isobar_anno_blanks] <- isobar_anno$values[isobar_anno_blanks - 1]
  lipid_table[, 2] <- inverse.rle(isobar_anno)
  #add fatty acid columns
  comp <- stri_extract_all(regex = "[0-9]+\\:[0-9]", str = lipid_table[, 3]) #get fatty acid chains
  nonsphingo_idx <- grep("SM", lipid_table[, 1], invert = TRUE)
  comp[nonsphingo_idx] <- lapply(comp[nonsphingo_idx], 
                                 function(charvec) 
                                   charvec[
                                     order(
                                       as.numeric(
                                         sub(":", ".", charvec, fixed = TRUE)))]) #sort by chain length
  lipid_table$FAshort <- sapply(comp, `[[`, 1)
  lipid_table$FAlong <- sapply(comp, `[[`, 2)
  dup_FAs_idx <- c(which(duplicated(lipid_table, MARGIN = c(5, 6), fromLast = T)), which(duplicated(lipid_table, MARGIN = c(5, 6))))
  return(lipid_table)
}