#' Read the clinical data set on sepsis cases
#'
#' @return the data.frame with data, column names are identical to file content (non-standard for R)
#' @export
#'
#' @examples
get_human_sepsis_data <- function(){
  data <- read.csv(file = "../../data/measurements/Summary human sample data.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, dec = ",", check.names = FALSE)
  data <- data[, -(which(colSds(as.matrix(data[, -1:-5])) == 0) + 5)] #remove columns with fixed values
  return(data)
}

#' Read the legend to the clinical data set on sepsis cases. The legend has two group assignments of different coarseness.
#'
#' @return data.frame, first column stores the column names in the measurement data table, second the fine grained grouping, third the coarse grained grouping
#' @export
#'
#' @examples
get_human_sepsis_legend <- function(){
  data <- read.csv(file = "../../data/measurements/Legend human sample data.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, dec = ",", check.names = FALSE)
  return(data)
}

#' Read the data set on provoked sepsis and controls in rat
#'
#' @return the data.frame with data, column names are identical to file content (non-standard for R)
#' @export
#'
#' @examples
get_rat_sepsis_data <- function(){
  data <- read.csv(file = "../../data/measurements/Summary rat sample data.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, dec = ",", check.names = FALSE, blank.lines.skip = TRUE, strip.white = TRUE)
  data <- data[apply(data, 1, function(x){ sum(is.na(x)) }) < ncol(data) - 20, ] #Strip rows without any non-NA value
  data <- data[, apply(data, 2, function(x) { sum(is.na(x)) }) < nrow(data)] #Strip columns without any non-NA value
  data <- data[, apply(data, 2, function(x){ length(unique(x))}) > 1] #Strip columns with constant values
  data$BE <- data$BE - min(data$BE, na.rm = TRUE)
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
  return(data)
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
max_norm <- function(x, subset = 1:ncol(data)){
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
  res[["tpr"]] <- sum(pred > 0 & ref > 0)/sum(ref > 0)
  res[["tnr"]] <- sum(pred <= 0 & ref <= 0)/sum(ref <= 0)
  res[["fpr"]] <- sum(pred > 0 & ref <= 0)/sum(ref <= 0)
  res[["fnr"]] <- sum(pred <= 0 & ref > 0)/sum(ref > 0)
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
  res <- lapply(1:5, function(x){ c(class1_indizes[[x]], class2_indizes[[x]]) })
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
#'
#' @return a list with members day_sig_u_diff and day_sig_t_diff. These are data.frames of p-values from U-test and t-test respectively.
#' @export
#'
#' @examples
human_sig_diffs_along_days <- function(data, corr_fdr = TRUE){
  day_survival_tab <- table(data[c("Day", "Survival")])
  ##Keep days with large enough sample count
  day_survival_tab <- day_survival_tab[rowMins(day_survival_tab) >= 3,]
  ##Get significant differences for each day
  day_sig_u_diff <- data.frame(Day = as.numeric(rownames(day_survival_tab)), matrix(0, nrow = nrow(day_survival_tab), ncol = ncol(data)-5))
  colnames(day_sig_u_diff)[-1] <- colnames(data)[-1:-5]
  day_sig_t_diff <- day_sig_u_diff
  non_short_stay_pats <- unique(subset(data, Day >= 0, Patient))[[1]]
  for (d in seq_along(day_sig_u_diff$Day)){
    #Reduce to significant differences
    day <- day_sig_u_diff$Day[d]
    for (n in seq_along(day_sig_u_diff[1,-1])){
      g1 <- subset(data, subset = Survival == "S" & Day == day & Patient %in% non_short_stay_pats, select = n + 5)
      g2 <- subset(data, subset = Survival == "NS" & Day == day & Patient %in% non_short_stay_pats, select = n + 5)
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
  time_sig_u_diff <- data.frame(Time = rownames(time_survival_tab), matrix(0, nrow = nrow(time_survival_tab), ncol = ncol(data)-4))
  colnames(time_sig_u_diff)[-1] <- colnames(data)[-1:-4]
  time_sig_t_diff <- time_sig_u_diff
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
        t_res <- t.test(x = g1, y = g2, var.equal = FALSE)
        time_sig_t_diff[d, n + 1] <- t_res$p.value
      }
      else{
        time_sig_u_diff[d, n + 1] <- NA
        time_sig_t_diff[d, n + 1] <- NA
      }
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
  res <- list(time_sig_u_diff = time_sig_u_diff, time_sig_t_diff = time_sig_t_diff)
}

col.na.omit <- function(data){
  return(t(na.omit(t(data))))
}