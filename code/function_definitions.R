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
  data <- data[rowSums(is.na(data[,-1:-4])) < ncol(data) - 4, ] #Strip rows without any non-NA value
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