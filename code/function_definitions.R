#' Read the clinical data set on sepsis cases
#'
#' @return the data.frame with data, column names are identical to file content (non-standard for R)
#' @export
#'
#' @examples
get_human_sepsis_data <- function(){
  data <- read.csv(file = "../../data/measurements/Summary human sample data.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, dec = ",", check.names = FALSE)
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