#Load libraries
library(reshape2)
library(data.table)
library(ggplot2)
library(gridExtra)
library(matrixStats)
library(heatmaply)
library(missRanger)
library(multcomp)
#library(multcompView)
library(parallel)
library(nlme)
#library(lme4)
library(car)
library(corpcor)
library(pracma)
library(psych)
library(VennDiagram)
library(igraph)
library(RColorBrewer)
library(tictoc)

#Import central functions
source("../function_definitions.R")

#Set path
out_dir <- "../../results/data_stats/"

#Make sure paths exist
if (!dir.exists(out_dir))
  dir.create(out_dir)

###########################
#Import data
###########################

##Import clinical data
human_data <- get_human_sepsis_data()

##Import corresponding group assignment
human_sepsis_legend <- get_human_sepsis_legend()
human_sepsis_legend$group[human_sepsis_legend$group == ""] <- human_sepsis_legend[human_sepsis_legend$group == "", 1]
human_sepsis_legend <- human_sepsis_legend[-1:-5, ]
human_sepsis_pheno_var_groups <- get_human_pheno_var_groups()
pheno_start <- which(human_sepsis_legend[,1] == "Urea")
human_sepsis_legend$group[match(human_sepsis_pheno_var_groups[,1], human_sepsis_legend[pheno_start:nrow(human_sepsis_legend), 1]) + pheno_start - 1] <- human_sepsis_pheno_var_groups$Mervyn_group

##Import experiment data
rat_sepsis_data <- get_rat_sepsis_data()

##Import corresponding group assignment
rat_sepsis_legend <- get_rat_sepsis_legend()

###########################
#Process data
###########################

#Impute missing values
colns <- colnames(human_data)
colnames(human_data) <- make.names(colnames(human_data))
human_data[,-1:-5] <- missRanger(data = human_data[,-1:-5])
colnames(human_data) <- colns

#Seperate septic and nonseptic patients
human_nonsepsis_data <- human_data[human_data$`CAP / FP` == "-", ]
human_sepsis_data <- human_data[human_data$`CAP / FP` != "-", ]

#Scale mesaurement values by standardizatio
human_data_normal <- human_data
human_data_normal[,-1:-5] <- scale(x = human_data_normal[,-1:-5])
human_sepsis_data_normal <- human_data_normal[human_data_normal$`CAP / FP` != "-", ]
human_sepsis_data_normal$Patient <- as.factor(human_sepsis_data_normal$Patient)
human_sepsis_data_normal$Day <- as.factor(human_sepsis_data_normal$Day)
human_nonsepsis_data_normal <- human_data_normal[human_data_normal$`CAP / FP` == "-", ]
human_nonsepsis_data_normal$Patient <- as.factor(human_nonsepsis_data_normal$Patient)
human_nonsepsis_data_normal$Day <- as.factor(human_nonsepsis_data_normal$Day)

#Find outlier sample
##Group metabolites
coarse_group_list <- human_sepsis_legend[human_sepsis_legend[,1] %in% colnames(human_data)[-1:-5], 3] #lucky ... no col in data without match in legend
human_data_grouped <- cbind(human_data[,1:5], matrix(0, nrow = nrow(human_data_normal), ncol=length(unique(coarse_group_list))))
colnames(human_data_grouped)[-1:-5] <- unique(coarse_group_list)
human_data_grouped$Patient <- as.factor(human_data_grouped$Patient)
human_data_grouped$Day <- as.factor(human_data_grouped$Day)
###Split for metabolites and "phenotypical" factors
split_start <- which(colnames(human_data) == "Urea")
pheno_sel <- split_start:ncol(human_data)
metab_sel <- 6:(split_start-1)
group_pheno_sel <- which(colnames(human_data_grouped) %in% unique(coarse_group_list[pheno_sel - 5]))
group_metab_sel <- which(colnames(human_data_grouped) %in% unique(coarse_group_list[metab_sel - 5]))
hd <- human_data
hd[, pheno_sel] <- scale(hd[, pheno_sel])
for (n in 1:nrow(human_data)){
  m_agg <- tapply(X = t(hd[n, metab_sel]), INDEX = factor(coarse_group_list[metab_sel - 5]), FUN = sum)
  p_agg <- tapply(X = t(hd[n, pheno_sel]), INDEX = factor(coarse_group_list[pheno_sel - 5]), FUN = mean)
  b_agg <- c(m_agg, p_agg)
  human_data_grouped[n, -1:-5] <- b_agg[match(colnames(human_data_grouped)[-1:-5], names(b_agg))]
}
human_sepsis_data_grouped <- human_data_grouped[human_data_grouped$`CAP / FP` != "-", ]
human_sepsis_data_grouped[, c("CAP / FP", "Patient", "Day")] <- lapply(human_sepsis_data_grouped[, c("CAP / FP", "Patient", "Day")], as.factor)
human_sepsis_data_normal_grouped <- human_sepsis_data_grouped
human_sepsis_data_normal_grouped[, -1:-5] <- scale(human_sepsis_data_normal_grouped[, -1:-5])
human_data_normal_grouped <- human_data_grouped
human_data_normal_grouped[, -1:-5] <- scale(human_data_normal_grouped[, -1:-5])

#Whiten each patient
pat_list <- human_data$Patient
pat_list <- table(pat_list)[table(pat_list) > 1]
max_pat_samples <- which(pat_list == max(pat_list))
pat_list <- names(pat_list)
target_pat <- max_pat_samples[which(human_data$`CAP / FP`[match(pat_list, human_data$Patient)][max_pat_samples] != "-")]
human_data <- subset(human_data, Patient %in% pat_list, c(1:5, metab_sel))
excl_cols <- list()
for (pat in pat_list){
  source_dat <- human_data[human_data$Patient == pat, -1:-5]
  source_dat <- as.matrix(scale(source_dat))
  C_s <- cov(source_dat)
  y <- which.xy(is.na(C_s))[, 2]
  excl_cols[[pat]] <- which(table(y) == ncol(C_s))
}
ecol <- unique(unlist(excl_cols)) + 5
human_data <- human_data[, -ecol]
human_whitened_data <- human_data
pat_list <- pat_list[-target_pat]
target_dat <- human_data[human_data$Patient == target_pat, -1:-5]
target_dat <- as.matrix(scale(target_dat))
C_t <- cov(target_dat) + diag(nrow = ncol(target_dat))
svd_C_t <- svd(C_t)
with(svd_C_t, sqrt_C_t <- u %*% diag(1/d) %*% t(v))
for (pat in pat_list){
  source_dat <- human_data[human_data$Patient == pat, -1:-5]
  source_dat <- as.matrix(scale(source_dat))
  C_s <- cov(source_dat) + diag(nrow = ncol(source_dat))
  svd_C_s <- svd(C_s)
  with(svd_C_s, isqrt_C_s <- u %*% diag(sqrt(1/d)) %*% t(v))
  whitened_dat <- source_dat %*% isqrt_C_s %*% sqrt_C_t
  human_whitened_data[human_data$Patient == pat, -1:-5] <- whitened_dat
}
p <- prcomp(human_whitened_data[, -1:-5])
plot(p$x[, 1:2], col = factor(human_whitened_data$`CAP / FP`))