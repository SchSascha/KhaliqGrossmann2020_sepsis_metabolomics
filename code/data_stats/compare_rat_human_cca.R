library(matrixStats)
library(data.table)
library(missRanger)
library(nscancor)
library(glmnet)
library(VennDiagram)

source("../function_definitions.R")

out_dir <- "../../results/data_stats/"

#Data input
human_sepsis_data <- get_human_sepsis_data()
rat_sepsis_data <- get_rat_sepsis_data()

#Seperate septic and nonseptic samples
human_nonsepsis_data <- human_sepsis_data[human_sepsis_data$`CAP / FP` == "-", ]
human_sepsis_data <- human_sepsis_data[human_sepsis_data$`CAP / FP` != "-", ]

rat_nonsepis_data <- subset(rat_sepsis_data, group == "control" & material == "plasma")
rat_sepis_data <- subset(rat_sepsis_data, group != "control" & material == "plasma")

#Clean and scale data
human_sepsis_data_normal <- human_sepsis_data
human_sepsis_data_normal[,-1:-5] <- scale(human_sepsis_data_normal[,-1:-5])

human_split_col <- which(colnames(human_sepsis_data) == "LCA")

col <- colnames(human_sepsis_data_normal)
colnames(human_sepsis_data_normal) <- make.names(col)
human_sepsis_data_normal[, 6:ncol(human_sepsis_data_normal)] <- missRanger(human_sepsis_data_normal[,6:ncol(human_sepsis_data_normal)])
colnames(human_sepsis_data_normal) <- col

human_sepsis_data_metab <- subset(human_sepsis_data_normal, Day < 2, 6:human_split_col)
human_sepsis_data_pheno <- subset(human_sepsis_data_normal, Day < 2, (human_split_col + 1):ncol(human_sepsis_data))

rat_sepsis_data_normal <- rat_sepsis_data
rat_sepsis_data_normal[,-1:-4] <- scale(rat_sepsis_data_normal[,-1:-4])

rat_split_col <- which(colnames(rat_sepsis_data) == "LCA")

col <- colnames(rat_sepsis_data_normal)
colnames(rat_sepsis_data_normal) <- make.names(col)
rat_sepsis_data_normal[, 5:ncol(rat_sepsis_data_normal)] <- missRanger(rat_sepsis_data_normal[,5:ncol(rat_sepsis_data_normal)])
colnames(rat_sepsis_data_normal) <- col

rat_sepsis_data_metab <- subset(rat_sepsis_data_normal, material == "plasma", select = 5:rat_split_col)
rat_sepsis_data_pheno <- subset(rat_sepsis_data_normal, material == "plasma", select = (rat_split_col + 1):ncol(rat_sepsis_data))

rat_sepsis_data_metab <- rat_sepsis_data_metab[, !colAnyNAs(as.matrix(rat_sepsis_data_metab))]
rat_sepsis_data_metab <- rat_sepsis_data_metab[, colSds(as.matrix(rat_sepsis_data_metab)) >= 10^-10]

#Reduce data to same vars
metab_var_isect <- intersect(colnames(human_sepsis_data_metab), colnames(rat_sepsis_data_metab))
pheno_var_isect <- intersect(colnames(human_sepsis_data_pheno), colnames(rat_sepsis_data_pheno))
rat_sepsis_data_metab <- subset(rat_sepsis_data_metab, TRUE, select = metab_var_isect)
human_sepsis_data_metab <- subset(human_sepsis_data_metab, TRUE, select = metab_var_isect)
rat_sepsis_data_pheno <- subset(rat_sepsis_data_pheno, TRUE, select = pheno_var_isect)
human_sepsis_data_pheno <- subset(human_sepsis_data_pheno, TRUE, select = pheno_var_isect)

#Run CCA
dfmax_w <- c(10, 8, 5)
ypredict <- function(x, yc, cc) {
  en <- glmnet::glmnet(x, yc, alpha = 0.5, intercept = FALSE,
                       dfmax = dfmax_w[cc], lower.limits = 0)
  W <- coef(en)
  return(W[2:nrow(W), ncol(W)])
}
dfmax_v <- c(20, 16, 10)
xpredict <- function(y, xc, cc) {
  en <- glmnet::glmnet(y, xc, alpha = 0.5, intercept = FALSE,
                       dfmax = dfmax_v[cc], lower.limits = 0)
  V <- coef(en)
  return(V[2:nrow(V), ncol(V)])
}

human_nscc <- nscancor(x = human_sepsis_data_metab, y = human_sepsis_data_pheno, xpredict = xpredict, ypredict = ypredict, nvar = 3, xscale = T, yscale = T)
rat_nscc <- nscancor(x = rat_sepsis_data_metab, y = rat_sepsis_data_pheno, xpredict = xpredict, ypredict = ypredict, nvar = 3, xscale = T, yscale = T)

all_nscc <- nscancor(x = rbind(human_sepsis_data_metab, rat_sepsis_data_metab), 
                     y = rbind(human_sepsis_data_pheno, rat_sepsis_data_pheno), 
                     xpredict = xpredict,
                     ypredict = ypredict,
                     nvar = 3, 
                     xscale = T,
                     yscale = T)

nx_all <- names(all_nscc$ycoef[all_nscc$ycoef[,1] > 0, 1])
nx_human <- names(human_nscc$ycoef[human_nscc$ycoef[,1] > 0, 1])
nx_rat <- names(rat_nscc$ycoef[rat_nscc$ycoef[,1] > 0, 1])

venn.diagram(x = list(Combined = nx_all, Human = nx_human, Rat = nx_rat), 
             filename = paste0(out_dir,"CCA_overlap.tiff"))
