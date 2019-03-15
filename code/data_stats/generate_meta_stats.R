#Load libraries
library(reshape2)
library(data.table)
library(stringi)
library(scales)
library(ggplot2)
library(gridExtra)
library(pheatmap)
library(matrixStats)
library(heatmaply)
library(missRanger)
library(multcomp)
library(vegan)
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
human_sepsis_legend <- human_sepsis_legend[-1:-6, ]
human_sepsis_legend <- human_sepsis_legend[human_sepsis_legend[, 1] %in% colnames(human_data), ] #make sure the legend content matches the measurement variables
human_sepsis_pheno_var_groups <- get_human_pheno_var_groups()
pheno_start <- which(human_sepsis_legend[,1] == "Urea")
human_sepsis_legend$group[match(human_sepsis_pheno_var_groups[,1], human_sepsis_legend[pheno_start:nrow(human_sepsis_legend), 1]) + pheno_start - 1] <- human_sepsis_pheno_var_groups$Mervyn_group

##Import meta data
human_meta_data <- get_human_meta_data()

###########################
#Process data
###########################

patient_numbers <- lapply(unique(human_data$Group), function(gr) unique(subset(human_data, Group == gr)$Patient))
names(patient_numbers) <- unique(human_data$Group)
patient_numbers <- lapply(patient_numbers, function(x) paste0("N", x))
pub_meta_mean_table <- lapply(patient_numbers, 
                         function(gr){ 
                           d <- human_meta_data[human_meta_data$`Patient ID` %in% gr, ]; 
                           df <- data.frame(
                             n = format(nrow(d), trim = T),
                             fp_cap = paste0(sum(d$Diagnosis == "fp"), "/", sum(d$Diagnosis == "pn")),
                             format(t(colMeans(d[c("Age", "Weight", "SOFA score", "APACHE II")])), digits = 2, nsmall = 1, trim = T),
                             male_sex = sum(d$Gender == "M"),
                             format(t(colSums(d[-1:-8])), digits = 1, nsmall = 0, trim = T),
                             stringsAsFactors = FALSE
                           ) })
pub_meta_add_table <- lapply(patient_numbers, 
                              function(gr){ 
                                d <- human_meta_data[human_meta_data$`Patient ID` %in% gr, ]; 
                                df <- data.frame(
                                  n = format(nrow(d), trim = T),
                                  fp_cap = "",
                                  format(t(sapply(d[c("Age", "Weight", "SOFA score", "APACHE II")], sd)), digits = 2, nsmall = 1, trim = T),
                                  male_sex = format(sum(d$Gender == "M") / nrow(d) * 100, digits = 1, nsmall = 1),
                                  format(t(colMeans(d[-1:-8])) * 100, digits = 2, nsmall = 1, trim = T),
                                  stringsAsFactors = FALSE
                                ) })
pub_meta_table <- Reduce("rbind", pub_meta_mean_table)
pub_meta_a_table <- Reduce("rbind", pub_meta_add_table)
rownames(pub_meta_table) <- names(patient_numbers)

pub_meta_table[, 3:6] <- t(as.data.frame(lapply(1:nrow(pub_meta_table), function(row) paste0(pub_meta_table[row, 3:6], " \U00B1 ", pub_meta_a_table[row, 3:6]))))
nc <- ncol(pub_meta_table)
pub_meta_table[, 7:nc] <- t(as.data.frame(lapply(1:nrow(pub_meta_table), function(row) paste0(pub_meta_table[row, 7:nc], " (", pub_meta_a_table[row, 7:nc], ")"))))
pub_meta_table <- t(pub_meta_table)
rownames(pub_meta_table) <- c("n =", "FP/CAP", gsub(x = rownames(pub_meta_table)[-1:-2], pattern = "\\.", replacement = " "))
rownames(pub_meta_table)[7] <- "Male sex - n (%)"
rownames(pub_meta_table)[c(3, 4, 6)] <- paste0(rownames(pub_meta_table)[c(3, 4, 6)], c(" - yr", " - kg", " score"))
pub_meta_table["FP/CAP", c("Nonsep-S", "Nonsep-NS")] <- "-"
pub_meta_table <- pub_meta_table[, c(2, 1, 3, 4)]
fwrite(x = as.data.frame(pub_meta_table), file = paste0(out_dir, "patient_meta_table.csv"), sep = "\t", row.names = TRUE, col.names = TRUE)
