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

##Import patient critical illness data
human_illness_data <- get_human_illness_data()

###########################
#Process data
###########################

#Add illness data to meta data
diag_to_add <- unique(human_illness_data$Diagnoses)
human_meta_data[diag_to_add] <- lapply(diag_to_add, function(diag) human_meta_data$`Patient ID` %in% human_illness_data$Number[human_illness_data$Diagnoses == diag])
human_meta_data[diag_to_add] <- lapply(human_meta_data[diag_to_add], as.numeric)

#Calculate statistics for differences in meta data between groups
hmd <- human_meta_data
hmd$Group <- human_data$Group[match(substring(hmd$`Patient ID`, 2), human_data$Patient)]
##Set comparisons
cmp_pairs <- list(list("Septic-S", "Septic-NS"), 
                  list("n.Septic-S", "Septic-S"), 
                  list("n.Septic-S", "n.Septic-NS"),
                  list(c("Septic-S", "Septic-NS"), c("n.Septic-S", "n.Septic-NS")))
##Set discrete variables to compare
cmp_vars_disc <- colnames(human_meta_data)[which(colnames(human_meta_data) == "Diabetes"):which(colnames(human_meta_data) == "Other hospital")]
cmp_vars_disc <- setdiff(cmp_vars_disc, "Homelessness")
cmp_vars_disc <- c("Gender", cmp_vars_disc)
##Calculate statistics for discrete variables
cmp_disc_res <- lapply(cmp_pairs, 
                       function(cmp) 
                         lapply(cmp_vars_disc, 
                                function(varb){
                                  vtab <- table(hmd[c("Group", varb)])
                                  if (length(cmp[[1]]) > 1){
                                    xvec <- colSums(vtab[cmp[[1]], ])
                                  }else{
                                    xvec <- as.numeric(vtab[cmp[[1]], ])
                                  }
                                  if (length(cmp[[2]]) > 1){
                                    yvec <- colSums(vtab[cmp[[2]], ])
                                  }else{
                                    yvec <- as.numeric(vtab[cmp[[2]], ])
                                  }
                                  mat <- matrix(c(xvec, yvec), ncol = 2, byrow = FALSE)
                                  return(chisq.test(x = mat, rescale.p = TRUE))
                                }))
cmp_disc_ps <- lapply(cmp_disc_res, sapply, `[[`, "p.value")
for (n in seq_along(cmp_disc_ps))
  names(cmp_disc_ps[[n]]) <- cmp_vars_disc
cmp_disc_sig <- lapply(cmp_disc_ps, sapply, `<`, 0.05)
cmp_disc_sig <- lapply(cmp_disc_sig, which)
cmp_disc_sig <- lapply(cmp_disc_sig, names)
##Set continuous variables
cmp_vars_cont <- colnames(human_meta_data)[c(4, 6, 7, 8)]
##Calculate statistics for continuous variables
cmp_cont_res <- lapply(cmp_pairs, 
                       function(cmp) 
                         lapply(cmp_vars_cont, 
                                function(varb){
                                  xvec <- subset(hmd, Group %in% cmp[[1]], varb)
                                  yvec <- subset(hmd, Group %in% cmp[[2]], varb)
                                  return(t.test(x = xvec, y = yvec))
                                }))
cmp_cont_ps <- lapply(cmp_cont_res, sapply, `[[`, "p.value")
for (n in seq_along(cmp_cont_ps))
  names(cmp_cont_ps[[n]]) <- cmp_vars_cont
cmp_cont_sig <- lapply(cmp_cont_ps, sapply, `<`, 0.05)
cmp_cont_sig <- lapply(cmp_cont_sig, which)
cmp_cont_sig <- lapply(cmp_cont_sig, names)
##Combine discrete and continuous results
cmp_sig <- list()
for (n in seq_along(cmp_pairs)){
  cmp_sig[[n]] <- c(cmp_cont_sig[[n]], cmp_disc_sig[[n]])
}

#Calculate means, SDs, counts and percentages
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

#Put together mean/count and sd/percentage
pub_meta_table[, 3:6] <- t(as.data.frame(lapply(1:nrow(pub_meta_table), function(row) paste0(pub_meta_table[row, 3:6], " \U00B1 ", pub_meta_a_table[row, 3:6]))))
nc <- ncol(pub_meta_table)
pub_meta_table[, 7:nc] <- t(as.data.frame(lapply(1:nrow(pub_meta_table), function(row) paste0(pub_meta_table[row, 7:nc], " (", pub_meta_a_table[row, 7:nc], ")"))))
pub_meta_table <- t(pub_meta_table)
#Prettify row names and table content
rownames(pub_meta_table) <- c("n =", "FP/CAP", gsub(x = rownames(pub_meta_table)[-1:-2], pattern = "\\.", replacement = " "))
rownames(pub_meta_table)[7] <- "Male sex - n (%)"
rownames(pub_meta_table)[c(3, 4, 6)] <- paste0(rownames(pub_meta_table)[c(3, 4, 6)], c(" - yr", " - kg", " score"))
pub_meta_table["FP/CAP", c("n.Septic-S", "n.Septic-NS")] <- "-"
pub_meta_table[(nrow(pub_meta_table) - length(diag_to_add) + 3):nrow(pub_meta_table), c("Septic-NS", "Septic-S")] <- "-"
pub_meta_table[nrow(pub_meta_table) - length(diag_to_add) + 1:2, c("n.Septic-NS", "n.Septic-S")] <- "-"
pub_meta_table <- pub_meta_table[, c(2, 1, 3, 4)]
pub_meta_table <- pub_meta_table[setdiff(rownames(pub_meta_table), "FP/CAP"), ]
if (unique(pub_meta_table["Homelessness", ]) == "0 (0.0)"){
  pub_meta_table <- pub_meta_table[setdiff(rownames(pub_meta_table), "Homelessness"), ]
}
pub_meta_table_u <- cbind(rownames(pub_meta_table), pub_meta_table)
new_rownames <- rep("", nrow(pub_meta_table_u))
new_rownames[match(c("Diabetes", "Beta blockers", "Tobacco use", "Employed", "Elective", "ED", "Fracture femur"), pub_meta_table_u[, 1])] <- 
  c("Co-morbidities - n (%)", "Medication use prior to submission - n (%)", "Social history - n (%)", "Occupational status - n (%)", "Type of admission - n (%)", "Admission source - n (%)", "Other critical illness - n (%)")
pub_meta_table_u <- cbind(new_rownames, pub_meta_table_u)
#Mark significant differences between groups in the rownames of the meta data table output
mark_idx <- lapply(lapply(cmp_sig, paste0, collapse = "|"), grep, x = pub_meta_table_u[, 2])
mark_idx_inv <- melt(mark_idx)
marks <- lapply(unique(mark_idx_inv$value), function(s) mark_idx_inv$L1[mark_idx_inv$value == s])
marks <- lapply(marks, function(m) paste0(letters[m], collapse = ", "))
marks_add <- paste0(pub_meta_table_u[unique(mark_idx_inv$value), 2], unlist(marks))
pub_meta_table_u[unique(mark_idx_inv$value), 2] <- marks_add
fwrite(x = as.data.frame(pub_meta_table_u), file = paste0(out_dir, "patient_meta_table.csv"), sep = "\t", row.names = FALSE, col.names = TRUE)