#Load libraries
library(reshape2)
library(data.table)
library(stringi)
library(scales)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(scatterpie)
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
library(xlsx)

#Import central functions
source("../function_definitions.R")

#Set path
out_dir <- "../../results/data_stats/"

#Make sure paths exist
if (!dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE)

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

##Import experiment data
rat_sepsis_data <- get_rat_sepsis_data()

##Import corresponding group assignment
rat_sepsis_legend <- get_rat_sepsis_legend()

##Import normal plasma ranges
human_healthy_ranges <- get_french_normal_data()

##Import validation data
human_sepsis_val_data <- get_Ferrario_validation_data()
human_sepsis_val_data <- human_sepsis_val_data[, -5] #column 5 has a lysoPC where all values at Day 6 are missing

###########################
#Process data
###########################

#Match metabolite names in healthy reference
met_mat_name <- match(human_healthy_ranges$Metabolite, human_sepsis_legend$name)
met_mat_id <- match(human_healthy_ranges$Metabolite, human_sepsis_legend[, 1])
met_mat_name[is.na(met_mat_name)] <- met_mat_id[is.na(met_mat_name)]
stopifnot(!is.na(met_mat_name))
human_healthy_ranges$ID <- human_sepsis_legend[met_mat_name, 1]

#Impute missing values
colns <- colnames(human_data)
colnames(human_data) <- make.names(colnames(human_data))
imputables <- 7:which(colnames(human_data) == "H1")
human_data[, imputables] <- missRanger(data = human_data[, imputables], seed = 105431)
colnames(human_data) <- colns

colns <- colnames(human_sepsis_val_data)
colnames(human_sepsis_val_data) <- make.names(colnames(human_sepsis_val_data))
human_sepsis_val_data[, -1:-4] <- missRanger(data = human_sepsis_val_data[, -1:-4])
colnames(human_sepsis_val_data) <- colns

#Seperate septic and non-Septic patients
human_nonsepsis_data <- human_data[human_data$`CAP / FP` == "-", ]
human_sepsis_data <- human_data[human_data$`CAP / FP` != "-", ]

#Scale mesaurement values by standardization
human_data_normal <- human_data
human_data_normal[,-1:-6] <- scale(x = human_data_normal[,-1:-6])
human_sepsis_data_normal <- human_data_normal[human_data_normal$`CAP / FP` != "-", ]
human_sepsis_data_normal$Patient <- as.factor(human_sepsis_data_normal$Patient)
human_sepsis_data_normal$Day <- as.factor(human_sepsis_data_normal$Day)
human_nonsepsis_data_normal <- human_data_normal[human_data_normal$`CAP / FP` == "-", ]
human_nonsepsis_data_normal$Patient <- as.factor(human_nonsepsis_data_normal$Patient)
human_nonsepsis_data_normal$Day <- as.factor(human_nonsepsis_data_normal$Day)

#Find outlier sample
##Group metabolites
coarse_group_list <- human_sepsis_legend[human_sepsis_legend[,1] %in% colnames(human_data)[-1:-6], 3] #lucky ... no col in data without match in legend
human_data_grouped <- cbind(human_data[,1:6], matrix(0, nrow = nrow(human_data_normal), ncol=length(unique(coarse_group_list))))
colnames(human_data_grouped)[-1:-6] <- unique(coarse_group_list)
human_data_grouped$Patient <- as.factor(human_data_grouped$Patient)
human_data_grouped$Day <- as.factor(human_data_grouped$Day)
###Split for metabolites and "phenotypical" factors
split_start <- which(colnames(human_data) == "Urea")
pheno_sel <- split_start:ncol(human_data)
metab_sel <- 7:(split_start-1)
group_pheno_sel <- which(colnames(human_data_grouped) %in% unique(coarse_group_list[pheno_sel - 6]))
group_metab_sel <- which(colnames(human_data_grouped) %in% unique(coarse_group_list[metab_sel - 6]))
hd <- human_data
hd[, pheno_sel] <- scale(hd[, pheno_sel])
for (n in 1:nrow(human_data)){
  m_agg <- tapply(X = t(hd[n, metab_sel]), INDEX = factor(coarse_group_list[metab_sel - 6]), FUN = sum)
  p_agg <- tapply(X = t(hd[n, pheno_sel]), INDEX = factor(coarse_group_list[pheno_sel - 6]), FUN = mean)
  b_agg <- c(m_agg, p_agg)
  human_data_grouped[n, -1:-6] <- b_agg[match(colnames(human_data_grouped)[-1:-6], names(b_agg))]
}
human_sepsis_data_grouped <- human_data_grouped[human_data_grouped$`CAP / FP` != "-", ]
human_sepsis_data_grouped[, c("CAP / FP", "Patient", "Day")] <- lapply(human_sepsis_data_grouped[, c("CAP / FP", "Patient", "Day")], as.factor)
human_sepsis_data_normal_grouped <- human_sepsis_data_grouped
human_sepsis_data_normal_grouped[, -1:-6] <- scale(human_sepsis_data_normal_grouped[, -1:-6])
human_data_normal_grouped <- human_data_grouped
human_data_normal_grouped[, -1:-6] <- scale(human_data_normal_grouped[, -1:-6])

#Build measurement characteristics table
met_group_table <- as.data.frame(table(human_sepsis_legend$group[metab_sel - 6]))

#Get data overview
##Get overview of sample distribution along days
human_sig_diff_res <- human_sig_diffs_along_days(human_sepsis_data, corr_fdr = TRUE)
day_sig_u_diff <- human_sig_diff_res$day_sig_u_diff
day_sig_t_diff <- human_sig_diff_res$day_sig_t_diff

#Get sig vars
sig_u_class <- na.omit(colnames(day_sig_u_diff[,-1])[colAnys(day_sig_u_diff[, -1] <= 0.05)])
sig_t_class <- na.omit(colnames(day_sig_t_diff[,-1])[colAnys(day_sig_t_diff[, -1] <= 0.05)])

#Get sig var pos
day_sig_t_diff_pos_long <- get_sig_var_pos(day_sig_t_diff, time_var = "Day")
day_sig_u_diff_pos_long <- get_sig_var_pos(day_sig_u_diff, time_var = "Day")

#Find time points in long time course data where changes are significant
tanova_day_set <- c(0,1,2,3)
names(tanova_day_set) <- tanova_day_set
full_tanova_data <- subset(human_data, Day %in% tanova_day_set)
full_tanova_data$Survival[full_tanova_data$`CAP / FP` == "-"] <- "non-Septic"

##Run repeated measures ANOVA as depicted in the R Companion at http://rcompanion.org/handbook/I_09.html
fml <- concentration ~ Day*Survival
met_set <- colnames(full_tanova_data)[metab_sel]
names(met_set) <- met_set
ftd_c <- full_tanova_data
ftd_c$Survival <- sub(pattern = "NS|S", x = as.character(ftd_c$Survival), replacement = "Seps")
ftd_c <- ftd_c[, -pheno_sel]
ftd_s <- subset(full_tanova_data, Survival != "non-Septic")
anova.car.c <- t3ANOVA(data = ftd_c, random = ~1|Patient, formula = fml, use.corAR = TRUE, col.set = met_set, id.vars = c("Survival", "Day", "Patient"), control = lmeControl(msMaxIter = 100))
anova.car.s <- t3ANOVA(data = ftd_s, random = ~1|Patient, formula = fml, use.corAR = TRUE, col.set = met_set, id.vars = c("Survival", "Day", "Patient"), control = lmeControl(msMaxIter = 100))
anova.car.c.pre.ps <- anova.car.c$ps
anova.car.s.pre.ps <- anova.car.s$ps
sig.anova.car.c.pre.class <- colnames(anova.car.c.pre.ps)[colAnys(anova.car.c.pre.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
sig.anova.car.s.pre.class <- colnames(anova.car.s.pre.ps)[colAnys(anova.car.s.pre.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
anova.car.c.ps <- (apply(anova.car.c.pre.ps, 2, function(row){ p.adjust(p = row, method = "fdr") }))
anova.car.s.ps <- (apply(anova.car.s.pre.ps, 2, function(row){ p.adjust(p = row, method = "fdr") }))
#anova.car.c.ps <- matrix(p.adjust(anova.car.c.ps[3:4, ], method = "fdr"), nrow = 2)
#anova.car.s.ps <- matrix(p.adjust(anova.car.s.ps[3:4, ], method = "fdr"), nrow = 2)
sig.anova.car.c.class <- colnames(anova.car.c.ps)[colAnys(anova.car.c.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
sig.anova.car.c.class <- intersect(sig.anova.car.c.class, sig.anova.car.c.pre.class)
sig.anova.car.s.class <- colnames(anova.car.s.ps)[colAnys(anova.car.s.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
sig.anova.car.s.class <- intersect(sig.anova.car.s.class, sig.anova.car.s.pre.class)
insig.anova.car.c.pre.class <- colnames(anova.car.c.pre.ps)[colAlls(anova.car.c.pre.ps[c("Survival", "Day:Survival"), ] > 0.05)]
insig.anova.car.s.pre.class <- colnames(anova.car.s.pre.ps)[colAlls(anova.car.s.pre.ps[c("Survival", "Day:Survival"), ] > 0.05)]
insig.anova.car.c.class <- colnames(anova.car.c.ps)[colAlls(anova.car.c.ps[c("Survival", "Day:Survival"), ] >= 0.95)]
insig.anova.car.c.class <- intersect(insig.anova.car.c.class, insig.anova.car.c.pre.class)
insig.anova.car.s.class <- colnames(anova.car.s.ps)[colAlls(anova.car.s.ps[c("Survival", "Day:Survival"), ] >= 0.95)]
insig.anova.car.s.class <- intersect(insig.anova.car.s.class, insig.anova.car.s.pre.class)

###Post hoc analysis of when the change happened
####Build models
fml <- concentration ~ DaySurv - 1 #without intercept because Anova(..., type = 3) gives an error with intercept
ftd_ph <- full_tanova_data
ftd_ph$DaySurv <- interaction(ftd_ph$Day, ftd_ph$Survival, drop = TRUE)
anova.car.ph.models <- mclapply(met_set, fit_lin_mod_lme, data = ftd_ph, formula = fml, id.vars = c("Day", "Patient", "DaySurv"), random = ~1|Patient, use.corAR = TRUE, control = lmeControl(msMaxIter = 100))
anova.car.ph.models <- mclapply(anova.car.ph.models, eval)
####Construct constrast matrix
contr.m <- get_S_NS_C_contrmat(tanova_day_set)
daySurv <- interaction(factor(full_tanova_data$Day), factor(full_tanova_data$Survival), drop = TRUE)
rownames(contr.m) <- levels(daySurv)
####Actual post hoc test
anova.car.ph.res <- mclapply(anova.car.ph.models, function(e) summary(glht(e, linfct = mcp(DaySurv = t(contr.m)))))
anova.car.ph.ps <- sapply(anova.car.ph.res, function(e) e$test$pvalues)
anova.car.ph.sig.contr <- lapply(data.frame(anova.car.ph.ps), function(col) which(col <= 0.05))
names(anova.car.ph.ps) <- names(anova.car.ph.models)
names(anova.car.ph.sig.contr) <- names(anova.car.ph.models)

###Also on phenomenological variables
fml <- concentration ~ Day*Survival
met_set <- colnames(full_tanova_data)[pheno_sel]
names(met_set) <- met_set
ftd_s <- subset(full_tanova_data, Survival != "non-Septic" & Day %in% tanova_day_set)
anova.car.s.pheno <- t3ANOVA(data = ftd_s, random = ~1|Patient, formula = fml, use.corAR = TRUE, col.set = met_set, id.vars = c("Survival", "Day", "Patient"))
anova.car.s.pheno.pre.ps <- anova.car.s.pheno$ps
sig.anova.car.s.pheno.pre.class <- colnames(anova.car.s.pheno.pre.ps)[colAnys(anova.car.s.pheno.pre.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
insig.anova.car.s.pheno.pre.class <- colnames(anova.car.s.pheno.pre.ps)[colAlls(anova.car.s.pheno.pre.ps[c("Survival", "Day:Survival"), ] > 0.05)]
anova.car.s.pheno.ps <- (apply(anova.car.s.pheno.pre.ps, 2, p.adjust, method = "fdr"))
sig.anova.car.s.pheno.class <- colnames(anova.car.s.pheno.ps)[colAnys(anova.car.s.pheno.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
sig.anova.car.s.pheno.class <- intersect(sig.anova.car.s.pheno.class, sig.anova.car.s.pheno.pre.class)
insig.anova.car.s.pheno.class <- colnames(anova.car.s.pheno.ps)[colAlls(anova.car.s.pheno.ps[c("Survival", "Day:Survival"), ] > 0.95)]
insig.anova.car.s.pheno.class <- intersect(insig.anova.car.s.pheno.class, insig.anova.car.s.pheno.pre.class)

###Post hoc analysis of when the change happened
####Build models, reuse met_set
fml <- concentration ~ DaySurv - 1
ftd_ph <- subset(full_tanova_data[full_tanova_data$`CAP / FP` != "-", ], Day %in% tanova_day_set)
ftd_ph$DaySurv <- interaction(ftd_ph$Day, ftd_ph$Survival, drop = TRUE)
anova.car.ph.pheno.models <- mclapply(met_set, fit_lin_mod_lme, data = ftd_ph, formula = fml, id.vars = c("Day", "Patient", "DaySurv"), random = ~1|Patient, use.corAR = TRUE, control = lmeControl(msMaxIter = 100))
anova.car.ph.pheno.models <- mclapply(anova.car.ph.pheno.models, function(e) try(eval(e)))
####Construct constrast matrix
contr.m <- get_S_NS_C_contrmat(tanova_day_set)[-1:-length(tanova_day_set), 5:8] #pheno vars only for septic patients, so skip non-Septic comparisons
daySurv <- interaction(factor(full_tanova_data$Day), factor(full_tanova_data$Survival), drop = TRUE)
rownames(contr.m) <- levels(daySurv)[-1:-4]
####Actual post hoc test
anova.car.ph.pheno.res <- mclapply(anova.car.ph.pheno.models, function(e) summary(glht(e, linfct = mcp(DaySurv = t(contr.m)))))
anova.car.ph.pheno.ps <- sapply(anova.car.ph.pheno.res, function(e) e$test$pvalues)
anova.car.ph.pheno.sig.contr <- lapply(data.frame(anova.car.ph.pheno.ps), function(col) which(col <= 0.05))
names(anova.car.ph.pheno.ps) <- names(anova.car.ph.pheno.models)
names(anova.car.ph.pheno.sig.contr) <- names(anova.car.ph.pheno.models)

#ANOVA to find where non-sepsis overlaps with either septic survivors or nonsurvivors
fml <- concentration ~ Day*Survival
met_set <- colnames(full_tanova_data)[metab_sel]
names(met_set) <- met_set
ftd_cvs <- full_tanova_data
ftd_cvns <- full_tanova_data
ftd_csvns <- full_tanova_data
ftd_cnsvs <- full_tanova_data
ftd_cvs <- ftd_cvs[ftd_cvs$`CAP / FP` == "-" | ftd_cvs$Survival == "S", ]
ftd_cvs$Survival[ftd_cvs$`CAP / FP` == "-"] <- "non-Septic"
ftd_cvns <- ftd_cvns[ftd_cvns$`CAP / FP` == "-" | ftd_cvns$Survival == "NS", ]
ftd_cvns$Survival[ftd_cvns$`CAP / FP` == "-"] <- "non-Septic"
ftd_csvns$Survival[ftd_csvns$Survival == "S" | ftd_csvns$`CAP / FP` == "-"] <- "C_or_S"
ftd_cnsvs$Survival[ftd_cnsvs$Survival == "NS" | ftd_csvns$`CAP / FP` == "-"] <- "C_or_NS"
anova.car.cvs <- t3ANOVA(data = ftd_cvs, random = ~1|Patient, formula = fml, use.corAR = TRUE, col.set = met_set, id.vars = c("Survival", "Day", "Patient"), control = lmeControl(msMaxIter = 100))
anova.car.cvns <- t3ANOVA(data = ftd_cvns, random = ~1|Patient, formula = fml, use.corAR = TRUE, col.set = met_set, id.vars = c("Survival", "Day", "Patient"), control = lmeControl(msMaxIter = 100))
anova.car.csvns <- t3ANOVA(data = ftd_csvns, random = ~1|Patient, formula = fml, use.corAR = TRUE, col.set = met_set, id.vars = c("Survival", "Day", "Patient"), control = lmeControl(msMaxIter = 100))
anova.car.cnsvs <- t3ANOVA(data = ftd_cnsvs, random = ~1|Patient, formula = fml, use.corAR = TRUE, col.set = met_set, id.vars = c("Survival", "Day", "Patient"), control = lmeControl(msMaxIter = 100))
anova.car.cvs.pre.ps <- anova.car.cvs$ps
anova.car.cvns.pre.ps <- anova.car.cvns$ps
anova.car.csvns.pre.ps <- anova.car.csvns$ps
anova.car.cnsvs.pre.ps <- anova.car.cnsvs$ps
sig.anova.car.cvs.pre.class <- colnames(anova.car.cvs.pre.ps)[colAnys(anova.car.cvs.pre.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
sig.anova.car.cvns.pre.class <- colnames(anova.car.cvns.pre.ps)[colAnys(anova.car.cvns.pre.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
sig.anova.car.csvns.pre.class <- colnames(anova.car.csvns.pre.ps)[colAnys(anova.car.csvns.pre.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
sig.anova.car.cnsvs.pre.class <- colnames(anova.car.cnsvs.pre.ps)[colAnys(anova.car.cnsvs.pre.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
anova.car.cvs.ps <- (apply(anova.car.cvs.pre.ps, 2, function(row){ p.adjust(p = row, method = "fdr") }))
anova.car.cvns.ps <- (apply(anova.car.cvns.pre.ps, 2, function(row){ p.adjust(p = row, method = "fdr") }))
anova.car.csvns.ps <- (apply(anova.car.csvns.pre.ps, 2, function(row){ p.adjust(p = row, method = "fdr") }))
anova.car.cnsvs.ps <- (apply(anova.car.cnsvs.pre.ps, 2, function(row){ p.adjust(p = row, method = "fdr") }))
sig.anova.car.cvs.class <- colnames(anova.car.cvs.ps)[colAnys(anova.car.cvs.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
sig.anova.car.cvs.class <- intersect(sig.anova.car.cvs.class, sig.anova.car.cvs.pre.class)
sig.anova.car.cvns.class <- colnames(anova.car.cvns.ps)[colAnys(anova.car.cvns.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
sig.anova.car.cvns.class <- intersect(sig.anova.car.cvns.class, sig.anova.car.cvns.pre.class)
sig.anova.car.csvns.class <- colnames(anova.car.csvns.ps)[colAnys(anova.car.csvns.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
sig.anova.car.csvns.class <- intersect(sig.anova.car.csvns.class, sig.anova.car.csvns.pre.class)
sig.anova.car.cnsvs.class <- colnames(anova.car.cnsvs.ps)[colAnys(anova.car.cnsvs.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
sig.anova.car.cnsvs.class <- intersect(sig.anova.car.cnsvs.class, sig.anova.car.cnsvs.pre.class)

c_is_s_sig <- intersect(sig.anova.car.cvs.class, sig.anova.car.cnsvs.class)
c_is_ns_sig <- intersect(sig.anova.car.cvns.class, sig.anova.car.csvns.class)

#ANOVA to find where non-septic Nonsurvivors differ from septic Nonsurvivors
fml <- concentration ~ Day*Survival
met_set <- colnames(full_tanova_data)[metab_sel]
names(met_set) <- met_set
ftd_s <- subset(full_tanova_data, Group %in% c("Septic-NS", "non-Septic-NS"))
ftd_s <- ftd_s[, -pheno_sel]
anova.car.nssepnon <- t3ANOVA(data = ftd_s, random = ~1|Patient, formula = fml, use.corAR = TRUE, col.set = met_set, id.vars = c("Survival", "Day", "Patient"), control = lmeControl(msMaxIter = 100))
anova.car.nssepnon.pre.ps <- anova.car.nssepnon$ps
sig.anova.car.nssepnon.pre.class <- colnames(anova.car.nssepnon.pre.ps)[colAnys(anova.car.nssepnon.pre.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
anova.car.nssepnon.ps <- (apply(anova.car.nssepnon.pre.ps, 2, function(row){ p.adjust(p = row, method = "fdr") }))
sig.anova.car.nssepnon.class <- colnames(anova.car.nssepnon.ps)[colAnys(anova.car.nssepnon.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
sig.anova.car.nssepnon.class <- intersect(sig.anova.car.nssepnon.class, sig.anova.car.nssepnon.pre.class)

#ANOVA to find where non-septic Nonsurvivors differ from non-septic Survivors
fml <- concentration ~ Day*Group
met_set <- colnames(full_tanova_data)[metab_sel]
names(met_set) <- met_set
ftd_s <- subset(full_tanova_data, Group %in% c("non-Septic-S", "non-Septic-NS"))
ftd_s <- ftd_s[, -pheno_sel]
anova.car.cnscs <- t3ANOVA(data = ftd_s, random = ~1|Patient, formula = fml, use.corAR = TRUE, col.set = met_set, id.vars = c("Group", "Day", "Patient"), control = lmeControl(msMaxIter = 100))
anova.car.cnscs.pre.ps <- anova.car.cnscs$ps
sig.anova.car.cnscs.pre.class <- colnames(anova.car.cnscs.pre.ps)[colAnys(anova.car.cnscs.pre.ps[c("Group", "Day:Group"), ] <= 0.05)]
anova.car.cnscs.ps <- (apply(anova.car.cnscs.pre.ps, 2, function(row){ p.adjust(p = row, method = "fdr") }))
sig.anova.car.cnscs.class <- colnames(anova.car.cnscs.ps)[colAnys(anova.car.cnscs.ps[c("Group", "Day:Group"), ] <= 0.05)]
sig.anova.car.cnscs.class <- intersect(sig.anova.car.cnscs.class, sig.anova.car.cnscs.pre.class)

##Human, metab concentration time course, only metabolites significant in control vs sepsis and with control overlapping S, statistical approach
#keep_set <- intersect(sig.anova.car.s.class, sig.anova.car.cvns.class)
keep_set <- intersect(colnames(anova.car.cvns.ps)[anova.car.cvns.ps["Survival", ] < 0.05], colnames(anova.car.s.ps)[anova.car.s.ps["Survival", ] < 0.05])
keep_set <- intersect(keep_set,colnames(anova.car.cvns.pre.ps)[anova.car.cvns.pre.ps["Survival", ] < 0.05])
keep_set <- intersect(keep_set,colnames(anova.car.s.pre.ps)[anova.car.s.pre.ps["Survival", ] < 0.05])
#hd <- human_data[, c(1:6, intersect(metab_sel, which(!colnames(human_data) %in% sig.anova.car.s.class)))]
hd <- human_data[, c(1:6, which(colnames(human_data) %in% keep_set))]
n_mets <- ncol(hd) - 6
ncols <- 5
p <- ggplot(data = subset(melt(hd, id.vars = 1:6), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Group, colour = Group)) + 
  facet_wrap(facets = ~ variable, ncol = ncols, nrow = ceiling(n_mets/ncols), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration, µM") +
  xlab("Day") +
  human_col_scale() +
  theme_bw()
ggsave(plot = p, filename = "human_metab_survival_nonsig_C_overlap_S_sig.png", path = out_dir, width = 12/4*ncols, height = 0.3 + 1.5 * ceiling(n_mets/ncols), units = "in")

#ANOVA for Ferrario data, sS vs sNS
val_anova <- t3ANOVA(data = human_sepsis_val_data, random = ~1|Patient, formula = concentration ~ Day + Day*Survival28 + Survival28, use.corAR = TRUE, col.set = colnames(human_sepsis_val_data[, c(-1:-4)]), id.vars = c("Survival28", "Day", "Patient"), control = lmeControl(msMaxIter = 100))
colnames(val_anova$ps) <- colnames(human_sepsis_val_data[, c(-1:-4)])
val_anova_ps <- val_anova$ps
val_anova_fdr <- apply(val_anova_ps, 2, p.adjust, method = "fdr")
sig.anova.car.val.s.pre.class <- colnames(val_anova$ps)[colAnys(val_anova_ps[3:4, ] < 0.05)]
sig.anova.car.val.s.class <- colnames(val_anova$ps)[colAnys(val_anova_fdr[3:4, ] < 0.05)]

#Save everything ANOVA to file
anova_complete_res <- grep(pattern = "anova.car", x = names(environment()), value = TRUE)
anova_complete_res <- grep(pattern = "models|res|terms|normality", x = anova_complete_res, invert = TRUE, value = TRUE)
save(list = anova_complete_res, file = paste0(out_dir, "ANOVA_complete_res.RData"))

#Build long format tables for plotting
##Scale measurement values at maximum concentration
human_sepsis_data_max_norm <- max_norm(x = human_sepsis_data, subset = -1:-6)
rat_sepsis_data_max_norm <- max_norm(x = rat_sepsis_data, subset = -1:-4)
##Melt scaled data into long form
human_sepsis_data_long_form <- melt(data = human_sepsis_data_max_norm, id.vars = c("Sample ID", "Patient", "Day", "Survival", "CAP / FP", "Group"))
rat_sepsis_data_long_form <- melt(data = rat_sepsis_data_max_norm, id.vars = c("Sample Identification", "material", "group", "time point"))

#Build tables for cluster-heatmaps
##Build covariance matrix with metabolite groups
human_sepsis_data_normal_grouped_metab_cov <- cbind(human_sepsis_data_normal_grouped[, 1:6], cov(t(human_sepsis_data_normal_grouped[,group_metab_sel]), use = "pairwise.complete.obs"))
human_sepsis_data_normal_grouped_pheno_cov <- cbind(subset(human_sepsis_data_normal_grouped, Day %in% tanova_day_set)[, 1:6], cov(t(subset(human_sepsis_data_normal_grouped, Day %in% tanova_day_set)[,group_pheno_sel]), use = "pairwise.complete.obs"))
human_sepsis_data_normal_grouped_pheno_cov <- human_sepsis_data_normal_grouped_pheno_cov[, !apply(human_sepsis_data_normal_grouped_pheno_cov, 2, function(x){all(is.na(x))})]
##Build covariance matrix with original metabolites
human_sepsis_data_normal_metab_cov <- cbind(human_sepsis_data_normal[, 1:6], cov(t(human_sepsis_data_normal[,metab_sel]), use = "pairwise.complete.obs"))
human_sepsis_data_normal_pheno_cov <- cbind(subset(human_sepsis_data_normal, Day %in% tanova_day_set)[, 1:6], cov(t(subset(human_sepsis_data_normal, Day %in% tanova_day_set)[,pheno_sel]), use = "pairwise.complete.obs"))
human_sepsis_data_normal_pheno_cov <- human_sepsis_data_normal_pheno_cov[, !apply(human_sepsis_data_normal_pheno_cov, 2, function(x){all(is.na(x))})]
##Build correlation matrix with metabolite groups
human_sepsis_data_normal_grouped_metab_cor <- cbind(human_sepsis_data_normal_grouped[, 1:6], cor(t(human_sepsis_data_normal_grouped[,group_metab_sel]), use = "pairwise.complete.obs"))
human_sepsis_data_normal_grouped_pheno_cor <- cbind(subset(human_sepsis_data_normal_grouped, Day %in% tanova_day_set)[, 1:6], cor(t(subset(human_sepsis_data_normal_grouped, Day %in% tanova_day_set)[, group_pheno_sel]), use = "pairwise.complete.obs"))
human_sepsis_data_normal_grouped_pheno_cor <- human_sepsis_data_normal_grouped_pheno_cor[, !apply(human_sepsis_data_normal_grouped_pheno_cor, 2, function(x){all(is.na(x))})]
##Build covariance matrix with original metabolites
human_sepsis_data_normal_metab_cor <- cbind(human_sepsis_data_normal[, 1:6], cor(t(human_sepsis_data_normal[, metab_sel]), use = "pairwise.complete.obs"))
human_sepsis_data_normal_pheno_cor <- cbind(subset(human_sepsis_data_normal, Day %in% tanova_day_set)[, 1:6], cor(t(subset(human_sepsis_data_normal, Day %in% tanova_day_set)[, pheno_sel]), use = "pairwise.complete.obs"))
human_sepsis_data_normal_pheno_cor <- human_sepsis_data_normal_pheno_cor[, !apply(human_sepsis_data_normal_pheno_cor, 2, function(x){all(is.na(x))})]
##Build covariance matrix with metabolite and phenom. var groups
human_sepsis_data_normal_grouped_metab_pheno_cov <- cbind(subset(human_sepsis_data_normal_grouped, Day %in% tanova_day_set)[, 1:6], cov(t(subset(human_sepsis_data_normal_grouped, Day %in% tanova_day_set)[, -1:-6]), use = "pairwise.complete.obs"))
##Build correlation matrix with metabolite and phenom. var groups
human_sepsis_data_normal_grouped_metab_pheno_cor <- cbind(subset(human_sepsis_data_normal_grouped, Day %in% tanova_day_set)[, 1:6], cor(t(subset(human_sepsis_data_normal_grouped, Day %in% tanova_day_set)[, -1:-6]), use = "pairwise.complete.obs"))
human_sepsis_data_normal_grouped_metab_pheno_cor <- human_sepsis_data_normal_grouped_metab_pheno_cor[, !apply(human_sepsis_data_normal_grouped_metab_pheno_cor, 2, function(x){all(is.na(x))})]
##Build covariance matrix with original metabolites and phenom. vars
human_sepsis_data_normal_metab_pheno_cov <- cbind(subset(human_sepsis_data_normal, Day %in% tanova_day_set)[, 1:6], cov(t(subset(human_sepsis_data_normal, Day %in% tanova_day_set)[, -1:-6]), use = "pairwise.complete.obs"))
##Build correlation matrix with original metabolites and phenom. vars
human_sepsis_data_normal_metab_pheno_cor <- cbind(subset(human_sepsis_data_normal, Day %in% tanova_day_set)[, 1:6], cor(t(subset(human_sepsis_data_normal, Day %in% tanova_day_set)[, -1:-6]), use = "pairwise.complete.obs"))
human_sepsis_data_normal_metab_pheno_cor <- human_sepsis_data_normal_metab_pheno_cor[, !apply(human_sepsis_data_normal_metab_pheno_cor, 2, function(x){all(is.na(x))})]

##Build covariance matrix of metabolites with original metabolites
###Survival-ignorant
human_sepsis_data_normal_conc_metab_cov <- cov(human_sepsis_data_normal[, metab_sel], use = "pairwise.complete.obs")
human_sepsis_data_normal_conc_pheno_cov <- cov(subset(human_sepsis_data_normal, Day %in% tanova_day_set)[, pheno_sel], use = "pairwise.complete.obs")
###Survival-regardent
human_sepsis_data_normal_NS_conc_metab_cov <- cov(subset(human_sepsis_data_normal, Survival == "NS", metab_sel), use = "pairwise.complete.obs")
human_sepsis_data_normal_S_conc_metab_cov <- cov(subset(human_sepsis_data_normal, Survival == "S", metab_sel), use = "pairwise.complete.obs")
human_sepsis_data_normal_NS_conc_pheno_cov <- cov(subset(human_sepsis_data_normal, Survival == "NS" & Day %in% tanova_day_set, pheno_sel), use = "pairwise.complete.obs")
human_sepsis_data_normal_S_conc_pheno_cov <- cov(subset(human_sepsis_data_normal, Survival == "S" & Day %in% tanova_day_set, pheno_sel), use = "pairwise.complete.obs")
##Build covariance matrix of metabolites with metabolites groups
###Survival-ignorant
human_sepsis_data_normal_grouped_conc_metab_cov <- cov(human_sepsis_data_normal_grouped[, group_metab_sel], use = "pairwise.complete.obs")
human_sepsis_data_normal_grouped_conc_pheno_cov <- cov(subset(human_sepsis_data_normal_grouped, Day %in% tanova_day_set)[, group_pheno_sel], use = "pairwise.complete.obs")
###Survival-regardent
human_sepsis_data_normal_NS_grouped_conc_metab_cov <- cov(subset(human_sepsis_data_normal_grouped, Survival == "NS", group_metab_sel), use = "pairwise.complete.obs")
human_sepsis_data_normal_S_grouped_conc_metab_cov <- cov(subset(human_sepsis_data_normal_grouped, Survival == "S", group_metab_sel), use = "pairwise.complete.obs")
human_sepsis_data_normal_NS_grouped_conc_pheno_cov <- cov(subset(human_sepsis_data_normal_grouped, Survival == "NS" & Day %in% tanova_day_set, group_pheno_sel), use = "pairwise.complete.obs")
human_sepsis_data_normal_S_grouped_conc_pheno_cov <- cov(subset(human_sepsis_data_normal_grouped, Survival == "S" & Day %in% tanova_day_set, group_pheno_sel), use = "pairwise.complete.obs")

##Build variance vectors of metabolites from normalized concentrations
###One method: variance of metabolite seperately for each patient
pat_first <- match(unique(human_sepsis_data_normal$Patient), human_sepsis_data_normal$Patient)
human_sepsis_data_normal_conc_var <- lapply(human_sepsis_data_normal[, -1:-6], aggregate, list(Patient = human_sepsis_data_normal$Patient), var, na.rm = T)
human_sepsis_data_normal_conc_var <- data.frame(Patient = human_sepsis_data_normal_conc_var[[1]]$Patient, 
                                                Survival = human_sepsis_data_normal$Survival[pat_first],
                                                lapply(human_sepsis_data_normal_conc_var, function(e) e$x))
human_sepsis_data_normal_grouped_conc_var <- lapply(human_sepsis_data_normal_grouped[, -1:-6], aggregate, list(Patient = human_sepsis_data_normal_grouped$Patient), var, na.rm = T)
human_sepsis_data_normal_grouped_conc_var <- data.frame(Patient = human_sepsis_data_normal_grouped_conc_var[[1]]$Patient, 
                                                        Survival = human_sepsis_data_normal_grouped$Survival[pat_first],
                                                        lapply(human_sepsis_data_normal_grouped_conc_var, function(e) e$x))
hsd <- subset(human_sepsis_data, Day %in% tanova_day_set)
human_sepsis_data_conc_var <- lapply(hsd[, -1:-6], aggregate, list(Patient = hsd$Patient), var, na.rm = T)
human_sepsis_data_conc_var <- data.frame(Patient = human_sepsis_data_conc_var[[1]]$Patient, 
                                                Survival = hsd$Survival[pat_first],
                                                lapply(human_sepsis_data_conc_var, function(e) e$x))
human_sepsis_data_conc_mean <- lapply(hsd[, -1:-6], aggregate, list(Patient = hsd$Patient), mean, na.rm = T)
human_sepsis_data_conc_mean <- data.frame(Patient = human_sepsis_data_conc_mean[[1]]$Patient, 
                                         Survival = hsd$Survival[pat_first],
                                         lapply(human_sepsis_data_conc_mean, function(e) e$x))

###Another method: variance of metabolite over all patients and days
human_sepsis_data_normal_patcentered <- human_sepsis_data_normal
human_sepsis_data_normal_patcentered[,-1:-6] <- Reduce("rbind", lapply(as.character(unique(human_sepsis_data_normal_patcentered$`Patient`)), function(x){ scale(subset(human_sepsis_data_normal_patcentered, Patient == x, -1:-6), center = T, scale = F) }))
human_sepsis_data_normal_NS_conc_metab_var <- colVars(as.matrix(subset(human_sepsis_data_normal_patcentered, Day %in% tanova_day_set & Survival == "NS", metab_sel)))
human_sepsis_data_normal_S_conc_metab_var <- colVars(as.matrix(subset(human_sepsis_data_normal_patcentered, Day %in% tanova_day_set & Survival == "S", metab_sel)))
human_sepsis_data_normal_NS_conc_pheno_var <- colVars(as.matrix(subset(human_sepsis_data_normal_patcentered, Day %in% tanova_day_set & Survival == "NS", pheno_sel)))
human_sepsis_data_normal_S_conc_pheno_var <- colVars(as.matrix(subset(human_sepsis_data_normal_patcentered, Day %in% tanova_day_set & Survival == "S", pheno_sel)))

###Another method: variance of metabolite over all patients seperately for each day
human_sepsis_data_normal_NS_conc_metab_day_var <- t(sapply(tanova_day_set, function(d){ colVars(as.matrix(subset(human_sepsis_data_normal_patcentered, Survival == "NS" & Day == d, metab_sel))) }))
human_sepsis_data_normal_S_conc_metab_day_var <- t(sapply(tanova_day_set, function(d){ colVars(as.matrix(subset(human_sepsis_data_normal_patcentered, Survival == "S" & Day == d, metab_sel))) }))
human_sepsis_data_normal_NS_conc_pheno_day_var <- t(sapply(tanova_day_set, function(d){ colVars(as.matrix(subset(human_sepsis_data_normal_patcentered, Survival == "NS" & Day == d, pheno_sel))) }))
human_sepsis_data_normal_S_conc_pheno_day_var <- t(sapply(tanova_day_set, function(d){ colVars(as.matrix(subset(human_sepsis_data_normal_patcentered, Survival == "S" & Day == d, pheno_sel))) }))

##Human, variance of pheno vars seperate by day
pheno_day_var_df <- rbind(cbind(data.frame(Day = tanova_day_set, Survival = "NS"), human_sepsis_data_normal_NS_conc_pheno_day_var), 
                          cbind(data.frame(Day = tanova_day_set, Survival = "S"), human_sepsis_data_normal_S_conc_pheno_day_var))
colnames(pheno_day_var_df) <- c("Day", "Survival", colnames(human_sepsis_data)[pheno_sel])
pheno_day_var_long_df <- melt(pheno_day_var_df, id.vars = c("Day", "Survival"))
pheno_day_var_long_df$group <- human_sepsis_legend$group[match(pheno_day_var_long_df$variable, human_sepsis_legend[[1]])]
metab_normal_day_var_df <- rbind(cbind(data.frame(Day = tanova_day_set, Survival = "NS"), human_sepsis_data_normal_NS_conc_metab_day_var), 
                          cbind(data.frame(Day = tanova_day_set, Survival = "S"), human_sepsis_data_normal_S_conc_metab_day_var))
colnames(metab_normal_day_var_df) <- c("Day", "Survival", colnames(human_sepsis_data)[metab_sel])
metab_day_var_long_df <- melt(metab_normal_day_var_df, id.vars = c("Day", "Survival"))
metab_day_var_long_df$group <- human_sepsis_legend$group[match(metab_day_var_long_df$variable, human_sepsis_legend[[1]])]

##Calculate post hoc p-values for variance time course
###Post hoc analysis of when the change happened
####Build models
fml <- value ~ DaySurv - 1
ftd_ph <- subset(metab_day_var_long_df, Day %in% tanova_day_set)
ftd_ph$DaySurv <- interaction(ftd_ph$Day, ftd_ph$Survival, drop = TRUE)
group_set <- unique(ftd_ph$group)
group_set <- group_set[group_set != "sugar"]
ftd_ph <- lapply(group_set, function(gr) subset(ftd_ph, group == gr))
names(ftd_ph) <- group_set
anova.car.ph.metab.variance.models <- mclapply(ftd_ph, function(dat) fit_lin_mod_lme(dvar = "value", data = dat, formula = fml, id.vars = c("Survival", "DaySurv"), random = ~1|Survival, use.corAR = TRUE, lmeControl(msMaxIter = 100)))
anova.car.ph.metab.variance.models <- mclapply(anova.car.ph.metab.variance.models, function(e) try(eval(e)))
####Construct constrast matrix
contr.m <- get_S_NS_C_contrmat(tanova_day_set)
daySurv <- interaction(factor(ftd_ph$Day), factor(ftd_ph$Survival), drop = TRUE)
rownames(contr.m) <- levels(daySurv)
contr.m <- contr.m[-1:-4, -1:-4]
####Actual post hoc test
anova.car.ph.metab.variance.res <- mclapply(anova.car.ph.metab.variance.models, function(e) summary(glht(e, linfct = mcp(DaySurv = t(contr.m)))))
anova.car.ph.metab.variance.ps <- sapply(anova.car.ph.metab.variance.res, function(e) e$test$pvalues)
anova.car.ph.metab.variance.sig.contr <- lapply(data.frame(anova.car.ph.metab.variance.ps), function(col) which(col <= 0.05))
names(anova.car.ph.metab.variance.ps) <- names(anova.car.ph.metab.variance.models)
names(anova.car.ph.metab.variance.sig.contr) <- names(anova.car.ph.metab.variance.models)

fml <- value ~ DaySurv - 1
ftd_ph <- subset(pheno_day_var_long_df, Day %in% tanova_day_set)
ftd_ph$DaySurv <- interaction(ftd_ph$Day, ftd_ph$Survival, drop = TRUE)
group_set <- unique(ftd_ph$group)
group_set <- group_set[group_set != "sugar"]
ftd_ph <- lapply(group_set, function(gr) subset(ftd_ph, group == gr))
names(ftd_ph) <- group_set
anova.car.ph.pheno.variance.models <- mclapply(ftd_ph, function(dat) fit_lin_mod_lme(dvar = "value", data = dat, formula = fml, id.vars = c("Survival", "DaySurv"), random = ~1|Survival, use.corAR = TRUE, lmeControl(msMaxIter = 100)))
anova.car.ph.pheno.variance.models <- mclapply(anova.car.ph.pheno.variance.models, function(e) try(eval(e)))
####Construct constrast matrix
contr.m <- get_S_NS_C_contrmat(tanova_day_set)
daySurv <- interaction(factor(ftd_ph$Day), factor(ftd_ph$Survival), drop = TRUE)
rownames(contr.m) <- levels(daySurv)
contr.m <- contr.m[-1:-4, -1:-4]
####Actual post hoc test
anova.car.ph.pheno.variance.res <- mclapply(anova.car.ph.pheno.variance.models, function(e) summary(glht(e, linfct = mcp(DaySurv = t(contr.m)))))
anova.car.ph.pheno.variance.ps <- sapply(anova.car.ph.pheno.variance.res, function(e) e$test$pvalues)
anova.car.ph.pheno.variance.sig.contr <- lapply(data.frame(anova.car.ph.pheno.variance.ps), function(col) which(col <= 0.05))
names(anova.car.ph.pheno.variance.ps) <- names(anova.car.ph.pheno.variance.models)
names(anova.car.ph.pheno.variance.sig.contr) <- names(anova.car.ph.pheno.variance.models)

##Find significantly correlating concentrations
if (file.exists("generate_stats_bootstrap1000_result.RData")){
  load("generate_stats_bootstrap1000_result.RData")
}else{
  human_sepsis_data_normal_conc_metab_corr <- list()
  cols_grouped_metab <- colnames(human_sepsis_data)[-1:-6]
  cols_metab <- colnames(human_sepsis_data)[-1:-6]
  bootstrap_repeats <- 1000
  tic()
  for (d in 0:2){
    corr_dat <- list()
    n_NS <- sum(human_sepsis_data$Survival == "NS" & human_sepsis_data$Day == d)
    n_S <- sum(human_sepsis_data$Survival == "S" & human_sepsis_data$Day == d)
    NS_corr <- my.corr.test(x = subset(human_sepsis_data_normal, Survival == "NS" & Day == d, select = -1:-6), adjust = "fdr")
    NS_grouped_corr <- my.corr.test(x = subset(human_sepsis_data_normal_grouped, Survival == "NS" & Day == d, select = -1:-6), adjust = "fdr")
    ##Get significant metabolite pairs in nonsurvivors
    xy <- which.xy(NS_grouped_corr$p <= 0.5)
    xy <- subset(xy, x < y) # x is row, y is column, so x < y means "upper triangle only"
    if (nrow(xy) > 0) {
      coeffs <- apply(xy, c(1), function(cc){ NS_grouped_corr$r[cc[1], cc[2]] })
      cdir <- ifelse(coeffs > 0, "+", "-")
      corr_dat[["NS_grouped_sig_pairs"]] <- paste0(cols_grouped_metab[xy$x], " ~(", cdir,") ", cols_grouped_metab[xy$y])
      corr_dat[["NS_grouped_sig_coeff"]] <- coeffs
    }
    else {
      corr_dat[["NS_grouped_sig_pairs"]] <- c()
      corr_dat[["NS_grouped_sig_coeff"]] <- c()
    }
    xy <- which.xy(NS_corr$p <= 0.05)
    xy <- subset(xy, x < y)
    if (nrow(xy) > 0){
      coeffs <- apply(xy, c(1), function(cc){ NS_corr$r[cc[1], cc[2]] })
      cdir <- ifelse(coeffs > 0, "+", "-")
      corr_dat[["NS_sig_pairs"]] <- paste0(cols_metab[xy$x], " ~(", cdir, ") ", cols_metab[xy$y])
      corr_dat[["NS_sig_coeff"]] <- coeffs
    }
    else {
      corr_dat[["NS_sig_pairs"]] <- c()
      corr_dat[["NS_sig_coeff"]] <- c()
    }
    
    #bootstrap_list <- list()
    bootstrap_list <- mclapply(X = 1:bootstrap_repeats, FUN = bootstrap_S_corr_fun, n_S = n_S, corr_dat = corr_dat)
    human_sepsis_data_normal_conc_metab_corr[[paste0("Day",d)]] <- bootstrap_list
  }
  toc()
  save(human_sepsis_data_normal_conc_metab_corr, file = "generate_stats_bootstrap1000_result.RData")
}
S_sig_pairs_day0_cutoff <- mean(unlist(lapply(human_sepsis_data_normal_conc_metab_corr$Day0, function(x){ length(x$S_sig_pairs) })))
S_sig_pairs_day1_cutoff <- mean(unlist(lapply(human_sepsis_data_normal_conc_metab_corr$Day1, function(x){ length(x$S_sig_pairs) })))
S_sig_pairs_day2_cutoff <- mean(unlist(lapply(human_sepsis_data_normal_conc_metab_corr$Day2, function(x){ length(x$S_sig_pairs) })))
S_sig_pairs_day0_intersect_cutoff <- mean(unlist(lapply(human_sepsis_data_normal_conc_metab_corr$Day0, function(x){ length(intersect(x$S_sig_pairs, x$NS_sig_pairs)) })))
S_sig_pairs_day1_intersect_cutoff <- mean(unlist(lapply(human_sepsis_data_normal_conc_metab_corr$Day1, function(x){ length(intersect(x$S_sig_pairs, x$NS_sig_pairs)) })))
S_sig_pairs_day2_intersect_cutoff <- mean(unlist(lapply(human_sepsis_data_normal_conc_metab_corr$Day2, function(x){ length(intersect(x$S_sig_pairs, x$NS_sig_pairs)) })))
S_sig_pairs_day0_tab <- table(unlist(lapply(human_sepsis_data_normal_conc_metab_corr$Day0, function(x){ x$S_sig_pairs })))
S_sig_pairs_day1_tab <- table(unlist(lapply(human_sepsis_data_normal_conc_metab_corr$Day1, function(x){ x$S_sig_pairs })))
S_sig_pairs_day2_tab <- table(unlist(lapply(human_sepsis_data_normal_conc_metab_corr$Day2, function(x){ x$S_sig_pairs })))
NS_sig_pairs_day0_tab <- table(unlist(lapply(human_sepsis_data_normal_conc_metab_corr$Day0, function(x){ x$NS_sig_pairs })))
NS_sig_pairs_day1_tab <- table(unlist(lapply(human_sepsis_data_normal_conc_metab_corr$Day1, function(x){ x$NS_sig_pairs })))
NS_sig_pairs_day2_tab <- table(unlist(lapply(human_sepsis_data_normal_conc_metab_corr$Day2, function(x){ x$NS_sig_pairs })))

##Build inverse covariance matrices
inv_human_sepsis_data_normal_conc_metab_cov <- pcor.shrink(na.omit(human_sepsis_data_normal[, metab_sel]))
class(inv_human_sepsis_data_normal_conc_metab_cov) <- "matrix"
inv_human_sepsis_data_normal_conc_pheno_cov <- pcor.shrink(na.omit(human_sepsis_data_normal[, pheno_sel]))
class(inv_human_sepsis_data_normal_conc_pheno_cov) <- "matrix"
inv_human_sepsis_data_normal_NS_conc_metab_cov <- pcor.shrink(na.omit(subset(human_sepsis_data_normal, Survival == "NS", metab_sel)))
class(inv_human_sepsis_data_normal_NS_conc_metab_cov) <- "matrix"
inv_human_sepsis_data_normal_S_conc_metab_cov <- pcor.shrink(na.omit(subset(human_sepsis_data_normal, Survival == "S", metab_sel)))
class(inv_human_sepsis_data_normal_S_conc_metab_cov) <- "matrix"
inv_human_sepsis_data_normal_NS_conc_pheno_cov <- pcor.shrink(na.omit(subset(human_sepsis_data_normal, Survival == "NS", pheno_sel)))
class(inv_human_sepsis_data_normal_NS_conc_pheno_cov) <- "matrix"
inv_human_sepsis_data_normal_S_conc_pheno_cov <- pcor.shrink(na.omit(subset(human_sepsis_data_normal, Survival == "S", pheno_sel)))
class(inv_human_sepsis_data_normal_S_conc_pheno_cov) <- "matrix"

####################
#CSV output
####################

#human significant metabolites per metabolite group, including pehno vars as one group
table_groups <- get_human_sepsis_legend()
sigs <- c(sig.anova.car.s.class, sig.anova.car.s.pheno.class)
groups <- table_groups[match(sigs, table_groups[, 1]), 3]
s_tab_list <- tapply(X = sigs, INDEX = list(Group = factor(groups)), FUN = c)
max_l <- max(sapply(s_tab_list, length))
s_tab_list <- lapply(s_tab_list, function(l) c(l, rep("", max_l - length(l))))
s_tab <- as.data.frame(s_tab_list)
s_tab <- s_tab[, match(make.names(unique(groups)), colnames(s_tab))]
colnames(s_tab) <- unique(groups)
write.table(x = s_tab, file = paste0(out_dir, "table_Survival_anova_sig_features.csv"), sep = "\t", quote = FALSE, row.names = FALSE)

#human significant metabolites per metabolite group, including pehno vars as one group
table_groups <- get_human_sepsis_legend()
sigs <- sig.anova.car.c.class
groups <- table_groups[match(sigs, table_groups[, 1]), 3]
c_tab_list <- tapply(X = sigs, INDEX = list(Group = factor(groups)), FUN = c)
max_l <- max(sapply(c_tab_list, length))
c_tab_list <- lapply(c_tab_list, function(l) c(l, rep("", max_l - length(l))))
c_tab <- as.data.frame(c_tab_list)
c_tab <- c_tab[, match(make.names(unique(groups)), colnames(c_tab))]
colnames(c_tab) <- unique(groups)
write.table(x = c_tab, file = paste0(out_dir, "table_control_anova_sig_features.csv"), sep = "\t", quote = FALSE, row.names = FALSE)

#human fold changes at day 0, septic-NS vs septic-S
SS <- subset(human_data, Day == 0 & Group == "Septic-S", sig.anova.car.s.class)
NS <- subset(human_data, Day == 0 & Group == "Septic-NS", sig.anova.car.s.class)
fold_change <- log2(colMeans(SS) / colMeans(NS))
sink(file = paste0(out_dir, "fold_change_geq_1_day0_sS_vs_sNS.txt"))
print(fold_change[which(abs(fold_change) > 1)])
sink()

#also for all metabolites, day 0, Septic-NS vs Septic-S
SS <- subset(human_data, Day == 0 & Group == "Septic-S", metab_sel)
NS <- subset(human_data, Day == 0 & Group == "Septic-NS", metab_sel)
fold_change <- log2(colMeans(SS) / colMeans(NS))
sink(file = paste0(out_dir, "fold_change_day0_sS_vs_sNS.txt"))
print(fold_change)
sink()

#human fold changes at day 0, septic-NS vs non-Septic-NS
sNS <- subset(human_data, Day == 0 & Group == "Septic-NS", sig.anova.car.nssepnon.class)
nNS <- subset(human_data, Day == 0 & Group == "non-Septic-NS", sig.anova.car.nssepnon.class)
fold_change <- log2(colMeans(sNS) / colMeans(nNS))
sink(file = paste0(out_dir, "fold_change_day0_sNS_vs_nNS.txt"))
print(fold_change)#[which(abs(fold_change) > 1)])
sink()

##XLSX table with all p and q values
anova.ptabc <- list()
anova.ptabc[[1]] <- as.data.frame(t(anova.car.c.pre.ps))
colnames(anova.ptabc[[1]]) <- paste0("p of ", colnames(anova.ptabc[[1]]))
anova.ptabs <- list()
anova.ptabs[[1]] <- as.data.frame(t(anova.car.s.pre.ps))
colnames(anova.ptabs[[1]]) <- paste0("p of ", colnames(anova.ptabs[[1]]))
anova.ptabs[[length(anova.ptabs) + 1]] <- as.data.frame(t(anova.car.s.pheno.pre.ps))
colnames(anova.ptabs[[length(anova.ptabs)]]) <- paste0("p of ", colnames(anova.ptabs[[length(anova.ptabs)]]))
anova.qtabc <- list()
anova.qtabc[[1]] <- as.data.frame(t(anova.car.c.ps))
colnames(anova.qtabc[[1]]) <- paste0("q of ", colnames(anova.qtabc[[1]]))
anova.qtabs <- list()
anova.qtabs[[1]] <- as.data.frame(t(anova.car.s.ps))
colnames(anova.qtabs[[1]]) <- paste0("q of ", colnames(anova.qtabs[[1]]))
anova.qtabs[[length(anova.qtabs) + 1]] <- as.data.frame(t(anova.car.s.pheno.ps))
colnames(anova.qtabs[[length(anova.qtabs)]]) <- paste0("q of ", colnames(anova.qtabs[[length(anova.qtabs)]]))
anova.pqtabc <- lapply(seq_along(anova.ptabc), function(n) cbind(anova.ptabc[[n]], anova.qtabc[[n]]))
anova.pqtabs <- lapply(seq_along(anova.ptabs), function(n) cbind(anova.ptabs[[n]], anova.qtabs[[n]]))
anova.pqtaball <- c(anova.pqtabc, anova.pqtabs)
sheetnames <- c("metab.", "biochem.")
sheetnames <- paste0(rep(sheetnames, times = 2), ", ", rep(c("non-Septic vs Septic", "Septic-S vs Septic-NS"), each = 2))
sheetnames <- sheetnames[-2]
pqtab_xlsx_file <- paste0(out_dir, "patient_ANOVA_pq_first_suppl.xlsx")
pqtab_workbook <- createWorkbook(type = "xlsx")
for (n in seq_along(anova.pqtaball)){
  sheet <- createSheet(wb = pqtab_workbook, sheetName = sheetnames[n])
  for (col in 1:10){
    setColumnWidth(sheet = sheet, colIndex = col, colWidth = 20)
  }
  sheet <- addDataFrame(x = anova.pqtaball[[n]], 
                        sheet = sheet, 
                        col.names = TRUE, 
                        row.names = TRUE, 
                        colnamesStyle = CellStyle(wb = pqtab_workbook, 
                                                  font = Font(wb = pqtab_workbook, isBold = TRUE), 
                                                  alignment = Alignment(horizontal = "ALIGN_CENTER")))
}
sheet <- createSheet(wb = pqtab_workbook, sheetName = "legend")
for (col in 1:10){
  setColumnWidth(sheet = sheet, colIndex = col, colWidth = 30)
}
legend_df <- read.xlsx(file = "../../data/measurements/Summary human sample data.xlsx", sheetIndex = 2)
colnames(legend_df)[1] <- ""
sheet <- addDataFrame(x = legend_df,
                      sheet = sheet, 
                      row.names = FALSE,
                      colnamesStyle = CellStyle(wb = pqtab_workbook, 
                                                font = Font(wb = pqtab_workbook, isBold = TRUE), 
                                                alignment = Alignment(horizontal = "ALIGN_CENTER")))
saveWorkbook(wb = pqtab_workbook, file = pqtab_xlsx_file)

##Another XLSX table with all the ANOVA p- and q-values but for non-Septic-NS vs non-Septic-S and Septic-NS
anova.ptabc <- list()
anova.ptabc[[1]] <- as.data.frame(t(anova.car.cnscs.pre.ps))
colnames(anova.ptabc[[1]]) <- paste0("p of ", colnames(anova.ptabc[[1]]))
anova.ptabs <- list()
anova.ptabs[[1]] <- as.data.frame(t(anova.car.nssepnon.pre.ps))
colnames(anova.ptabs[[1]]) <- paste0("p of ", colnames(anova.ptabs[[1]]))
anova.qtabc <- list()
anova.qtabc[[1]] <- as.data.frame(t(anova.car.cnscs.pre.ps))
colnames(anova.qtabc[[1]]) <- paste0("q of ", colnames(anova.qtabc[[1]]))
anova.qtabs <- list()
anova.qtabs[[1]] <- as.data.frame(t(anova.car.nssepnon.ps))
colnames(anova.qtabs[[1]]) <- paste0("q of ", colnames(anova.qtabs[[1]]))
anova.pqtabc <- lapply(seq_along(anova.ptabc), function(n) cbind(anova.ptabc[[n]], anova.qtabc[[n]]))
anova.pqtabs <- lapply(seq_along(anova.ptabs), function(n) cbind(anova.ptabs[[n]], anova.qtabs[[n]]))
anova.pqtaball <- c(anova.pqtabc, anova.pqtabs)
sheetnames <- c("metab.")
sheetnames <- paste0(rep(sheetnames, times = 2), ", ", c("non-Septic-S vs non-Septic-NS", "non-Septic-NS vs Septic-NS"))
pqtab_xlsx_file <- paste0(out_dir, "patient_ANOVA_pq_second_suppl.xlsx")
pqtab_workbook <- createWorkbook(type = "xlsx")
for (n in seq_along(anova.pqtaball)){
  sheet <- createSheet(wb = pqtab_workbook, sheetName = sheetnames[n])
  for (col in 1:10){
    setColumnWidth(sheet = sheet, colIndex = col, colWidth = 20)
  }
  sheet <- addDataFrame(x = anova.pqtaball[[n]], 
                        sheet = sheet, 
                        col.names = TRUE, 
                        row.names = TRUE, 
                        colnamesStyle = CellStyle(wb = pqtab_workbook, 
                                                  font = Font(wb = pqtab_workbook, isBold = TRUE), 
                                                  alignment = Alignment(horizontal = "ALIGN_CENTER")))
}
sheet <- createSheet(wb = pqtab_workbook, sheetName = "legend")
for (col in 1:10){
  setColumnWidth(sheet = sheet, colIndex = col, colWidth = 30)
}
legend_df <- read.xlsx(file = "../../data/measurements/Summary human sample data.xlsx", sheetIndex = 2)
colnames(legend_df)[1] <- ""
sheet <- addDataFrame(x = legend_df,
                      sheet = sheet, 
                      row.names = FALSE,
                      colnamesStyle = CellStyle(wb = pqtab_workbook, 
                                                font = Font(wb = pqtab_workbook, isBold = TRUE), 
                                                alignment = Alignment(horizontal = "ALIGN_CENTER")))
saveWorkbook(wb = pqtab_workbook, file = pqtab_xlsx_file)


####################
#Plot data
####################

##Human, graphs of significant correlations; NS, S and overlap
###Build data structures
S_sig_pairs_day0_tab <- sort(S_sig_pairs_day0_tab, decreasing = TRUE)
S_sig_pairs_day1_tab <- sort(S_sig_pairs_day1_tab, decreasing = TRUE)
S_sig_pairs_day2_tab <- sort(S_sig_pairs_day2_tab, decreasing = TRUE)
S_sig_pairs_day0_tab <- S_sig_pairs_day0_tab[1:ceiling(S_sig_pairs_day0_cutoff)]
S_sig_pairs_day1_tab <- S_sig_pairs_day1_tab[1:ceiling(S_sig_pairs_day1_cutoff)]
S_sig_pairs_day2_tab <- S_sig_pairs_day2_tab[1:ceiling(S_sig_pairs_day2_cutoff)]

split_pat <- " ~\\(.\\) "

S_sig_pairs_day0_mets <- strsplit(x = names(S_sig_pairs_day0_tab), split = split_pat)
S_sig_pairs_day1_mets <- strsplit(x = names(S_sig_pairs_day1_tab), split = split_pat)
S_sig_pairs_day2_mets <- strsplit(x = names(S_sig_pairs_day2_tab), split = split_pat)
NS_sig_pairs_day0_mets <- strsplit(x = names(NS_sig_pairs_day0_tab), split = split_pat)
NS_sig_pairs_day1_mets <- strsplit(x = names(NS_sig_pairs_day1_tab), split = split_pat)
# NS_sig_pairs_day2_mets <- strsplit(x = names(NS_sig_pairs_day2_tab), split = split_pat) # no sig pairs at day 2

S_d0_from <- unlist(lapply(S_sig_pairs_day0_mets, function(x) x[1]))
NS_d0_from <- unlist(lapply(NS_sig_pairs_day0_mets, function(x) x[1]))
S_d0_to <- unlist(lapply(S_sig_pairs_day0_mets, function(x) x[2]))
NS_d0_to <- unlist(lapply(NS_sig_pairs_day0_mets, function(x) x[2]))
S_d1_from <- unlist(lapply(S_sig_pairs_day1_mets, function(x) x[1]))
NS_d1_from <- unlist(lapply(NS_sig_pairs_day1_mets, function(x) x[1]))
S_d1_to <- unlist(lapply(S_sig_pairs_day1_mets, function(x) x[2]))
NS_d1_to <- unlist(lapply(NS_sig_pairs_day1_mets, function(x) x[2]))
S_d2_from <- unlist(lapply(S_sig_pairs_day2_mets, function(x) x[1]))
#NS_d2_from <- unlist(lapply(NS_sig_pairs_day2_mets, function(x) x[1]))
S_d2_to <- unlist(lapply(S_sig_pairs_day2_mets, function(x) x[2]))
#NS_d2_to <- unlist(lapply(NS_sig_pairs_day2_mets, function(x) x[2]))

S_d0_sign <- grepl(x = names(S_sig_pairs_day0_tab), pattern = "(+)", fixed = TRUE)
NS_d0_sign <- grepl(x = names(NS_sig_pairs_day0_tab), pattern = "(+)", fixed = TRUE)
S_d1_sign <- grepl(x = names(S_sig_pairs_day1_tab), pattern = "(+)", fixed = TRUE)
NS_d1_sign <- grepl(x = names(NS_sig_pairs_day1_tab), pattern = "(+)", fixed = TRUE)
S_d2_sign <- grepl(x = names(S_sig_pairs_day2_tab), pattern = "(+)", fixed = TRUE)
NS_d2_sign <- grepl(x = names(NS_sig_pairs_day2_tab), pattern = "(+)", fixed = TRUE)

##Human, graph of significant concentration correlations at day 0
graph_day0_nodes <- data.frame(metabolite = human_sepsis_legend[,1], connectivity = 0, class = "", stringsAsFactors = FALSE)
combined_mets <- c(unlist(S_sig_pairs_day0_mets), unlist(NS_sig_pairs_day0_mets))
cmt <- table(combined_mets)
graph_day0_nodes$connectivity[match(names(cmt), graph_day0_nodes$metabolite)] <- cmt
graph_day0_nodes <- subset(graph_day0_nodes, connectivity > 0)

from <- c(S_d0_from, NS_d0_from)
to <- c(S_d0_to, NS_d0_to)
graph_day0_edges <- data.frame(from = from, to = to, class = "", stringsAsFactors = FALSE)
graph_day0_edges <- subset(graph_day0_edges, from %in% graph_day0_nodes$metabolite & to %in% graph_day0_nodes$metabolite)
graph_day0_edges$class[1:length(S_d0_from)] <- "Survivor"
graph_day0_edges$class[(length(S_d0_from)+1):nrow(graph_day0_edges)] <- "Nonsurvivor"
graph_day0_edges$class[duplicated(graph_day0_edges[c("from", "to")])] <- "Overlap"
graph_day0_edges$class[duplicated(graph_day0_edges[c("from", "to")], fromLast = TRUE)] <- "Overlap"
graph_day0_edges$lty <- c(1 + as.numeric(grepl(x = names(S_sig_pairs_day0_tab), pattern = "(-)", fixed = TRUE)), 1 + as.numeric(grepl(x = names(NS_sig_pairs_day0_tab), pattern = "(-)", fixed = TRUE)))

igraph_day0_df <- graph.data.frame(graph_day0_edges, directed = FALSE, vertices = graph_day0_nodes)
colorset <- colorspace::rainbow_hcl(3)
E(igraph_day0_df)$color <- colorset[as.numeric(as.factor(E(igraph_day0_df)$class))]
E(igraph_day0_df)$width <- 2
V(igraph_day0_df)$color <- "grey"
V(igraph_day0_df)$size <- 1.2 * sqrt(V(igraph_day0_df)$connectivity)
V(igraph_day0_df)$label.cex <- 0.8

png(filename = paste0(out_dir, "metab_sig_corr_graph_day0.png"), width = 1024, height = 1024, units = "px")
plot(igraph_day0_df, ylim = c(-1.2, 1))
legend(x=0, y=-1.1, c("Nonsurvivor", "Overlap", "Survivor"), pch=22, col="#777777", pt.bg=colorset, pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()

##Human, graph of significant overlapping concentration correlations at day 0
graph_shared_day0_nodes <- data.frame(metabolite = human_sepsis_legend[,1], connectivity = 0, class = "", stringsAsFactors = FALSE)
shared_sig_pairs_day0_sel <- paste0(S_d0_from, S_d0_to) %in% paste0(NS_d0_from, NS_d0_to)
shared_sig_pairs_day0_mets <- c(S_d0_from[shared_sig_pairs_day0_sel], S_d0_to[shared_sig_pairs_day0_sel])
shared_sig_pairs_day0 <- names(S_sig_pairs_day0_tab)[shared_sig_pairs_day0_sel]
cmt <- table(shared_sig_pairs_day0_mets)
graph_shared_day0_nodes$connectivity[match(names(cmt), graph_shared_day0_nodes$metabolite)] <- cmt
graph_shared_day0_nodes <- subset(graph_shared_day0_nodes, connectivity > 0)

graph_shared_day0_edges <- data.frame(from = S_d0_from[shared_sig_pairs_day0_sel], to = S_d0_to[shared_sig_pairs_day0_sel], class = "", stringsAsFactors = FALSE)
graph_shared_day0_edges$class <- "Overlap"
graph_shared_day0_edges$lty <- 1 + as.numeric(grepl(x = shared_sig_pairs_day0, pattern = "-"))

igraph_shared_day0_df <- graph.data.frame(d = graph_shared_day0_edges, directed = FALSE, vertices = graph_shared_day0_nodes)
colorset <- colorspace::rainbow_hcl(3)
edge_type <- shared_sig_pairs_day0 %in% names(S_sig_pairs_day1_tab) + 3 * shared_sig_pairs_day0 %in% names(NS_sig_pairs_day1_tab)
edge_type[edge_type == 0] <- 2
E(igraph_shared_day0_df)$color <- colorset[edge_type]
E(igraph_shared_day0_df)$width <- 4
legend_groups <- human_sepsis_legend
legend_groups$group[pheno_sel - 6] <- "phenomenological variable"
legend_groups <- legend_groups[match(V(igraph_shared_day0_df)$name, legend_groups[,1]), ]
vertex_classes <- unique(legend_groups$group)
vertex_membership <- as.numeric(as.factor(legend_groups$group))
num_groups <- max(vertex_membership)
V(igraph_shared_day0_df)$color <- brewer.pal(num_groups, "Set3")[vertex_membership]
V(igraph_shared_day0_df)$frame.color <- "#ffffff"

V(igraph_shared_day0_df)$size <- 1.2 * sqrt(V(igraph_shared_day0_df)$connectivity) + 2
V(igraph_shared_day0_df)$label.cex <- 1

png(filename = paste0(out_dir, "metab_sig_corr_graph_shared_day0.png"), width = 1024, height = 1024, units = "px")
plot(igraph_shared_day0_df, ylim = c(-1.2, 1), main = "Significantly correlating metabolite concentrations at day 0")
legend(x=0.5, y=-1.1, c("Pair present in day 1 survivors", "Not present at day 1", "Pair present in day 1 nonsurvivors"), pch=22, col="#777777", pt.bg=colorset, pt.cex=2, cex=1, bty="n", ncol=1)
legend(x=-0.5, y=-1.1, vertex_classes, pch=21, col="#777777", pt.bg=brewer.pal(num_groups, "Set3")[as.numeric(as.factor(vertex_classes))], pt.cex=2, cex=1, bty="n", ncol=1)
dev.off()

##Human, graph of significant concentration correlations at day 1
graph_day1_nodes <- data.frame(metabolite = human_sepsis_legend[,1], connectivity = 0, class = "", stringsAsFactors = FALSE)
combined_mets <- c(unlist(S_sig_pairs_day1_mets), unlist(NS_sig_pairs_day1_mets))
cmt <- table(combined_mets)
graph_day1_nodes$connectivity[match(names(cmt), graph_day1_nodes$metabolite)] <- cmt
graph_day1_nodes <- subset(graph_day1_nodes, connectivity > 0)

from <- unlist(lapply(c(S_sig_pairs_day1_mets, NS_sig_pairs_day1_mets), function(x) x[1]))
to <- unlist(lapply(c(S_sig_pairs_day1_mets, NS_sig_pairs_day1_mets), function(x) x[2]))
graph_day1_edges <- data.frame(from = from, to = to, class = "", stringsAsFactors = FALSE)
graph_day1_edges$class[1:length(S_sig_pairs_day1_mets)] <- "Survivor"
graph_day1_edges$class[(length(S_sig_pairs_day1_mets)+1):nrow(graph_day1_edges)] <- "Nonsurvivor"

igraph_day1_df <- graph.data.frame(d = graph_day1_edges, directed = FALSE, vertices = graph_day1_nodes)
colorset <- colorspace::rainbow_hcl(3)
edge_type <- (graph_day1_edges$class == "Survivor") + 3 * (graph_day1_edges$class == "Nonsurvivor")
edge_type[edge_type == 0] <- 2
E(igraph_day1_df)$color <- colorset[edge_type]
E(igraph_day1_df)$width <- 4
legend_groups <- human_sepsis_legend
legend_groups$group[pheno_sel - 6] <- "phenomenological variable"
legend_groups <- legend_groups[match(V(igraph_day1_df)$name, legend_groups[,1]), ]
vertex_classes <- unique(legend_groups$group)
vertex_membership <- as.numeric(as.factor(legend_groups$group))
num_groups <- max(vertex_membership)
V(igraph_day1_df)$color <- brewer.pal(num_groups, "Set3")[vertex_membership]
V(igraph_day1_df)$frame.color <- "#ffffff"
V(igraph_day1_df)$size <- 1.2 * sqrt(V(igraph_day1_df)$connectivity) + 2
V(igraph_day1_df)$label.cex <- 1
#igraph_day1_community <- create.communities(as.numeric(as.factor(human_sepsis_legend$group[match(V(igraph_day1_df)$name, human_sepsis_legend[,1])])))

png(filename = paste0(out_dir, "metab_sig_corr_graph_day1.png"), width = 1024, height = 1024, units = "px")
plot(igraph_day1_df, ylim = c(-1.2, 1), main = "Significantly correlating metabolite concentrations at day 1")
legend(x=0.5, y=-1.1, c("Pair present in day 1 survivors", "Pair present in day 1 nonsurvivors"), pch=22, col="#777777", pt.bg=colorset[-2], pt.cex=2, cex=1, bty="n", ncol=1)
legend(x=-0.5, y=-1.1, vertex_classes, pch=21, col="#777777", pt.bg=brewer.pal(num_groups, "Set3")[as.numeric(as.factor(vertex_classes))], pt.cex=2, cex=1, bty="n", ncol=1)
dev.off()

##Human, data for table with number of metabolite correlations per day
col_metab <- colnames(human_sepsis_data_normal)
S_d0_sig_metab_pheno <- S_d0_from %in% col_metab[metab_sel] & S_d0_to %in% col_metab[pheno_sel]
NS_d0_sig_metab_pheno <- NS_d0_from %in% col_metab[metab_sel] & NS_d0_to %in% col_metab[pheno_sel]
Overlap_d0_sig_metab_pheno <- paste0(S_d0_from, S_d0_to)[S_d0_sig_metab_pheno] %in% paste0(NS_d0_from, NS_d0_to)[NS_d0_sig_metab_pheno]
sum(S_d0_sig_metab_pheno)
sum(NS_d0_sig_metab_pheno)
sum(Overlap_d0_sig_metab_pheno)
S_d1_sig_metab_pheno <- S_d1_from %in% col_metab[metab_sel] & S_d1_to %in% col_metab[pheno_sel]
NS_d1_sig_metab_pheno <- NS_d1_from %in% col_metab[metab_sel] & NS_d1_to %in% col_metab[pheno_sel]
Overlap_d1_sig_metab_pheno <- paste0(S_d1_from, S_d1_to)[S_d1_sig_metab_pheno] %in% paste0(NS_d1_from, NS_d1_to)[NS_d1_sig_metab_pheno]
sum(S_d1_sig_metab_pheno)
sum(NS_d1_sig_metab_pheno)
sum(Overlap_d1_sig_metab_pheno)
S_d2_sig_metab_pheno <- S_d2_from %in% col_metab[metab_sel] & S_d2_to %in% col_metab[pheno_sel]
sum(S_d2_sig_metab_pheno)

S_d0_sig_metab_metab <- S_d0_from %in% col_metab[metab_sel] & S_d0_to %in% col_metab[metab_sel]
NS_d0_sig_metab_metab <- NS_d0_from %in% col_metab[metab_sel] & NS_d0_to %in% col_metab[metab_sel]
Overlap_d0_sig_metab_metab <- paste0(S_d0_from, S_d0_to)[S_d0_sig_metab_metab] %in% paste0(NS_d0_from, NS_d0_to)[NS_d0_sig_metab_metab]
sum(S_d0_sig_metab_metab)
sum(NS_d0_sig_metab_metab)
sum(Overlap_d0_sig_metab_metab)
S_d1_sig_metab_metab <- S_d1_from %in% col_metab[metab_sel] & S_d1_to %in% col_metab[metab_sel]
NS_d1_sig_metab_metab <- NS_d1_from %in% col_metab[metab_sel] & NS_d1_to %in% col_metab[metab_sel]
Overlap_d1_sig_metab_metab <- paste0(S_d1_from, S_d1_to)[S_d1_sig_metab_metab] %in% paste0(NS_d1_from, NS_d1_to)[NS_d1_sig_metab_metab]
sum(S_d1_sig_metab_metab)
sum(NS_d1_sig_metab_metab)
sum(Overlap_d1_sig_metab_metab)
S_d2_sig_metab_metab <- S_d2_from %in% col_metab[metab_sel] & S_d2_to %in% col_metab[metab_sel]
sum(S_d2_sig_metab_metab)

##Human, connectivity statistics
mean(table(c(S_d0_from, S_d0_to)))
mean(table(c(NS_d0_from, NS_d0_to)))

mean(table(c(S_d1_from, S_d1_to)))
mean(table(c(NS_d1_from, NS_d1_to)))

mean(table(c(S_d2_from, S_d2_to)))

S_d0_igraph <- graph.data.frame(d = data.frame(from = S_d0_from, to = S_d0_to, stringsAsFactors = TRUE), directed = FALSE, vertices = union(S_d0_from, S_d0_to))
V(S_d0_igraph)$size <- 8
S_d0_community <- cluster_fast_greedy(S_d0_igraph)
sum(crossing(S_d0_community, S_d0_igraph))
modularity(S_d0_community)
mean(sizes(S_d0_community))
length(S_d0_community)

NS_d0_igraph <- graph.data.frame(d = data.frame(from = NS_d0_from, to = NS_d0_to, stringsAsFactors = TRUE), directed = FALSE, vertices = union(NS_d0_from, NS_d0_to))
V(NS_d0_igraph)$size <- 8
NS_d0_community <- cluster_fast_greedy(NS_d0_igraph)
sum(crossing(NS_d0_community, NS_d0_igraph))
modularity(NS_d0_community)
mean(sizes(NS_d0_community))
length(NS_d0_community)

S_d1_igraph <- graph.data.frame(d = data.frame(from = S_d1_from, to = S_d1_to, stringsAsFactors = TRUE), directed = FALSE, vertices = union(S_d1_from, S_d1_to))
V(S_d1_igraph)$size <- 8
S_d1_community <- cluster_fast_greedy(S_d1_igraph)
sum(crossing(S_d1_community, S_d1_igraph))
modularity(S_d1_community)
mean(sizes(S_d1_community))
length(S_d1_community)

NS_d1_igraph <- graph.data.frame(d = data.frame(from = NS_d1_from, to = NS_d1_to, stringsAsFactors = TRUE), directed = FALSE, vertices = union(NS_d1_from, NS_d1_to))
V(NS_d1_igraph)$size <- 8
NS_d1_community <- cluster_fast_greedy(NS_d1_igraph)
sum(crossing(NS_d1_community, NS_d1_igraph))
modularity(NS_d1_community)
mean(sizes(NS_d1_community))
length(NS_d1_community)

S_d2_igraph <- graph.data.frame(d = data.frame(from = S_d2_from, to = S_d2_to, stringsAsFactors = TRUE), directed = FALSE, vertices = union(S_d2_from, S_d2_to))
V(S_d2_igraph)$size <- 8
S_d2_community <- cluster_fast_greedy(S_d2_igraph)
sum(crossing(S_d2_community, S_d2_igraph))
modularity(S_d2_community)
mean(sizes(S_d2_community))
length(S_d2_community)

png(filename = paste0(out_dir, "human_sepsis_S_day0_community_graph.png"), width = 800, height = 800)
plot(S_d0_community, S_d0_igraph)
dev.off()
png(filename = paste0(out_dir, "human_sepsis_NS_day0_community_graph.png"), width = 800, height = 800)
plot(NS_d0_community, NS_d0_igraph)
dev.off()
png(filename = paste0(out_dir, "human_sepsis_S_day1_community_graph.png"), width = 800, height = 800)
plot(S_d1_community, S_d1_igraph)
dev.off()
png(filename = paste0(out_dir, "human_sepsis_NS_day1_community_graph.png"), width = 800, height = 800)
plot(NS_d1_community, NS_d1_igraph)
dev.off()
png(filename = paste0(out_dir, "human_sepsis_S_day2_community_graph.png"), width = 800, height = 800)
plot(S_d2_community, S_d2_igraph)
dev.off()

##Plot old or by now uninteresting plots - deactivated to speed up sourcing
#source("generate_stats_old_plot_container.R")

##Human, cluster heatmap, significant metabs only
###Try with pheatmap
x <- human_data
# x[, -1:-6] <- (tanh(scale(x[, -1:-6])) + 1) / 2 #tanh(z-score)
x[, -1:-6] <- apply(scale(x[, -1:-6]), 2, 
                    function(col){ 
                      q <- quantile(x = col, probs = c(0.05, 0.95), na.rm = TRUE) #cut at 5th and 95th percentile
                      col[col > q[2]] <- q[2]
                      col[col < q[1]] <- q[1]
                      return(col)
                    })
x$Day <- as.numeric(as.character(x$Day))
x$Group <- as.character(x$Group)
x$Group <- reorder(x$Group, (1 * (x$Group == "non-Septic-S")) + (2 * (x$Group == "non-Septic-NS")) + (3 * (x$Group == "Septic-S")) + (4 * (x$Group == "Septic-NS")))
x$Day <- reorder(x$Day, x$Day)
x$`CAP / FP` <- reorder(x$`CAP / FP`, (x$`CAP / FP` == "-") + (2 * (x$`CAP / FP` == "FP")) + (3 * (x$`CAP / FP` == "CAP")))
x <- x[order(x$`CAP / FP`),]
x <- x[order(x$Group),]
x <- x[order(x$Day),]
#xm <- x[, c(1:6, metab_sel)]
xm <- x
xmt <- data.frame(t(xm[, -1:-6]))
rownames(xmt) <- colnames(xm[,-1:-6])
colnames(xmt) <- xm$`Sample ID`

xmts <- xmt[rownames(xmt) %in% union(sig.anova.car.s.class, sig.anova.car.s.pheno.class), ]
xmtc <- xmt[rownames(xmt) %in% sig.anova.car.c.class, ]

anno_col <- xm[c("Day", "CAP / FP", "Group")]
rownames(anno_col) <- xm$`Sample ID`
anno_col$Day <- (anno_col$Day)
anno_row <- data.frame(metab_group = human_sepsis_legend$group[match(rownames(xmts), human_sepsis_legend[, 1])])
anno_row$metab_group <- c(as.character(anno_row$metab_group[1:(which(rownames(xmts) == "Albumin") - 1)]), rep("clinical parameter", nrow(anno_row) - which(rownames(xmts) == "Albumin") + 1))
anno_row$metab_group <- factor(anno_row$metab_group)
colnames(anno_row) <- "Metab. group"
rownames(anno_row) <- rownames(xmts)
pat_order <- order(xm$Group, xm$`CAP / FP`, xm$Patient, xm$Day)
y_labels <- rownames(xmts)
y_overlap <- which(y_labels %in% sig.anova.car.nssepnon.class)
y_expr <- sapply(y_labels, function(x) eval(bquote(expression(.(x)))))
y_expr[y_overlap] <- sapply(y_labels[y_overlap], function(x) eval(bquote(expression(bold(.(x))))))
phm <- pheatmap(mat = xmts[, pat_order], 
         border_color = NA,
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = anno_col[pat_order, ],
         annotation_row = anno_row,
         annotation_colors = list(Group = c(`Septic-S` = hue_pal()(4)[3], `Septic-NS` = hue_pal()(4)[1], `non-Septic-S` = hue_pal()(4)[2], `non-Septic-NS` = hue_pal()(4)[4]),
                                  `CAP / FP` = c(`-` = "Grey70", `CAP` = brewer_pal()(2)[1], `FP` = brewer_pal()(2)[2]), 
                                  Day = setNames(seq_gradient_pal(low = "grey97", high = colors()[472])(seq(0, 1, length.out = length(unique(anno_col$Day)))), unique(anno_col$Day)),
                                  `Metab. group` = setNames(hue_pal()(length(unique(anno_row$`Metab. group`))), unique(anno_row$`Metab. group`))),
         gaps_col = which(diff(as.numeric(anno_col$Group[pat_order])) != 0),
         gaps_row = which(diff(as.numeric(anno_row$`Metab. group`)) != 0),
         labels_row = y_expr,
         labels_col = rep("", ncol(xmts)),
         filename = paste0(out_dir, "human_heatmap_metab_sig_s_pheat_thresh005.png"),
         silent = TRUE,
         width = 9,
         height = 10)
##Save for combination with ROC plots in different skript
saveRDS(object = phm, file = "sig_feat_all_pats_heatmap.RData")

##Human, metabolite concentration time course, all metabolites, all groups
hd <- human_data[, -pheno_sel]
n_mets <- ncol(hd) - 6
ncols <- 8
p <- ggplot(data = subset(melt(hd, id.vars = 1:6), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Group, colour = Group)) + 
  facet_wrap(facets = ~ variable, ncol = ncols, nrow = ceiling(n_mets/ncols), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration, µM") +
  xlab("Day") +
  human_col_scale() +
  theme_bw()
ggsave(plot = p, filename = paste0("human_metab_control_time_course_all.png"), path = out_dir, width = 12/4*ncols, height = 0.3 + 1.5 * ceiling(n_mets/ncols), units = "in")

##Human, metabolite concentration time course, all metabolites, all groups, with normal ranges
hd <- human_data[, -pheno_sel]
hd$Group <- factor(hd$Group, levels = c(unique(hd$Group), "Healthy_Fr"))
hd <- subset(melt(hd, id.vars = 1:6), Day %in% 0:3)
hd$variable <- factor(hd$variable)
hhr <- human_healthy_ranges
orig_nrow_hhr <- nrow(hhr)
hhr <- hhr[rep(1:orig_nrow_hhr, each = length(0:3)), ]
hhr$Day <- rep(0:3, times = orig_nrow_hhr)
hhr$variable <- hhr$ID
hhr$Group <- "Healthy_Fr"
hhr$value <- hhr$Median
hhr <- hhr[order(match(hhr$ID, hd$variable)), ]
hhr <- subset(hhr, variable %in% hd$variable)
hhr$variable <- factor(hhr$variable)
n_mets <- length(unique(hd$variable))
ncols <- 8
p <- ggplot(data = hd, mapping = aes(x = Day, y = value, group = Group, colour = Group)) + 
  facet_wrap(facets = ~ variable, ncol = ncols, nrow = ceiling(n_mets/ncols), scales = "free_y") +
  #geom_ribbon(data = hhr, mapping = aes(x = Day, ymin = LQ, ymax = UQ, group = Survival, colour = Survival), fill = "grey70", linetype = 0, alpha = 0.5, inherit.aes = FALSE) +
  geom_rect(data = hhr, mapping = aes(xmin = -Inf, xmax = Inf, ymin = LQ, ymax = UQ, colour = Group), fill = "grey70", linetype = 0, alpha = 0.2, inherit.aes = FALSE) +
  geom_hline(data = hhr, mapping = aes(yintercept = Median), colour = "grey40") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration, µM") +
  xlab("Day") +
  human_col_scale(black_color = "grey40") +
  theme_bw()
ggsave(plot = p, filename = "human_metab_control_time_course_all_w_normal_ranges.png", path = out_dir, width = 12/4*ncols, height = 0.3 + 1.5 * ceiling(n_mets/ncols), units = "in")

##Human, metabolite concentration time course, all metabolites, Survival
hsd <- human_sepsis_data[, -pheno_sel]
n_mets <- ncol(hsd) - 6
ncols <- 8
p <- ggplot(data = subset(melt(hsd, id.vars = 1:6), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Group, colour = Group)) + 
  facet_wrap(facets = ~ variable, ncol = ncols, nrow = ceiling(n_mets/ncols), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration, µM") +
  xlab("Day") +
  human_col_scale() +
  theme_bw()
ggsave(plot = p, filename = paste0("human_metab_survival_time_course_all.png"), path = out_dir, width = 12/4*ncols, height = 0.3 + 1.5 * ceiling(n_mets/ncols), units = "in")

##Human, pheno var concentration time course, all variables
hd <- human_data[, -metab_sel]
n_mets <- ncol(hd) - 6
ncols <- 8
p <- ggplot(data = subset(melt(hd, id.vars = 1:6), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Group, colour = Group)) + 
  facet_wrap(facets = ~ variable, ncol = ncols, nrow = ceiling(n_mets/ncols), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration, µM") +
  xlab("Day") +
  human_col_scale() +
  theme_bw()
ggsave(plot = p, filename = paste0("human_pheno_control_time_course.png"), path = out_dir, width = 12/4*ncols, height = 0.3 + 1.5 * ceiling(n_mets/ncols), units = "in")

##Human, metab var where nonsurvivors are either higher or lower than survivors
hsd <- human_sepsis_data[, c(1:6, which(colnames(human_sepsis_data) %in% c("lysoPC a C24:0", "C14", "PC aa C36:6")))]
ncols <- 3
n_mets <- ncol(hsd) - 6
p <- ggplot(data = subset(melt(hsd, id.vars = 1:6), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Group, colour = Group)) + 
  facet_wrap(facets = ~ variable, ncol = ncols, nrow = ceiling(n_mets/ncols), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration, µM") +
  xlab("Day") +
  ggtitle("Metabolites where Survivors are inbetween Nonsurvivors") +
  human_col_scale() +
  theme_bw()
ggsave(plot = p, filename = paste0("human_metab_survival_time_course_NS_S_NS.png"), path = out_dir, width = 12/4*ncols, height = 0.3 + 1.5 * ceiling(n_mets/ncols), units = "in")

##Human, Creatinie vs Blood Creatinine
hsd <- subset(human_sepsis_data, Day %in% 0:3)
hsd$Patient <- factor(hsd$Patient)
p <- ggplot(data = hsd, mapping = aes_string(x = "Creatinine", y = "`Blood Creatinine`", colour = "Group", shape = "`CAP / FP`")) +
  geom_point() +
  human_col_scale() +
  xlab("Biocrates Creatinine") + 
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(plot = p, filename = "human_Creatinine_reliability.png", path = out_dir, width = 5, height = 4, units = "in")

##Human, metabolite concentration time course, only metabolites with p-val < 0.05 and FDR < 0.05 in Survival or Day:Survival from type III repeated measures ANOVA
fs <- c("s", "c")
gs <- c("by Sepsis survival", "between Sepsis and Control")
ts <- list((length(tanova_day_set)+1):(length(tanova_day_set)*3), seq_along(tanova_day_set))
ss <- list(sig.anova.car.s.class, sig.anova.car.c.class)
for (n in seq_along(ss)[sapply(ss, length) > 0]){
  #h_time_course_sig_diff_dat <- melt(max_norm(full_tanova_data, subset = -1:-6), id.vars = c("Patient", "Survival", "Day", "Sample ID", "CAP / FP"))
  h_time_course_sig_diff_dat <- melt(full_tanova_data, id.vars = c("Patient", "Survival", "Day", "Sample ID", "CAP / FP", "Group"))
  if (n == 1){
    h_time_course_sig_diff_dat <- subset(h_time_course_sig_diff_dat, subset = Survival != "non-Septic")
    cms <- colMaxs(as.matrix(subset(full_tanova_data, !grepl(pattern = "non-Septic", x = Survival), -1:-6)))
  }else{
    h_time_course_sig_diff_dat$Group <- factor(c("Septic", "non-Septic")[1 + (h_time_course_sig_diff_dat$Survival == "non-Septic")])
    cms <- colMaxs(as.matrix(full_tanova_data[, -1:-6]))
  }
  h_time_course_sig_diff_dat <- subset(h_time_course_sig_diff_dat, subset = as.character(variable) %in% ss[[n]])
  h_time_course_sig_diff_dat$variable <- factor(as.character(h_time_course_sig_diff_dat$variable), levels = ss[[n]], ordered = TRUE)
  n_mets <- length(unique(h_time_course_sig_diff_dat$variable))
  h_time_course_group <- data.frame(variable = sort(unique(h_time_course_sig_diff_dat$variable)), value = 1.3, Day = 3, text = coarse_group_list[match(sort(unique(h_time_course_sig_diff_dat$variable)), colnames(human_sepsis_data)[-1:-6])])
  h_time_course_group$value <- cms[match(h_time_course_group$variable, colnames(full_tanova_data[, -1:-6]))]
  h_time_course_group$value <- h_time_course_group$value * 1.3
  lv <- sapply(anova.car.ph.sig.contr, length)
  h_time_course_sig_times <- data.frame(t = "", variable = rep("", sum(lv)), text = "*", value = 1.1, stringsAsFactors = FALSE)
  h_time_course_sig_times$t <- Reduce("c", anova.car.ph.sig.contr)
  h_time_course_sig_times$variable <- factor(rep(names(anova.car.ph.sig.contr), times = lv))
  h_time_course_sig_times$Day <- rep(seq_along(tanova_day_set), times = 2)[h_time_course_sig_times$t]
  h_time_course_sig_times <- subset(h_time_course_sig_times, variable %in% h_time_course_sig_diff_dat$variable & t %in% ts[[n]])
  h_time_course_sig_times$value <- cms[match(h_time_course_sig_times$variable, colnames(full_tanova_data)[-1:-6])]
  h_time_course_sig_times$value <- h_time_course_sig_times$value * 1.1
  hhr <- human_healthy_ranges
  orig_nrow_hhr <- nrow(hhr)
  hhr <- hhr[rep(1:orig_nrow_hhr, each = length(0:3)), ]
  hhr$Day <- rep(0:3, times = orig_nrow_hhr)
  hhr$variable <- hhr$ID
  hhr$Group <- "Healthy_Fr"
  hhr$value <- hhr$Median
  hhr <- hhr[hhr$variable %in% h_time_course_sig_diff_dat$variable, ]
  hhr <- hhr[order(match(hhr$ID, h_time_course_sig_diff_dat$variable)), ]
  hhr$variable <- factor(hhr$variable)
  n_cols <- 5
  h_time_course_sig_diff_plot <- ggplot(h_time_course_sig_diff_dat, aes(x = factor(Day), y = value, group = Group, color = Group)) +
    facet_wrap(facets = ~ variable, ncol = n_cols, nrow = ceiling(n_mets/n_cols), scales = "free_y") +
    stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
    geom_text(data = h_time_course_group, mapping = aes(x = Day, y = value, label = text), inherit.aes = FALSE, size = 2.8) +
    geom_rect(data = hhr, mapping = aes(xmin = -Inf, xmax = Inf, ymin = LQ, ymax = UQ, colour = Group), fill = "grey70", linetype = 0, alpha = 0.2, inherit.aes = FALSE) +
    geom_hline(data = hhr, mapping = aes(yintercept = Median), colour = "grey40") +
    geom_point(data = h_time_course_sig_times, mapping = aes(x = Day, y = value), inherit.aes = FALSE, shape = 8, size = 2.5) +
    geom_point(position = position_dodge(width = 0.2))
  if (n == 1){
    h_time_course_sig_diff_plot <- h_time_course_sig_diff_plot +
      human_col_scale(black_color = "grey40", levels = c("Septic-NS", "", "Septic-S", "", "Healthy_Fr"))
  }else{
    h_time_course_sig_diff_plot <- h_time_course_sig_diff_plot + 
      human_col_scale(black_color = "grey40", levels = c("Septic", "non-Septic", "", "", "Healthy_Fr"))
  }
  h_time_course_sig_diff_plot <- h_time_course_sig_diff_plot +
    ylab("Concentration, µM") + 
    xlab("Day") + 
    ggtitle(paste0("Metabolites significantly differing ", gs[n], " by repeated measures ANOVA")) + 
    theme_bw()
  ggsave(plot = h_time_course_sig_diff_plot, filename = paste0("human_metab_time_course_car_rm_anova_", fs[n], "_sig_diff.png"), path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets/5), units = "in")
}

##Human, pheno var concentration time course, only metabolites with p-val < 0.05 at day1 from type III repeated measures ANOVA, Survival
fs <- c("s", "c")
gs <- c("by Sepsis survival", "between Sepsis and Control")
ts <- list((length(tanova_day_set)+1):(length(tanova_day_set)*3), seq_along(tanova_day_set))
ss <- list(sig.anova.car.s.pheno.class)
for (n in seq_along(ss)[sapply(ss, length) > 0]){ #only runs for sepsis
  #h_time_course_sig_diff_dat <- melt(max_norm(full_tanova_data, subset = -1:-6), id.vars = c("Patient", "Survival", "Day", "Sample ID", "CAP / FP"))
  h_time_course_sig_diff_dat <- melt(full_tanova_data, id.vars = c("Patient", "Survival", "Day", "Sample ID", "CAP / FP", "Group"))
  if (n == 1){
    h_time_course_sig_diff_dat <- subset(h_time_course_sig_diff_dat, subset = Survival != "non-Septic" & Day %in% tanova_day_set)
    cms <- colMaxs(as.matrix(subset(full_tanova_data, !grepl(pattern = "non-Septic", x = Survival), -1:-6)))
  }else{
    h_time_course_sig_diff_dat$Group <- factor(c("Septic", "non-Septic")[1 + (h_time_course_sig_diff_dat$Survival == "non-Septic")])
    cms <- colMaxs(as.matrix(full_tanova_data[, -1:-6]))
  }
  h_time_course_sig_diff_dat <- subset(h_time_course_sig_diff_dat, subset = as.character(variable) %in% ss[[n]])
  h_time_course_sig_diff_dat$variable <- factor(as.character(h_time_course_sig_diff_dat$variable), levels = ss[[n]], ordered = TRUE)
  n_mets <- length(unique(h_time_course_sig_diff_dat$variable))
  h_time_course_group <- data.frame(variable = sort(unique(h_time_course_sig_diff_dat$variable)), value = 1.3, Day = 3, text = coarse_group_list[match(sort(unique(h_time_course_sig_diff_dat$variable)), colnames(human_sepsis_data)[-1:-6])])
  h_time_course_group$value <- cms[match(h_time_course_group$variable, colnames(full_tanova_data[, -1:-6]))]
  h_time_course_group$value <- h_time_course_group$value * 1.3
  lv <- sapply(anova.car.ph.pheno.sig.contr, length)
  h_time_course_sig_times <- data.frame(t = "", variable = rep("", sum(lv)), text = "*", value = 1.1, stringsAsFactors = FALSE)
  h_time_course_sig_times$t <- Reduce("c", anova.car.ph.pheno.sig.contr)
  h_time_course_sig_times$variable <- factor(rep(names(anova.car.ph.pheno.sig.contr), times = lv))
  h_time_course_sig_times$Day <- rep(seq_along(tanova_day_set), times = 2)[h_time_course_sig_times$t]
  h_time_course_sig_times <- subset(h_time_course_sig_times, variable %in% h_time_course_sig_diff_dat$variable & t %in% ts[[n]])
  h_time_course_sig_times$value <- cms[match(h_time_course_sig_times$variable, colnames(full_tanova_data)[-1:-6])]
  h_time_course_sig_times$value <- h_time_course_sig_times$value * 1.1
  n_cols <- 5
  h_time_course_sig_diff_plot <- ggplot(h_time_course_sig_diff_dat, aes(x = factor(Day), y = value, group = Group, color = Group)) +
    facet_wrap(facets = ~ variable, ncol = n_cols, nrow = ceiling(n_mets/n_cols), scales = "free_y") +
    geom_point(position = position_dodge(width = 0.2)) +
    #geom_boxplot(mapping = aes(x = factor(Day), color = Survival, y = value), inherit.aes = FALSE) + 
    stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
    geom_text(data = h_time_course_group, mapping = aes(x = Day, y = value, label = text), inherit.aes = FALSE, size = 2.8) +
    geom_point(data = h_time_course_sig_times, mapping = aes(x = Day, y = value), inherit.aes = FALSE, shape = 8, size = 2.5) +
    # ylim(c(-0.05, 1.4)) +
    # ylab("Concentration relative to max value") +
    human_col_scale() +
    ylab("Concentration, µM") + 
    xlab("Day") + 
    ggtitle(paste0("Phenomenological variables significantly differing ", gs[n], " by repeated measures ANOVA")) + 
    theme_bw()
  ggsave(plot = h_time_course_sig_diff_plot, filename = paste0("human_pheno_time_course_car_rm_anova_", fs[[n]], "_sig_diff.png"), path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets/5), units = "in")
}

##Human, metabolite concentrations signif. diff. between Septic-NS and non-Septic-NS in type 3 ANOVA
h_time_course_sig_diff_dat <- melt(full_tanova_data, id.vars = c("Patient", "Survival", "Day", "Sample ID", "CAP / FP", "Group"))
h_time_course_sig_diff_dat <- subset(h_time_course_sig_diff_dat, variable %in% sig.anova.car.nssepnon.class & Group %in% c("non-Septic-NS", "Septic-NS"))
h_time_course_sig_diff_dat$Group <- factor(h_time_course_sig_diff_dat$Group)
cms <- colMaxs(as.matrix(subset(full_tanova_data, grepl(pattern = "-NS", x = Group), -1:-6)))
h_time_course_group <- data.frame(variable = sort(unique(h_time_course_sig_diff_dat$variable)), value = 1, Day = 3, text = coarse_group_list[match(sort(unique(h_time_course_sig_diff_dat$variable)), colnames(human_sepsis_data)[-1:-6])])
h_time_course_group$value <- cms[match(h_time_course_group$variable, colnames(full_tanova_data[, -1:-6]))]
hhr <- human_healthy_ranges
orig_nrow_hhr <- nrow(hhr)
hhr <- hhr[rep(1:orig_nrow_hhr, each = length(0:3)), ]
hhr$Day <- rep(0:3, times = orig_nrow_hhr)
hhr$variable <- hhr$ID
hhr$Group <- "Healthy_Fr"
hhr$value <- hhr$Median
hhr <- hhr[hhr$variable %in% h_time_course_sig_diff_dat$variable, ]
hhr <- hhr[order(match(hhr$ID, h_time_course_sig_diff_dat$variable)), ]
hhr$variable <- factor(hhr$variable)
h_time_course_group$value <- pmax(h_time_course_group$value * 1.3, hhr[match(h_time_course_group$variable, hhr$Metabolite), ]$UQ * 1.3, na.rm = TRUE)
n_cols <- 5
n_mets <- length(unique(h_time_course_sig_diff_dat$variable))
h_time_course_sig_diff_plot <- ggplot(h_time_course_sig_diff_dat, aes(x = factor(Day), y = value, group = Group, color = Group)) +
  facet_wrap(facets = ~ variable, ncol = n_cols, nrow = ceiling(n_mets/n_cols), scales = "free_y") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  geom_text(data = h_time_course_group, mapping = aes(x = Day, y = value, label = text), inherit.aes = FALSE, size = 2.8) +
  geom_rect(data = hhr, mapping = aes(xmin = -Inf, xmax = Inf, ymin = LQ, ymax = UQ, colour = Group), fill = "grey70", linetype = 0, alpha = 0.2, inherit.aes = FALSE) +
  geom_hline(data = hhr, mapping = aes(yintercept = Median), colour = "grey40") +
  geom_point(position = position_dodge(width = 0.2)) +
  human_col_scale(black_color = "grey40", levels = c("Septic-NS", "", "", "non-Septic-NS", "Healthy_Fr")) +
  ylab("Concentration, µM") + 
  xlab("Day") + 
  ggtitle(paste0("Metabolites significantly differing between Septic-NS and non-Septic-NS by repeated measures ANOVA")) + 
  theme_bw()
ggsave(plot = h_time_course_sig_diff_plot, filename = paste0("human_metab_time_course_car_rm_anova_septicns_nonsepns_sig_diff.png"), path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets/5), units = "in")

##Human, metabolite concentrations signif. diff. between non-Septic-S and non-Septic-NS in type 3 ANOVA
h_time_course_sig_diff_dat <- melt(full_tanova_data, id.vars = c("Patient", "Survival", "Day", "Sample ID", "CAP / FP", "Group"))
h_time_course_sig_diff_dat <- subset(h_time_course_sig_diff_dat, variable %in% sig.anova.car.cnscs.class & Group %in% c("non-Septic-NS", "non-Septic-S"))
h_time_course_sig_diff_dat$Group <- factor(h_time_course_sig_diff_dat$Group)
cms <- colMaxs(as.matrix(subset(full_tanova_data, grepl(pattern = "-NS", x = Group), -1:-6)))
h_time_course_group <- data.frame(variable = sort(unique(h_time_course_sig_diff_dat$variable)), value = 1, Day = 3, text = coarse_group_list[match(sort(unique(h_time_course_sig_diff_dat$variable)), colnames(human_sepsis_data)[-1:-6])])
h_time_course_group$value <- cms[match(h_time_course_group$variable, colnames(full_tanova_data[, -1:-6]))]
hhr <- human_healthy_ranges
orig_nrow_hhr <- nrow(hhr)
hhr <- hhr[rep(1:orig_nrow_hhr, each = length(0:3)), ]
hhr$Day <- rep(0:3, times = orig_nrow_hhr)
hhr$variable <- hhr$ID
hhr$Group <- "Healthy_Fr"
hhr$value <- hhr$Median
hhr <- hhr[hhr$variable %in% h_time_course_sig_diff_dat$variable, ]
hhr <- hhr[order(match(hhr$ID, h_time_course_sig_diff_dat$variable)), ]
hhr$variable <- factor(hhr$variable)
h_time_course_group$value <- pmax(h_time_course_group$value * 1.3, hhr[match(h_time_course_group$variable, hhr$Metabolite), ]$UQ * 1.3, na.rm = TRUE)
n_cols <- 5
n_mets <- length(unique(h_time_course_sig_diff_dat$variable))
h_time_course_sig_diff_plot <- ggplot(h_time_course_sig_diff_dat, aes(x = factor(Day), y = value, group = Group, color = Group)) +
  facet_wrap(facets = ~ variable, ncol = n_cols, nrow = ceiling(n_mets/n_cols), scales = "free_y") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  geom_text(data = h_time_course_group, mapping = aes(x = Day, y = value, label = text), inherit.aes = FALSE, size = 2.8) +
  geom_rect(data = hhr, mapping = aes(xmin = -Inf, xmax = Inf, ymin = LQ, ymax = UQ, colour = Group), fill = "grey70", linetype = 0, alpha = 0.2, inherit.aes = FALSE) +
  geom_hline(data = hhr, mapping = aes(yintercept = Median), colour = "grey40") +
  geom_point(position = position_dodge(width = 0.2)) +
  human_col_scale(black_color = "grey40", levels = c("", "non-Septic-S", "", "non-Septic-NS", "Healthy_Fr")) +
  ylab("Concentration, µM") + 
  xlab("Day") + 
  ggtitle(paste0("Metabolites significantly differing between non-Septic-S and non-Septic-NS by repeated measures ANOVA")) + 
  theme_bw()
ggsave(plot = h_time_course_sig_diff_plot, filename = paste0("human_metab_time_course_car_rm_anova_cns_cs_sig_diff.png"), path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets/5), units = "in")

##Safety corridor concept with Survivor min-max as corridor
pat_dev_score <- subset(human_data, Day %in% tanova_day_set, select = setdiff(which(!(colnames(human_data) %in% sig.anova.car.s.class)), pheno_sel)) # select already includes cols 1:6
colnames(pat_dev_score)[na.omit(match(human_sepsis_legend[which(human_sepsis_legend$group == "amino acid"), 1], colnames(pat_dev_score)))] <- human_sepsis_legend$name[human_sepsis_legend$group == "amino acid" & human_sepsis_legend[, 1] %in% colnames(pat_dev_score)] #assignment OK, checked manually
pat_dev_max <- colMaxs(as.matrix(pat_dev_score[pat_dev_score$Survival == "S", -1:-6]))
pat_dev_min <- colMins(as.matrix(pat_dev_score[pat_dev_score$Survival == "S", -1:-6]))
udev <- pat_dev_score[, -1:-6] > matrix(pat_dev_max, ncol = ncol(pat_dev_score) - 6, nrow = nrow(pat_dev_score), byrow = TRUE)
ldev <- pat_dev_score[, -1:-6] < matrix(pat_dev_min, ncol = ncol(pat_dev_score) - 6, nrow = nrow(pat_dev_score), byrow = TRUE)
sdev <- aggregate(udev | ldev, by = list(Patient = pat_dev_score$Patient), FUN = max) #count same metabolite at different time points as one deviation
met_pat_score <- data.frame(Metabolite = colnames(sdev)[-1], score = colSums(sdev[, -1]), stringsAsFactors = FALSE)
met_pat_score_sep <- data.frame(Metabolite = colnames(sdev[sdev$Patient %in% human_sepsis_data$Patient, ])[-1], score = colSums(sdev[sdev$Patient %in% human_sepsis_data$Patient, -1]), stringsAsFactors = FALSE)
met_pat_score_nonsep <- data.frame(Metabolite = colnames(sdev[sdev$Patient %in% human_nonsepsis_data$Patient, ])[-1], score = colSums(sdev[sdev$Patient %in% human_nonsepsis_data$Patient, -1]), stringsAsFactors = FALSE)
subset(met_pat_score_sep, score >= 4 & !(Metabolite %in% met_pat_score_nonsep[met_pat_score_nonsep$score > 0, 1])) # Metabolite devs specific for septic-NS
sdev_melted <- subset(melt(sdev, id.vars = 1), value == 1)
score_thresh <- 3
minmax_corridor_met_sel <- subset(met_pat_score, score > score_thresh)
minmax_corridor_met_sel$color <- rev(grey_pal()(max(minmax_corridor_met_sel$score) - score_thresh))[minmax_corridor_met_sel$score - score_thresh]
hsd <- subset(pat_dev_score, Day %in% tanova_day_set, select = c(colnames(human_data)[1:6], minmax_corridor_met_sel$Metabolite))
hsd <- cbind(hsd[1:6], hsd[-1:-6][order(colMeans(as.matrix(subset(hsd, Survival == "S", -1:-6))), decreasing = TRUE)])
hsd <- melt(hsd, id.vars = 1:6)
pat_path <- melt(subset(pat_dev_score, Survival == "NS" & Day %in% tanova_day_set), id.vars = 1:6)
pat_path <- subset(pat_path, variable %in% minmax_corridor_met_sel$Metabolite)
pat_path <- subset(pat_path, interaction(variable, Patient) %in% interaction(sdev_melted$variable, sdev_melted$Patient))
pat_path$x <- match(pat_path$variable, unique(hsd$variable))
pat_path$x <- pat_path$x + (scale(seq_along(unique(pat_path$Day)))/4)[unlist(lapply(rle(pat_path$Patient)[["lengths"]], function(to) 1:to))]
pat_path$metabolite_group <- human_sepsis_legend$group[match(pat_path$variable, human_sepsis_legend[, 1])]
pat_path$metabolite_group[is.na(pat_path$metabolite_group)] <- human_sepsis_legend$group[match(pat_path$variable, human_sepsis_legend$name)][is.na(pat_path$metabolite_group)]
rsize <- rel(4.5)
p <- ggplot(data = subset(hsd, Survival == "S"), mapping = aes(y = value, x = variable, color = Group)) +
  stat_summary(fun.y = "mean", fun.ymin = "min", fun.ymax = "max", position = position_dodge(width = 0.7), geom = "errorbar", size = rel(0.8)) +
  #stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax = "mean", position = position_dodge(width = 0.7), geom = "errorbar") +
  geom_line(data = pat_path, mapping = aes(y = value, x = x, color = Group, group = interaction(variable, Patient), linetype = factor(Patient)), size = rel(0.8), inherit.aes = FALSE) +
  geom_tile(mapping = aes(x = variable, y = 2e-3, fill = metabolite_group), width = 1, height = 0.15, data = pat_path, inherit.aes = FALSE) +
  geom_point(mapping = aes(x = Metabolite, y = 1.45e-3, shape = factor(score)), data = minmax_corridor_met_sel, inherit.aes = FALSE) +
  geom_tile(mapping = aes(x = Metabolite, y = 1.45e-3), width = 1, height = 0.15, fill = minmax_corridor_met_sel$color, data = minmax_corridor_met_sel, inherit.aes = FALSE) +
  guides(color = guide_legend(title = "Patient Group", order = 1, override.aes = list(size = 2.5), direction = "vertical"),
        fill = guide_legend(title = "Metabolite\nGroup", order = 3, nrow = 2, override.aes = list(size = 8)),
        shape = guide_legend(title = "Number of patients\nwith deviation", override.aes = list(shape = 15, size = 8, colour = sort(unique(minmax_corridor_met_sel$color), decreasing = TRUE)), order = 4, keywidth = rel(1.1), keyheight = rel(1.1), ncol = 2),
        linetype = guide_legend(title = "NS Patient", override.aes = list(color = hue_pal()(4)[c(4, 1)[1 + (subset(hsd, !duplicated(Patient) & Survival == "NS", "Group")[[1]] == "Septic-NS")]], size = 1.25), order = 2, nrow = 3, keywidth = rel(6), direction = "vertical")) +
  scale_color_discrete(drop = FALSE) +
  scale_y_log10(expand = c(0,0)) +
  ylab("Concentration, µM") +
  xlab("") +
  human_col_scale() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = rsize), panel.grid = element_line(colour = 0), text = element_text(size = rsize), axis.text.y = element_text(size = rsize), legend.text = element_text(size = rsize), legend.box = "vertical", legend.position = "bottom")
ggsave(filename = paste0("human_metab_nonsig_single_plot_minmax_all_pats.png"), path = out_dir, plot = p, width = 16, height = 18, units = "in")
p_main <- p

hsd_lpc24 <- subset(human_data, Day %in% tanova_day_set)
hsd_lpc24 <- melt(hsd_lpc24, id.vars = 1:6)
hsd_lpc24 <- subset(hsd_lpc24, as.character(variable) == "lysoPC a C24:0" & Survival == "S")
hsd_lpc24_max <- aggregate(x = hsd_lpc24$value, by = list(Group = hsd_lpc24$Group), FUN = max)
hsd_lpc24_min <- aggregate(x = hsd_lpc24$value, by = list(Group = hsd_lpc24$Group), FUN = min)
pp_lpc24 <- subset(melt(subset(human_data, Survival == "NS" & Day %in% tanova_day_set), id.vars = 1:6), as.character(variable) == "lysoPC a C24:0")
rsize <- rel(4.5)
p <- ggplot(data = pp_lpc24, mapping = aes(y = value, x = Day, color = Group, group = Patient, linetype = factor(Patient))) +
  facet_wrap(~ variable, ncol = 1, nrow = 1) +
  geom_line(size = rel(0.8)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  geom_hline(data = hsd_lpc24_min, mapping = aes(color = Group, yintercept = x), size = rel(0.8)) +
  geom_hline(data = hsd_lpc24_max, mapping = aes(color = Group, yintercept = x), size = rel(0.8)) +
  # guides(color = guide_legend(order = 1),
  #        linetype = guide_legend(title = "NS Patient", override.aes = list(color = hue_pal()(4)[c(4, 1)[1 + (subset(hsd, !duplicated(Patient) & Survival == "NS", "Group")[[1]] == "Septic-NS")]]), order = 4)) +
  guides(color = "none", linetype = "none", group = "none") +
  scale_color_discrete(drop = FALSE) +
  #scale_y_continuous(trans = pseudo_log_trans(sigma = 0.25, base = 2), expand = c(0, 0), limits = c(-10, 1000)) +
  #ylab("Concentration, µM") +
  ylab("") +
  xlab("Day") +
  human_col_scale() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, size = 16), axis.text.y = element_text(size = 16), panel.grid = element_line(colour = 0), text = element_text(size = rsize), strip.text = element_text(size = rsize), strip.background = element_blank())
ggsave(filename = paste0("human_metab_nonsig_single_plot_minmax_all_pats_lysoPCaC24.png"), path = out_dir, plot = p, width = 5, height = 5, units = "in")
p_insert <- p

p_combined <- ggdraw(p_main) + 
  draw_plot(p_insert, 0.71, 0.715, 0.28, 0.28) +
  draw_plot_label(label = c("A", "B"), x = c(0, 0.71), y = c(1, 0.995), size = 20)
ggsave(filename = "human_metab_nonsig_single_plot_minmax_all_pats_panel.png", path = out_dir, plot = p_combined, width = 16, height = 18, units = "in")
ggsave(filename = "human_metab_nonsig_single_plot_minmax_all_pats_panel.svg", path = out_dir, plot = p_combined, width = 16, height = 18, units = "in")

##Generalized corridor defined by min and max of survivor (septic and non-Septic)
pat_dev_score <- subset(human_data, Day %in% 0:3, select = setdiff(which(!(colnames(human_data) %in% sig.anova.car.s.class)), pheno_sel)) # select already includes cols 1:6
pat_dev_max <- colMaxs(as.matrix(pat_dev_score[pat_dev_score$Survival == "S", -1:-6]))
pat_dev_min <- colMins(as.matrix(pat_dev_score[pat_dev_score$Survival == "S", -1:-6]))
udev <- pat_dev_score[, -1:-6] > matrix(pat_dev_max, ncol = ncol(pat_dev_score) - 6, nrow = nrow(pat_dev_score), byrow = TRUE)
ldev <- pat_dev_score[, -1:-6] < matrix(pat_dev_min, ncol = ncol(pat_dev_score) - 6, nrow = nrow(pat_dev_score), byrow = TRUE)
sdev <- aggregate(udev | ldev, by = list(Patient = pat_dev_score$Patient), FUN = max) #count same metabolite at different time points as one deviation
dev_score <- data.frame(Patient = sdev$Patient, score = rowSums(sdev[, -1]))
dev_score$Survival <- pat_dev_score$Survival[match(dev_score$Patient, pat_dev_score$Patient)]
dev_score$Group <- pat_dev_score$Group[match(dev_score$Patient, pat_dev_score$Patient)]
{
  set.seed(2003)
  m <- subset(human_data, Day %in% tanova_day_set)
  m <- m[, -pheno_sel]
  m <- m[, !(colnames(m) %in% sig.anova.car.s.class)]
  di <- which(m$Survival == "NS")
  dev_rep_res <- sapply(X = 1:1000, FUN = sim_dev, n = nrow(m), d = ncol(m) - 6, dev_idx = di, sample_groups = m$Patient)
}
p_thresh <- apply(dev_rep_res, 1, quantile, p = 0.95, names = FALSE)
p_thresh_sig <- rbind(dev_score$score[match(as.numeric(names(p_thresh)), dev_score$Patient)], 
                      dev_score$score[match(as.numeric(names(p_thresh)), dev_score$Patient)] > p_thresh)
p_thresh_nonu <- table(p_thresh_sig[1, ])
p_thresh_sig <- p_thresh_sig[, order(p_thresh_sig[1, ])]
p_thresh_sig <- rbind(p_thresh_sig, unlist(lapply(names(p_thresh_nonu), function(nonu) 1:p_thresh_nonu[nonu] - 0.5)))
p_thresh_sig <- data.frame(t(p_thresh_sig))
colnames(p_thresh_sig) <- c("score", "include", "y")
p_thresh_sig_UK_minmax <- p_thresh_sig
p_thresh_UK_minmax <- p_thresh
p <- ggplot(data = dev_score, mapping = aes(fill = Group, x = score)) +
  geom_histogram(position = position_stack(), bins = max(dev_score$score) + 1) +
  geom_point(data = subset(p_thresh_sig, include == 1), mapping = aes(x = score, y = y), shape = 8, size = 1.2, inherit.aes = FALSE) +
  human_col_scale(aesthetics = "fill") +
  ylab("Number of Patients") +
  xlab("Number of metabolites outside of the safe corridor at Days 0-3") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.direction = "horizontal", legend.position = "bottom")
ggsave(plot = p, filename = "generalized_safe_corridor_minmax.png", path = out_dir, width = 8, height = 4, units = "in")
dev_score$x <- factor("Our study")
pbox <- ggplot(data = subset(dev_score, Group %in% c("Septic-NS", "non-Septic-NS")), mapping = aes(fill = Group, x = score, y = x, colour = Group)) +
  geom_dotplot(stackdir = "center", binaxis = "x") + 
  human_col_scale(aesthetics = c("fill", "colour")) +
  xlab("Number of metabolites outside of the safe corridor at Days 0-3") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.direction = "horizontal", legend.position = "bottom")
ggsave(plot = pbox, filename = "generalized_safe_corridor_minmax_box.png", path = out_dir, width = 8, height = 4, units = "in")
dev_score_UK_minmax <- dev_score
print(paste0("Contribution of high deviation and low deviation to score is ", sum(udev), " and ", sum(ldev), " respectively."))
w <- which.xy(udev | ldev) # tell me which variables make a difference
mtab <- sort(table(w[, 2]), decreasing = TRUE)
names(mtab) <- colnames(pat_dev_score)[-1:-6][as.numeric(names(mtab))]
mtab
dev_mets_out_UK <- data.frame(Name = names(mtab), stringsAsFactors = FALSE)
dev_mets_out_UK$Group <- human_sepsis_legend$group[match(names(mtab), human_sepsis_legend[, 1])]
dev_mets_out_UK$PatCount <- colSums(sdev[, names(mtab)])
uk_minmax_mtab <- mtab
###Find out p-values for number of patients per deviation
{
  set.seed(2005)
  m <- subset(human_data, Day %in% tanova_day_set)
  m <- m[, -pheno_sel]
  m <- m[, !(colnames(m) %in% sig.anova.car.s.class)]
  di <- which(m$Group == "Septic-NS")
  dev_rep_res <- sapply(X = 1:1000, FUN = sim_dev_met, n = nrow(m), d = ncol(m) - 6, dev_idx = di, sample_groups = m$Patient)
  set.seed(2005)
  di_nns <- which(m$Group == "non-Septic-NS")
  dev_rep_res_nns <- sapply(X = 1:1000, FUN = sim_dev_met, n = nrow(m), d = ncol(m) - 6, dev_idx = di_nns, sample_groups = m$Patient)
  set.seed(2005)
  m <- m[, c(1:6, which(colnames(m) %in% human_sepsis_legend[human_sepsis_legend$group %in% c("acylcarnitine", "lysophosphatidylcholine"), 1]))]
  dev_rep_res_ac_lpc <- sapply(X = 1:1000, FUN = sim_dev, n = nrow(m), d = ncol(m) - 6, dev_idx = di, sample_groups = m$Patient)
}
sink(file = paste0(out_dir, "safety_corridor_deviation_count_p_values.txt"))
###Find out how likely a number of deviations in Septic-NS is
dev_met_p_template <- rev(cumsum(rev(table(as.numeric(dev_rep_res)))) / sum(table(as.numeric(dev_rep_res)))) #rev()'ed to sum up from the extreme end
print(paste0("p of number of Septic-NS patients with as many deviations:"))
print(dev_met_p_template)
dev_met_p_tab <- colSums(subset(sdev, Patient %in% human_data$Patient[human_data$Group == "Septic-NS"], -1))
dev_met_p_tab <- dev_met_p_template[match(dev_met_p_tab, names(dev_met_p_template))]
dev_met_q_tab <- p.adjust(dev_met_p_tab, method = "fdr")
names(dev_met_q_tab) <- colnames(subset(sdev, Patient %in% human_data$Patient[human_data$Group == "Septic-NS"], -1))[-1]
print("metabolites where at least number of Septic-NS patients with deviation has q <= 0.05")
print(dev_met_q_tab[names(which(dev_met_q_tab <= 0.05))])
###Find out how likely a number of deviations in Septic-NS is when no non-Septic-NS patient deviates
dev_rep_res_sns_nns <- dev_rep_res
dev_rep_res_sns_nns[dev_rep_res_nns > 0] <- 0
dev_met_p_template_sns_nns <- rev(cumsum(rev(table(as.numeric(dev_rep_res_sns_nns)))) / sum(table(as.numeric(dev_rep_res_sns_nns)))) #rev()'ed to sum up from the extreme end
print("p of number of Septic-NS patients deviating with as at least as many deviations while no non-Septic-NS patient deviating:")
print(dev_met_p_template_sns_nns)
dev_met_sns_nns_p_tab <- colSums(subset(sdev, Patient %in% human_data$Patient[human_data$Group == "Septic-NS"], -1))
dev_met_sns_nns_p_tab[colSums(subset(sdev, Patient %in% human_data$Patient[human_data$Group == "non-Septic-NS"], -1)) >= 1] <- 0
dev_met_sns_nns_p_tab <- dev_met_p_template_sns_nns[match(dev_met_sns_nns_p_tab, names(dev_met_p_template_sns_nns))]
dev_met_sns_nns_q_tab <- p.adjust(dev_met_sns_nns_p_tab, method = "fdr")
names(dev_met_sns_nns_p_tab) <- colnames(subset(sdev, Patient %in% human_data$Patient[human_data$Group == "Septic-NS"], -1))
names(dev_met_sns_nns_q_tab) <- names(dev_met_sns_nns_p_tab)
print("metabolites where the number of Septic-NS patients with deviation is above this and no non-Septic-NS patient deviates:")
print(dev_met_sns_nns_q_tab[names(which(dev_met_sns_nns_q_tab <= 0.05))])
###Find out how likely you get 7/8 patients deviating for at least one AC or PC #jump
dev_rep_res_7_8 <- colSums(dev_rep_res_ac_lpc != 0)
dev_met_p_template_7_8 <- rev(cumsum(rev(table(as.numeric(dev_rep_res_7_8))))) / sum(table(as.numeric(dev_rep_res_7_8)))
print("p of number of Septic-NS patients deviating in at least this many ACs or PCs:")
print(dev_met_p_template_7_8)
sink()

###Find minimal combination of features to seperate sNS from sS
min_met_set_for_dev <- subset(sdev, Patient %in% human_sepsis_data$Patient[human_sepsis_data$Group == "Septic-NS"])
min_met_set_for_dev <- min_met_set_for_dev[, c(1, 1 + which(colAnys(as.matrix(min_met_set_for_dev[, -1]))))]
pat_dev_idx_per_met <- lapply(lapply(min_met_set_for_dev[, -1], as.logical), which)
pat_dev_idx_met_crossover <- lapply(pat_dev_idx_per_met, function(e) lapply(pat_dev_idx_per_met, union, e))
full_idx <- Reduce("rbind", lapply(lapply(pat_dev_idx_met_crossover, lapply, length), sapply, `==`, nrow(min_met_set_for_dev)))
xy <- which.xy(full_idx)
colnames(min_met_set_for_dev[, -1])[c(xy$x[1], xy$y[1])]
colnames(min_met_set_for_dev[, -1])[c(xy$x[2], xy$y[2])]
###Simulate how likely this is by chance
{
  set.seed(2145)
  m <- subset(human_data, Day %in% tanova_day_set)
  m <- m[, -pheno_sel]
  m <- m[, !(colnames(m) %in% sig.anova.car.s.class)]
  di <- which(m$Survival == "NS")
  dev_rep_res <- sapply(X = 1:1000, FUN = sim_dev_triple, n = nrow(m), d = ncol(m) - 6, dev_idx = di, sample_groups = m$Patient)
}
dev_rep_tab <- table(dev_rep_res)
print("P-value of having three metabolites whose deviations cover all Septic-NS patients")
print(sum(dev_rep_tab[-1]) / sum(dev_rep_tab))
###

sink(file = paste0(out_dir, "devs_UK_sNS_met_counts.csv"))
print("Test if lysoPC a C24:0, PC aa C36:0 and lysoPC a C18:2 provide enough deviations to cover all sNS patients")
min_met_set_for_dev$`lysoPC a C24:0` + min_met_set_for_dev$`PC aa C36:0` + min_met_set_for_dev$`lysoPC a C18:2`
print("ColSums and RowSums of deviations per metabolite group")
print("lysoPC")
colSums(min_met_set_for_dev[, 1 + grep(pattern = "lysoPC", x = colnames(min_met_set_for_dev)[-1])])
rowSums(min_met_set_for_dev[, 1 + grep(pattern = "lysoPC", x = colnames(min_met_set_for_dev)[-1])])
print("Acylcarnitines")
colSums(min_met_set_for_dev[, 1 + grep(pattern = "^C", x = colnames(min_met_set_for_dev)[-1])])
rowSums(min_met_set_for_dev[, 1 + grep(pattern = "^C", x = colnames(min_met_set_for_dev)[-1])])
print("Sphingolipids")
colSums(min_met_set_for_dev[, 1 + grep(pattern = "SM", x = colnames(min_met_set_for_dev)[-1])])
rowSums(min_met_set_for_dev[, 1 + grep(pattern = "SM", x = colnames(min_met_set_for_dev)[-1])])
sink()

nNS_min_met_set <- subset(sdev, Patient %in% human_data$Patient[human_data$Group == "non-Septic-NS"])
rowSums(nNS_min_met_set[c("lysoPC a C24:0", "PC aa C36:0", "lysoPC a C18:2")])

##Generalized corridor validation on Ferrario data
vd <- human_sepsis_val_data[, c(-3:-4, -4 + union(-which(val_anova_fdr[3, ] < 0.05), -which(val_anova_fdr[4, ] < 0.05)))]
###Survivor mean +- SD corridor
vd_mean <- colMeans(vd[vd$Survival28 == "S", -1:-2])
vd_sd <- colSds(as.matrix(vd[vd$Survival28 == "S", -1:-2]))
sd_mul <- 4.5
udev <- vd[, -1:-2] > matrix(vd_mean + sd_mul * vd_sd, ncol = ncol(vd) - 2, nrow = nrow(vd), byrow = TRUE)
ldev <- vd[, -1:-2] < matrix(vd_mean - sd_mul * vd_sd, ncol = ncol(vd) - 2, nrow = nrow(vd), byrow = TRUE)
sdev <- aggregate(udev | ldev, by = list(Patient = vd$Patient), FUN = max) #count same metabolite at two time points as one deviation
dev_score <- data.frame(Patient = sdev$Patient, score = rowSums(sdev[, -1]))
dev_score$Survival <- vd$Survival28[match(dev_score$Patient, vd$Patient)]
p <- ggplot(data = dev_score, mapping = aes(fill = Survival, x = score)) +
  geom_histogram(position = position_stack(), bins = max(dev_score$score) + 1) +
  ylab("Number of Patients") +
  xlab("Number of metabolites outside of the safe corridor at Days 0 & 7") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(plot = p, filename = "generalized_safe_corridor_Ferrario_SD4.5.png", path = out_dir, width = 6, height = 3, units = "in")
print(paste0("Contribution of high deviation and low deviation to score is ", sum(udev), " and ", sum(ldev), " respectively."))
w <- which.xy(udev | ldev) # tell me which variables make a difference
mtab <- sort(table(w[, 2]), decreasing = TRUE)
names(mtab) <- colnames(vd)[-1:-2][as.numeric(names(mtab))]
mtab
###Survivor min-max corridor
vd_max <- colMaxs(as.matrix(vd[vd$Survival28 == "S", -1:-2]))
vd_min <- colMins(as.matrix(vd[vd$Survival28 == "S", -1:-2]))
udev <- vd[, -1:-2] > matrix(vd_max, ncol = ncol(vd) - 2, nrow = nrow(vd), byrow = TRUE)
ldev <- vd[, -1:-2] < matrix(vd_min, ncol = ncol(vd) - 2, nrow = nrow(vd), byrow = TRUE)
sdev <- aggregate(udev | ldev, by = list(Patient = vd$Patient), FUN = max) #count same metabolite at two time points as one deviation
dev_score <- data.frame(Patient = sdev$Patient, score = rowSums(sdev[, -1]))
dev_score$Survival <- vd$Survival28[match(dev_score$Patient, vd$Patient)]
{
  set.seed(1003)
  m <- vd
  m <- m[, !(colnames(m) %in% sig.anova.car.val.s.class)]
  di <- which(m$Survival == "NS")
  dev_rep_res <- sapply(X = 1:1000, FUN = sim_dev, n = nrow(m), d = ncol(m) - 1, dev_idx = di, sample_groups = m$Patient)
}
p_thresh <- apply(dev_rep_res, 1, quantile, p = 0.95, names = FALSE)
p <- ggplot(data = dev_score, mapping = aes(fill = Survival, x = score)) +
  geom_histogram(position = position_stack(), bins = max(dev_score$score) + 1) +
  ylab("Number of Patients") +
  xlab("Number of metabolites outside of the safe corridor at Days 0 & 7") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(plot = p, filename = "generalized_safe_corridor_Ferrario_minmax.png", path = out_dir, width = 6, height = 3, units = "in")
dev_score_vd_minmax <- dev_score
dev_score_vd_minmax$x <- factor("Ferrario et al., 2016")
dev_score_vd_minmax$Group <- paste0("Septic-", dev_score_vd_minmax$Survival)
dev_score_cb_minmax <- rbind(dev_score_UK_minmax, dev_score_vd_minmax)
dev_score_cb_minmax$x <- factor(dev_score_cb_minmax$x)
if (t.test(x = subset(dev_score_UK_minmax, Group == "Septic-NS", "score"), y = subset(dev_score_UK_minmax, Group == "non-Septic-NS", "score"))$p.value > 0.05){
  dev_score_cb_minmax <- subset(dev_score_cb_minmax, Group != "non-Septic-NS")
  xmax <- 0.9
  legx <- 0.12
} else {
  xmax <- 1.3
  legx <- 0.2
}
pbox <- ggplot(data = subset(dev_score_cb_minmax, Survival == "NS"), mapping = aes(fill = Group, x = x, y = score, colour = Group)) +
  geom_dotplot(stackdir = "center", binaxis = "y", dotsize = 0.7) + 
  human_col_scale(aesthetics = c("fill", "colour")) +
  scale_x_discrete(limits = levels(dev_score_cb_minmax$x)[2:1]) +
  #geom_path(inherit.aes = FALSE, data = data.frame(y = min(subset(p_thresh_sig_UK_minmax, include == 1)$score) - 1, x = 1 + c(0.6, 1.4)), mapping = aes(x, y), linetype = 2) +
  #geom_path(inherit.aes = FALSE, data = data.frame(y = mean(p_thresh), x = c(0.6, 1.4)), mapping = aes(x, y), linetype = 2) +
  geom_rect(data = data.frame(xmin = 0.4, xmax = xmax, ymin = 56, ymax = 70), mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white", colour = "black", inherit.aes = FALSE) +
  coord_flip() +
  ylab("Number of metabolites outside of the safe corridor per patient") +
  xlab("") +
  scale_y_continuous(limits = c(0, 70), expand = c(0,0)) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.grid.major.y = element_line(color = "grey80"), legend.direction = "vertical", legend.position = c(0.9, legx), legend.title = element_blank(), legend.background = element_blank())
ggsave(plot = pbox, filename = "generalized_safe_corridor_minmax_box_UK_Ferrario.png", path = out_dir, width = 7, height = 2, units = "in")
ggsave(plot = pbox, filename = "generalized_safe_corridor_minmax_box_UK_Ferrario.svg", path = out_dir, width = 7, height = 2, units = "in")
print(paste0("Contribution of high deviation and low deviation to score is ", sum(udev), " and ", sum(ldev), " respectively."))
w <- which.xy(udev | ldev) # tell me which variables make a difference
mtab <- sort(table(w[, 2]), decreasing = TRUE)
names(mtab) <- colnames(vd)[-1:-2][as.numeric(names(mtab))]
mtab
length(intersect(names(uk_minmax_mtab), names(mtab)))
sort(intersect(names(uk_minmax_mtab), names(mtab))) # metabolites that deviate in UK and Ferrario data
intersect(setdiff(colnames(human_sepsis_val_data)[-1:-4], sig.anova.car.val.s.class), # metabolites that could deviate in both data sets
          setdiff(colnames(human_sepsis_data)[-1:-6], sig.anova.car.val.s.class))
dev_mets_out_Ferrario <- data.frame(Name = names(mtab), stringsAsFactors = FALSE)
dev_mets_out_Ferrario$Group <- human_sepsis_legend$group[match(names(mtab), human_sepsis_legend[, 1])]
dev_mets_out_Ferrario$PatCount <- colSums(sdev[, names(mtab)])
dev_mets_out_all <- merge(x = dev_mets_out_UK, y = dev_mets_out_Ferrario, by = c("Name", "Group"), all = TRUE, suffixes = c("_UK", "_Ferrario"))
dev_mets_out_all <- dev_mets_out_all[order(dev_mets_out_all$Group, dev_mets_out_all$Name), ]
dev_mets_out_all$PatCount_UK[is.na(dev_mets_out_all$PatCount_UK)] <- "-" #not 0 because it was significant in sS vs sNS
dev_mets_out_all$PatCount_Ferrario[is.na(dev_mets_out_all$PatCount_Ferrario) & !(dev_mets_out_all$Name %in% sig.anova.car.val.s.class)] <- 0
dev_mets_out_all$PatCount_Ferrario[dev_mets_out_all$Name %in% sig.anova.car.val.s.class] <- "-" #not 0 because significant in sS vs sNS
dev_mets_out_all$PatCount_Ferrario[!(dev_mets_out_all$Name %in% colnames(human_sepsis_val_data))] <- "-" #not 0 because metabolite not in data set
fwrite(x = dev_mets_out_all, file = paste0(out_dir, "generalized_safe_corridor_minmax_dev_mets.csv"))

###Find number of deviations by chance per metabolite and overlap by chance between Ferrario et al and our study
n_possible_overlap <- length(intersect(colnames(vd)[-1:-2], colnames(pat_dev_score)[-1:-6]))
n_overlap <- length(intersect(names(uk_minmax_mtab), names(mtab)))
chisq.test(x = c(n_overlap, n_possible_overlap - n_overlap), y = c(3, 1)) #3 out of 4 SMs
chisq.test(x = c(n_overlap, n_possible_overlap - n_overlap), y = c(5, 2)) #5 out of 7 LysoPCs
chisq.test(x = c(n_overlap, n_possible_overlap - n_overlap), y = c(2, 1)) #2 out of 3 ACs
#All p-values = 1 ...

###Plot all metabolites we mention in the text
dev_mets_out_all_pie <- dev_mets_out_all
mets_mentioned <- c("C18:1-OH", paste0("PC ae ", c("C34:2", "C38:0", "C38:5", "C38:6", "C42:0", "C42:2", "C42:3", "C40:1", "C40:2")), "PC aa C36:0", "lysoPC a C18:2", "lysoPC a C24:0")
dev_mets_out_all_pie <- subset(dev_mets_out_all_pie, Name %in% mets_mentioned)
dev_mets_out_all_pie$IPatCount_UK <- sum(human_data$Survival[!duplicated(human_data$Patient)] == "NS") - as.numeric(dev_mets_out_all_pie$PatCount_UK)
dev_mets_out_all_pie$IPatCount_Ferrario <- sum(vd$Survival28[!duplicated(vd$Patient)] == "NS") - as.numeric(dev_mets_out_all_pie$PatCount_Ferrario)
dev_mets_out_all_pie$rFerrario <- 1
dev_mets_out_all_pie$rFerrario[dev_mets_out_all_pie$PatCount_Ferrario == "-"] <- 0
dev_mets_out_all_pie[is.na(dev_mets_out_all_pie)] <- 0
dev_mets_out_all_pie[dev_mets_out_all_pie == "-"] <- 1
dev_mets_out_all_pie$PatCount_UK <- as.numeric(dev_mets_out_all_pie$PatCount_UK)
dev_mets_out_all_pie$PatCount_Ferrario <- as.numeric(dev_mets_out_all_pie$PatCount_Ferrario)
dev_mets_out_all_pie$Pos <- 1:nrow(dev_mets_out_all_pie)
dev_mets_out_all_pie_m <- melt(dev_mets_out_all_pie, id.vars = c("Name", "Pos", "Group", "rFerrario"))
dev_mets_out_all_pie_m$source <- stri_extract(str = dev_mets_out_all_pie_m$variable, regex = "UK|Ferrario")
dev_mets_out_all_pie_m$variable <- as.character(dev_mets_out_all_pie_m$variable)
dev_mets_out_all_pie_m$variable[stri_startswith(str = dev_mets_out_all_pie_m$variable, fixed = "I")] <- "Nondeviating"
dev_mets_out_all_pie_m$variable[stri_startswith(str = dev_mets_out_all_pie_m$variable, fixed = "P")] <- "Deviating"
dev_mets_out_all_pie_m <- dcast(data = dev_mets_out_all_pie_m, formula = Pos + Name + Group + rFerrario + source ~ variable)
dev_mets_out_all_pie_m$y <- (1:2)[1 + stri_detect(str = dev_mets_out_all_pie_m$source, fixed = "UK")]
dev_mets_out_all_pie_m$rFerrario[dev_mets_out_all_pie_m$source == "UK"] <- 1
p_dev_metabs_pats <- ggplot(data = dev_mets_out_all_pie_m) +
  geom_scatterpie(mapping = aes(x = Pos, y = y, r = rFerrario/4), data = dev_mets_out_all_pie_m, cols = c("Deviating", "Nondeviating"), legend_name = "Fraction") + 
  coord_fixed() +
  xlab("") +
  ylab("") + 
  guides(type = guide_legend(title = "Fraction of Patients", ncol = 2)) +
  scale_y_continuous(labels = c("Our study", "Ferrario et al., 2016"), breaks = 2:1) +
  scale_x_continuous(labels = dev_mets_out_all_pie$Name, breaks = dev_mets_out_all_pie$Pos) +
  scale_fill_manual(breaks = c("Deviating", "Nondeviating"), values = c("black", "white")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.grid = element_blank(), panel.grid.major.y = element_line(color = "grey80"), legend.direction = "vertical", legend.position = "right", legend.background = element_blank())
ggsave(plot = p_dev_metabs_pats, filename = "safe_corridor_minmax_metabs_UK_Ferrario.png", path = out_dir, width = 9, height = 2.5, units = "in")
ggsave(plot = p_dev_metabs_pats, filename = "safe_corridor_minmax_metabs_UK_Ferrario.svg", path = out_dir, width = 9, height = 2.5, units = "in")

###Find out what number of deviating metabolites we find by chance
{
  set.seed(1014)
  m <- human_data[, -pheno_sel]
  m <- m[, !(colnames(m) %in% sig.anova.car.s.class)]
  di <- which(m$Survival == "NS")
  dev_rep_res <- sapply(X = 1:100, FUN = sim_dev_met, n = nrow(m), d = ncol(m) - 1, dev_idx = di, sample_groups = m$Survival)
}
sapply(dev_rep_res, `[[`)

##Compare C4, lysoPC a C28:0, -:1 and SM C22:3 S-vs-NS in Ferrario et al. visually
hsvd <- subset(human_sepsis_val_data, C4 < 300, select = c("Survival28", "Day", "C4", "lysoPC a C28:1", "SM C22:3"))
hsvd$Day <- paste0("Day ", hsvd$Day)
hsvd <- na.omit(hsvd)
hsvd_C4_tsig <- t.test(x = hsvd$C4[hsvd$Survival28 == "S" & hsvd$Day == "Day 0"], y = hsvd$C4[hsvd$Survival28 == "NS" & hsvd$Day == "Day 0"])
hsvd_lpc281_tsig <- t.test(x = hsvd$`lysoPC a C28:1`[hsvd$Survival28 == "S" & hsvd$Day == "Day 6"], y = hsvd$`lysoPC a C28:1`[hsvd$Survival28 == "NS" & hsvd$Day == "Day 6"])
hsvdp <- ggplot(data = melt(hsvd, id.vars = 1:2), mapping = aes(y = value, x = variable, colour = Survival28)) +
  facet_wrap( ~ Day) +
  geom_boxplot() + 
  ylab("Concentration, µM") +
  xlab("Metabolite") +
  scale_y_log10() +
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(plot = hsvdp, filename = "best_4_UK_feats_in_Ferrario.png", path = out_dir, width = 6, height = 6, units = "in")

##Human, metab concentration time course, only metabolites significant in control vs sepsis and with control overlapping S
keep_set <- c("ADMA", "Kynurenine", 
              paste0("lysoPC a C", c("16:0", "16:1", "17:0", "18:1", "18:2", "20:3", "20:4", "21:0")),
              paste0("PC aa C", c("40:5", "42:6")),
              paste0("PC ae C", c("36:4", "40:1", "40:4", "44:3")))
#hd <- human_data[, c(1:6, intersect(metab_sel, which(!colnames(human_data) %in% sig.anova.car.s.class)))]
hd <- human_data[, c(1:6, which(colnames(human_data) %in% keep_set))]
n_mets <- ncol(hd) - 6
ncols <- 5
p <- ggplot(data = subset(melt(hd, id.vars = 1:6), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Group, colour = Group)) + 
  facet_wrap(facets = ~ variable, ncol = ncols, nrow = ceiling(n_mets/ncols), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration, µM") +
  xlab("Day") +
  human_col_scale() +
  theme_bw()
ggsave(plot = p, filename = "human_metab_survival_nonsig_C_overlap_S.png", path = out_dir, width = 12/4*ncols, height = 0.3 + 1.5 * ceiling(n_mets/ncols), units = "in")

##Human, metab concentration time course, only metabolites significant in control vs sepsis and with control overlapping NS
keep_set <- c(paste0("PC aa C", c("26:0", "28:1", "36:5")),
              paste0("PC ae C", c("32:1", "34:0", "34:1", "40:2", "42:5")),
              "SM C24:1", "SM C26:1")
#hd <- human_data[, c(1:6, intersect(metab_sel, which(!colnames(human_data) %in% sig.anova.car.s.class)))]
hd <- human_data[, c(1:6, which(colnames(human_data) %in% keep_set))]
n_mets <- ncol(hd) - 6
ncols <- 5
p <- ggplot(data = subset(melt(hd, id.vars = 1:6), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Group, colour = Group)) + 
  facet_wrap(facets = ~ variable, ncol = ncols, nrow = ceiling(n_mets/ncols), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration, µM") +
  xlab("Day") +
  human_col_scale() +
  theme_bw()
ggsave(plot = p, filename = "human_metab_survival_nonsig_C_overlap_NS.png", path = out_dir, width = 12/4*ncols, height = 0.3 + 1.5 * ceiling(n_mets/ncols), units = "in")

##Human, metab concentration time course, only metabolites significant in control vs sepsis and with control overlapping NS, statistical approach
keep_set <- c_is_ns_sig
#hd <- human_data[, c(1:6, intersect(metab_sel, which(!colnames(human_data) %in% sig.anova.car.s.class)))]
hd <- human_data[, c(1:6, which(colnames(human_data) %in% keep_set))]
n_mets <- ncol(hd) - 6
ncols <- 5
p <- ggplot(data = subset(melt(hd, id.vars = 1:6), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Group, colour = Group)) + 
  facet_wrap(facets = ~ variable, ncol = ncols, nrow = ceiling(n_mets/ncols), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration, µM") +
  xlab("Day") +
  human_col_scale() +
  theme_bw()
ggsave(plot = p, filename = "human_metab_survival_nonsig_C_overlap_NS_sig.png", path = out_dir, width = 12/4*ncols, height = 0.3 + 1.5 * ceiling(n_mets/ncols), units = "in")

##Human, metabolite groups
###Given groups, Survival
n_mets <- ncol(human_sepsis_data_grouped[, -group_pheno_sel]) - 6
p <- ggplot(data = subset(melt(human_sepsis_data_grouped[, -group_pheno_sel], id.vars = 1:6), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Group, colour = Group)) + 
  facet_wrap(facets = ~ variable, ncol = 4, nrow = ceiling(n_mets/4), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration, µM") +
  xlab("Day") +
  human_col_scale() +
  theme_bw()
ggsave(plot = p, filename = paste0("human_metab_group_survival_time_course.png"), path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets/4), units = "in")

###Given groups, Contral vs Sepsis
hdg <- human_data_grouped[, -group_pheno_sel]
n_mets <- ncol(hdg) - 6
ncols <- 4
p <- ggplot(data = subset(melt(hdg, id.vars = 1:6), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Group, colour = Group)) + 
  facet_wrap(facets = ~ variable, ncol = ncols, nrow = ceiling(n_mets/ncols), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.3)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration, µM") +
  xlab("Day") +
  human_col_scale() +
  theme_bw()
ggsave(plot = p, filename = paste0("human_metab_group_control_time_course.png"), path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets/ncols), units = "in")

###Acylcarnitine made up groups
ac_transferase_grouped <- human_data[, c(1:6, which(human_sepsis_legend$group == "acylcarnitine") + 6)]
ac_transferase_grouped[, -1:-6] <- lapply(ac_transferase_grouped[, -1:-6], `/`, ac_transferase_grouped$C0)
ac_transferase_grouped <- transform(ac_transferase_grouped, 
                                    ShortChainAC = C2 + C3 + C4, 
                                    MediumChainAC = C5 + `C6 (C4:1-DC)` + C8 + C9 + C10 + C12, 
                                    LongChainAC = C14 + C16 + C18)
ac_transferase_grouped <- ac_transferase_grouped[, -which(colnames(ac_transferase_grouped) == "C0")]
n_mets <- ncol(ac_transferase_grouped) - 6
n_cols <- 6
p <- ggplot(data = subset(melt(ac_transferase_grouped, id.vars = 1:6), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Group, colour = Group)) + 
  facet_wrap(facets = ~ variable, ncol = n_cols, nrow = ceiling(n_mets/n_cols), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  human_col_scale() +
  xlab("Day") +
  theme_bw()
ggsave(plot = p, filename = paste0("human_metab_AC_group_time_course.png"), path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets/n_cols), units = "in")

###PCaa UFA grouped
ufa_grouped <- human_data[, c(1:6, which(human_sepsis_legend$group == "phosphatidylcholine") + 1)]
ufa_grouped
ufas <- unique(sub(pattern = "C[0-9]{2}", replacement = "C..", x = colnames(ufa_grouped)[-1:-6]))
ufa_sum <- as.data.frame(sapply(lapply(ufas, grep, x = colnames(ufa_grouped)), function(cols) rowSums(ufa_grouped[, cols])))
colnames(ufa_sum) <- sub(pattern = "C..:", replacement = "CXX:", x = ufas, fixed = TRUE)
ufa_sum <- cbind(ufa_grouped[, 1:6], ufa_sum)
n_mets <- ncol(ufa_sum) - 6
n_col <- 7
p <- ggplot(data = subset(melt(ufa_sum, id.vars = 1:6), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Group, colour = Group)) + 
  facet_wrap(facets = ~ variable, ncol = n_col, nrow = ceiling(n_mets/n_col), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.3)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration, µM") +
  xlab("Day") +
  human_col_scale() +
  theme_bw()
ggsave(plot = p, filename = paste0("human_metab_PC_group_time_course.png"), path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets/n_col), units = "in")

###PC UFA grouped, ratio of aa to ae
ufa_grouped <- human_data[, c(1:6, which(human_sepsis_legend$group == "phosphatidylcholine") + 1)]
ufas <- unique(sub(pattern = "C[0-9]{2}", replacement = "C..", x = colnames(ufa_grouped)[-1:-6]))
ufa_sum <- as.data.frame(sapply(lapply(ufas, grep, x = colnames(ufa_grouped)), function(cols) rowSums(ufa_grouped[, cols])))
ufa_ratio <- as.data.frame(sapply(1:7, function(col) ufa_sum[, col] / ufa_sum[, col + 7]))
colnames(ufa_ratio) <- sub(pattern = "PC aa C..:", replacement = "CXX:", x = ufas[1:7], fixed = TRUE)
ufa_ratio <- cbind(ufa_grouped[, 1:6], ufa_ratio)
n_mets <- ncol(ufa_ratio) - 6
n_col <- 7
p <- ggplot(data = subset(melt(ufa_ratio, id.vars = 1:6), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Group, colour = Group)) + 
  facet_wrap(facets = ~ variable, ncol = n_col, nrow = ceiling(n_mets/n_col), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Ratio PC aa/PCae") +
  xlab("Day") +
  human_col_scale() +
  theme_bw()
ggsave(plot = p, filename = paste0("human_metab_PC_rel_group_time_course.png"), path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets/n_col), units = "in")

###Fatty acid carrier ratios
fa_type_grouped <- human_data[, c(1:6, grep(pattern = "PC|SM", x = human_sepsis_legend[,1]) + 6)]
fa_type_grouped$TotalPCaa <- rowSums(fa_type_grouped[, grep(pattern = "PC aa", x = colnames(fa_type_grouped))])
fa_type_grouped$TotalPCae <- rowSums(fa_type_grouped[, grep(pattern = "PC ae", x = colnames(fa_type_grouped))])
fa_type_grouped$TotalSM <- rowSums(fa_type_grouped[, grep(pattern = "SM ", x = colnames(fa_type_grouped))])
fa_type_grouped$TotalLysoPC <- rowSums(fa_type_grouped[, grep(pattern = "lysoPC", x = colnames(fa_type_grouped))])
fa_type_grouped$PCaaToSM <- fa_type_grouped$TotalPCaa / fa_type_grouped$TotalSM
fa_type_grouped$PCaeToSM <- fa_type_grouped$TotalPCae / fa_type_grouped$TotalSM
fa_type_grouped$PCaeToPCaa <- fa_type_grouped$TotalPCaa / fa_type_grouped$TotalPCae
fa_type_grouped$lysoPCToPCaa <- fa_type_grouped$TotalLysoPC / fa_type_grouped$TotalPCaa
fa_type_grouped$lysoPCToPCae <- fa_type_grouped$TotalLysoPC / fa_type_grouped$TotalPCae
fa_type_grouped$lysoPCToSM <- fa_type_grouped$TotalLysoPC / fa_type_grouped$TotalSM
fa_type_grouped <- fa_type_grouped[, c(1:6, which(colnames(fa_type_grouped) == "PCaaToSM"):ncol(fa_type_grouped))]
n_mets <- ncol(fa_type_grouped)
p <- ggplot(data = subset(melt(fa_type_grouped, id.vars = 1:6), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Group, colour = Group)) + 
  facet_wrap(facets = ~ variable, ncol = 4, nrow = ceiling(n_mets/4), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration ratio") +
  xlab("Day") +
  human_col_scale() +
  theme_bw()
ggsave(plot = p, filename = paste0("human_metab_FA_carrier_group_time_course.png"), path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets/4), units = "in")

##Human, pheno vars time course, Control vs Sepsis
n_mets <- ncol(human_sepsis_data[, -metab_sel]) - 6
hd <- human_sepsis_data[, -metab_sel]
p <- ggplot(data = subset(melt(hd, id.vars = 1:6), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Group, colour = Group)) + 
  facet_wrap(facets = ~ variable, ncol = 4, nrow = ceiling(n_mets/4), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration, µM") +
  xlab("Day") +
  human_col_scale(drop = TRUE) +
  theme_bw()
ggsave(plot = p, filename = paste0("human_pheno_group_control_time_course.png"), path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets/4), units = "in")
