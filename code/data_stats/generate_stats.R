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
coarse_group_list <- human_sepsis_legend[human_sepsis_legend[,1] %in% colnames(human_sepsis_data)[-1:-5], 3] #lucky ... no col in data without match in legend
human_sepsis_data_grouped <- cbind(human_sepsis_data[,1:5], matrix(0, nrow = nrow(human_sepsis_data_normal), ncol=length(unique(coarse_group_list))))
colnames(human_sepsis_data_grouped)[-1:-5] <- unique(coarse_group_list)
human_sepsis_data_grouped$Patient <- as.factor(human_sepsis_data_grouped$Patient)
human_sepsis_data_grouped$Day <- as.factor(human_sepsis_data_grouped$Day)
###Split for metabolites and "phenotypical" factors
split_start <- which(colnames(human_sepsis_data) == "Urea")
pheno_sel <- split_start:ncol(human_sepsis_data)
metab_sel <- 6:(split_start-1)
group_pheno_sel <- which(colnames(human_sepsis_data_grouped) %in% unique(coarse_group_list[pheno_sel - 5]))
group_metab_sel <- which(colnames(human_sepsis_data_grouped) %in% unique(coarse_group_list[metab_sel - 5]))
hsd <- human_sepsis_data
hsd[, pheno_sel] <- scale(hsd[, pheno_sel])
for (n in 1:nrow(human_sepsis_data)){
  m_agg <- tapply(X = t(hsd[n, metab_sel]), INDEX = factor(coarse_group_list[metab_sel - 5]), FUN = sum)
  p_agg <- tapply(X = t(hsd[n, pheno_sel]), INDEX = factor(coarse_group_list[pheno_sel - 5]), FUN = mean)
  b_agg <- c(m_agg, p_agg)
  human_sepsis_data_grouped[n, -1:-5] <- b_agg[match(colnames(human_sepsis_data_grouped)[-1:-5], names(b_agg))]
}
human_sepsis_data_normal_grouped <- human_sepsis_data_grouped
human_sepsis_data_normal_grouped[, -1:-5] <- scale(human_sepsis_data_grouped[, -1:-5])
##Cluster plot, see outlier?
X11();plot(hclust(dist(human_sepsis_data_normal[, metab_sel])))
#human_sepsis_data <- human_sepsis_data[-11, ]

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
full_tanova_data$Survival[full_tanova_data$`CAP / FP` == "-"] <- "Nonsep"

##Run repeated measures ANOVA as depicted in the R Companion at http://rcompanion.org/handbook/I_09.html
fml <- concentration ~ Day*Survival
met_set <- colnames(full_tanova_data)[metab_sel]
names(met_set) <- met_set
ftd_c <- full_tanova_data
ftd_c$Survival <- sub(pattern = "NS|S", x = as.character(ftd_c$Survival), replacement = "Seps")
ftd_s <- subset(full_tanova_data, Survival != "Nonsep")
anova.car.c <- t3ANOVA(data = ftd_c, random = ~1|Patient, formula = fml, use.corAR = TRUE, col.set = met_set, id.vars = c("Survival", "Day", "Patient"))
anova.car.s <- t3ANOVA(data = ftd_s, random = ~1|Patient, formula = fml, use.corAR = TRUE, col.set = met_set, id.vars = c("Survival", "Day", "Patient"))
anova.car.c.pre.ps <- anova.car.c$ps
anova.car.s.pre.ps <- anova.car.s$ps
sig.anova.car.c.pre.class <- colnames(anova.car.c.pre.ps)[colAnys(anova.car.c.pre.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
sig.anova.car.s.pre.class <- colnames(anova.car.s.pre.ps)[colAnys(anova.car.s.pre.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
insig.anova.car.c.pre.class <- colnames(anova.car.c.pre.ps)[colAlls(anova.car.c.pre.ps[c("Survival", "Day:Survival"), ] > 0.05)]
insig.anova.car.s.pre.class <- colnames(anova.car.s.pre.ps)[colAlls(anova.car.s.pre.ps[c("Survival", "Day:Survival"), ] > 0.05)]
anova.car.c.ps <- (apply(anova.car.c.pre.ps, 2, function(row){ p.adjust(p = row, method = "fdr") }))
anova.car.s.ps <- (apply(anova.car.s.pre.ps, 2, function(row){ p.adjust(p = row, method = "fdr") }))
#anova.car.c.ps <- matrix(p.adjust(anova.car.c.ps[3:4, ], method = "fdr"), nrow = 2)
#anova.car.s.ps <- matrix(p.adjust(anova.car.s.ps[3:4, ], method = "fdr"), nrow = 2)
sig.anova.car.c.class <- colnames(anova.car.c.ps)[colAnys(anova.car.c.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
sig.anova.car.c.class <- intersect(sig.anova.car.c.class, sig.anova.car.c.pre.class)
sig.anova.car.s.class <- colnames(anova.car.s.ps)[colAnys(anova.car.s.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
sig.anova.car.s.class <- intersect(sig.anova.car.s.class, sig.anova.car.s.pre.class)
insig.anova.car.c.class <- colnames(anova.car.c.ps)[colAlls(anova.car.c.ps[c("Survival", "Day:Survival"), ] >= 0.95)]
insig.anova.car.c.class <- intersect(insig.anova.car.c.class, insig.anova.car.c.pre.class)
insig.anova.car.s.class <- colnames(anova.car.s.ps)[colAlls(anova.car.s.ps[c("Survival", "Day:Survival"), ] >= 0.95)]
insig.anova.car.s.class <- intersect(insig.anova.car.s.class, insig.anova.car.s.pre.class)

###Post hoc analysis of when the change happened
####Build models
fml <- concentration ~ DaySurv - 1 #without intercept because Anova(..., type = 3) gives an error with intercept
ftd_ph <- full_tanova_data
ftd_ph$DaySurv <- interaction(ftd_ph$Day, ftd_ph$Survival, drop = TRUE)
anova.car.ph.models <- lapply(met_set, fit_lin_mod_lme, data = ftd_ph, formula = fml, id.vars = c("Day", "Patient", "DaySurv"), random = ~1|Patient, use.corAR = TRUE)
anova.car.ph.models <- lapply(anova.car.ph.models, eval)
####Construct constrast matrix
ld <- length(tanova_day_set)
sd <- seq_along(tanova_day_set)
daySurv <- interaction(factor(full_tanova_data$Day), factor(full_tanova_data$Survival), drop = TRUE)
contr.m <- matrix(0, ncol = ld * 2, nrow = ld * 3)
colnames(contr.m) <- c(paste0("Seps vs Comp, D", tanova_day_set), paste0("S vs NS, D", tanova_day_set))
rownames(contr.m) <- levels(daySurv)
#####Sepsis vs Nonsepsis
contr.m[sd, sd] <- diag(-2, nrow = ld, ncol = ld)
contr.m[sd + ld, sd] <- diag(1, nrow = ld, ncol = ld)
contr.m[sd + 2 * ld, sd] <- diag(1, nrow = ld, ncol = ld)
#####Septic survivors vs septic nonsurvivors
contr.m[sd + ld, sd + ld] <- diag(-1, nrow = ld, ncol = ld)
contr.m[sd + 2 * ld, sd + ld] <- diag(1, nrow = ld, ncol = ld)
####Actual post hoc test
anova.car.ph.res <- lapply(anova.car.ph.models, function(e) summary(glht(e, linfct = mcp(DaySurv = t(contr.m)))))
anova.car.ph.ps <- sapply(anova.car.ph.res, function(e) e$test$pvalues)
anova.car.ph.sig.contr <- lapply(data.frame(anova.car.ph.ps), function(col) which(col <= 0.05))
names(anova.car.ph.ps) <- names(anova.car.ph.models)
names(anova.car.ph.sig.contr) <- names(anova.car.ph.models)

###Also on phenomenological variables
fml <- concentration ~ Day*Survival
met_set <- colnames(full_tanova_data)[pheno_sel]
names(met_set) <- met_set
ftd_c <- full_tanova_data
ftd_c$Survival <- sub(pattern = "NS|S", x = as.character(ftd_c$Survival), replacement = "Seps")
ftd_s <- subset(full_tanova_data, Survival != "Nonsep" & Day %in% tanova_day_set)
anova.car.s.pheno <- t3ANOVA(data = ftd_s, random = ~1|Patient, formula = fml, use.corAR = TRUE, col.set = met_set, id.vars = c("Survival", "Day", "Patient"))
anova.car.s.pheno.pre.ps <- anova.car.s.pheno$ps
sig.anova.car.s.pheno.pre.class <- colnames(anova.car.s.pheno.pre.ps)[colAnys(anova.car.s.pheno.pre.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
insig.anova.car.s.pheno.pre.class <- colnames(anova.car.s.pheno.pre.ps)[colAlls(anova.car.s.pheno.pre.ps[c("Survival", "Day:Survival"), ] > 0.05)]
anova.car.s.pheno.ps <- (apply(anova.car.s.pheno.pre.ps, 2, p.adjust, method = "fdr"))
sig.anova.car.s.pheno.class <- colnames(anova.car.s.pheno.ps)[colAnys(anova.car.s.pheno.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
sig.anova.car.s.pheno.class <- intersect(sig.anova.car.s.pheno.class, sig.anova.car.s.pheno.pre.class)
insig.anova.car.s.pheno.class <- colnames(anova.car.s.pheno.ps)[colAlls(anova.car.s.pheno.ps[c("Survival", "Day:Survival"), ] > 0.95)]
insig.anova.car.s.pheno.class <- intersect(insig.anova.car.s.pheno.class, insig.anova.car.s.pheno.pre.class)
anova.car.c.pheno <- t3ANOVA(data = ftd_c, random = ~1|Patient, formula = fml, use.corAR = TRUE, col.set = met_set, id.vars = c("Survival", "Day", "Patient"))
anova.car.c.pheno.pre.ps <- anova.car.c.pheno$ps
sig.anova.car.c.pheno.pre.class <- colnames(anova.car.c.pheno.pre.ps)[colAnys(anova.car.c.pheno.pre.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
anova.car.c.pheno.ps <- (apply(anova.car.c.pheno.pre.ps, 2, p.adjust, method = "fdr"))
sig.anova.car.c.pheno.class <- colnames(anova.car.c.pheno.ps)[colAnys(anova.car.c.pheno.ps[c("Survival", "Day:Survival"), ] <= 0.05)]
sig.anova.car.c.pheno.class <- intersect(sig.anova.car.c.pheno.class, sig.anova.car.c.pheno.pre.class)

###Post hoc analysis of when the change happened
####Build models, reuse met_set
fml <- concentration ~ DaySurv - 1
ftd_ph <- subset(full_tanova_data, Day %in% tanova_day_set)
ftd_ph$DaySurv <- interaction(ftd_ph$Day, ftd_ph$Survival, drop = TRUE)
anova.car.ph.pheno.models <- lapply(met_set, fit_lin_mod_lme, data = ftd_ph, formula = fml, id.vars = c("Day", "Patient", "DaySurv"), random = ~1|Patient, use.corAR = TRUE)
anova.car.ph.pheno.models <- lapply(anova.car.ph.pheno.models, eval)
####Construct constrast matrix
ld <- length(tanova_day_set)
sd <- seq_along(tanova_day_set)
daySurv <- interaction(factor(full_tanova_data$Day), factor(full_tanova_data$Survival), drop = TRUE)
contr.m <- matrix(0, ncol = ld * 2, nrow = ld * 3)
colnames(contr.m) <- c(paste0("Seps vs Comp, D", tanova_day_set), paste0("S vs NS, D", tanova_day_set))
rownames(contr.m) <- levels(daySurv)
#####Sepsis vs Nonsepsis
contr.m[sd, sd] <- diag(-2, nrow = ld, ncol = ld)
contr.m[sd + ld, sd] <- diag(1, nrow = ld, ncol = ld)
contr.m[sd + 2 * ld, sd] <- diag(1, nrow = ld, ncol = ld)
#####Septic survivors vs septic nonsurvivors
contr.m[sd + ld, sd + ld] <- diag(-1, nrow = ld, ncol = ld)
contr.m[sd + 2 * ld, sd + ld] <- diag(1, nrow = ld, ncol = ld)
####Actual post hoc test
anova.car.ph.pheno.res <- lapply(anova.car.ph.pheno.models, function(e) summary(glht(e, linfct = mcp(DaySurv = t(contr.m)))))
anova.car.ph.pheno.ps <- sapply(anova.car.ph.pheno.res, function(e) e$test$pvalues)
anova.car.ph.pheno.sig.contr <- lapply(data.frame(anova.car.ph.pheno.ps), function(col) which(col <= 0.05))
names(anova.car.ph.pheno.ps) <- names(anova.car.ph.pheno.models)
names(anova.car.ph.pheno.sig.contr) <- names(anova.car.ph.pheno.models)
#Save everything ANOVA to file
anova_complete_res <- grep(pattern = "anova.car", x = names(environment()), value = TRUE)
anova_complete_res <- grep(pattern = "models|res|terms|normality", x = anova_complete_res, invert = TRUE, value = TRUE)
save(list = anova_complete_res, file = paste0(out_dir, "ANOVA_complete_res.RData"))

#Build long format tables for plotting
##Scale measurement values at maximum concentration
human_sepsis_data_max_norm <- max_norm(x = human_sepsis_data, subset = -1:-5)
rat_sepsis_data_max_norm <- max_norm(x = rat_sepsis_data, subset = -1:-4)
##Melt scaled data into long form
human_sepsis_data_long_form <- melt(data = human_sepsis_data_max_norm, id.vars = c("Sample ID", "Patient", "Day", "Survival", "CAP / FP"))
rat_sepsis_data_long_form <- melt(data = rat_sepsis_data_max_norm, id.vars = c("Sample Identification", "material", "group", "time point"))

#Build tables for cluster-heatmaps
##Build covariance matrix with metabolite groups
human_sepsis_data_normal_grouped_metab_cov <- cbind(human_sepsis_data_normal_grouped[, 1:5], cov(t(human_sepsis_data_normal_grouped[,group_metab_sel]), use = "pairwise.complete.obs"))
human_sepsis_data_normal_grouped_pheno_cov <- cbind(human_sepsis_data_normal_grouped[, 1:5], cov(t(human_sepsis_data_normal_grouped[,group_pheno_sel]), use = "pairwise.complete.obs"))
human_sepsis_data_normal_grouped_pheno_cov <- human_sepsis_data_normal_grouped_pheno_cov[, !apply(human_sepsis_data_normal_grouped_pheno_cov, 2, function(x){all(is.na(x))})]
##Build covariance matrix with original metabolites
human_sepsis_data_normal_metab_cov <- cbind(human_sepsis_data_normal[, 1:5], cov(t(human_sepsis_data_normal[,metab_sel]), use = "pairwise.complete.obs"))
human_sepsis_data_normal_pheno_cov <- cbind(human_sepsis_data_normal[, 1:5], cov(t(human_sepsis_data_normal[,pheno_sel]), use = "pairwise.complete.obs"))
human_sepsis_data_normal_pheno_cov <- human_sepsis_data_normal_pheno_cov[, !apply(human_sepsis_data_normal_pheno_cov, 2, function(x){all(is.na(x))})]
##Build correlation matrix with metabolite groups
human_sepsis_data_normal_grouped_metab_cor <- cbind(human_sepsis_data_normal_grouped[, 1:5], cor(t(human_sepsis_data_normal_grouped[,group_metab_sel]), use = "pairwise.complete.obs"))
human_sepsis_data_normal_grouped_pheno_cor <- cbind(human_sepsis_data_normal_grouped[, 1:5], cor(t(human_sepsis_data_normal_grouped[,group_pheno_sel]), use = "pairwise.complete.obs"))
human_sepsis_data_normal_grouped_pheno_cor <- human_sepsis_data_normal_grouped_pheno_cor[, !apply(human_sepsis_data_normal_grouped_pheno_cor, 2, function(x){all(is.na(x))})]
##Build covariance matrix with original metabolites
human_sepsis_data_normal_metab_cor <- cbind(human_sepsis_data_normal[, 1:5], cor(t(human_sepsis_data_normal[, metab_sel]), use = "pairwise.complete.obs"))
human_sepsis_data_normal_pheno_cor <- cbind(human_sepsis_data_normal[, 1:5], cor(t(human_sepsis_data_normal[, pheno_sel]), use = "pairwise.complete.obs"))
human_sepsis_data_normal_pheno_cor <- human_sepsis_data_normal_pheno_cor[, !apply(human_sepsis_data_normal_pheno_cor, 2, function(x){all(is.na(x))})]
##Build covariance matrix with metabolite and phenom. var groups
human_sepsis_data_normal_grouped_metab_pheno_cov <- cbind(human_sepsis_data_normal_grouped[, 1:5], cov(t(human_sepsis_data_normal_grouped[, -1:-5]), use = "pairwise.complete.obs"))
##Build correlation matrix with metabolite and phenom. var groups
human_sepsis_data_normal_grouped_metab_pheno_cor <- cbind(human_sepsis_data_normal_grouped[, 1:5], cor(t(human_sepsis_data_normal_grouped[, -1:-5]), use = "pairwise.complete.obs"))
human_sepsis_data_normal_grouped_metab_pheno_cor <- human_sepsis_data_normal_grouped_metab_pheno_cor[, !apply(human_sepsis_data_normal_grouped_metab_pheno_cor, 2, function(x){all(is.na(x))})]
##Build covariance matrix with original metabolites and phenom. vars
human_sepsis_data_normal_metab_pheno_cov <- cbind(human_sepsis_data_normal[, 1:5], cov(t(human_sepsis_data_normal[, -1:-5]), use = "pairwise.complete.obs"))
##Build correlation matrix with original metabolites and phenom. vars
human_sepsis_data_normal_metab_pheno_cor <- cbind(human_sepsis_data_normal[, 1:5], cor(t(human_sepsis_data_normal[, -1:-5]), use = "pairwise.complete.obs"))
human_sepsis_data_normal_metab_pheno_cor <- human_sepsis_data_normal_metab_pheno_cor[, !apply(human_sepsis_data_normal_metab_pheno_cor, 2, function(x){all(is.na(x))})]

##Build covariance matrix of metabolites with original metabolites
###Survival-ignorant
human_sepsis_data_normal_conc_metab_cov <- cov(human_sepsis_data_normal[, metab_sel], use = "pairwise.complete.obs")
human_sepsis_data_normal_conc_pheno_cov <- cov(human_sepsis_data_normal[, pheno_sel], use = "pairwise.complete.obs")
###Survival-regardent
human_sepsis_data_normal_NS_conc_metab_cov <- cov(subset(human_sepsis_data_normal, Survival == "NS", metab_sel), use = "pairwise.complete.obs")
human_sepsis_data_normal_S_conc_metab_cov <- cov(subset(human_sepsis_data_normal, Survival == "S", metab_sel), use = "pairwise.complete.obs")
human_sepsis_data_normal_NS_conc_pheno_cov <- cov(subset(human_sepsis_data_normal, Survival == "NS", pheno_sel), use = "pairwise.complete.obs")
human_sepsis_data_normal_S_conc_pheno_cov <- cov(subset(human_sepsis_data_normal, Survival == "S", pheno_sel), use = "pairwise.complete.obs")
##Build covariance matrix of metabolites with metabolites groups
###Survival-ignorant
human_sepsis_data_normal_grouped_conc_metab_cov <- cov(human_sepsis_data_normal_grouped[, group_metab_sel], use = "pairwise.complete.obs")
human_sepsis_data_normal_grouped_conc_pheno_cov <- cov(human_sepsis_data_normal_grouped[, group_pheno_sel], use = "pairwise.complete.obs")
###Survival-regardent
human_sepsis_data_normal_NS_grouped_conc_metab_cov <- cov(subset(human_sepsis_data_normal_grouped, Survival == "NS", group_metab_sel), use = "pairwise.complete.obs")
human_sepsis_data_normal_S_grouped_conc_metab_cov <- cov(subset(human_sepsis_data_normal_grouped, Survival == "S", group_metab_sel), use = "pairwise.complete.obs")
human_sepsis_data_normal_NS_grouped_conc_pheno_cov <- cov(subset(human_sepsis_data_normal_grouped, Survival == "NS", group_pheno_sel), use = "pairwise.complete.obs")
human_sepsis_data_normal_S_grouped_conc_pheno_cov <- cov(subset(human_sepsis_data_normal_grouped, Survival == "S", group_pheno_sel), use = "pairwise.complete.obs")

##Build variance vectors of metabolites from normalized concentrations
###One method: variance of metabolite seperately for each patient
pat_first <- match(unique(human_sepsis_data_normal$Patient), human_sepsis_data_normal$Patient)
human_sepsis_data_normal_conc_var <- lapply(human_sepsis_data_normal[, -1:-5], aggregate, list(Patient = human_sepsis_data_normal$Patient), var, na.rm = T)
human_sepsis_data_normal_conc_var <- data.frame(Patient = human_sepsis_data_normal_conc_var[[1]]$Patient, 
                                                Survival = human_sepsis_data_normal$Survival[pat_first],
                                                lapply(human_sepsis_data_normal_conc_var, function(e) e$x))
human_sepsis_data_normal_grouped_conc_var <- lapply(human_sepsis_data_normal_grouped[, -1:-5], aggregate, list(Patient = human_sepsis_data_normal_grouped$Patient), var, na.rm = T)
human_sepsis_data_normal_grouped_conc_var <- data.frame(Patient = human_sepsis_data_normal_grouped_conc_var[[1]]$Patient, 
                                                        Survival = human_sepsis_data_normal_grouped$Survival[pat_first],
                                                        lapply(human_sepsis_data_normal_grouped_conc_var, function(e) e$x))

###Another method: variance of metabolite over all patients and days
human_sepsis_data_normal_patcentered <- human_sepsis_data_normal
human_sepsis_data_normal_patcentered[,-1:-5] <- Reduce("rbind", lapply(as.character(unique(human_sepsis_data_normal_patcentered$`Patient`)), function(x){ scale(subset(human_sepsis_data_normal_patcentered, Patient == x, -1:-5), center = T, scale = F) }))
human_sepsis_data_normal_NS_conc_metab_var <- colVars(as.matrix(subset(human_sepsis_data_normal_patcentered, Day %in% tanova_day_set & Survival == "NS", metab_sel)))
human_sepsis_data_normal_S_conc_metab_var <- colVars(as.matrix(subset(human_sepsis_data_normal_patcentered, Day %in% tanova_day_set & Survival == "S", metab_sel)))
human_sepsis_data_normal_NS_conc_pheno_var <- colVars(as.matrix(subset(human_sepsis_data_normal_patcentered, Day %in% tanova_day_set & Survival == "NS", pheno_sel)))
human_sepsis_data_normal_S_conc_pheno_var <- colVars(as.matrix(subset(human_sepsis_data_normal_patcentered, Day %in% tanova_day_set & Survival == "S", pheno_sel)))

###Another method: variance of metabolite over all patients seperately for each day
human_sepsis_data_normal_NS_conc_metab_day_var <- t(sapply(c(tanova_day_set, 5), function(d){ colVars(as.matrix(subset(human_sepsis_data_normal_patcentered, Survival == "NS" & Day == d, metab_sel))) }))
human_sepsis_data_normal_S_conc_metab_day_var <- t(sapply(c(tanova_day_set, 5), function(d){ colVars(as.matrix(subset(human_sepsis_data_normal_patcentered, Survival == "S" & Day == d, metab_sel))) }))
human_sepsis_data_normal_NS_conc_pheno_day_var <- t(sapply(c(tanova_day_set, 5), function(d){ colVars(as.matrix(subset(human_sepsis_data_normal_patcentered, Survival == "NS" & Day == d, pheno_sel))) }))
human_sepsis_data_normal_S_conc_pheno_day_var <- t(sapply(c(tanova_day_set, 5), function(d){ colVars(as.matrix(subset(human_sepsis_data_normal_patcentered, Survival == "S" & Day == d, pheno_sel))) }))

##Find significantly correlating concentrations
if (file.exists("generate_stats_bootstrap1000_result.RData")){
  load("generate_stats_bootstrap1000_result.RData")
}else{
  human_sepsis_data_normal_conc_metab_corr <- list()
  cols_grouped_metab <- colnames(human_sepsis_data)[-1:-5]
  cols_metab <- colnames(human_sepsis_data)[-1:-5]
  bootstrap_repeats <- 1000
  tic()
  for (d in 0:2){
    corr_dat <- list()
    n_NS <- sum(human_sepsis_data$Survival == "NS" & human_sepsis_data$Day == d)
    n_S <- sum(human_sepsis_data$Survival == "S" & human_sepsis_data$Day == d)
    NS_corr <- my.corr.test(x = subset(human_sepsis_data_normal, Survival == "NS" & Day == d, select = -1:-5), adjust = "fdr")
    NS_grouped_corr <- my.corr.test(x = subset(human_sepsis_data_normal_grouped, Survival == "NS" & Day == d, select = -1:-5), adjust = "fdr")
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
    bootstrap_list <- mclapply(X = 1:bootstrap_repeats, FUN = bootstrap_S_corr_fun, mc.cores = 7, n_S = n_S, corr_dat = corr_dat)
    human_sepsis_data_normal_conc_metab_corr[[paste0("Day",d)]] <- bootstrap_list
  }
  toc()
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
#Plot data
####################

##Human, Venn diagram of significant metabolite correlations for ungrouped metabolites
for (d in seq_along(human_sepsis_data_normal_conc_metab_corr)){
  png(filename = paste0(out_dir, "human_normal_metab_sig_pairs_Venn_day", d - 1, ".png"), width = 500, height = 500, units = "px")
  corr_dat <- human_sepsis_data_normal_conc_metab_corr[[d]]
  v1 <- length(corr_dat$NS_sig_pairs)
  v2 <- length(corr_dat$S_sig_pairs)
  vc <- length(intersect(corr_dat$NS_sig_pairs, corr_dat$S_sig_pairs))
  day_surv_table <- table(human_sepsis_data[c("Day", "Survival")])
  grid.newpage()
  g <- draw.pairwise.venn(area1 = v1, area2 = v2, cross.area = vc, category = c(paste0("Nonsurvivors, n=", day_surv_table[d, 1]) , paste0("Survivors, n=", day_surv_table[d, 2])), fill = c("blue", "red"), lwd = 0, cat.pos = c(0,0), scaled = F, cex = 2, cat.cex = 2, fontfamily = "arial", cat.fontfamily = "arial")
  grid.arrange(gTree(children = g), top = textGrob(paste0("Number of significantly correlating\nmetabolites at day ", d-1), gp = gpar(cex = 2, font = "arial")))
  dev.off()
}

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
legend_groups$group[pheno_sel - 1] <- "phenomenological variable"
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
legend_groups$group[pheno_sel - 1] <- "phenomenological variable"
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

##Human, cluster-heatmap, all metabolites, but groups at the side
x <- na.omit(data.frame(group = coarse_group_list[metab_sel - 5], human_sepsis_data_normal_conc_metab_cov));
heatmaply(x = x[, -1], row_side_colors = x[c("group")], key.title = "Cov", dendrogram = FALSE, showticklabels = FALSE, file = paste0(out_dir, "human_normal_metab_nonclust.png"))
##Human, cluster-heatmap, phenomenological vars, but groups at the side
x <- na.omit(data.frame(group = coarse_group_list[pheno_sel - 5], inv_human_sepsis_data_normal_conc_pheno_cov));
heatmaply(x = x[, -1], row_side_colors = x[c("group")], key.title = "pcor", dendrogram = FALSE, showticklabels = FALSE, file = paste0(out_dir, "human_normal_pheno_nonclust_pcor.png"))

x <- na.omit(data.frame(group = coarse_group_list[metab_sel - 5], inv_human_sepsis_data_normal_S_conc_metab_cov));
heatmaply(x = x[, -1], row_side_colors = x[c("group")], key.title = "pcor", dendrogram = FALSE, showticklabels = FALSE, file = paste0(out_dir, "human_normal_metab_S_nonclust_pcor.png"))
x <- na.omit(data.frame(group = coarse_group_list[metab_sel - 5], inv_human_sepsis_data_normal_NS_conc_metab_cov));
heatmaply(x = x[, -1], row_side_colors = x[c("group")], key.title = "pcor", dendrogram = FALSE, showticklabels = FALSE, file = paste0(out_dir, "human_normal_metab_NS_nonclust_pcor.png"))
x <- na.omit(data.frame(group = coarse_group_list[pheno_sel - 5], inv_human_sepsis_data_normal_S_conc_pheno_cov));
heatmaply(x = x[, -1], row_side_colors = x[c("group")], key.title = "pcor", dendrogram = FALSE, showticklabels = FALSE, file = paste0(out_dir, "human_normal_pheno_S_nonclust_pcor.png"))
x <- na.omit(data.frame(group = coarse_group_list[pheno_sel - 5], inv_human_sepsis_data_normal_NS_conc_pheno_cov));
heatmaply(x = x[, -1], row_side_colors = x[c("group")], key.title = "pcor", dendrogram = FALSE, showticklabels = FALSE, file = paste0(out_dir, "human_normal_pheno_NS_nonclust_pcor.png"))
rm("x")

##Human, cluster-heatmap, all phenom. vars, groups at the top
x <- na.omit(human_sepsis_data_normal[, c(1:5, pheno_sel)]);
heatmaply(x = x[, -1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = data.frame(Group = coarse_group_list[pheno_sel - 5]), plot_method = "plotly", margins = c(100,50,0,150))
rm("x")

##Human, cluster-heatmap, phenom. var groups
x <- na.omit(human_sepsis_data_normal_grouped[, c(1:5, group_pheno_sel)])
h <- heatmaply(x = x[, -1:-5], row_side_colors = x[c("Survival", "CAP / FP")], plot_method = "ggplot", margins = c(150,100,30,0), key.title = "Normalized\nConcentration", main = "Phenom. groups, survivors and nonsurvivors", subplot_heights = c(0.1, 0.9))
h$width <- 800
h$height <- 1100
export(p = h, file = paste0(out_dir, "human_normal_pheno_grouped.png"))
x <- na.omit(subset(human_sepsis_data_normal_grouped[, c(1:5, group_pheno_sel)], Survival == "S"))
h <- heatmaply(x = x[, -1:-5], row_side_colors = x[c("CAP / FP")], plot_method = "ggplot", margins = c(150,50,30,150), key.title = "Normalized\nConcentration", main = "Phenom. groups, survivors", subplot_heights = c(0.1, 0.9))
h$width <- 800
h$height <- 800
export(p = h, file = paste0(out_dir, "human_normal_pheno_S_grouped.png"))
x <- na.omit(subset(human_sepsis_data_normal_grouped[, c(1:5, group_pheno_sel)], Survival == "NS"))
h <- heatmaply(x = x[, -1:-5], row_side_colors = x[c("CAP / FP")], plot_method = "ggplot", margins = c(150,100,30,150), key.title = "Normalized\nConcentration", main = "Phenom. groups, nonsurvivors", subplot_heights = c(0.1, 0.9))
h$width <- 800
h$height <- 600
export(p = h, file = paste0(out_dir, "human_normal_pheno_NS_grouped.png"))
rm("x")

##Human, cluster-heatmap, metab groups
x <- na.omit(human_sepsis_data_normal_grouped[, c(1:5, group_metab_sel)])
h <- heatmaply(x = x[, -1:-5], row_side_colors = x[c("Survival", "CAP / FP")], plot_method = "ggplot", margins = c(150,100,30,0), key.title = "Normalized\nConcentration", main = "Metab. groups, survivors and nonsurvivors", subplot_heights = c(0.1, 0.9))
h$width <- 800
h$height <- 1100
export(p = h, file = paste0(out_dir, "human_normal_metab_grouped.png"))
x <- na.omit(subset(human_sepsis_data_normal_grouped[, c(1:5, group_metab_sel)], Survival == "S"))
h <- heatmaply(x = x[, -1:-5], row_side_colors = x[c("CAP / FP")], plot_method = "ggplot", margins = c(150,100,30,150), key.title = "Normalized\nConcentration", main = "Metab. groups, survivors", subplot_heights = c(0.1, 0.9))
h$width <- 800
h$height <- 800
export(p = h, file = paste0(out_dir, "human_normal_metab_S_grouped.png"))
x <- na.omit(subset(human_sepsis_data_normal_grouped[, c(1:5, group_metab_sel)], Survival == "NS"))
h <- heatmaply(x = x[, -1:-5], row_side_colors = x[c("CAP / FP")], plot_method = "ggplot", margins = c(150,100,30,150), key.title = "Normalized\nConcentration", main = "Metab. groups, nonsurvivors", subplot_heights = c(0.1, 0.9))
h$width <- 800
h$height <- 600
export(p = h, file = paste0(out_dir, "human_normal_metab_NS_grouped.png"))
rm("x")

##Human, cluster-heatmap, coarse grouped metabolites
x <- na.omit(human_sepsis_data_normal[,c(1:5, metab_sel)])
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], file = paste0(out_dir, "human_normal_metab.png"), main = "Metabolite profile does not cluster survival well", key.title = "Concentration\n(normalized)", showticklabels = TRUE)
x <- na.omit(human_sepsis_data_normal[,c(1:5, pheno_sel)])
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], file = paste0(out_dir, "human_normal_pheno.png"), main = "Phenomenological profile has survival clusters", key.title = "Concentration\n(normalized)", showticklabels = TRUE)
x <- na.omit(human_sepsis_data_normal_grouped[,c(1:5, group_metab_sel)])
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], file = paste0(out_dir, "human_normal_metab_grouped.png"), main = "Metablite group profiles somewhat cluster survival", key.title = "Concentration\n(normalized)", showticklabels = TRUE)
x <- na.omit(human_sepsis_data_normal_grouped[,c(1:5, group_pheno_sel)])
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], file = paste0(out_dir, "human_normal_pheno_grouped.png"), main = "Pheno var group profiles [?] survival", key.title = "Concentration\n(normalized)", showticklabels = TRUE)
rm("x")

##Human, cluster-heatmap, one per metabolite group
x <- subset(rbind(human_sepsis_data, human_nonsepsis_data), Day %in% tanova_day_set)
x <- max_norm(x, -1:-5)
x$Day <- as.numeric(as.character(x$Day))
x$Survival <- as.character(x$Survival)
x$Survival[x$`CAP / FP` == "-"] <- "Control"
x$Survival <- reorder(x$Survival, (x$Survival == "Control") + (2 * (x$Survival == "S")) + (3 * (x$Survival == "NS")))
x$Day <- reorder(x$Day, x$Day)
x <- x[order(x$Survival),]
x <- x[order(x$Day),]
xm <- x[, c(1:5, metab_sel)]
xmt <- data.frame(t(xm[, -1:-5]))
rownames(xmt) <- colnames(xm[,-1:-5])
colnames(xmt) <- xm$`Sample ID`

control.sig <- rownames(xmt[, -1:-5]) %in% sig.anova.car.c.class
survival.sig <- rownames(xmt[, -1:-5]) %in% sig.anova.car.s.class
mat_sigs <- data.frame(control.sig = control.sig, survival.sig = survival.sig, stringsAsFactors = FALSE)
mat_sigs <- lapply(mat_sigs, function(x){ c("nonsignif.", "p < 0.05")[x + 1] })
mat_sigs <- data.frame(mat_sigs)
colnames(mat_sigs) <- c("Control vs Sepsis", "S vs NS")

lower_margin <- 85
for (met_group in unique(coarse_group_list[metab_sel - 5])){
  group_sel <- coarse_group_list[metab_sel - 5] %in% met_group
  #xfplotdat <- t(max_norm(t(xmt[group_sel, xm$material == mat])))
  xfplotdat <- xmt[group_sel, ]
  sel <- !rowAlls(is.na(xfplotdat))
  top_row_h <- 0.03 * 76/sum(sel)
  subplot_h <- c(top_row_h, 1 - top_row_h)
  if (sum(sel) > 1){
    h <- heatmaply(x = xfplotdat[sel, ],
                   dendrogram = "row", 
                   plot_method = "plotly", 
                   col_side_colors = subset(xm, TRUE , c("Day", "Survival")), 
                   row_side_colors = mat_sigs[group_sel, ][sel, ], 
                   key.title = "concentration", 
                   margins = c(lower_margin,100,NA,50), 
                   height = lower_margin + round(sum(sel) * 1000/76),
                   subplot_heights = subplot_h, 
                   subplot_widths = c(0.91, 0.03, 0.06))
    h$width <- 2200
    #h$height <- 1200
    h$height <- lower_margin + round(sum(sel) * 1100/76) #1200 is a good height for 76 rows of metabolites
    export(p = h, file = paste0(out_dir, "human_heatmap_", met_group, ".png"))
  }
  else
    print(paste0(met_group, " has too few metabolites, possibly after filtering all-NA rows"))
}

#Also for pheno vars
xm <- x[, c(1:5, pheno_sel)]
xmt <- data.frame(t(xm[, -1:-5]))
rownames(xmt) <- colnames(xm[,-1:-5])
colnames(xmt) <- xm$`Sample ID`
control.sig <- rownames(xmt[, -1:-5]) %in% sig.anova.car.c.pheno.class
survival.sig <- rownames(xmt[, -1:-5]) %in% sig.anova.car.s.pheno.class
mat_sigs <- data.frame(control.sig = control.sig, survival.sig = survival.sig, stringsAsFactors = FALSE)
mat_sigs <- lapply(mat_sigs, function(x){ c("nonsignif.", "p < 0.05")[x + 1] })
mat_sigs <- data.frame(mat_sigs)
colnames(mat_sigs) <- c("Control vs Sepsis", "S vs NS")

lower_margin <- 85
xfplotdat <- xmt
sel <- !rowAlls(is.na(xfplotdat))
top_row_h <- 0.03 * 76/sum(sel)
subplot_h <- c(top_row_h, 1 - top_row_h)
h <- heatmaply(x = xfplotdat[sel, ],
               dendrogram = "row", 
               plot_method = "plotly", 
               col_side_colors = subset(xm, TRUE , c("Day", "Survival")), 
               row_side_colors = mat_sigs[sel, ], 
               key.title = "concentration", 
               margins = c(lower_margin,100,NA,50), 
               height = lower_margin + round(sum(sel) * 1000/76),
               subplot_heights = subplot_h, 
               subplot_widths = c(0.91, 0.03, 0.06))
h$width <- 2200
h$height <- lower_margin + round(sum(sel) * 1100/76) #1200 is a good height for 76 rows of metabolites
export(p = h, file = paste0(out_dir, "human_heatmap_pheno.png"))

##Human, variance of metabolites, NS vs. S
metab_var_df <- human_sepsis_data_normal_conc_var[,c(2, metab_sel - 3)]
colnames(metab_var_df)[-1] <- colnames(human_sepsis_data)[metab_sel]
metab_var_long_df <- melt(metab_var_df, id.vars = "Survival")
metab_var_long_df$group <- human_sepsis_legend$group[match(metab_var_long_df$variable, human_sepsis_legend[,1])]
metab_var_long_df$value <- log10(metab_var_long_df$value)
metab_var_long_df <- metab_var_long_df[is.finite(metab_var_long_df$value),]
metab_var_long_df <- subset(metab_var_long_df, variable != "sugar")

##Do F-test for variance difference
metab_var_sig_df <- data.frame(group = setdiff(unique(metab_var_long_df$group), "sugar"), sig = 1, value = 0, stringsAsFactors = FALSE)
for (gr in metab_var_sig_df$group){
  b1 <- subset(metab_var_long_df, group == gr & Survival == "S", value)
  b2 <- subset(metab_var_long_df, group == gr & Survival == "NS", value)
  if (nrow(b1) > 1 & nrow(b2) > 1)
    metab_var_sig_df$sig[metab_var_sig_df$group == gr] <- t.test(x = b1, y = b2)$p.value
}
metab_var_sig_df$sig <- p.adjust(p = metab_var_sig_df$sig, method = "fdr")
metab_var_sig_df$value = 2
metab_var_sig_df <- subset(metab_var_sig_df, sig <= 0.05)
  
metab_var_plot <- ggplot(data = metab_var_long_df, mapping = aes(x = group, y = value, fill = Survival)) +
  geom_boxplot() +
  geom_point(data = metab_var_sig_df, mapping = aes(x = group, y = value), inherit.aes = FALSE, shape = 8, size = 1) +
  ylab("variance (log)") + 
  #scale_y_log10() +
  ggtitle("Concentration variance over all days\ndiffers with regard to survival") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = metab_var_plot, filename = "human_all_days_grouped_metab_var_patientcentered.png", path = out_dir, width = 7, height = 4, units = "in")

##Human, variance of phenomenological vars, NS vs. S
pheno_var_df <- na.omit(human_sepsis_data_normal_conc_var[,c(2, pheno_sel- 3)])
colnames(pheno_var_df)[-1] <- colnames(human_sepsis_data)[pheno_sel]
colnames(pheno_var_df)[colnames(pheno_var_df) == "Creatinine.1"] <- "Creatinine"
pheno_var_long_df <- melt(pheno_var_df, id.vars = "Survival")
pheno_var_long_df$group <- human_sepsis_pheno_var_groups$group[match(pheno_var_long_df$variable, human_sepsis_pheno_var_groups[,1])]
pheno_var_long_df$value <- log10(pheno_var_long_df$value)
pheno_var_long_df <- pheno_var_long_df[is.finite(pheno_var_long_df$value),]

##Do F-test for variance difference
pheno_var_sig_df <- data.frame(group = unique(pheno_var_long_df$group), sig = 1, value = 0, stringsAsFactors = FALSE)
for (gr in pheno_var_sig_df$group){
  b1 <- subset(pheno_var_long_df, group == gr & Survival == "S", value)
  b2 <- subset(pheno_var_long_df, group == gr & Survival == "NS", value)
  if (nrow(b1) > 1 & nrow(b2) > 1)
    pheno_var_sig_df$sig[pheno_var_sig_df$group == gr] <- t.test(x = b1, y = b2)$p.value
}
pheno_var_sig_df$sig <- p.adjust(p = pheno_var_sig_df$sig, method = "fdr")
pheno_var_sig_df$value = 4
pheno_var_sig_df <- subset(pheno_var_sig_df, sig <= 0.05)

pheno_var_plot <- ggplot(data = pheno_var_long_df, mapping = aes(x = group, y = value, fill = Survival)) +
  geom_boxplot() +
  geom_point(data = pheno_var_sig_df, mapping = aes(x = group, y = value), inherit.aes = FALSE, shape = 8, size = 1) +
  #scale_y_log10() +
  ylab("variance (log)") + 
  ggtitle("Concentration variance over all days\ndiffers with regard to survival") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = pheno_var_plot, filename = "human_all_days_grouped_pheno_var_patientcentered.png", path = out_dir, width = 9, height = 4, units = "in")

##Human, variance of metab vars over all days but metabolite over all patients
pheno_day_var_df <- rbind(cbind(data.frame(Survival = "NS"), t(human_sepsis_data_normal_NS_conc_pheno_var)), 
                          cbind(data.frame(Survival = "S"), t(human_sepsis_data_normal_S_conc_pheno_var)))
colnames(pheno_day_var_df) <- c("Survival", colnames(human_sepsis_data)[pheno_sel])
pheno_day_var_long_df <- melt(pheno_day_var_df, id.vars = c("Survival"))
pheno_day_var_long_df$group <- human_sepsis_legend$group[match(pheno_day_var_long_df$variable, human_sepsis_legend[[1]])]
pheno_day_var_plot <- ggplot(data = pheno_day_var_long_df, mapping = aes(x = as.factor(group), y = value, fill = Survival)) +
  geom_boxplot(outlier.size = 0.7) +
  ylab("variance") +
  xlab("clinical parameter") +
  #ggtitle("Patient-wise concentration variances differ") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "human_all_days_grouped_pheno_var.png", path = out_dir, plot = pheno_day_var_plot, width = 9, height = 3.5, units = "in")

##Human, variance of metab vars seperate by day
metab_day_var_df <- rbind(cbind(data.frame(Survival = "NS"), t(human_sepsis_data_normal_NS_conc_metab_var)), 
                          cbind(data.frame(Survival = "S"), t(human_sepsis_data_normal_S_conc_metab_var)))
colnames(metab_day_var_df) <- c("Survival", colnames(human_sepsis_data)[metab_sel])
metab_day_var_long_df <- melt(metab_day_var_df, id.vars = c("Survival"))
metab_day_var_long_df$group <- human_sepsis_legend$group[match(metab_day_var_long_df$variable, human_sepsis_legend[[1]])]
metab_day_var_plot <- ggplot(data = metab_day_var_long_df, mapping = aes(x = group, y = value, fill = Survival)) +
  geom_boxplot(outlier.size = 0.7) +
  ylab("variance") +
  xlab("metabolite group") +
  #ggtitle("Patient-wise concentration variances differ") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "human_all_days_grouped_metab_var.png", path = out_dir, plot = metab_day_var_plot, width = 8, height = 3.5, units = "in")

##Human, variance of metab vars seperate by day
pheno_day_var_df <- rbind(cbind(data.frame(Day = c(tanova_day_set, 5), Survival = "NS"), human_sepsis_data_normal_NS_conc_pheno_day_var), 
                          cbind(data.frame(Day = c(tanova_day_set, 5), Survival = "S"), human_sepsis_data_normal_S_conc_pheno_day_var))
colnames(pheno_day_var_df) <- c("Day", "Survival", colnames(human_sepsis_data)[pheno_sel])
pheno_day_var_long_df <- melt(pheno_day_var_df, id.vars = c("Day", "Survival"))
pheno_day_var_long_df$group <- human_sepsis_legend$group[match(pheno_day_var_long_df$variable, human_sepsis_legend[[1]])]
pheno_day_var_plot <- ggplot(data = pheno_day_var_long_df, mapping = aes(x = factor(Day), y = value, fill = Survival)) +
  facet_wrap( ~ group, ncol = 5, nrow = 2, scales = "free_y") +
  geom_boxplot(outlier.size = 0.7) +
  ylab("variance") +
  xlab("Day") +
  #ggtitle("Patient-wise concentration variances differ between days") +
  theme_bw()
ggsave(filename = "human_seperate_days_grouped_pheno_var.png", path = out_dir, plot = pheno_day_var_plot, width = 9, height = 3.5, units = "in")

##Human, variance of metab vars seperate by day
metab_day_var_df <- rbind(cbind(data.frame(Day = c(tanova_day_set, 5), Survival = "NS"), human_sepsis_data_normal_NS_conc_metab_day_var), 
                          cbind(data.frame(Day = c(tanova_day_set, 5), Survival = "S"), human_sepsis_data_normal_S_conc_metab_day_var))
colnames(metab_day_var_df) <- c("Day", "Survival", colnames(human_sepsis_data)[metab_sel])
metab_day_var_long_df <- melt(metab_day_var_df, id.vars = c("Day", "Survival"))
metab_day_var_long_df$group <- human_sepsis_legend$group[match(metab_day_var_long_df$variable, human_sepsis_legend[[1]])]
metab_day_var_plot <- ggplot(data = metab_day_var_long_df, mapping = aes(x = factor(Day), y = value, fill = Survival)) +
  facet_wrap( ~ group, ncol = 4, nrow = 2, scales = "free_y") +
  geom_boxplot(outlier.size = 0.7) +
  ylab("variance") +
  xlab("Day") +
  #ggtitle("Metabolite concentration variances differ between days") +
  theme_bw()
ggsave(filename = "human_seperate_days_grouped_metab_var.png", path = out_dir, plot = metab_day_var_plot, width = 8, height = 3.5, units = "in")


##Human, covariance cluster-heatmap, metabolites
x <- human_sepsis_data_normal_metab_cov
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_metab_cov.png"), main = "Metabolite profile covariance has mainly patient clusters", key.title = "Cov", showticklabels = FALSE)
x <- na.omit(human_sepsis_data_normal_pheno_cov)
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_pheno_cov.png"), main = "Phenomenological profile covariance has patient\n and survival clusters", key.title = "Cov", showticklabels = FALSE)
x <- na.omit(human_sepsis_data_normal_metab_pheno_cov)
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_metab_and_pheno_cov.png"), main = "", key.title = "Cov", showticklabels = FALSE, margins = c(100,50,0,150), plot_method = "ggplot")
x <- na.omit(human_sepsis_data_normal_grouped_metab_pheno_cov)
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_grouped_metab_and_pheno_cov.png"), main = "", key.title = "Cov", showticklabels = FALSE, margins = c(100,50,0,150), plot_method = "ggplot")
x <- na.omit(human_sepsis_data_normal_metab_pheno_cor)
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_metab_and_pheno_cor.png"), main = "", key.title = "Cor", showticklabels = FALSE, margins = c(100,50,0,150), plot_method = "ggplot")
x <- na.omit(human_sepsis_data_normal_grouped_metab_pheno_cor)
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_grouped_metab_and_pheno_cor.png"), main = "", key.title = "Cor", showticklabels = FALSE, margins = c(100,50,0,150), plot_method = "ggplot")
rm("x")

#TODO: switch to similarity for patient signature visualization
x <- as.matrix(na.omit(subset(human_sepsis_data_normal, TRUE, select = pheno_sel)))
x <- cbind(human_sepsis_data_normal[,1:5], tcrossprod(x))
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_metab_and_pheno_cov.png"), main = "", key.title = "Dist", showticklabels = FALSE, margins = c(100,50,0,150), plot_method = "ggplot")

##Human, cluster-heatmap of non-CAP and non-FP patients, coarse grouped metabolites; error because nothing es left
#heatmaply(x = subset(human_sepsis_data_normal_grouped, subset = human_sepsis_data_normal_grouped$`CAP / FP` %in% c("-"), select = c(-1,-2,-3,-5)))

##Human, cluster-heatmap of patient covariance matrix, coarse grouped metabolites, survival marked
x <- human_sepsis_data_normal_grouped_metab_cov
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_grouped_metab_cov.png"))
x <- na.omit(human_sepsis_data_normal_grouped_pheno_cov)
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_grouped_pheno_cov.png"))
rm("x")

##Human, cluster-heatmap of patient correlation matrix, ungrouped metabolites, survival marked
x <- human_sepsis_data_normal_metab_cor
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_metab_cor.png"), margin = c(100,50,0,150), key.title = "Cor", showticklabels = FALSE)
x <- na.omit(human_sepsis_data_normal_pheno_cor)
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_pheno_cor.png"), main = "Phenomenological profile correlation gives patient\n and survival clusters")
rm("x")

##Human, cluster-heatmap of sample distance matrix, ungrouped metabolites, then ungrouped pheno vars, survival marked
x <- human_sepsis_data_normal
heatmaply(x = as.matrix(dist(x[, metab_sel])), row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_metab_dist.png"), margin = c(100,50,0,150), key.title = "Euclidean\ndistance", showticklabels = FALSE)
heatmaply(x = as.matrix(dist(x[, pheno_sel])), row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_pheno_dist.png"), margin = c(100,50,0,150), key.title = "Euclidean\ndistance", showticklabels = FALSE)
rm(x)

##Human, cluster-heatmat of patient correlation matrix, coarse grouped metabolites, survival and CAP/FP marked
x <- human_sepsis_data_normal_grouped_metab_cor
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_grouped_metab_cor.png"), margins = c(100, 50, 0, 150), showticklabels = F, k_row = 2, key.title = "Cor")
x <- na.omit(human_sepsis_data_normal_grouped_pheno_cor)
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_grouped_pheno_cor.png"), main = "Profiles of grouped phenomenological variables cluster nothing")
rm("x")

##Human, cluster-heatmap of patient covariance matrix, coarse grouped everything, survival-regardent
x <- human_sepsis_data_normal
for (group in setdiff(unique(coarse_group_list), "sugar")){
  col_sel <- which(coarse_group_list == group)
  h <- heatmaply(x = x[, col_sel + 5], row_side_colors = x[c("Survival", "CAP / FP")], key.title = "Normlized\nConcentrations", margins = c(100,80,0,250), k_row = 2, subplot_heights = c(0.08, 0.92))
  h$width <- 1600
  h$height <- 1200
  export(p = h, file = paste0(out_dir, "human_normal_conc_", group, ".png"))
}
rm("x")

##Human, cluster-heatmap of metabolite covariance matrix, survival-ignorant
x <- human_sepsis_data_normal_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_conc_metab_cov.png"), key.title = "Cov")
x <- human_sepsis_data_normal_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_conc_pheno_cov.png"), key.title = "Cov")
rm("x")

##Human, cluster-heatmap of metabolite covariance matrix, survival-regardent
x <- human_sepsis_data_normal_NS_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_NS_conc_metab_cov.png"), key.title = "Cov")
x <- human_sepsis_data_normal_S_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_S_conc_metab_cov.png"), key.title = "Cov")
x <- human_sepsis_data_normal_NS_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_NS_conc_pheno_cov.png"), key.title = "Cov")
x <- human_sepsis_data_normal_S_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_S_conc_pheno_cov.png"), key.title = "Cov")
rm("x")

##Human, cluster-heatmap of covariance matrix of grouped metabolites, survival-ignorant
x <- human_sepsis_data_normal_grouped_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_grouped_conc_metab_cov.png"))
x <- human_sepsis_data_normal_grouped_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_grouped_conc_pheno_cov.png"))
rm("x")

##Human, cluster-heatmap of covariance matrix of grouped metabolites, survival-regardent
x <- human_sepsis_data_normal_NS_grouped_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_NS_grouped_conc_metab_cov.png"), key.title = "Cov", margins = c(200,200,0,150))
x <- human_sepsis_data_normal_S_grouped_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_S_grouped_conc_metab_cov.png"), key.title = "Cov", margins = c(200,200,0,150))
x <- human_sepsis_data_normal_NS_grouped_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_NS_grouped_conc_pheno_cov.png"), key.title = "Cov", margins = c(200,200,0,150))
x <- human_sepsis_data_normal_S_grouped_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_S_grouped_conc_pheno_cov.png"), key.title = "Cov", margins = c(200,200,0,150))
rm("x")

##Human, p-val (t-test) plot of differences between non-survivors and survivors, grouped by day
day_sig_t_diff$method <- "t-test"
day_sig_u_diff$method <- "U-test"
day_sig_ut_dat <- rbind(melt(day_sig_t_diff, id.vars = c("Day", "method")), melt(day_sig_u_diff, id.vars = c("Day", "method")))
h_day_sig_u_t_diff_plot <- ggplot(day_sig_ut_dat, aes(x = Day, y = value)) +
  facet_grid(method ~ .) +
  geom_point(position = position_jitter(width = 0.1), size = 0.7) +
  geom_hline(yintercept = 0.05) +
  scale_y_log10() +
  ylab("p value, FDR-corrected") +
  ggtitle("U-test and t-test give similar results\n- metabolite diffs (non-)survivors") +
  theme_bw()
ggsave(plot = h_day_sig_u_t_diff_plot, filename = "human_all_days_survival_sig_diff.png", path = out_dir, width = 4, height = 4, units = "in")

##Human, metabolite concentration time course, only metabolites with p-val < 0.05 at any day
###Prepare significant diff data
##Reduce to significantly different metabolites
human_sepsis_data_long_form_sig <- subset(human_sepsis_data_long_form, subset = as.character(variable) %in% sig_t_class)
human_sepsis_data_long_form_sig$variable <- factor(as.character(human_sepsis_data_long_form_sig$variable), levels = sig_t_class, ordered = TRUE)
h_time_course_sig_diff_dat <- subset(na.omit(human_sepsis_data_long_form_sig), Day %in% c(0,1,2,3,5,7,14,21,28))
h_time_course_sigs <- subset(day_sig_t_diff_pos_long, variable %in% h_time_course_sig_diff_dat$variable & Day %in% c(0,1,2,3,5,7,14,21,28))
h_time_course_sigs$value <- 0.9
h_time_course_group <- data.frame(variable = sort(unique(h_time_course_sig_diff_dat$variable)), value = 0.9, Day = 20, text = coarse_group_list[match(sort(unique(h_time_course_sig_diff_dat$variable)), colnames(human_sepsis_data)[-1:-5])])
h_time_course_sig_diff_plot <- ggplot(h_time_course_sig_diff_dat, aes(x = Day, y = value, group = Survival, color = Survival)) +
  facet_wrap(facets = ~ variable, ncol = 7, nrow = ceiling(length(sig_t_class)/7)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", size = 0.5) +
  geom_point(data = h_time_course_sigs, mapping = aes(x = Day, y = value), shape = 8, inherit.aes = FALSE, size = 0.8) +
  geom_text(data = h_time_course_group, mapping = aes(x = Day, y = value, label = text), inherit.aes = FALSE, size = 2.8) +
  ylab("Concentration relative to max value") +
  ggtitle("Metabolites significantly differing for survival by t-test at any time point") + 
  theme_bw()
ggsave(plot = h_time_course_sig_diff_plot, filename = "human_metab_time_course_t_sig_diff.png", path = out_dir, width = 14, height = 10, units = "in")

##Human, metabolite concentration time course, only metabolites with p-val < 0.05 and FDR < 0.05 in Survival or Day:Survival from type III repeated measures ANOVA
fs <- c("s", "c")
gs <- c("by Sepsis survival", "between Sepsis and Control")
ts <- list((length(tanova_day_set)+1):(length(tanova_day_set)*3), seq_along(tanova_day_set))
ss <- list(sig.anova.car.s.class, sig.anova.car.c.class)
for (n in seq_along(ss)[sapply(ss, length) > 0]){
  #h_time_course_sig_diff_dat <- melt(max_norm(full_tanova_data, subset = -1:-5), id.vars = c("Patient", "Survival", "Day", "Sample ID", "CAP / FP"))
  h_time_course_sig_diff_dat <- melt(full_tanova_data, id.vars = c("Patient", "Survival", "Day", "Sample ID", "CAP / FP"))
  if (n == 1){
    h_time_course_sig_diff_dat <- subset(h_time_course_sig_diff_dat, subset = Survival != "Nonsep")
    cms <- colMaxs(as.matrix(subset(full_tanova_data, !grepl(pattern = "Nonsep", x = Survival), -1:-5)))
  }
  else{
    h_time_course_sig_diff_dat$Survival <- factor(c("Sepsis", "Nonsepsis")[1 + (h_time_course_sig_diff_dat$Survival == "Nonsep")])
    cms <- colMaxs(as.matrix(full_tanova_data[, -1:-5]))
  }
  h_time_course_sig_diff_dat <- subset(h_time_course_sig_diff_dat, subset = as.character(variable) %in% ss[[n]])
  h_time_course_sig_diff_dat$variable <- factor(as.character(h_time_course_sig_diff_dat$variable), levels = ss[[n]], ordered = TRUE)
  n_mets <- length(unique(h_time_course_sig_diff_dat$variable))
  h_time_course_group <- data.frame(variable = sort(unique(h_time_course_sig_diff_dat$variable)), value = 1.3, Day = 3, text = coarse_group_list[match(sort(unique(h_time_course_sig_diff_dat$variable)), colnames(human_sepsis_data)[-1:-5])])
  h_time_course_group$value <- cms[match(h_time_course_group$variable, colnames(full_tanova_data[, -1:-5]))]
  h_time_course_group$value <- h_time_course_group$value * 1.3
  lv <- sapply(anova.car.ph.sig.contr, length)
  h_time_course_sig_times <- data.frame(t = "", variable = rep("", sum(lv)), text = "*", value = 1.1, stringsAsFactors = FALSE)
  h_time_course_sig_times$t <- Reduce("c", anova.car.ph.sig.contr)
  h_time_course_sig_times$variable <- factor(rep(names(anova.car.ph.sig.contr), times = lv))
  h_time_course_sig_times$Day <- rep(seq_along(tanova_day_set), times = 2)[h_time_course_sig_times$t]
  h_time_course_sig_times <- subset(h_time_course_sig_times, variable %in% h_time_course_sig_diff_dat$variable & t %in% ts[[n]])
  h_time_course_sig_times$value <- cms[match(h_time_course_sig_times$variable, colnames(full_tanova_data)[-1:-5])]
  h_time_course_sig_times$value <- h_time_course_sig_times$value * 1.1
  h_time_course_sig_diff_plot <- ggplot(h_time_course_sig_diff_dat, aes(x = factor(Day), y = value, group = Survival, color = Survival)) +
    facet_wrap(facets = ~ variable, ncol = 5, nrow = ceiling(n_mets/5), scales = "free_y") +
    geom_point(position = position_dodge(width = 0.2)) +
    #geom_boxplot(mapping = aes(x = factor(Day), color = Survival, y = value), inherit.aes = FALSE) + 
    stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
    geom_text(data = h_time_course_group, mapping = aes(x = Day, y = value, label = text), inherit.aes = FALSE, size = 2.8) +
    geom_point(data = h_time_course_sig_times, mapping = aes(x = Day, y = value), inherit.aes = FALSE, shape = 8, size = 2.5) +
    #ylim(c(-0.05, 1.4)) +
    #ylab("Concentration relative to max value") +
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
ss <- list(sig.anova.car.s.pheno.class, sig.anova.car.c.pheno.class)
for (n in seq_along(ss)[sapply(ss, length) > 0]){
  #h_time_course_sig_diff_dat <- melt(max_norm(full_tanova_data, subset = -1:-5), id.vars = c("Patient", "Survival", "Day", "Sample ID", "CAP / FP"))
  h_time_course_sig_diff_dat <- melt(full_tanova_data, id.vars = c("Patient", "Survival", "Day", "Sample ID", "CAP / FP"))
  if (n == 1){
    h_time_course_sig_diff_dat <- subset(h_time_course_sig_diff_dat, subset = Survival != "Nonsep" & Day %in% tanova_day_set)
    cms <- colMaxs(as.matrix(subset(full_tanova_data, !grepl(pattern = "Nonsep", x = Survival), -1:-5)))
  }else{
    h_time_course_sig_diff_dat$Survival <- factor(c("Sepsis", "Nonsepsis")[1 + (h_time_course_sig_diff_dat$Survival == "Nonsep")])
    cms <- colMaxs(as.matrix(full_tanova_data[, -1:-5]))
  }
  h_time_course_sig_diff_dat <- subset(h_time_course_sig_diff_dat, subset = as.character(variable) %in% ss[[n]])
  h_time_course_sig_diff_dat$variable <- factor(as.character(h_time_course_sig_diff_dat$variable), levels = ss[[n]], ordered = TRUE)
  n_mets <- length(unique(h_time_course_sig_diff_dat$variable))
  h_time_course_group <- data.frame(variable = sort(unique(h_time_course_sig_diff_dat$variable)), value = 1.3, Day = 3, text = coarse_group_list[match(sort(unique(h_time_course_sig_diff_dat$variable)), colnames(human_sepsis_data)[-1:-5])])
  h_time_course_group$value <- cms[match(h_time_course_group$variable, colnames(full_tanova_data[, -1:-5]))]
  h_time_course_group$value <- h_time_course_group$value * 1.3
  lv <- sapply(anova.car.ph.pheno.sig.contr, length)
  h_time_course_sig_times <- data.frame(t = "", variable = rep("", sum(lv)), text = "*", value = 1.1, stringsAsFactors = FALSE)
  h_time_course_sig_times$t <- Reduce("c", anova.car.ph.pheno.sig.contr)
  h_time_course_sig_times$variable <- factor(rep(names(anova.car.ph.pheno.sig.contr), times = lv))
  h_time_course_sig_times$Day <- rep(seq_along(tanova_day_set), times = 2)[h_time_course_sig_times$t]
  h_time_course_sig_times <- subset(h_time_course_sig_times, variable %in% h_time_course_sig_diff_dat$variable & t %in% ts[[n]])
  h_time_course_sig_times$value <- cms[match(h_time_course_sig_times$variable, colnames(full_tanova_data)[-1:-5])]
  h_time_course_sig_times$value <- h_time_course_sig_times$value * 1.1
  h_time_course_sig_diff_plot <- ggplot(h_time_course_sig_diff_dat, aes(x = factor(Day), y = value, group = Survival, color = Survival)) +
    facet_wrap(facets = ~ variable, ncol = 5, nrow = ceiling(n_mets/5), scales = "free_y") +
    geom_point(position = position_dodge(width = 0.2)) +
    #geom_boxplot(mapping = aes(x = factor(Day), color = Survival, y = value), inherit.aes = FALSE) + 
    stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
    geom_text(data = h_time_course_group, mapping = aes(x = Day, y = value, label = text), inherit.aes = FALSE, size = 2.8) +
    geom_point(data = h_time_course_sig_times, mapping = aes(x = Day, y = value), inherit.aes = FALSE, shape = 8, size = 2.5) +
    # ylim(c(-0.05, 1.4)) +
    # ylab("Concentration relative to max value") +
    ylab("Concentration, µM") + 
    xlab("Day") + 
    ggtitle(paste0("Phenomenological variables significantly differing ", gs[n], " by repeated measures ANOVA")) + 
    theme_bw()
  ggsave(plot = h_time_course_sig_diff_plot, filename = paste0("human_pheno_time_course_car_rm_anova_", fs[[n]], "_sig_diff.png"), path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets/5), units = "in")
}






##Human, metab concentration time course, only metabolites with p-val > 0.05 and FDR > 0.95 in Survival or Day:Survival from type III repeated measures ANOVA
n_mets <- length(insig.anova.car.s.class)
p <- ggplot(data = subset(melt(human_sepsis_data[, c(1:5, which(colnames(human_sepsis_data) %in% insig.anova.car.s.class))], id.vars = 1:5), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Survival, colour = Survival)) + 
  facet_wrap(facets = ~ variable, ncol = 4, nrow = ceiling(n_mets/4), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration, µM") +
  xlab("Day") +
  theme_bw()
plot(p)

##Human, metab concentration time course, only metabolites with p-val > 0.05 and FDR > 0.95 in Survival or Day:Survival from type III repeated measures ANOVA
n_mets <- length(insig.anova.car.c.class)
p <- ggplot(data = subset(melt(human_data[, c(1:5, which(colnames(human_data) %in% insig.anova.car.c.class))], id.vars = 1:5), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Survival, colour = Survival)) + 
  facet_wrap(facets = ~ variable, ncol = 4, nrow = ceiling(n_mets/4), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration, µM") +
  xlab("Day") +
  theme_bw()
plot(p)

##Human, metab concentration time course, only metabolites with p-val > 0.05 and FDR > 0.95 in Survival or Day:Survival from type III repeated measures ANOVA
n_mets <- length(insig.anova.car.s.pheno.pre.class)
p <- ggplot(data = subset(melt(human_sepsis_data[, c(1:5, which(colnames(human_sepsis_data) %in% insig.anova.car.s.pheno.pre.class))], id.vars = 1:5), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Survival, colour = Survival)) + 
  facet_wrap(facets = ~ variable, ncol = 4, nrow = ceiling(n_mets/4), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration, µM") +
  xlab("Day") +
  theme_bw()
plot(p)

##Human, metab concentration time course, only metabolites with p-val > 0.05 and FDR > 0.95 in Survival or Day:Survival from type III repeated measures ANOVA
n_mets <- length(insig.anova.car.c.pheno.pre.class)
p <- ggplot(data = subset(melt(human_data[, c(1:5, which(colnames(human_data) %in% insig.anova.car.c.pheno.pre.class))], id.vars = 1:5), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Survival, colour = Survival)) + 
  facet_wrap(facets = ~ variable, ncol = 4, nrow = ceiling(n_mets/4), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration, µM") +
  xlab("Day") +
  theme_bw()
plot(p)





##Human, metabolite groups
n_mets <- ncol(human_sepsis_data_grouped[, -group_pheno_sel]) - 5
p <- ggplot(data = subset(melt(human_sepsis_data_grouped[, -group_pheno_sel], id.vars = 1:5), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Survival, colour = Survival)) + 
  facet_wrap(facets = ~ variable, ncol = 4, nrow = ceiling(n_mets/4), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration, µM") +
  xlab("Day") +
  theme_bw()
ggsave(plot = p, filename = paste0("human_metab_group_time_course.png"), path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets/4), units = "in")

##Human, metabolites vs survival, p < 0.05, second day only, 
hp2 <- ggplot(data = subset(x = human_sepsis_data_long_form_sig, subset = Day == 0), mapping = aes(x = Survival, y = value)) + 
  facet_wrap(facets = ~ variable, nrow = 5, ncol = 10) +
  geom_boxplot(outlier.size = 0.5) + 
  ylab("Scaled concentration") +
  ggtitle("Human data from second day; metabolites differing by U-test with p < 0.05 after FDR correction") +
  theme_bw()
ggsave(plot = hp2, filename = "human_day0_metab_conc_vs_survival_sig.png", path = out_dir, width = 13, height = 7)

##Rat, metabolites vs survival, first measurement only
rp1 <- ggplot(data = subset(x = rat_sepsis_data_long_form, subset = `time point` == "6h"), mapping = aes(x = group, y = value)) +
  facet_wrap(facets = ~ variable, nrow = 17, ncol = 16) +
  geom_boxplot(outlier.size = 0.5) +
  ylab("Scaled concentration") +
  ggtitle("Rat data at 6h") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
ggsave(plot = rp1, filename = "rat_meas0_metab_conc_vs_survival.png", path = out_dir, width = 17, height = 17)

##Rat, metabolites vs tissue type, first measurement only, no control samples
rp2 <- ggplot(data = subset(x = rat_sepsis_data_long_form, subset = `time point` == "6h" & group != "control"), mapping = aes(x = material, y = value)) +
  facet_wrap(facets = ~ variable, nrow = 17, ncol = 16) +
  geom_boxplot(outlier.size = 0.5) + 
  ylab("scaled concentration") +
  ggtitle("Rat data at 6h without control samples") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
ggsave(plot = rp2, filename = "rat_meas0_noncontrol_metab_conc_vs_tissue_type.png", path = out_dir, width = 17, height = 17)
