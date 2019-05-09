#Load libraries
library(reshape2)
library(data.table)
library(ggplot2)
library(ggrepel)
library(matrixStats)
library(kernlab)
library(pairwiseCI)
library(heatmaply)
library(missRanger)
library(webshot)
library(htmlwidgets)
library(forcats)
library(car)
library(nlme)
library(multcomp)

#Import central functions
source("../function_definitions.R")

#Set path
out_dir <- "../../results/data_stats_rat_surv_vs_nonsurv/"

#Make sure paths exist
if (!dir.exists(out_dir))
  dir.create(out_dir)

###########################
#Import data
###########################

##Import clinical data
human_sepsis_data <- get_human_sepsis_data()

##Import corresponding group assignment
human_sepsis_legend <- get_human_sepsis_legend()
human_sepsis_legend$group[human_sepsis_legend$group == ""] <- human_sepsis_legend[human_sepsis_legend$group == "", 1]
human_sepsis_legend <- human_sepsis_legend[-1:-6, ]

##Import experiment data
rat_sepsis_data <- get_rat_sepsis_data()

###Remove outlier rats
rat_sepsis_data <- rat_sepsis_data[!rat_sepsis_data$`Sample Identification` %in% c("060H", "039L"),]

##Import corresponding group assignment
rat_sepsis_legend <- get_rat_sepsis_legend()
rat_sepsis_legend$group[rat_sepsis_legend$group == ""] <- rat_sepsis_legend[rat_sepsis_legend$group == "", 1]

##Import the metabolite names that Anna found significantly different
annas_mets <- get_annas_rat_sig_diffs()

###########################
#Process data
###########################

#Impute missing data
pheno_sel <- (which(colnames(rat_sepsis_data) == "H1")+1):ncol(rat_sepsis_data)
metab_sel <- 5:which(colnames(rat_sepsis_data) == "H1")
cols <- colnames(rat_sepsis_data)
colnames(rat_sepsis_data) <- make.names(cols)
for (mat in setdiff(unique(rat_sepsis_data$material), "plasma")){
  idat <- subset(rat_sepsis_data, material == mat, metab_sel)
  rat_sepsis_data[rat_sepsis_data$material == mat, metab_sel] <- missRanger(idat)
}
idat <- subset(rat_sepsis_data, material == "plasma", c(metab_sel, pheno_sel))
rat_sepsis_data[rat_sepsis_data$material == "plasma", c(metab_sel, pheno_sel)] <- missRanger(idat)
colnames(rat_sepsis_data) <- cols

#Remove dynamic physiological properties like heart rate
rat_sepsis_data <- rat_sepsis_data[, -which(colnames(rat_sepsis_data) %in% c("HR", "SV", "CO", "EF", "Resp Rate", "Temperature"))]#
pheno_sel <- (which(colnames(rat_sepsis_data) == "H1")+1):ncol(rat_sepsis_data)

#Calculate concentration ratios for survivors/non-survivors
rat_only_sepsis_data <- subset(rat_sepsis_data, group != "control" & `time point` != "72h", c(1:9, 11:30, 32:44))
rosdm <- na.omit(melt(rat_only_sepsis_data, id.vars = 1:4))
#rosdm$variable <- as.character(rosdm$variable)
rat_ratio_ci <- pairwiseCI(formula = value ~ group, by = c("material", "time point", "variable"), data = rosdm, method = "Param.ratio")
rrci <- melt(rat_ratio_ci$byout)
rrci <- subset(rrci, !L2 %in% c("method", "compnames"))
rrci[c("material", "time", "metabolite")] <- tstrsplit(rrci$L1, split = ".", fixed = TRUE)
rrci <- dcast(data = rrci, metabolite ~ material + time + L2)
#rrci <- rrci[, -which(colnames(rrci) == "L1")]
rat_ratio_ci <- rrci

#Get data overview
##Get overview of sample distribution along days
rat_plasma_sig_diff_res <- rat_sig_diffs_along_time(subset(rat_sepsis_data, material == "plasma" & group %in% c("septic survivor", "septic non-survivor")), corr_fdr = T)
rat_plasma_time_sig_u_diff <- rat_plasma_sig_diff_res$time_sig_u_diff
rat_plasma_time_sig_t_diff <- rat_plasma_sig_diff_res$time_sig_t_diff

rat_heart_sig_diff_res <- rat_sig_diffs_along_time(subset(rat_sepsis_data, material == "heart" & group %in% c("septic survivor", "septic non-survivor")), corr_fdr = T)
rat_heart_time_sig_u_diff <- rat_heart_sig_diff_res$time_sig_u_diff
rat_heart_time_sig_t_diff <- rat_heart_sig_diff_res$time_sig_t_diff

rat_liver_sig_diff_res <- rat_sig_diffs_along_time(subset(rat_sepsis_data, material == "liver" & group %in% c("septic survivor", "septic non-survivor")), corr_fdr = T)
rat_liver_time_sig_u_diff <- rat_liver_sig_diff_res$time_sig_u_diff
rat_liver_time_sig_t_diff <- rat_liver_sig_diff_res$time_sig_t_diff

#Get significant differences at all time points
cols <- colnames(rat_sepsis_data)[-1:-4]
rat_liver_sig_t_class <- na.omit(cols[colAnys(as.matrix(rat_liver_time_sig_t_diff[, -1]) <= 0.05)])
rat_plasma_sig_t_class <- na.omit(cols[colAnys(as.matrix(rat_plasma_time_sig_t_diff[, -1]) <= 0.05)])
rat_heart_sig_t_class <- na.omit(cols[colAnys(as.matrix(rat_heart_time_sig_t_diff[, -1]) <= 0.05)])

#Combine all differences
rat_sig_t_class <- union(union(rat_liver_sig_t_class, rat_plasma_sig_t_class), rat_heart_sig_t_class)

#Get sig var pos in time
rat_plasma_time_sig_t_diff_pos_long <- get_sig_var_pos(rat_plasma_time_sig_t_diff)
rat_liver_time_sig_t_diff_pos_long <- get_sig_var_pos(rat_liver_time_sig_t_diff)

#Get fold change with respect to control rats
rat_sepsis_data_fc <- rat_sepsis_data
for (mat in unique(rat_sepsis_data$material)){
  ctrl_subset <- rat_sepsis_data$group == "control" & rat_sepsis_data$material == mat
  surv_subset <- rat_sepsis_data$group == "septic survivor" & rat_sepsis_data$material == mat
  nonsurv_subset <- rat_sepsis_data$group == "septic non-survivor" & rat_sepsis_data$material == mat
  for (time in unique(rat_sepsis_data[surv_subset, ]$`time point`)){
    ctrl_col_mean <- colMeans(as.matrix(subset(rat_sepsis_data, ctrl_subset & `time point` == time, select = -1:-4)))
    rat_sepsis_data_fc[surv_subset & rat_sepsis_data$`time point` == time, -1:-4] <- t(log(t(subset(rat_sepsis_data, surv_subset & `time point` == time, select = -1:-4)) / ctrl_col_mean))
  }
  for (time in unique(rat_sepsis_data[nonsurv_subset, ]$`time point`)){
    ctrl_col_mean <- colMeans(as.matrix(subset(rat_sepsis_data, ctrl_subset & `time point` == time, select = -1:-4)))
    rat_sepsis_data_fc[nonsurv_subset & rat_sepsis_data_fc$`time point` == time, -1:-4] <- t(log(t(subset(rat_sepsis_data, nonsurv_subset & `time point` == time, select = -1:-4)) / ctrl_col_mean))
  }
}
rat_sepsis_data_fc <- subset(rat_sepsis_data_fc, group != "control")

rat_plasma_sig_diff_ctrl_surv_res <- rat_sig_diffs_along_time(subset(rat_sepsis_data, material == "plasma" & group %in% c("septic survivor", "control")), corr_fdr = T)
rat_plasma_time_fc_diff_ctrl_surv <- rat_plasma_sig_diff_ctrl_surv_res$time_fold_change
rat_plasma_sig_diff_ctrl_nonsurv_res <- rat_sig_diffs_along_time(subset(rat_sepsis_data, material == "plasma" & group %in% c("septic non-survivor", "control")), corr_fdr = T)
rat_plasma_time_fc_diff_ctrl_nonsurv <- rat_plasma_sig_diff_ctrl_nonsurv_res$time_fold_change

rat_liver_sig_diff_ctrl_surv_res <- rat_sig_diffs_along_time(subset(rat_sepsis_data, material == "liver" & group %in% c("septic survivor", "control")), corr_fdr = T)
rat_liver_time_fc_diff_ctrl_surv <- rat_liver_sig_diff_ctrl_surv_res$time_fold_change
rat_liver_sig_diff_ctrl_nonsurv_res <- rat_sig_diffs_along_time(subset(rat_sepsis_data, material == "liver" & group %in% c("septic non-survivor", "control")), corr_fdr = T)
rat_liver_time_fc_diff_ctrl_nonsurv <- rat_liver_sig_diff_ctrl_nonsurv_res$time_fold_change

rat_heart_sig_diff_ctrl_surv_res <- rat_sig_diffs_along_time(subset(rat_sepsis_data, material == "heart" & group %in% c("septic survivor", "control")), corr_fdr = T)
rat_heart_time_fc_diff_ctrl_surv <- rat_heart_sig_diff_ctrl_surv_res$time_fold_change
rat_heart_sig_diff_ctrl_nonsurv_res <- rat_sig_diffs_along_time(subset(rat_sepsis_data, material == "heart" & group %in% c("septic non-survivor", "control")), corr_fdr = T)
rat_heart_time_fc_diff_ctrl_nonsurv <- rat_heart_sig_diff_ctrl_nonsurv_res$time_fold_change

#Find differences with ANOVA
met_set <- colnames(rat_sepsis_data)
met_set <- met_set[colSums(is.na(rat_sepsis_data)) < 2]
met_set <- met_set[-1:-4]
names(met_set) <- met_set
fml <- concentration ~ group*time
##Split analysis into one for control vs sepsis and one for survival vs nonsurvival
anova.car.c <- list()
anova.car.s <- list()
for (mat in unique(rat_sepsis_data$material)){
  ftd_c <- subset(rat_sepsis_data, material == mat)
  ftd_c$group <- sub(pattern = "septic.+", x = as.character(ftd_c$group), replacement = "septic")
  ftd_c$time <- ftd_c$`time point`
  ftd_s <- subset(rat_sepsis_data, group != "control" & material == mat)
  ftd_s$time <- ftd_s$`time point`
  ftd_s <- subset(ftd_s, time != "72h")
  anova.car.c[[mat]] <- t3ANOVA(col.set = met_set, data = ftd_c, formula = fml, id.vars = c("group", "time"), random = ~1|group, use.corAR = FALSE)
  anova.car.s[[mat]] <- t3ANOVA(col.set = met_set, data = ftd_s, formula = fml, id.vars = c("group", "time"), random = ~1|group, use.corAR = FALSE)
}
rat.sig.anova.car.c.pre.class <- lapply(anova.car.c, function(e) colnames(e$ps)[colAnys(e$ps[c("group", "group:time"), ] <= 0.05)])
rat.sig.anova.car.s.pre.class <- lapply(anova.car.s, function(e) colnames(e$ps)[colAnys(e$ps[c("group", "group:time"), ] <= 0.05)])
anova.car.c.pre.ps <- lapply(anova.car.c, function(e) e$ps)
anova.car.s.pre.ps <- lapply(anova.car.s, function(e) e$ps)
anova.car.c.ps <- lapply(anova.car.c.pre.ps, apply, 2, p.adjust, method = "fdr") #Note: does *not* qualitatively differ to FDR correction for .c. and .s. seperately
anova.car.s.ps <- lapply(anova.car.s.pre.ps, apply, 2, p.adjust, method = "fdr") #Note: does *not* qualitatively differ to FDR correction for .c. and .s. seperately
#anova.car.c.ps <- lapply(anova.car.c.ps, t) #apply above returns a p-by-n matrix from an n-by-p input
#anova.car.s.ps <- lapply(anova.car.s.ps, t) #apply above returns a p-by-n matrix from an n-by-p input
rat.sig.anova.car.c.class <- lapply(anova.car.c.ps, function(e) colnames(e)[colAnys(e[c("group", "group:time"), ] <= 0.05)])
rat.sig.anova.car.s.class <- lapply(anova.car.s.ps, function(e) colnames(e)[colAnys(e[c("group", "group:time"), ] <= 0.05)])
for (mat in names(rat.sig.anova.car.c.class)){
  #rat.sig.anova.car.c.class[[mat]] <- intersect(rat.sig.anova.car.c.class[[mat]], names(which(anova.car.c.var.homog.p[[mat]] > 0.05)))
  #rat.sig.anova.car.c.class[[mat]] <- intersect(rat.sig.anova.car.c.class[[mat]], names(which(anova.car.c.normality.p[[mat]] > 0.05)))
  rat.sig.anova.car.c.class[[mat]] <- intersect(rat.sig.anova.car.c.class[[mat]], rat.sig.anova.car.c.pre.class[[mat]])
}
for (mat in names(rat.sig.anova.car.s.class)){
  #rat.sig.anova.car.s.class[[mat]] <- intersect(rat.sig.anova.car.s.class[[mat]], names(which(anova.car.s.var.homog.p[[mat]] > 0.05)))
  #rat.sig.anova.car.s.class[[mat]] <- intersect(rat.sig.anova.car.s.class[[mat]], names(which(anova.car.s.normality.p[[mat]] > 0.05)))
  rat.sig.anova.car.s.class[[mat]] <- intersect(rat.sig.anova.car.s.class[[mat]], rat.sig.anova.car.s.pre.class[[mat]])
}

#Get time point of differences with post hoc test
fml <- concentration ~ timeGroup - 1
anova.car.ph.models <- list()
for (mat in unique(rat_sepsis_data$material)){
  ftd_ph <- subset(rat_sepsis_data, material == mat)
  ftd_ph$time <- ftd_ph$`time point`
  ftd_ph$timeGroup <- interaction(ftd_ph$time, ftd_ph$group, drop = TRUE)
  anova.car.ph.models[[mat]] <- lapply(met_set, fit_lin_mod_lme, data = ftd_ph, formula = fml, id.vars = c("time", "group", "timeGroup"), random = ~1|group, use.corAR = FALSE)
}
anova.car.ph.models <- lapply(anova.car.ph.models, lapply, function(e) try(eval(e)) )
anova.car.ph.models <- lapply(anova.car.ph.models, function(l) l[!sapply(l, is.character)] )
for (mat in unique(rat_sepsis_data$material)){
  ret <- intersect(names(anova.car.ph.models[[mat]]), union(names(anova.car.c[[mat]]$models), names(anova.car.s[[mat]]$models)))
  anova.car.ph.models[[mat]] <- anova.car.ph.models[[mat]][ret]
}
timeGroup <- interaction(factor(rat_sepsis_data$`time point`), factor(rat_sepsis_data$group), drop = TRUE)
contr.m <- matrix(0, nrow = 8, ncol = 5) #8 levels ({6h,24h}X{S,NS,C} + {72}X{S,C}); 5 contrasts (C vs {NS,S} at {6h,24h,72h], and NS vs S at {6h,24h})
colnames(contr.m) <- c("Seps vs Comp, 24h", "Seps vs Comp, 6h", "Seps vs Comp, 72h", "S vs NS, 24h", "S vs NS, 6h")
rownames(contr.m) <- levels(timeGroup)
#Control vs Sepsis
contr.m[1:2, 1:2] <- diag(x = -2, nrow = 2, ncol = 2)
contr.m[3, 3] <- -1
contr.m[4:5, 1:2] <- diag(x = 1, nrow = 2, ncol = 2)
contr.m[6:7, 1:2] <- diag(x = 1, nrow = 2, ncol = 2)
contr.m[8, 3] <- 1
#Survivors vs Nonsurvivors
contr.m[4:5, 4:5] <- diag(x = -1, nrow = 2, ncol = 2)
contr.m[6:7, 4:5] <- diag(x = 1, nrow = 2, ncol = 2)
#Actual post hoc test
anova.car.ph.res <- lapply(anova.car.ph.models, lapply, function(e) summary(glht(e, linfct = mcp(timeGroup = t(contr.m)))))
anova.car.ph.ps <- lapply(anova.car.ph.res, sapply, function(e) e$test$pvalues)
anova.car.ph.sig.contr <- lapply(anova.car.ph.ps, function(e) lapply(data.frame(e), function(col) which(col <= 0.05)))
for (mat in unique(rat_sepsis_data$material)){
  rownames(anova.car.ph.ps[[mat]]) <- rownames(anova.car.ph.res[[mat]][[1]]$linfct)
  names(anova.car.ph.sig.contr[[mat]]) <- names(anova.car.ph.res[[mat]])
}

#Find differences in pheno vars with ANOVA
met_set <- colnames(rat_sepsis_data)
met_set <- met_set[(which(met_set == "H1") + 1):length(met_set)]
names(met_set) <- met_set
fml <- concentration ~ group*time
##Split analysis into one for control vs sepsis and one for survival vs nonsurvival
mat <- "plasma"
ftd_c <- subset(rat_sepsis_data, material == mat)
ftd_c$group <- sub(pattern = "septic.+", x = as.character(ftd_c$group), replacement = "septic")
ftd_c$time <- ftd_c$`time point`
ftd_s <- subset(rat_sepsis_data, group != "control" & material == mat)
ftd_s$time <- ftd_s$`time point`
ftd_s <- subset(ftd_s, time != "72h")
anova.car.c.pheno <- t3ANOVA(col.set = met_set, data = ftd_c, formula = fml, id.vars = c("group", "time"), random = ~1|group, use.corAR = FALSE)
anova.car.s.pheno <- t3ANOVA(col.set = met_set, data = ftd_s, formula = fml, id.vars = c("group", "time"), random = ~1|group, use.corAR = FALSE)
anova.car.s.pheno.pre.ps <- anova.car.s.pheno$ps
anova.car.c.pheno.pre.ps <- anova.car.c.pheno$ps
rat.sig.anova.c.pheno.car.pre.class <- colnames(anova.car.c.pheno.pre.ps)[colAnys(anova.car.c.pheno.pre.ps[c("group", "group:time"), ] <= 0.05)]
rat.sig.anova.s.pheno.car.pre.class <- colnames(anova.car.s.pheno.pre.ps)[colAnys(anova.car.s.pheno.pre.ps[c("group", "group:time"), ] <= 0.05)]
anova.car.c.pheno.ps <- apply(anova.car.c.pheno.pre.ps, 2, p.adjust, method = "fdr") #Note: does *not* qualitatively differ to FDR correction for .c. and .s. seperately
#anova.car.c.pheno.ps <- t(anova.car.c.pheno.ps) #apply above returns a p-by-n matrix from an n-by-p input
anova.car.s.pheno.ps <- apply(anova.car.s.pheno.pre.ps, 2, p.adjust, method = "fdr") #Note: does *not* qualitatively differ to FDR correction for .s. and .s. seperately
#anova.car.s.pheno.ps <- t(anova.car.s.pheno.ps) #apply above returns a p-by-n matrix from an n-by-p input
rat.sig.anova.car.c.pheno.class <- colnames(anova.car.c.pheno.ps)[colAnys(anova.car.c.pheno.ps[c("group", "group:time"), ] <= 0.05)]
rat.sig.anova.car.s.pheno.class <- colnames(anova.car.s.pheno.ps)[colAnys(anova.car.s.pheno.ps[c("group", "group:time"), ] <= 0.05)]
#rat.sig.anova.car.c.pheno.class <- intersect(rat.sig.anova.car.c.pheno.class, names(which(anova.car.c.pheno.var.homog.p > 0.05)))
#rat.sig.anova.car.c.pheno.class <- intersect(rat.sig.anova.car.c.pheno.class, names(which(anova.car.c.pheno.normality.p > 0.05)))
rat.sig.anova.car.c.pheno.class <- intersect(rat.sig.anova.car.c.pheno.class, rat.sig.anova.c.pheno.car.pre.class)
#rat.sig.anova.car.s.pheno.class <- intersect(rat.sig.anova.car.s.pheno.class, names(which(anova.car.s.pheno.var.homog.p > 0.05)))
#rat.sig.anova.car.s.pheno.class <- intersect(rat.sig.anova.car.s.pheno.class, names(which(anova.car.s.pheno.normality.p > 0.05)))
rat.sig.anova.car.s.pheno.class <- intersect(rat.sig.anova.car.s.pheno.class, rat.sig.anova.s.pheno.car.pre.class)

#Get time point of differences with post hoc test
fml <- concentration ~ timeGroup - 1
anova.car.pheno.ph.levene <- list()
anova.car.pheno.ph.models <- list()
ftd_ph <- subset(rat_sepsis_data, material == "plasma")
ftd_ph$time <- ftd_ph$`time point`
ftd_ph$timeGroup <- interaction(ftd_ph$time, ftd_ph$group, drop = TRUE)
anova.car.pheno.ph.models <- lapply(met_set, fit_lin_mod_lme, data = ftd_ph, formula = fml, id.vars = c("time", "group", "timeGroup"), random = ~1|group, use.corAR = FALSE)
anova.car.pheno.ph.models <- lapply(anova.car.pheno.ph.models, eval)
anova.car.pheno.ph.models <- anova.car.pheno.ph.models[sapply(anova.car.pheno.ph.models, is.list)]
ret <- intersect(names(anova.car.pheno.ph.models), union(names(anova.car.c.pheno$models), names(anova.car.s.pheno$models)))
anova.car.pheno.ph.models <- anova.car.pheno.ph.models[ret]
timeGroup <- interaction(factor(rat_sepsis_data$`time point`), factor(rat_sepsis_data$group), drop = TRUE)
contr.m <- matrix(0, nrow = 8, ncol = 5) #8 levels ({6h,24h}X{S,NS,C} + {72}X{S,C}); 5 contrasts (C vs {NS,S} at {6h,24h,72h], and NS vs S at {6h,24h})
colnames(contr.m) <- c("Seps vs Comp, 24h", "Seps vs Comp, 6h", "Seps vs Comp, 72h", "S vs NS, 24h", "S vs NS, 6h")
rownames(contr.m) <- levels(timeGroup)
#Control vs Sepsis
contr.m[1:2, 1:2] <- diag(x = -2, nrow = 2, ncol = 2)
contr.m[3, 3] <- -1
contr.m[4:5, 1:2] <- diag(x = 1, nrow = 2, ncol = 2)
contr.m[6:7, 1:2] <- diag(x = 1, nrow = 2, ncol = 2)
contr.m[8, 3] <- 1
#Survivors vs Nonsurvivors
contr.m[4:5, 4:5] <- diag(x = -1, nrow = 2, ncol = 2)
contr.m[6:7, 4:5] <- diag(x = 1, nrow = 2, ncol = 2)
#Actual post hoc test
anova.car.pheno.ph.res <- lapply(anova.car.pheno.ph.models, function(e) summary(glht(e, linfct = mcp(timeGroup = t(contr.m)))))
anova.car.pheno.ph.ps <- sapply(anova.car.pheno.ph.res, function(e) e$test$pvalues)
anova.car.pheno.ph.sig.contr <- lapply(data.frame(anova.car.pheno.ph.ps), function(col) which(col <= 0.05))
rownames(anova.car.pheno.ph.ps) <- rownames(anova.car.pheno.ph.res[[1]]$linfct)
names(anova.car.pheno.ph.sig.contr) <- names(anova.car.pheno.ph.res)

#Save everything ANOVA to file
anova_complete_res <- grep(pattern = "anova.car", x = names(environment()), value = TRUE)
anova_complete_res <- grep(pattern = "models|res|terms|normality", x = anova_complete_res, invert = TRUE, value = TRUE)
save(list = anova_complete_res, file = paste0(out_dir, "ANOVA_complete_res.RData"))

#Check correlations between tissue types on ungrouped metabolites
rat_sepsis_data_normal <- rat_sepsis_data
rat_sepsis_data_normal[, -1:-4] <- scale(rat_sepsis_data[-1:-4])
cross_mat_corr <- list()
for (comp in list(c("heart", "plasma"), c("heart", "liver"), c("plasma", "liver"))){
  comp_str <- paste(comp, collapse = " ~ ")
  cross_mat_corr[[comp_str]] <- list()
  met_set <- colnames(rat_sepsis_data_normal)
  met_set <- met_set[which(!colAnys(is.na(subset(rat_sepsis_data_normal, material %in% comp))))]
  met_set <- met_set[-1:-4]
  for (met in met_set){
    test_data <- subset(rat_sepsis_data_normal, material %in% comp, c("Sample Identification", "material", met))
    colnames(test_data)[1] <- "ID"
    test_data$ID <- substr(test_data$ID, start = 1, stop = 3)
    test_data <- dcast(test_data, ID ~ material, value.var = met)
    cross_mat_corr[[comp_str]][[met]] <- cor.test(x = test_data[,2], y = test_data[,3])
  }
}
cross_mat_coeff <- lapply(cross_mat_corr, sapply, function(r){ r$estimate })
cross_mat_p <- lapply(cross_mat_corr, sapply, function(r){ r$p.value })
cross_mat_fdr <- lapply(cross_mat_p, p.adjust, method = "fdr")

#Check correlations between tissue types on grouped metabolites
coarse_group_list <- rat_sepsis_legend[match(colnames(rat_sepsis_data)[-1:-4], rat_sepsis_legend[,1]), 3] #lucky ... no col in data without match in legend
coarse_group_list[is.na(coarse_group_list)] <- "Plasma Creatinine"
rat_sepsis_data_grouped <- cbind(rat_sepsis_data[,1:4], matrix(0, nrow = nrow(rat_sepsis_data), ncol=length(unique(coarse_group_list))))
colnames(rat_sepsis_data_grouped)[-1:-4] <- unique(coarse_group_list)
rsd <- rat_sepsis_data
rsd[, pheno_sel] <- scale(rsd[, pheno_sel])
for (n in 1:nrow(rat_sepsis_data)){
  m_agg <- tapply(X = t(rsd[n, metab_sel]), INDEX = factor(coarse_group_list[metab_sel - 4]), FUN = sum)
  p_agg <- tapply(X = t(rsd[n, pheno_sel]), INDEX = factor(coarse_group_list[pheno_sel - 4]), FUN = mean)
  b_agg <- c(m_agg, p_agg)
  rat_sepsis_data_grouped[n, -1:-4] <- b_agg[match(colnames(rat_sepsis_data_grouped)[-1:-4], names(b_agg))]
}
cross_mat_group_corr <- list()
for (comp in list(c("heart", "plasma"), c("heart", "liver"), c("plasma", "liver"))){
  comp_str <- paste(comp, collapse = " ~ ")
  cross_mat_group_corr[[comp_str]] <- list()
  for (met in unique(coarse_group_list[1:199])){
    test_data <- subset(rat_sepsis_data_grouped, material %in% comp, c("Sample Identification", "material", met))
    colnames(test_data)[1] <- "ID"
    test_data$ID <- substr(test_data$ID, start = 1, stop = 3)
    test_data <- dcast(test_data, ID ~ material, value.var = met)
    if (all(colSums(is.na(test_data)) < nrow(test_data)))
      cross_mat_group_corr[[comp_str]][[met]] <- cor.test(x = test_data[,2], y = test_data[,3])
  }
}
cross_mat_group_coeff <- lapply(cross_mat_group_corr, sapply, function(r){ r$estimate })
cross_mat_group_p <- lapply(cross_mat_group_corr, sapply, function(r){ r$p.value })
cross_mat_group_fdr <- lapply(cross_mat_group_p, p.adjust, method = "fdr")

cross_mat_df <- lapply(names(cross_mat_coeff), function(n) data.frame(Metabolite = names(cross_mat_p[[n]]), Corr = cross_mat_coeff[[n]], p = cross_mat_p[[n]], FDR = cross_mat_fdr[[n]]))
names(cross_mat_df) <- names(cross_mat_coeff)
cross_mat_group_df <- lapply(names(cross_mat_group_coeff), function(n) data.frame(Metabolite = names(cross_mat_group_p[[n]]), Corr = cross_mat_group_coeff[[n]], p = cross_mat_group_p[[n]], FDR = cross_mat_group_fdr[[n]]))
names(cross_mat_group_df) <- names(cross_mat_group_p)
for (n in names(cross_mat_df)){
  path <- paste0(out_dir, "cross_material_correlation_", sub(pattern = " ~ ", replacement = "_", x = n), ".csv")
  fwrite(x = cross_mat_df[[n]], file = path, sep = "\t")
}

#Build long format tables for plotting
##Scale measurement values at maximum concentration
rat_sepsis_data_max_norm <- max_norm(x = rat_sepsis_data, subset = -1:-4)
##Melt scaled data into long form
rat_sepsis_data_long_form <- melt(data = rat_sepsis_data_max_norm, id.vars = c("Sample Identification", "material", "group", "time point"))
##Reduce to significantly different metabolites
rat_sepsis_data_long_form_sig <- subset(rat_sepsis_data_long_form, subset = as.character(variable) %in% rat_sig_t_class)
rat_sepsis_data_long_form_sig$variable <- factor(as.character(rat_sepsis_data_long_form_sig$variable), levels = rat_sig_t_class, ordered = TRUE)

#Build tables for cluster-heatmaps
##Scale mesaurement values by standardization
rat_sepsis_data_normal <- rat_sepsis_data
rat_sepsis_data_normal[,-1:-4] <- scale(x = rat_sepsis_data[,-1:-4])
rat_sepsis_data_normal$`time point` <- as.factor(rat_sepsis_data_normal$`time point`)
##Group metabolites
coarse_group_list <- rat_sepsis_legend[rat_sepsis_legend[,1] %in% colnames(rat_sepsis_data[,-1:-4]), 3]
rat_sepsis_data_normal_grouped <- cbind(rat_sepsis_data[,1:4], matrix(0, nrow = nrow(rat_sepsis_data_normal), ncol=length(unique(coarse_group_list))))
colnames(rat_sepsis_data_normal_grouped)[-1:-4] <- unique(coarse_group_list)
rat_sepsis_data_normal_grouped$`time point` <- as.factor(rat_sepsis_data_normal_grouped$`time point`)
for (n in 1:nrow(rat_sepsis_data_normal)){
  rat_sepsis_data_normal_grouped[n, -1:-4] <- aggregate(t(rat_sepsis_data_normal[n,-1:-4]), by = list(coarse_group_list), FUN = mean, na.action = na.omit)[,2]
}
##Split for metabolites and "phenotypical" factors
if ("HR" %in% colnames(rat_sepsis_data)){
  split_start <- which(colnames(rat_sepsis_data) == "HR")
}else{
  split_start <- which(colnames(rat_sepsis_data) == "pH")
}

pheno_sel <- split_start:ncol(rat_sepsis_data)
metab_sel <- 5:(split_start-1)
group_pheno_sel <- which(colnames(rat_sepsis_data_normal_grouped) %in% colnames(rat_sepsis_data[,pheno_sel]))
group_metab_sel <- which(colnames(rat_sepsis_data_normal_grouped) %in% unique(coarse_group_list[metab_sel - 4]))
##Build covariance matrix with metabolite groups
rat_sepsis_data_normal_grouped_metab_cov <- cbind(rat_sepsis_data_normal_grouped[, 1:4], cov(t(rat_sepsis_data_normal_grouped[,group_metab_sel]), use = "pairwise.complete.obs"))
rat_sepsis_data_normal_grouped_metab_cov <- rat_sepsis_data_normal_grouped_metab_cov[, !apply(rat_sepsis_data_normal_grouped_metab_cov, 2, function(x){all(is.na(x))})]
##Build covariance matrix with original metabolites
rat_sepsis_data_normal_metab_cov <- cbind(rat_sepsis_data_normal[, 1:4], cov(t(rat_sepsis_data_normal[,metab_sel]), use = "pairwise.complete.obs"))
rat_sepsis_data_normal_metab_cov <- rat_sepsis_data_normal_metab_cov[, !apply(rat_sepsis_data_normal_metab_cov, 2, function(x){all(is.na(x))})]
rat_sepsis_data_normal_pheno_cov <- cbind(rat_sepsis_data_normal[, 1:4], cov(t(rat_sepsis_data_normal[,pheno_sel]), use = "pairwise.complete.obs"))
rat_sepsis_data_normal_pheno_cov <- rat_sepsis_data_normal_pheno_cov[, !apply(rat_sepsis_data_normal_pheno_cov, 2, function(x){all(is.na(x))})]
##Build correlation matrix with metabolite groups
rat_sepsis_data_normal_grouped_metab_cor <- cbind(rat_sepsis_data_normal_grouped[, 1:4], cor(t(rat_sepsis_data_normal_grouped[,group_metab_sel]), use = "pairwise.complete.obs"))
rat_sepsis_data_normal_grouped_metab_cor <- rat_sepsis_data_normal_grouped_metab_cor[, !apply(rat_sepsis_data_normal_grouped_metab_cor, 2, function(x){all(is.na(x))})]
##Build covariance matrix with original metabolites
rat_sepsis_data_normal_metab_cor <- cbind(rat_sepsis_data_normal[, 1:4], cor(t(rat_sepsis_data_normal[, metab_sel]), use = "pairwise.complete.obs"))
rat_sepsis_data_normal_metab_cor <- rat_sepsis_data_normal_metab_cor[, !apply(rat_sepsis_data_normal_metab_cor, 2, function(x){all(is.na(x))})]
rat_sepsis_data_normal_pheno_cor <- cbind(rat_sepsis_data_normal[, 1:4], cor(t(rat_sepsis_data_normal[, pheno_sel]), use = "pairwise.complete.obs"))
rat_sepsis_data_normal_pheno_cor <- rat_sepsis_data_normal_pheno_cor[, !apply(rat_sepsis_data_normal_pheno_cor, 2, function(x){all(is.na(x))})]

##Build covariance matrix of metabolites with original metabolites
###Survival-ignorant
rat_sepsis_data_normal_conc_metab_cov <- cov(rat_sepsis_data_normal[, metab_sel], use = "pairwise.complete.obs")
rat_sepsis_data_normal_conc_pheno_cov <- cov(rat_sepsis_data_normal[, pheno_sel], use = "pairwise.complete.obs")
###Survival-regardent
rat_sepsis_data_normal_NS_conc_metab_cov <- cov(subset(rat_sepsis_data_normal, group == "septic non-survivor", metab_sel), use = "pairwise.complete.obs")
rat_sepsis_data_normal_S_conc_metab_cov <- cov(subset(rat_sepsis_data_normal, group == "septic survivor", metab_sel), use = "pairwise.complete.obs")
rat_sepsis_data_normal_NS_conc_pheno_cov <- cov(subset(rat_sepsis_data_normal, group == "septic non-survivor", pheno_sel), use = "pairwise.complete.obs")
rat_sepsis_data_normal_S_conc_pheno_cov <- cov(subset(rat_sepsis_data_normal, group == "septic survivor", pheno_sel), use = "pairwise.complete.obs")
##Build covariance matrix of metabolites with metabolites groups
###Survival-ignorant
rat_sepsis_data_normal_grouped_conc_metab_cov <- cov(rat_sepsis_data_normal_grouped[, group_metab_sel], use = "pairwise.complete.obs")
###Survival-regardent
rat_sepsis_data_normal_NS_grouped_conc_metab_cov <- cov(subset(rat_sepsis_data_normal_grouped, group == "septic non-survivor", group_metab_sel), use = "pairwise.complete.obs")
rat_sepsis_data_normal_S_grouped_conc_metab_cov <- cov(subset(rat_sepsis_data_normal_grouped, group == "septic survivor", group_metab_sel), use = "pairwise.complete.obs")

####################
#Plot data
####################

##Rat, cluster-heatmap, metabolites, and phenomenological vars, fold change with respect to control
x <- rat_sepsis_data_fc
x$group <- reorder(x$group, (x$group == "septic survivor") + (2 * (x$group == "septic non-survivor")))
x$`time point` <- reorder(x$`time point`, (x$`time point` == "6h") + (2 * (x$`time point` == "24h")) + (3 * (x$`time point` == "72h")))
x <- x[order(x$group),]
x <- x[order(x$`time point`),]
xm <- x[,c(1:4, metab_sel)]
xm[, which(colAnys(is.na(xm)))] <- NULL
xm[, 4 + which(colAnys(is.infinite(as.matrix(xm[,5:ncol(xm)]))))] <- NULL
xmt <- data.frame(t(xm[, -1:-4]))
rownames(xmt) <- colnames(xm[,-1:-4])
colnames(xmt) <- xm$`Sample Identification`

mat_sigs <- list()
mat_sig_subset <- colnames(rat_plasma_sig_diff_ctrl_surv_res$time_sig_t_diff) %in% rownames(xmt)
mat_sigs[["plasma"]] <- as.matrix(colAnys(rat_plasma_sig_diff_ctrl_surv_res$time_sig_t_diff[,-1] <= 0.05) + 2*colAnys(rat_plasma_sig_diff_ctrl_nonsurv_res$time_sig_t_diff[,-1] <= 0.05))[mat_sig_subset]
mat_sigs[["liver"]] <- as.matrix(colAnys(rat_liver_sig_diff_ctrl_surv_res$time_sig_t_diff[,-1] <= 0.05) + 2*colAnys(rat_liver_sig_diff_ctrl_nonsurv_res$time_sig_t_diff[,-1] <= 0.05))[mat_sig_subset]
mat_sigs[["heart"]] <- as.matrix(colAnys(rat_heart_sig_diff_ctrl_surv_res$time_sig_t_diff[,-1] <= 0.05) + 2*colAnys(rat_heart_sig_diff_ctrl_nonsurv_res$time_sig_t_diff[,-1] <= 0.05))[mat_sig_subset]
mat_sigs <- lapply(mat_sigs, function(x){ c(NA, "Sig diff Surv", "Sig diff Nonsurv")[x + 1] })

for (mat in unique(xm$material)){
  h <- heatmaply(x = xmt[, xm$material == mat], dendrogram = "none", plot_method = "plotly", col_side_colors = subset(xm, material == mat, c("time point", "group")), row_side_colors = data.frame(sig_reg = mat_sigs[[mat]]), key.title = "fold change with\nrespect to control", margins = c(50,100,0,50), subplot_heights = c(0.03, 0.97))
  #png(filename = paste0(out_dir, "rat_fold_change_", mat, "_pheno.png"), width = 800, height = 1600)
  h$width <- 800
  h$height <- 1600
  export(p = h, file = paste0(out_dir, "rat_fold_change_", mat, "_metab.png"))
}

xp <- x[,c(1:4, pheno_sel)]
xp <- xp[, colnames(xp) != "Acetoacetate"] # Acetoacetate has NA at 6h measurements
xp <- na.omit(xp)
xp[, 4 + which(colAnys(is.infinite(as.matrix(xm[,5:ncol(xp)]))))] <- NULL
xpt <- data.frame(t(xp[, -1:-4]))
rownames(xpt) <- colnames(xp[,-1:-4])
colnames(xpt) <- xp$`Sample Identification`

mat <- "plasma"
mat_sig_subset <- colnames(rat_plasma_sig_diff_ctrl_surv_res$time_sig_t_diff) %in% rownames(xpt)
mat_sig_pl <- as.matrix(colAnys(rat_plasma_sig_diff_ctrl_surv_res$time_sig_t_diff[,-1] <= 0.05) + 2*colAnys(rat_plasma_sig_diff_ctrl_nonsurv_res$time_sig_t_diff[,-1] <= 0.05))[mat_sig_subset]
mat_sig_pheno <- c(NA, "Sig diff Surv", "Sig diff Nonsurv")[mat_sig_pl + 1]
h <- heatmaply(file = paste0(out_dir, "rat_fold_change_", mat, "_pheno.png"), x = xpt[, xp$material == mat], dendrogram = "row", plot_method = "plotly", col_side_colors = subset(xp, material == mat, c("time point", "group")), row_side_colors = data.frame(sig_reg = mat_sig_pheno), key.title = "fold change with\nrespect to control", margins = c(100,80,0,250), subplot_heights = c(0.05, 0.95))

rm(x, xm, xp, xpt, xmt, mat_sigs)

#The same plots broken down by metabolite group, and accompanied by ANOVA significance
x <- rat_sepsis_data
x$group <- reorder(x$group, (x$group == "control") + (2 * (x$group == "septic survivor")) + (3 * (x$group == "septic non-survivor")))
x$`time point` <- reorder(x$`time point`, (x$`time point` == "6h") + (2 * (x$`time point` == "24h")) + (3 * (x$`time point` == "72h")))
x <- x[order(x$group),]
x <- x[order(x$`time point`),]
xm <- x[,c(1:4, metab_sel)]
xmt <- data.frame(t(xm[, -1:-4]))
rownames(xmt) <- colnames(xm[,-1:-4])
colnames(xmt) <- xm$`Sample Identification`

mat_sigs <- list()
for (mat in unique(rat_sepsis_data$material)){
  control.sig <- rownames(xmt[,-1:-4]) %in% rat.sig.anova.car.c.class[[mat]]
  survival.sig <- rownames(xmt[,-1:-4]) %in% rat.sig.anova.car.s.class[[mat]]
  mat_sigs[[mat]] <- data.frame(control.sig = control.sig, survival.sig = survival.sig, stringsAsFactors = FALSE)
}
mat_sigs <- lapply(mat_sigs, lapply, function(x){ c("nonsignif.", "p < 0.05")[x + 1] })
mat_sigs <- lapply(mat_sigs, as.data.frame)
for (mat in names(mat_sigs)){
  colnames(mat_sigs[[mat]]) <- c("Control vs Sepsis", "S vs NS")
}

for (mat in unique(xm$material)){
  for (met_group in unique(coarse_group_list[metab_sel - 4])){
    group_sel <- coarse_group_list[metab_sel - 4] %in% met_group
    #xfplotdat <- t(max_norm(t(xmt[group_sel, xm$material == mat])))
    xfplotdat <- t(max_norm(t(xmt[group_sel, xm$material == mat])))
    sel <- !rowAlls(is.na(xfplotdat))
    lower_margin <- 85
    top_row_h <- 0.03 * 76/sum(sel)
    subplot_h <- c(top_row_h, 1 - top_row_h)
    if (sum(sel) > 1){
      h <- heatmaply(x = xfplotdat[sel, ], 
                     dendrogram = "row", 
                     plot_method = "plotly", 
                     col_side_colors = subset(xm, material == mat, c("time point", "group")), 
                     row_side_colors = mat_sigs[[mat]][group_sel, ][sel, ], 
                     key.title = "concentration", 
                     margins = c(lower_margin,100,NA,50), 
                     height = lower_margin + round(sum(sel) * 1100/76),
                     subplot_heights = subplot_h)
      h$width <- 1000
      #h$height <- 1200
      h$height <- lower_margin + round(sum(sel) * 1100/76) #1200 is a good height for 76 rows of metabolites
      export(p = h, file = paste0(out_dir, "rat_heatmap_", mat, "_", met_group, ".png"))
    }
    else
      print(paste0(met_group, " in ", mat, " has too few metabolites, possibly after filtering all-NA rows"))
  }
}

#Again for pheno vars
xm <- x[x$material == "plasma", c(1:4, pheno_sel)]
xmt <- data.frame(t(xm[, -1:-4]))
rownames(xmt) <- colnames(xm[,-1:-4])
colnames(xmt) <- xm$`Sample Identification`
control.sig <- rownames(xmt[, -1:-4]) %in% rat.sig.anova.car.c.pheno.class
survival.sig <- rownames(xmt[, -1:-4]) %in% rat.sig.anova.car.s.pheno.class
mat_sigs <- data.frame(control.sig = control.sig, survival.sig = survival.sig, stringsAsFactors = FALSE)
mat_sigs <- lapply(mat_sigs, function(x){ c("nonsignif.", "p < 0.05")[x + 1] })
mat_sigs <- data.frame(mat_sigs)
colnames(mat_sigs) <- c("Control vs Sepsis", "S vs NS")

lower_margin <- 85
xfplotdat <- xmt
sel <- !rowAlls(is.na(xfplotdat))
top_row_h <- 0.03 * 76/sum(sel)
subplot_h <- c(top_row_h, 1 - top_row_h)
h <- heatmaply(x = t(max_norm(t(xfplotdat[sel, ]))),
               dendrogram = "row", 
               plot_method = "plotly", 
               col_side_colors = subset(xm, TRUE , c("time point", "group")), 
               row_side_colors = mat_sigs[sel, ], 
               key.title = "concentration", 
               margins = c(lower_margin,100,NA,50), 
               height = lower_margin + round(sum(sel) * 1000/76),
               subplot_heights = subplot_h, 
               subplot_widths = c(0.86, 0.06, 0.08))
h$width <- 1000
h$height <- lower_margin + round(sum(sel) * 1100/76) #1200 is a good height for 76 rows of metabolites
export(p = h, file = paste0(out_dir, "rat_heatmap_pheno.png"))

##Rat, cluster-heatmap, metabolites
x <- rat_sepsis_data_normal[,c(1:4, metab_sel)]
x <- x[, apply(x, 2, function(x){ sum(is.na(x)) == 0 })]
heatmaply(x = x[,-1:-4], row_side_colors = x[c("group")], file = paste0(out_dir, "rat_normal_metab.png"), main = "Metabolite profile does not cluster survival well", key.title = "Concentration\n(normalized)")
x <- na.omit(rat_sepsis_data_normal[,c(1:4, pheno_sel)])
rownames(x) <- paste0(x$`Sample Identification`, ", ", x$`time point`)
heatmaply(x = x[,-1:-4], row_side_colors = x[c("group")], file = paste0(out_dir, "rat_normal_pheno.png"), main = "Phenomenological profile has survival clusters", key.title = "Concentration\n(normalized)", column_text_angle = 90, showticklabels = c(T,T))
x <- na.omit(rat_sepsis_data_normal_grouped[,c(1:4, group_metab_sel)])
rownames(x) <- paste0(x$`Sample Identification`, ", ", x$`time point`)
heatmaply(x = x[,-1:-4], row_side_colors = x[c("group")], file = paste0(out_dir, "rat_normal_grouped_metab.png"), main = "Metablite group profiles somewhat cluster survival/control", key.title = "Concentration\n(normalized)", showticklabels = c(T,T))
rm("x")

##Rat, covariance cluster-heatmap, metabolites
x <- na.omit(rat_sepsis_data_normal_metab_cov)
rownames(x) <- paste0(x$`Sample Identification`, ", ", x$`time point`)
colnames(x)[-1:-4] <- paste0(x$`Sample Identification`, ", ", x$`time point`)
heatmaply(x = x[,-1:-4], row_side_colors = x["group"], col_side_colors = x["material"], file = paste0(out_dir, "rat_normal_metab_cov.png"), main = "", key.title = "Cov", showticklabels = TRUE)
x <- na.omit(rat_sepsis_data_normal_pheno_cov)
rownames(x) <- paste0(x$`Sample Identification`, ", ", x$`time point`)
colnames(x)[-1:-4] <- paste0(x$`Sample Identification`, ", ", x$`time point`)
heatmaply(x = x[,-1:-4], row_side_colors = x[c("group")], file = paste0(out_dir, "rat_normal_pheno_cov.png"), main = "Metabolite profile covariance clusters survival", key.title = "Cov", showticklabels = TRUE)
rm("x")

##rat, cluster-heatmap of non-CAP and non-FP patients, coarse grouped metabolites; error because nothing es left
#heatmaply(x = subset(rat_sepsis_data_normal_grouped, subset = rat_sepsis_data_normal_grouped$`CAP / FP` %in% c("-"), select = c(-1,-2,-3,-5)))

##rat, cluster-heatmap of patient covariance matrix, coarse grouped metabolites, survival marked
x <- na.omit(rat_sepsis_data_normal_grouped_metab_cov)
rownames(x) <- paste0(x$`Sample Identification`, ", ", x$`time point`)
colnames(x)[-1:-4] <- paste0(x$`Sample Identification`, ", ", x$`time point`)
heatmaply(x = x[,-1:-4], row_side_colors = x[c("group")], file = paste0(out_dir, "rat_normal_grouped_metab_cov.png"))
rm("x")

##rat, cluster-heatmap of patient correlation matrix, ungrouped metabolites, survival marked
x <- rat_sepsis_data_normal_metab_cor
rownames(x) <- paste0(x$`Sample Identification`, ", ", x$`time point`)
colnames(x)[-1:-4] <- paste0(x$`Sample Identification`, ", ", x$`time point`)
heatmaply(x = x[,-1:-4], row_side_colors = x[c("group")], file = paste0(out_dir, "rat_normal_metab_cor.png"), main = "Metabolite profile correlation gives\nonly material clusters")
x <- na.omit(rat_sepsis_data_normal_pheno_cor)
rownames(x) <- paste0(x$`Sample Identification`, ", ", x$`time point`)
colnames(x)[-1:-4] <- paste0(x$`Sample Identification`, ", ", x$`time point`)
heatmaply(x = x[,-1:-4], row_side_colors = x[c("group")], file = paste0(out_dir, "rat_normal_pheno_cor.png"), main = "Phenomenological profile correlation gives survival/control clusters")
rm("x")

##Rat, cluster-heatmat of patient correlation matrix, coarse grouped metabolites, survival and CAP/FP marked
x <- na.omit(rat_sepsis_data_normal_grouped_metab_cor)
rownames(x) <- paste0(x$`Sample Identification`, ", ", x$`time point`)
colnames(x)[-1:-4] <- paste0(x$`Sample Identification`, ", ", x$`time point`)
heatmaply(x = x[,-1:-4], row_side_colors = x[c("group")], file = paste0(out_dir, "rat_normal_grouped_metab_cor.png"), main = "Profiles of grouped metabolites cluster nothing")
rm("x")

##rat, cluster-heatmap of metabolite covariance matrix, survival-ignorant
x <- rat_sepsis_data_normal_conc_metab_cov
x <- x[-which(colnames(x) == "TUDCA"), -which(colnames(x) == "TUDCA")]
heatmaply(x = x, file = paste0(out_dir, "rat_normal_conc_metab_cov.png"), key.title = "Cov")
x <- rat_sepsis_data_normal_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "rat_normal_conc_pheno_cov.png"), key.title = "Cov")
rm("x")

##rat, cluster-heatmap of metabolite covariance matrix, survival-regardent
x <- rat_sepsis_data_normal_NS_conc_metab_cov
x <- x[-which(colnames(x) == "TUDCA"), -which(colnames(x) == "TUDCA")]
heatmaply(x = x, file = paste0(out_dir, "rat_normal_NS_conc_metab_cov.png"), key.title = "Cov")
x <- rat_sepsis_data_normal_S_conc_metab_cov
x <- x[-which(colnames(x) == "TUDCA"), -which(colnames(x) == "TUDCA")]
heatmaply(x = x, file = paste0(out_dir, "rat_normal_S_conc_metab_cov.png"), key.title = "Cov")
x <- rat_sepsis_data_normal_NS_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "rat_normal_NS_conc_pheno_cov.png"), key.title = "Cov")
x <- rat_sepsis_data_normal_S_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "rat_normal_S_conc_pheno_cov.png"), key.title = "Cov")
rm("x")

##rat, cluster-heatmap of covariance matrix of grouped metabolites, survival-ignorant
x <- rat_sepsis_data_normal_grouped_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "rat_normal_grouped_conc_metab_cov.png"))
rm("x")

##rat, cluster-heatmap of covariance matrix of grouped metabolites, survival-regardent
x <- rat_sepsis_data_normal_NS_grouped_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "rat_normal_NS_grouped_conc_metab_cov.png"), key.title = "Cov")
x <- rat_sepsis_data_normal_S_grouped_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "rat_normal_S_grouped_conc_metab_cov.png"), key.title = "Cov")
rm("x")

##Rat, p-val (t-test) plot of differences between non-survivors and survivors, grouped by time point
rat_plasma_time_sig_t_diff$material <- "plasma"
rat_heart_time_sig_t_diff$material <- "heart"
rat_liver_time_sig_t_diff$material <- "liver"
rat_time_sig_t_diff <- rbind(rat_plasma_time_sig_t_diff, rat_heart_time_sig_t_diff, rat_liver_time_sig_t_diff)
r_time_sig_t_diff_plot <- ggplot(melt(rat_time_sig_t_diff, id.vars = c("Time", "material")), aes(x = reorder(Time, -as.numeric(Time)), y = value)) +
#r_time_sig_t_diff_plot <- ggplot(melt(rat_time_sig_t_diff, id.vars = c("Time", "material")), aes(x = Time, y = value)) +
  facet_grid(material ~ .) +
  geom_point(position = position_jitter(width = 0.1), size = 0.7) +
  geom_hline(yintercept = 0.05) +
  scale_y_log10() +
  xlab("time") +
  ylab("p value, FDR-corrected") +
  ggtitle("Rat metabolites differ at two time points\nin non-/survival, most in plasma") +
  theme_bw()
ggsave(plot = r_time_sig_t_diff_plot, filename = "rat_all_days_survival_sig_diff.png", path = out_dir, width = 4, height = 4, units = "in")

##Rat, annas significantly different metabolites
x <- subset(na.omit(rat_sepsis_data_long_form), variable %in% annas_mets & material %in% c("liver") & group %in% c("septic survivor", "septic non-survivor"))
colnames(x)[colnames(x) == "time point"] <- "time_point"
x <- subset(x, time_point %in% c("6h", "24h"))
ggplot(x, aes(x = group, y = value, group = time_point, color = time_point)) +
  facet_wrap(facets = ~ variable + time_point + material, ncol = 10, scales = "free_y") +
  geom_jitter(height = 0, width = 0.3) +
  theme_bw()

##rat, metabolite concentration time courses, ANOVA results, control vs. sepsis
for (mat in names(rat.sig.anova.car.c.class)){
  r_time_course_sig_diff_dat <- subset(rat_sepsis_data, material == mat)
  r_time_course_sig_diff_dat <- max_norm(r_time_course_sig_diff_dat, -1:-4)
  r_time_course_sig_diff_dat <- melt(r_time_course_sig_diff_dat, id.vars = c("Sample Identification", "material", "group", "time point"))
  r_time_course_sig_diff_dat$time <- r_time_course_sig_diff_dat$`time point`
  r_time_course_sig_diff_dat <- na.omit(subset(r_time_course_sig_diff_dat, variable %in% rat.sig.anova.car.c.class[[mat]]))
  r_time_course_sig_diff_dat$variable <- factor(r_time_course_sig_diff_dat$variable, levels = unique(r_time_course_sig_diff_dat$variable))
  r_time_course_sig_diff_dat$group[r_time_course_sig_diff_dat$group != "control"] <- "septic"
  n_mets <- length(unique(r_time_course_sig_diff_dat$variable))
  r_time_course_group <- unique(data.frame(variable = r_time_course_sig_diff_dat$variable, value = 1.4, Time = 2, text = coarse_group_list[match(r_time_course_sig_diff_dat$variable, colnames(rat_sepsis_data)[-1:-4])]))
  lv <- sapply(anova.car.ph.sig.contr[[mat]], length)
  r_time_course_sig_times <- data.frame(t = "", variable = rep("", sum(lv)), text = "*", value = 1.1, stringsAsFactors = FALSE)
  r_time_course_sig_times$t <- Reduce("c", anova.car.ph.sig.contr[[mat]])
  r_time_course_sig_times$variable <- factor(rep(names(anova.car.ph.sig.contr[[mat]]), times = lv))
  r_time_course_sig_times$Time <- c(2, 1, 3, 2, 1)[r_time_course_sig_times$t]
  r_time_course_sig_times <- subset(r_time_course_sig_times, variable %in% r_time_course_sig_diff_dat$variable & t %in% 1:3)
  r_time_course_sig_diff_plot <- ggplot(r_time_course_sig_diff_dat, aes(x = time, y = value, group = group, color = group, facets = variable)) +
    facet_wrap(facets = ~ variable, ncol = 6, nrow = ceiling(n_mets/6)) +
    #geom_boxplot(mapping = aes_string(x = time, color = group, y = value), inherit.aes = FALSE) +
    geom_point(position = position_dodge(width = 0.3)) +
    stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
    #stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", size = 0.5) +
    geom_text(data = r_time_course_group, mapping = aes(x = Time, y = value, label = text), inherit.aes = FALSE, size = 3.2) +
    geom_text(data = r_time_course_sig_times, mapping = aes(x = Time, y = value, label = text), inherit.aes = FALSE, size = 8) +
    scale_x_discrete(limits = c("6h", "24h", "72h")) +
    ylim(0, 1.5) +
    ylab("Concentration relative to max value") +
    ggtitle(paste0("Significantly different metabolites in Sepsis vs Control at any time point in ", mat)) + 
    theme_bw()
  ggsave(plot = r_time_course_sig_diff_plot, filename = paste0("rat_", mat, "_time_course_sig_c_anova_diff.png"), path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets/6), units = "in")
}

##rat, metabolite concentration time courses, ANOVA results, septic survivors vs septic nonsurvivors
for (mat in names(rat.sig.anova.car.s.class)){
  r_time_course_sig_diff_dat <- subset(rat_sepsis_data, material == mat)
  r_time_course_sig_diff_dat <- max_norm(r_time_course_sig_diff_dat, -1:-4)
  r_time_course_sig_diff_dat <- melt(r_time_course_sig_diff_dat, id.vars = c("Sample Identification", "material", "group", "time point"))
  r_time_course_sig_diff_dat <- subset(r_time_course_sig_diff_dat, variable %in% rat.sig.anova.car.s.class[[mat]] & group != "control")
  n_mets <- length(unique(r_time_course_sig_diff_dat$variable))
  r_time_course_group <- unique(data.frame(variable = r_time_course_sig_diff_dat$variable, value = 1.4, Time = 1.5, text = coarse_group_list[match(r_time_course_sig_diff_dat$variable, colnames(rat_sepsis_data)[-1:-4])]))
  lv <- sapply(anova.car.ph.sig.contr[[mat]], length)
  r_time_course_sig_times <- data.frame(t = "", variable = rep("", sum(lv)), text = "*", value = 1.1, stringsAsFactors = FALSE)
  r_time_course_sig_times$t <- Reduce("c", anova.car.ph.sig.contr[[mat]])
  r_time_course_sig_times$variable <- rep(names(anova.car.ph.sig.contr[[mat]]), times = lv)
  r_time_course_sig_times$Time <- c(2, 1, 3, 2, 1)[r_time_course_sig_times$t]
  r_time_course_sig_times <- subset(r_time_course_sig_times, variable %in% r_time_course_sig_diff_dat$variable & t %in% 4:5)
  r_time_course_sig_diff_plot <- ggplot(na.omit(r_time_course_sig_diff_dat), aes_string(x = "`time point`", y = "value", group = "group", color = "group")) +
    facet_wrap(facets = ~ variable, ncol = 6, nrow = ceiling(n_mets/6)) +
    #geom_boxplot(mapping = aes_string(x = "`time point`", color = "group", y = "value"), inherit.aes = FALSE) +
    geom_point(position = position_dodge(width = 0.3)) +
    stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
    #stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", size = 0.5) +
    geom_text(data = r_time_course_group, mapping = aes(x = Time, y = value, label = text), inherit.aes = FALSE, size = 3.2) +
    geom_text(data = r_time_course_sig_times, mapping = aes(x = Time, y = value, label = text), inherit.aes = FALSE, size = 8) +
    scale_x_discrete(limits = c("6h", "24h")) +
    ylim(0, 1.5) +
    ylab("Concentration relative to max value") +
    ggtitle(paste0("Significantly different metabolites in survival at any time point in ", mat)) + 
    theme_bw()
  ggsave(plot = r_time_course_sig_diff_plot, filename = paste0("rat_", mat, "_time_course_sig_s_anova_diff.png"), path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets/6), units = "in")
}

##rat, pheno var time courses, ANOVA results, control vs. sepsis
r_time_course_sig_diff_dat <- subset(rat_sepsis_data, material == "plasma")
r_time_course_sig_diff_dat <- max_norm(r_time_course_sig_diff_dat, -1:-4)
r_time_course_sig_diff_dat <- melt(r_time_course_sig_diff_dat, id.vars = c("Sample Identification", "material", "group", "time point"))
r_time_course_sig_diff_dat$time <- r_time_course_sig_diff_dat$`time point`
r_time_course_sig_diff_dat <- na.omit(subset(r_time_course_sig_diff_dat, variable %in% rat.sig.anova.car.c.pheno.class))
r_time_course_sig_diff_dat$variable <- factor(r_time_course_sig_diff_dat$variable, levels = unique(r_time_course_sig_diff_dat$variable))
r_time_course_sig_diff_dat$group[r_time_course_sig_diff_dat$group != "control"] <- "septic"
n_mets <- length(unique(r_time_course_sig_diff_dat$variable))
r_time_course_group <- unique(data.frame(variable = r_time_course_sig_diff_dat$variable, value = 1.4, Time = 2, text = coarse_group_list[match(r_time_course_sig_diff_dat$variable, colnames(rat_sepsis_data)[-1:-4])]))
lv <- sapply(anova.car.pheno.ph.sig.contr, length)
r_time_course_sig_times <- data.frame(t = "", variable = rep("", sum(lv)), text = "*", value = 1.1, stringsAsFactors = FALSE)
r_time_course_sig_times$t <- Reduce("c", anova.car.pheno.ph.sig.contr)
r_time_course_sig_times$variable <- factor(rep(names(anova.car.pheno.ph.sig.contr), times = lv))
r_time_course_sig_times$Time <- c(2, 1, 3, 2, 1)[r_time_course_sig_times$t]
r_time_course_sig_times <- subset(r_time_course_sig_times, variable %in% r_time_course_sig_diff_dat$variable & t %in% 1:3)
r_time_course_sig_diff_plot <- ggplot(r_time_course_sig_diff_dat, aes(x = time, y = value, group = group, color = group, facets = variable)) +
  facet_wrap(facets = ~ variable, ncol = 6, nrow = ceiling(n_mets/6)) +
  #geom_boxplot(mapping = aes_string(x = time, color = group, y = value), inherit.aes = FALSE) +
  geom_point(position = position_dodge(width = 0.3)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  #stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", size = 0.5) +
  geom_text(data = r_time_course_group, mapping = aes(x = Time, y = value, label = text), inherit.aes = FALSE, size = 3.2) +
  geom_text(data = r_time_course_sig_times, mapping = aes(x = Time, y = value, label = text), inherit.aes = FALSE, size = 8) +
  scale_x_discrete(limits = c("6h", "24h", "72h")) +
  ylim(0, 1.5) +
  ylab("Concentration relative to max value") +
  ggtitle("Significantly different clinical params in Sepsis vs Control at any time point") + 
  theme_bw()
ggsave(plot = r_time_course_sig_diff_plot, filename = paste0("rat_pheno_time_course_sig_c_anova_diff.png"), path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets/6), units = "in")

##rat, metabolite concentration time courses, ANOVA results, septic survivors vs septic nonsurvivors
r_time_course_sig_diff_dat <- subset(rat_sepsis_data, material == "plasma")
r_time_course_sig_diff_dat <- max_norm(r_time_course_sig_diff_dat, -1:-4)
r_time_course_sig_diff_dat <- melt(r_time_course_sig_diff_dat, id.vars = c("Sample Identification", "material", "group", "time point"))
r_time_course_sig_diff_dat <- subset(r_time_course_sig_diff_dat, variable %in% rat.sig.anova.car.s.pheno.class & group != "control")
n_mets <- length(unique(r_time_course_sig_diff_dat$variable))
r_time_course_group <- unique(data.frame(variable = r_time_course_sig_diff_dat$variable, value = 1.4, Time = 1.5, text = coarse_group_list[match(r_time_course_sig_diff_dat$variable, colnames(rat_sepsis_data)[-1:-4])]))
lv <- sapply(anova.car.pheno.ph.sig.contr, length)
r_time_course_sig_times <- data.frame(t = "", variable = rep("", sum(lv)), text = "*", value = 1.1, stringsAsFactors = FALSE)
r_time_course_sig_times$t <- Reduce("c", anova.car.pheno.ph.sig.contr)
r_time_course_sig_times$variable <- rep(names(anova.car.pheno.ph.sig.contr), times = lv)
r_time_course_sig_times$Time <- c(2, 1, 3, 2, 1)[r_time_course_sig_times$t]
r_time_course_sig_times <- subset(r_time_course_sig_times, variable %in% r_time_course_sig_diff_dat$variable & t %in% 4:5)
r_time_course_sig_diff_plot <- ggplot(na.omit(r_time_course_sig_diff_dat), aes_string(x = "`time point`", y = "value", group = "group", color = "group")) +
  facet_wrap(facets = ~ variable, ncol = 6, nrow = ceiling(n_mets/6)) +
  #geom_boxplot(mapping = aes_string(x = "`time point`", color = "group", y = "value"), inherit.aes = FALSE) +
  geom_point(position = position_dodge(width = 0.3)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  #stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", size = 0.5) +
  geom_text(data = r_time_course_group, mapping = aes(x = Time, y = value, label = text), inherit.aes = FALSE, size = 3.2) +
  geom_text(data = r_time_course_sig_times, mapping = aes(x = Time, y = value, label = text), inherit.aes = FALSE, size = 8) +
  scale_x_discrete(limits = c("6h", "24h")) +
  ylim(0, 1.5) +
  ylab("Concentration relative to max value") +
  ggtitle("Significantly different clinical params in survival at any time point") + 
  theme_bw()
ggsave(plot = r_time_course_sig_diff_plot, filename = paste0("rat_pheno_time_course_sig_s_anova_diff.png"), path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets/6), units = "in")

##rat, plasma, metabolite concentration time course, only metabolites with p-val < 0.05 at any time point
r_time_course_sig_diff_dat <- subset(rat_sepsis_data_long_form_sig, group %in% c("septic survivor", "septic non-survivor") & material == "plasma" & variable %in% rat_plasma_sig_t_class)
r_time_course_sigs <- subset(rat_plasma_time_sig_t_diff_pos_long, variable %in% r_time_course_sig_diff_dat$variable)
r_time_course_sigs$value <- 1.2
r_time_course_group <- data.frame(variable = sort(unique(r_time_course_sig_diff_dat$variable)), value = 1.5, Time = 2.4, text = coarse_group_list[match(sort(unique(r_time_course_sig_diff_dat$variable)), colnames(rat_sepsis_data)[-1:-4])])
r_time_course_sig_diff_plot <- ggplot(na.omit(r_time_course_sig_diff_dat), aes_string(x = "`time point`", y = "value", group = "group", color = "group")) +
  facet_wrap(facets = ~ variable, ncol = 5, nrow = ceiling(length(rat_plasma_time_sig_t_diff)/5)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", size = 0.5) +
  geom_point(data = r_time_course_sigs, mapping = aes(x = Time, y = value), shape = 8, inherit.aes = FALSE, size = 0.8) +
  geom_text(data = r_time_course_group, mapping = aes(x = Time, y = value, label = text), inherit.aes = FALSE, size = 3.2) +
  scale_x_discrete(limits = c("6h", "24h", "72h")) +
  ylim(0, 1.7) +
  ylab("Concentration relative to max value") +
  ggtitle("Metabolites significantly differing for survival at any time point in plasma") + 
  theme_bw()
ggsave(plot = r_time_course_sig_diff_plot, filename = "rat_metab_plasma_time_course_sig_diff.png", path = out_dir, width = 10, height = 4, units = "in")

##rat, liver, metabolite concentration time course, only metabolites with p-val < 0.05 at any time point
r_time_course_sig_diff_dat <- subset(rat_sepsis_data_long_form_sig, group %in% c("septic survivor", "septic non-survivor") & material == "liver" & variable %in% rat_liver_sig_t_class)
r_time_course_sigs <- subset(rat_liver_time_sig_t_diff_pos_long, variable %in% r_time_course_sig_diff_dat$variable)
m_val <- max(r_time_course_sig_diff_dat$value)
r_time_course_sigs$value <- m_val + 0.1
r_time_course_group <- data.frame(variable = sort(unique(r_time_course_sig_diff_dat$variable)), value = m_val + 0.3, Time = 2.4, text = coarse_group_list[match(sort(unique(r_time_course_sig_diff_dat$variable)), colnames(rat_sepsis_data)[-1:-4])])
r_time_course_sig_diff_plot <- ggplot(na.omit(r_time_course_sig_diff_dat), aes_string(x = "`time point`", y = "value", group = "group", color = "group")) +
  facet_wrap(facets = ~ variable, ncol = 5, nrow = ceiling(length(rat_liver_time_sig_t_diff)/5)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", size = 0.5) +
  geom_point(data = r_time_course_sigs, mapping = aes(x = Time, y = value), shape = 8, inherit.aes = FALSE, size = 0.8) +
  geom_text(data = r_time_course_group, mapping = aes(x = Time, y = value, label = text), inherit.aes = FALSE, size = 3.2) +
  scale_x_discrete(limits = c("6h", "24h", "72h")) +
  ylim(0, m_val + 0.4) +
  ylab("Concentr. relat. to max value") +
  ggtitle("Metabolites significantly differing for survival at any time point in liver") + 
  theme_bw()
ggsave(plot = r_time_course_sig_diff_plot, filename = "rat_metab_liver_time_course_sig_diff.png", path = out_dir, width = 8, height = 2.5, units = "in")

##rat, metabolite groups
###Given groups
for (mat in unique(rat_sepsis_data$material)){
  rsg <- rat_sepsis_data_grouped[, -group_pheno_sel]
  n_mets <- ncol(rsg) - 4
  rsg$group <- factor(rsg$group, levels = c("septic non-survivor", "control", "septic survivor"))
  p <- ggplot(data = subset(melt(rsg[, -group_pheno_sel], id.vars = 1:4), material %in% mat), mapping = aes_string(x = "`time point`", y = "value", group = "group", colour = "group")) + 
    facet_wrap(facets = ~ variable, ncol = 4, nrow = ceiling(n_mets/4), scales = "free_y") +
    geom_point(position = position_dodge(width = 0.2)) +
    stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
    scale_x_discrete(limits = c("6h", "24h", "72h")) +
    ylab("Concentration, µM") +
    xlab("Day") +
    theme_bw()
  ggsave(plot = p, filename = paste0("rat_metab_", mat, "_group_time_course.png"), path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets/4), units = "in")
}

###Acylcarnitine made up groups
ac_transferase_grouped <- rat_sepsis_data[, c(1:4, which(rat_sepsis_legend$group == "acylcarnitine") + 4)]
ac_transferase_grouped[, -1:-4] <- lapply(ac_transferase_grouped[, -1:-4], `/`, ac_transferase_grouped$C0)
actg <- colnames(ac_transferase_grouped)
ac_transferase_grouped <- transform(ac_transferase_grouped, 
                                    ShortChainAC = C2 + C3 + C4, 
                                    MediumChainAC = C5 + `C6 (C4:1-DC)` + C8 + C9 + C10 + C12, 
                                    LongChainAC = C14 + C16 + C18)
colnames(ac_transferase_grouped)[1:length(actg)] <- actg
ac_transferase_grouped$group <- factor(ac_transferase_grouped$group, levels = c("septic non-survivor", "control", "septic survivor"))
n_mets <- ncol(ac_transferase_grouped) - 1
p_ncol <- 4
for (mat in unique(rat_sepsis_data$material)){
  p <- ggplot(data = subset(melt(ac_transferase_grouped[, -5], id.vars = 1:4), material %in% mat), mapping = aes_string(x = "`time point`", y = "value", group = "group", colour = "group")) + 
    facet_wrap(facets = ~ variable, ncol = p_ncol, nrow = ceiling(n_mets/p_ncol), scales = "free_y") +
    geom_point(position = position_dodge(width = 0.2)) +
    stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
    scale_x_discrete(limits = c("6h", "24h", "72h")) +
    ylab("Concentration relative to L-Carnitine") +
    xlab("Time point") +
    theme_bw()
  ggsave(plot = p, filename = paste0("rat_metab_", mat, "_AC_group_time_course.png"), path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets/p_ncol), units = "in")
}

###Fatty acid carrier ratios
fa_type_grouped <- rat_sepsis_data[, c(1:4, grep(pattern = "PC|SM", x = rat_sepsis_legend[,1]) + 4)]
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
fa_type_grouped <- fa_type_grouped[, c(1:5, which(colnames(fa_type_grouped) == "PCaaToSM"):ncol(fa_type_grouped))]
fa_type_grouped$group <- factor(fa_type_grouped$group, levels = c("septic non-survivor", "control", "septic survivor"))
n_mets <- ncol(fa_type_grouped)
p_ncol <- 4
for (mat in unique(rat_sepsis_data$material)){
  p <- ggplot(data = subset(melt(fa_type_grouped, id.vars = 1:5), material %in% mat), mapping = aes_string(x = "`time point`", y = "value", group = "group", colour = "group")) + 
    facet_wrap(facets = ~ variable, ncol = p_ncol, nrow = ceiling(n_mets/p_ncol), scales = "free_y") +
    geom_point(position = position_dodge(width = 0.2)) +
    stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
    scale_x_discrete(limits = c("6h", "24h", "72h")) +
    ylab("Concentration ratio") +
    xlab("Time point") +
    theme_bw()
  ggsave(plot = p, filename = paste0("rat_metab_", mat, "_FA_carrier_group_time_course.png"), path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets/p_ncol), units = "in")
}

###PC UFA grouped
ufa_grouped <- rat_sepsis_data[, c(1:4, which(rat_sepsis_legend$group == "phosphatidylcholine"))]
ufas <- unique(sub(pattern = "C[0-9]{2}", replacement = "C..", x = colnames(ufa_grouped)[-1:-4]))
ufa_sum <- as.data.frame(sapply(lapply(ufas, grep, x = colnames(ufa_grouped)), function(cols) rowSums(ufa_grouped[, cols])))
colnames(ufa_sum) <- sub(pattern = "C..:", replacement = "CXX:", x = ufas, fixed = TRUE)
ufa_sum <- cbind(ufa_grouped[, 1:4], ufa_sum)
ufa_sum$group <- factor(ufa_sum$group, levels = c("septic non-survivor", "control", "septic survivor"))
n_mets <- ncol(ufa_sum) - 4
n_col <- 7
for (mat in unique(ufa_sum$material)){
  p <- ggplot(data = subset(melt(ufa_sum, id.vars = 1:4), material %in% mat), mapping = aes_string(x = "`time point`", y = "value", group = "group", colour = "group")) + 
    facet_wrap(facets = ~ variable, ncol = n_col, nrow = ceiling(n_mets/n_col), scales = "free_y") +
    geom_point(position = position_dodge(width = 0.2)) +
    stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
    ylab("Concentration, µM") +
    xlab("Time point") +
    scale_x_discrete(limits = c("6h", "24h", "72h")) +
    theme_bw()
  ggsave(plot = p, filename = paste0("rat_metab_PC_", mat, "_group_time_course.png"), path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets/n_col), units = "in")
}

###PC UFA grouped, ratio of aa to ae
ufa_grouped <- rat_sepsis_data[, c(1:4, which(rat_sepsis_legend$group == "phosphatidylcholine"))]
ufas <- unique(sub(pattern = "C[0-9]{2}", replacement = "C..", x = colnames(ufa_grouped)[-1:-4]))
ufa_sum <- as.data.frame(sapply(lapply(ufas, grep, x = colnames(ufa_grouped)), function(cols) rowSums(ufa_grouped[, cols])))
ufa_ratio <- as.data.frame(sapply(1:7, function(col) ufa_sum[, col] / ufa_sum[, col + 7]))
colnames(ufa_ratio) <- sub(pattern = "PC aa C..:", replacement = "CXX:", x = ufas[1:7], fixed = TRUE)
ufa_ratio <- cbind(ufa_grouped[, 1:4], ufa_ratio)
ufa_ratio$group <- factor(ufa_ratio$group, levels = c("septic non-survivor", "control", "septic survivor"))
n_mets <- ncol(ufa_ratio) - 4
n_col <- 7
for (mat in unique(ufa_ratio$material)){
  p <- ggplot(data = subset(melt(ufa_ratio, id.vars = 1:4), material %in% mat), mapping = aes_string(x = "`time point`", y = "value", group = "group", colour = "group")) + 
    facet_wrap(facets = ~ variable, ncol = n_col, nrow = ceiling(n_mets/n_col), scales = "free_y") +
    geom_point(position = position_dodge(width = 0.2)) +
    stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
    ylab("Ratio PC aa/PCae") +
    xlab("Time point") +
    scale_x_discrete(limits = c("6h", "24h", "72h")) +
    theme_bw()
  ggsave(plot = p, filename = paste0("rat_metab_PC_rel_", mat, "_group_time_course.png"), path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets/n_col), units = "in")
}

##rat, concentration ratios survivor vs nonsurvivor
# n_col <- 10
# n_mets <- length(unique(rat_ratio_ci$metabolite))
# p <- ggplot(data = rat_ratio_ci, mapping = aes(x = interaction(material, time), y = estimate)) + 
#   facet_wrap(facets = ~ metabolite, nrow = ceiling(n_mets / n_col), ncol = n_col) +
#   geom_point() +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90))
# ggsave(plot = p, file = "rat_ac_snsratios_all_mats.png", path = out_dir, width = 12, height = 0.3 + 1.5 * ceiling(n_mets / n_col), units = "in")

p <- ggplot(data = subset(rat_ratio_ci, metabolite %in% c("C0", "C2", "C4", "C4:1", "C3-DC (C4-OH)", "C6 (C4:1-DC)", "C5-DC (C6-OH)", "C6:1", "C8", "C10", "C10:1", "C12", "C12:1", "C12-OH", "C14", "C14:1", "C14-OH", "C16", "C16:1", "C16-OH")), 
            mapping = aes(x = liver_6h_estimate, 
                          xmin = liver_6h_lower, 
                          xmax = liver_6h_upper, 
                          y = plasma_6h_estimate, 
                          ymin = plasma_6h_lower, 
                          ymax = plasma_6h_upper)) +
  geom_abline(slope = 1, intercept = 0, size = 0.2) +
  geom_point(shape = 5, size = 0.1) +
  geom_errorbar(size = 0.1) +
  geom_errorbarh(size = 0.1) +
  geom_text(mapping = aes(label = metabolite), size = 2, hjust = 0, vjust = 0) +
  geom_text(x = 1.45, y = 1.465, label = "Ideal correspondence", angle = 46, size = 2) +
  xlab("AC ratio + CI, S vs NS in Liver at 6h") +
  ylab("AC ratio + CI, S vs NS in Plasma at 6h") + 
  coord_cartesian(xlim = c(0.75, 2), ylim = c(0.4, 1.5)) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(plot = p, filename = "rat_S_vs_NS_ac_ratios_liver_vs_plasma_6h.png", path = out_dir)

p <- ggplot(data = subset(rat_ratio_ci, metabolite %in% c("C0", "C2", "C4", "C4:1", "C3-DC (C4-OH)", "C6 (C4:1-DC)", "C5-DC (C6-OH)", "C6:1", "C8", "C10", "C10:1", "C12", "C12:1", "C12-OH", "C14", "C14:1", "C14-OH", "C16", "C16:1", "C16-OH")), 
            mapping = aes(x = liver_24h_estimate, 
                          xmin = liver_24h_lower, 
                          xmax = liver_24h_upper, 
                          y = plasma_24h_estimate, 
                          ymin = plasma_24h_lower, 
                          ymax = plasma_24h_upper)) +
  geom_abline(slope = 1, intercept = 0, size = 0.2) +
  geom_point(shape = 5, size = 0.1) +
  geom_errorbar(size = 0.1) +
  geom_errorbarh(size = 0.1) +
  geom_text(mapping = aes(label = metabolite), size = 2, hjust = 0, vjust = 0) +
  geom_text(x = 1.45, y = 1.465, label = "Ideal correspondence", angle = 46, size = 2) +
  xlab("AC ratio + CI, S vs NS in Liver at 24h") +
  ylab("AC ratio + CI, S vs NS in Plasma at 24h") + 
  coord_cartesian(xlim = c(0.75, 2), ylim = c(0.4, 1.5)) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(plot = p, filename = "rat_S_vs_NS_ac_ratios_liver_vs_plasma_24h.png", path = out_dir)

p <- ggplot(data = subset(rat_ratio_ci, metabolite %in% c("C0", "C2", "C4", "C4:1", "C3-DC (C4-OH)", "C6 (C4:1-DC)", "C5-DC (C6-OH)", "C6:1", "C8", "C10", "C10:1", "C12", "C12:1", "C12-OH", "C14", "C14:1", "C14-OH", "C16", "C16:1", "C16-OH")), 
            mapping = aes(x = heart_6h_estimate, 
                          xmin = heart_6h_lower, 
                          xmax = heart_6h_upper, 
                          y = plasma_6h_estimate, 
                          ymin = plasma_6h_lower, 
                          ymax = plasma_6h_upper)) +
  geom_abline(slope = 1, intercept = 0, size = 0.2) +
  geom_point(shape = 5, size = 0.1) +
  geom_errorbar(size = 0.1) +
  geom_errorbarh(size = 0.1) +
  geom_text(mapping = aes(label = metabolite), size = 2, hjust = 0, vjust = 0) +
  geom_text(x = 1.445, y = 1.465, label = "Ideal correspondence", angle = 68, size = 2) +
  xlab("AC ratio + CI, S vs NS in Heart at 6h") +
  ylab("AC ratio + CI, S vs NS in Plasma at 6h") + 
  coord_cartesian(xlim = c(0, 3), ylim = c(0.4, 1.5)) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(plot = p, filename = "rat_S_vs_NS_ac_ratios_heart_vs_plasma_6h.png", path = out_dir)

p <- ggplot(data = subset(rat_ratio_ci, metabolite %in% c("C0", "C2", "C4", "C4:1", "C3-DC (C4-OH)", "C6 (C4:1-DC)", "C5-DC (C6-OH)", "C6:1", "C8", "C10", "C10:1", "C12", "C12:1", "C12-OH", "C14", "C14:1", "C14-OH", "C16", "C16:1", "C16-OH")), 
            mapping = aes(x = heart_24h_estimate, 
                          xmin = heart_24h_lower, 
                          xmax = heart_24h_upper, 
                          y = plasma_24h_estimate, 
                          ymin = plasma_24h_lower, 
                          ymax = plasma_24h_upper)) +
  geom_abline(slope = 1, intercept = 0, size = 0.2) +
  geom_point(shape = 5, size = 0.1) +
  geom_errorbar(size = 0.1) +
  geom_errorbarh(size = 0.1) +
  geom_text(mapping = aes(label = metabolite), size = 2, hjust = 0, vjust = 0) +
  geom_text(x = 1.445, y = 1.465, label = "Ideal correspondence", angle = 68, size = 2) +
  xlab("AC ratio + CI, S vs NS in Heart at 24h") +
  ylab("AC ratio + CI, S vs NS in Plasma at 24h") + 
  coord_cartesian(xlim = c(0, 3), ylim = c(0.4, 1.5)) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(plot = p, filename = "rat_S_vs_NS_ac_ratios_heart_vs_plasma_24h.png", path = out_dir)

p <- ggplot(data = subset(rat_ratio_ci, metabolite %in% c("C0", "C2", "C4", "C4:1", "C3-DC (C4-OH)", "C6 (C4:1-DC)", "C5-DC (C6-OH)", "C6:1", "C8", "C10", "C10:1", "C12", "C12:1", "C12-OH", "C14", "C14:1", "C14-OH", "C16", "C16:1", "C16-OH")), 
            mapping = aes(x = heart_6h_estimate, 
                          xmin = heart_6h_lower, 
                          xmax = heart_6h_upper, 
                          y = liver_6h_estimate, 
                          ymin = liver_6h_lower, 
                          ymax = liver_6h_upper)) +
  geom_abline(slope = 1, intercept = 0, size = 0.2) +
  geom_point(shape = 5, size = 0.1) +
  geom_errorbar(size = 0.1) +
  geom_errorbarh(size = 0.1) +
  geom_text(mapping = aes(label = metabolite), size = 2, hjust = 0, vjust = 0) +
  geom_text(x = 0.445, y = 0.465, label = "Ideal correspondence", angle = 60, size = 2) +
  xlab("AC ratio + CI, S vs NS in Heart at 6h") +
  ylab("AC ratio + CI, S vs NS in Liver at 6h") + 
  coord_cartesian(xlim = c(0, 3), ylim = c(0.4, 2)) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(plot = p, filename = "rat_S_vs_NS_ac_ratios_heart_vs_liver_6h.png", path = out_dir)

p <- ggplot(data = subset(rat_ratio_ci, metabolite %in% c("C0", "C2", "C4", "C4:1", "C3-DC (C4-OH)", "C6 (C4:1-DC)", "C5-DC (C6-OH)", "C6:1", "C8", "C10", "C10:1", "C12", "C12:1", "C12-OH", "C14", "C14:1", "C14-OH", "C16", "C16:1", "C16-OH")), 
            mapping = aes(x = heart_24h_estimate, 
                          xmin = heart_24h_lower, 
                          xmax = heart_24h_upper, 
                          y = liver_24h_estimate, 
                          ymin = liver_24h_lower, 
                          ymax = liver_24h_upper)) +
  geom_abline(slope = 1, intercept = 0, size = 0.2) +
  geom_point(shape = 5, size = 0.1) +
  geom_errorbar(size = 0.1) +
  geom_errorbarh(size = 0.1) +
  geom_text(mapping = aes(label = metabolite), size = 2, hjust = 0, vjust = 0) +
  geom_text(x = 0.445, y = 0.465, label = "Ideal correspondence", angle = 60, size = 2) +
  xlab("AC ratio + CI, S vs NS in Heart at 24h") +
  ylab("AC ratio + CI, S vs NS in Liver at 24h") + 
  coord_cartesian(xlim = c(0, 3), ylim = c(0.4, 2)) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(plot = p, filename = "rat_S_vs_NS_ac_ratios_heart_vs_liver_24h.png", path = out_dir)

#TODO: add comparisons of sham rat vs septic S rat ratios, maybe also sham vs septic NS

ml_human_best_feat_set <- c("C4", "lysoPC a C28:0", "lysoPC a C28:1")
data_feat_set <- subset(rat_sepsis_data, material %in% c("liver", "plasma") & `time point` %in% c("6h", "24h"))
data_feat_set <- cbind(data_feat_set[, 1:4], data_feat_set[ml_human_best_feat_set])
data_feat_set <- melt(data_feat_set, id.vars = 1:4)
human_best_feat_plot <- ggplot(data = data_feat_set, mapping = aes(y = value, x = material, color = group)) +
  facet_grid(facets = `time point` ~ variable) + 
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = "dodge") +
  scale_discrete_manual(values = scales::hue_pal()(3), limits = unique(data_feat_set$group)[c(3, 1, 2)], aesthetics = c("fill", "color")) +
  scale_y_log10() +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.direction = "horizontal", legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = "human_survival_best_feats_plasma_liver.png", plot = human_best_feat_plot, path = out_dir, width = 6, height = 6, units = "in")

##Rat, metabolites vs survival, first measurement only
# rp1 <- ggplot(data = subset(x = rat_sepsis_data_long_form, subset = `time point` == "6h"), mapping = aes(x = group, y = value)) +
#   facet_wrap(facets = ~ variable, nrow = 17, ncol = 16) +
#   geom_boxplot(outlier.size = 0.5) +
#   ylab("Scaled concentration") +
#   ggtitle("Rat data at 6h") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90))
# ggsave(plot = rp1, filename = "rat_meas0_metab_conc_vs_survival.png", path = out_dir, width = 17, height = 17)
# 
# ##Rat, metabolites vs tissue type, first measurement only, no control samples
# rp2 <- ggplot(data = subset(x = rat_sepsis_data_long_form, subset = `time point` == "6h" & group != "control"), mapping = aes(x = material, y = value)) +
#   facet_wrap(facets = ~ variable, nrow = 17, ncol = 16) +
#   geom_boxplot(outlier.size = 0.5) + 
#   ylab("scaled concentration") +
#   ggtitle("Rat data at 6h without control samples") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90))
# ggsave(plot = rp2, filename = "rat_meas0_noncontrol_metab_conc_vs_tissue_type.png", path = out_dir, width = 17, height = 17)