#Load libraries
library(reshape2)
library(data.table)
library(ggplot2)
library(matrixStats)
library(kernlab)
library(heatmaply)
library(missRanger)
library(webshot)
library(htmlwidgets)
library(forcats)
library(car)
library(nlme)

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
human_sepsis_legend <- human_sepsis_legend[-1:-5, ]

##Import experiment data
rat_sepsis_data <- get_rat_sepsis_data()
rat_sepsis_data <- rat_sepsis_data[, -which(colnames(rat_sepsis_data) %in% c("HR", "SV", "CO", "EF", "Resp Rate", "Temperature"))]
###Remove outlier rats
rat_sepsis_data <- rat_sepsis_data[!rat_sepsis_data$`Sample Identification` %in% c("60H", "39L"),] 

##Import corresponding group assignment
rat_sepsis_legend <- get_rat_sepsis_legend()
rat_sepsis_legend$group[rat_sepsis_legend$group == ""] <- rat_sepsis_legend[rat_sepsis_legend$group == "", 1]

##Import the metabolite names that Anna found significantly different
annas_mets <- get_annas_rat_sig_diffs()

###########################
#Process data
###########################

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

#Find time points in long time course data where changes are significant
met_set <- colnames(rat_sepsis_data)
met_set <- met_set[colSums(is.na(rat_sepsis_data)) < 2]
met_set <- met_set[-1:-4]
fml <- concentration ~ group*time #TODO: work on formula and random effects specification
anova.car.models <- list()
anova.car.levene <- list()
for (mat in unique(rat_sepsis_data$material)){
  material.models <- list()
  material.levene <- list()
  for (met in met_set){
    test_data <- subset(rat_sepsis_data, material == mat & group != "control", c("group", "time point", met))
    test_data <- test_data[test_data["time point"] != "72h", ]
    colnames(test_data)[ncol(test_data)] <- "concentration"
    colnames(test_data)[2] <- "time"
    test_data$group <- as.factor(test_data$group)
    test_data$time <- as.factor(test_data$time)
    try(material.models[[met]] <- lme(fml, random = ~ 1|time, data = test_data, method = "REML"))
    try(material.levene[[met]] <- leveneTest(fml, data = test_data))
  }
  anova.car.models[[mat]] <- material.models
  anova.car.levene[[mat]] <- material.levene
}
anova.car.var.homog.p <- lapply(anova.car.levene, sapply, function(e) e[[3]][[1]])
anova.car.normality.p <- lapply(anova.car.models, sapply, function(e) shapiro.test(residuals(e))$p.value)
anova.car.terms <- list()
anova.car.ps <- list()
for (mat in unique(rat_sepsis_data$material)){
  anova.car.terms[[mat]] <- rownames(Anova(anova.car.models[[mat]][[1]], type = 3))
  anova.car.ps[[mat]] <- sapply(anova.car.models[[mat]], function(x){ Anova(x, type = 3)[[3]] }, USE.NAMES = TRUE)
  rownames(anova.car.ps[[mat]]) <- anova.car.terms[[mat]]
  anova.car.ps[[mat]] <- apply(anova.car.ps[[mat]], c(2), function(row){ p.adjust(p = row, method = "fdr") })
}
rat_sig_anova.car_class <- lapply(anova.car.ps, function(e) colnames(e)[colAnys(e[c("group", "group:time"), ] <= 0.05)])
for (mat in names(rat_sig_anova.car_class)){
  rat_sig_anova.car_class[[mat]] <- setdiff(rat_sig_anova.car_class[[mat]], names(which(anova.car.var.homog.p[[mat]] < 0.05)))
  rat_sig_anova.car_class[[mat]] <- setdiff(rat_sig_anova.car_class[[mat]], names(which(anova.car.normality.p[[mat]] < 0.05)))
}

#Check correlations between tissue types on ungrouped metabolites
rat_sepsis_data_normal <- rat_sepsis_data
rat_sepsis_data_normal[, -1:-4] <- scale(rat_sepsis_data[-1:-4])
cross_mat_corr <- list()
for (comp in list(c("heart", "plasma"), c("heart", "liver"), c("plasma", "liver"))){
  comp_str <- paste(comp, collapse = " ~ ")
  cross_mat_corr[[comp_str]] <- list()
  for (met in met_set){
    test_data <- subset(rat_sepsis_data_normal, material %in% comp & group != "control", c("Sample Identification", "material", met))
    colnames(test_data)[1] <- "ID"
    test_data$ID <- substr(test_data$ID, start = 1, stop = 3)
    test_data <- dcast(test_data, ID ~ material, value.var = met)
    cross_mat_corr[[comp_str]][[met]] <- cor.test(x = test_data[,2], y = test_data[,3])
  }
}
cross_mat_coeff <- lapply(cross_mat_corr, sapply, function(r){ r$estimate })
cross_mat_p <- lapply(cross_mat_corr, sapply, function(r){ r$p.value })
cross_mat_p <- lapply(cross_mat_p, p.adjust, method = "fdr")

#Check correlations between tissue types on grouped metabolites
coarse_group_list <- rat_sepsis_legend[match(colnames(rat_sepsis_data)[-1:-4], rat_sepsis_legend[,1]), 3] #lucky ... no col in data without match in legend
coarse_group_list[is.na(coarse_group_list)] <- "Plasma Creatinine"
rat_sepsis_data_grouped <- cbind(rat_sepsis_data[,1:4], matrix(0, nrow = nrow(rat_sepsis_data), ncol=length(unique(coarse_group_list))))
colnames(rat_sepsis_data_grouped)[-1:-4] <- unique(coarse_group_list)
for (n in 1:nrow(rat_sepsis_data)){
  agg <- aggregate(t(as.matrix(rat_sepsis_data_normal[n,-1:-4])), by = list(coarse_group_list), FUN = mean, na.action = na.omit, na.rm = T) #output is disordered!
  rat_sepsis_data_grouped[n, -1:-4] <- agg[match(colnames(rat_sepsis_data_grouped)[-1:-4], agg[,1]), 2]
}
cross_mat_group_corr <- list()
for (comp in list(c("heart", "plasma"), c("heart", "liver"), c("plasma", "liver"))){
  comp_str <- paste(comp, collapse = " ~ ")
  cross_mat_group_corr[[comp_str]] <- list()
  for (met in unique(coarse_group_list[1:199])){
    test_data <- subset(rat_sepsis_data_grouped, material %in% comp & group != "!control", c("Sample Identification", "material", met))
    colnames(test_data)[1] <- "ID"
    test_data$ID <- substr(test_data$ID, start = 1, stop = 3)
    test_data <- dcast(test_data, ID ~ material, value.var = met)
    if (all(colSums(is.na(test_data)) < nrow(test_data)))
      cross_mat_group_corr[[comp_str]][[met]] <- cor.test(x = test_data[,2], y = test_data[,3])
  }
}
cross_mat_group_coeff <- lapply(cross_mat_group_corr, sapply, function(r){ r$estimate })
cross_mat_group_p <- lapply(cross_mat_group_corr, sapply, function(r){ r$p.value })
cross_mat_group_p <- lapply(cross_mat_group_p, p.adjust, method = "fdr")


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
  h <- heatmaply(x = xmt[, xm$material == mat], dendrogram = "row", plot_method = "plotly", col_side_colors = subset(xm, material == mat, c("time point", "group")), row_side_colors = data.frame(sig_reg = mat_sigs[[mat]]), key.title = "fold change with\nrespect to control", margins = c(50,100,0,50), subplot_heights = c(0.03, 0.97))
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

##Human, cluster-heatmat of patient correlation matrix, coarse grouped metabolites, survival and CAP/FP marked
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


##rat, metabolite concentration time courses, ANOVA results
for (mat in names(rat_sig_anova.car_class)){
  r_time_course_sig_diff_dat <- subset(rat_sepsis_data, group %in% c("septic survivor", "septic non-survivor") & material == mat)
  r_time_course_sig_diff_dat <- max_norm(r_time_course_sig_diff_dat, -1:-4)
  r_time_course_sig_diff_dat <- melt(r_time_course_sig_diff_dat, id.vars = c("Sample Identification", "material", "group", "time point"))
  r_time_course_sig_diff_dat <- subset(r_time_course_sig_diff_dat, variable %in% rat_sig_anova.car_class[[mat]])
  r_time_course_group <- unique(data.frame(variable = r_time_course_sig_diff_dat$variable, value = 1.2, Time = 1.5, text = coarse_group_list[match(r_time_course_sig_diff_dat$variable, colnames(rat_sepsis_data)[-1:-4])]))
  r_time_course_sig_diff_plot <- ggplot(na.omit(r_time_course_sig_diff_dat), aes_string(x = "`time point`", y = "value", group = "group", color = "group")) +
    facet_wrap(facets = ~ variable, ncol = 6, nrow = ceiling(length(rat_plasma_time_sig_t_diff)/6)) +
    geom_boxplot(mapping = aes_string(x = "`time point`", color = "group", y = "value"), inherit.aes = FALSE) +
    #stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
    #stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", size = 0.5) +
    geom_text(data = r_time_course_group, mapping = aes(x = Time, y = value, label = text), inherit.aes = FALSE, size = 3.2) +
    scale_x_discrete(limits = c("6h", "24h")) +
    ylim(0, 1.4) +
    ylab("Concentration relative to max value") +
    ggtitle(paste0("Metabolites significantly differing for survival at any time point in ", mat)) + 
    theme_bw()
  ggsave(plot = r_time_course_sig_diff_plot, filename = paste0("rat_", mat, "_time_course_sig_anova_diff.png"), path = out_dir, width = 12, height = 0.1 + 0.4 * ceiling(length(rat_plasma_time_sig_t_diff)/6), units = "in")
}

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

##Human, metabolites vs survival, p < 0.05, second day only, 
# hp2 <- ggplot(data = subset(x = human_sepsis_data_long_form_sig, subset = Day == 0), mapping = aes(x = Survival, y = value)) + 
#   facet_wrap(facets = ~ variable, nrow = 5, ncol = 10) +
#   geom_boxplot(outlier.size = 0.5) + 
#   ylab("Scaled concentration") +
#   ggtitle("Human data from second day; metabolites differing by U-test with p < 0.05 after FDR correction") +
#   theme_bw()
# ggsave(plot = hp2, filename = "human_day0_metab_conc_vs_survival_sig.png", path = out_dir, width = 13, height = 7)
# 
# ##Rat, metabolites vs survival, first measurement only
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
