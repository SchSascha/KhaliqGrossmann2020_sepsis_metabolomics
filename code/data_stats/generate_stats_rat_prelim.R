#Load libraries
library(reshape2)
library(ggplot2)
library(matrixStats)
library(kernlab)
library(heatmaply)
library(missRanger)
library(TANOVA)

#Import central functions
source("../function_definitions.R")

#Set path
out_dir <- "../../results/data_stats_rat/"

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

##Import corresponding group assignment
rat_sepsis_legend <- get_rat_sepsis_legend()
rat_sepsis_legend$group[rat_sepsis_legend$group == ""] <- rat_sepsis_legend[rat_sepsis_legend$group == "", 1]

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

#Get day 1 stats
cols <- colnames(rat_sepsis_data)[-1:-4]
rat_liver_sig_u_class <- na.omit(cols[rat_liver_time_sig_u_diff[rat_liver_time_sig_u_diff$Time == "6h", -1] <= 0.05])
rat_liver_sig_t_class <- na.omit(cols[rat_liver_time_sig_t_diff[rat_liver_time_sig_t_diff$Time == "6h", -1] <= 0.05])
rat_plasma_sig_u_class <- na.omit(cols[rat_plasma_time_sig_u_diff[rat_plasma_time_sig_u_diff$Time == "6h", -1] <= 0.05])
rat_plasma_sig_t_class <- na.omit(cols[rat_plasma_time_sig_t_diff[rat_plasma_time_sig_t_diff$Time == "6h", -1] <= 0.05])
rat_heart_sig_u_class <- na.omit(cols[rat_heart_time_sig_u_diff[rat_heart_time_sig_u_diff$Time == "6h", -1] <= 0.05])
rat_heart_sig_t_class <- na.omit(cols[rat_heart_time_sig_t_diff[rat_heart_time_sig_t_diff$Time == "6h", -1] <= 0.05])

rat_sig_t_class <- union(union(rat_liver_sig_t_class, rat_plasma_sig_t_class), rat_heart_sig_t_class)

#Find time points in long time course data where changes are significant
##Construct tANOVA arguments; tANOVA requires equal number of replicates for each time point of a factor level
# full_tanova_data <- subset(human_sepsis_data, Day %in% rownames(table(Day))[table(Day) > 1])
# tanova_day_set <- c(0,1,2,3,5,7,14)
# #tanova_day_set <- c(0,1,2,3,5)
# tanova_patient_set <- unique(full_tanova_data$Patient)
# for (n in tanova_day_set){
#   tanova_patient_set <- intersect(tanova_patient_set, full_tanova_data$Patient[full_tanova_data$Day == n])
# }
# full_tanova_data <- subset(full_tanova_data, Patient %in% tanova_patient_set & Day %in% tanova_day_set)
# full_tanova_data <- full_tanova_data[order(full_tanova_data$Survival),]
# t_data <- t(scale(full_tanova_data[,-1:-5]))
# f1 <- as.numeric(as.factor(full_tanova_data$Survival))
# f2 <- 0
# tp <- as.numeric(as.factor(full_tanova_data$Day))
# ##Clean metabolite data
# t_data <- na.omit(t_data)
# ##Run tANOVA
# tanova_res <- tanova(data = t_data, f1 = f1, f2 = f2, tp = tp, test.type = 2, robustify = FALSE, eb = FALSE, B = 1000)
# tanova_sig_metabs <- tanova_res$obj$gene.order[tanova_res$obj$pvalue <= 0.05]

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
split_start <- which(colnames(rat_sepsis_data) == "HR")
pheno_sel <- split_start:ncol(rat_sepsis_data)
metab_sel <- 5:(split_start-1)
group_pheno_sel <- which(colnames(rat_sepsis_data_normal_grouped) %in% colnames(rat_sepsis_data[,pheno_sel]))
group_metab_sel <- which(colnames(rat_sepsis_data_normal_grouped) %in% unique(coarse_group_list[metab_sel - 4]))
##Build covariance matrix with metabolite groups
rat_sepsis_data_normal_grouped_metab_cov <- cbind(rat_sepsis_data_normal_grouped[, 1:4], cov(t(rat_sepsis_data_normal_grouped[,group_metab_sel]), use = "pairwise.complete.obs"))
rat_sepsis_data_normal_grouped_pheno_cov <- cbind(rat_sepsis_data_normal_grouped[, 1:4], cov(t(rat_sepsis_data_normal_grouped[,group_pheno_sel]), use = "pairwise.complete.obs"))
rat_sepsis_data_normal_grouped_pheno_cov <- rat_sepsis_data_normal_grouped_pheno_cov[, !apply(rat_sepsis_data_normal_grouped_pheno_cov, 2, function(x){all(is.na(x))})]
##Build covariance matrix with original metabolites
rat_sepsis_data_normal_metab_cov <- cbind(rat_sepsis_data_normal[, 1:4], cov(t(rat_sepsis_data_normal[,metab_sel]), use = "pairwise.complete.obs"))
rat_sepsis_data_normal_pheno_cov <- cbind(rat_sepsis_data_normal[, 1:4], cov(t(rat_sepsis_data_normal[,pheno_sel]), use = "pairwise.complete.obs"))
rat_sepsis_data_normal_pheno_cov <- rat_sepsis_data_normal_pheno_cov[, !apply(rat_sepsis_data_normal_pheno_cov, 2, function(x){all(is.na(x))})]
##Build correlation matrix with metabolite groups
rat_sepsis_data_normal_grouped_metab_cor <- cbind(rat_sepsis_data_normal_grouped[, 1:4], cor(t(rat_sepsis_data_normal_grouped[,group_metab_sel]), use = "pairwise.complete.obs"))
rat_sepsis_data_normal_grouped_pheno_cor <- cbind(rat_sepsis_data_normal_grouped[, 1:4], cor(t(rat_sepsis_data_normal_grouped[,group_pheno_sel]), use = "pairwise.complete.obs"))
rat_sepsis_data_normal_grouped_pheno_cor <- rat_sepsis_data_normal_grouped_pheno_cor[, !apply(rat_sepsis_data_normal_grouped_pheno_cor, 2, function(x){all(is.na(x))})]
##Build covariance matrix with original metabolites
rat_sepsis_data_normal_metab_cor <- cbind(rat_sepsis_data_normal[, 1:4], cov(t(rat_sepsis_data_normal[, metab_sel]), use = "pairwise.complete.obs"))
rat_sepsis_data_normal_pheno_cor <- cbind(rat_sepsis_data_normal[, 1:4], cov(t(rat_sepsis_data_normal[, pheno_sel]), use = "pairwise.complete.obs"))
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
rat_sepsis_data_normal_grouped_conc_pheno_cov <- cov(rat_sepsis_data_normal_grouped[, group_pheno_sel], use = "pairwise.complete.obs")
###Survival-regardent
rat_sepsis_data_normal_NS_grouped_conc_metab_cov <- cov(subset(rat_sepsis_data_normal_grouped, group == "septic non-survivor", group_metab_sel), use = "pairwise.complete.obs")
rat_sepsis_data_normal_S_grouped_conc_metab_cov <- cov(subset(rat_sepsis_data_normal_grouped, group == "septic survivor", group_metab_sel), use = "pairwise.complete.obs")
rat_sepsis_data_normal_NS_grouped_conc_pheno_cov <- cov(subset(rat_sepsis_data_normal_grouped, group == "septic non-survivor", group_pheno_sel), use = "pairwise.complete.obs")
rat_sepsis_data_normal_S_grouped_conc_pheno_cov <- cov(subset(rat_sepsis_data_normal_grouped, group == "septic survivor", group_pheno_sel), use = "pairwise.complete.obs")


####################
#Plot data
####################

##Human, cluster-heatmap, coarse grouped metabolites
x <- rat_sepsis_data_normal_metab_cov
heatmaply(x = x[,-1:-4], row_side_colors = x[c("group")], file = paste0(out_dir, "rat_normal_grouped_metab.png"), main = "Metabolite profile covariance clusters nothing", key.title = "Cov", showticklabels = FALSE)
x <- na.omit(rat_sepsis_data_normal_pheno_cov)
heatmaply(x = x[,-1:-4], row_side_colors = x[c("group")], file = paste0(out_dir, "rat_normal_grouped_pheno.png"), main = "Phenomenological profile covariance clusters survival/control better", key.title = "Cov", showticklabels = FALSE)
rm("x")

##rat, cluster-heatmap of non-CAP and non-FP patients, coarse grouped metabolites; error because nothing es left
#heatmaply(x = subset(rat_sepsis_data_normal_grouped, subset = rat_sepsis_data_normal_grouped$`CAP / FP` %in% c("-"), select = c(-1,-2,-3,-5)))

##rat, cluster-heatmap of patient covariance matrix, coarse grouped metabolites, survival marked
x <- rat_sepsis_data_normal_grouped_metab_cov
heatmaply(x = x[,-1:-4], row_side_colors = x[c("group")], file = paste0(out_dir, "rat_normal_grouped_metab_cov.png"))
x <- na.omit(rat_sepsis_data_normal_grouped_pheno_cov)
heatmaply(x = x[,-1:-4], row_side_colors = x[c("group")], file = paste0(out_dir, "rat_normal_grouped_pheno_cov.png"), main = "Profiles of grouped phenomenological variables cluster nothing")
rm("x")

##rat, cluster-heatmap of patient correlation matrix, ungrouped metabolites, survival marked
x <- rat_sepsis_data_normal_metab_cor
heatmaply(x = x[,-1:-4], row_side_colors = x[c("group")], file = paste0(out_dir, "rat_normal_metab_cor.png"), main = "Metabolite profile correlation gives mainly patient clusters")
x <- na.omit(rat_sepsis_data_normal_pheno_cor)
heatmaply(x = x[,-1:-4], row_side_colors = x[c("group")], file = paste0(out_dir, "rat_normal_pheno_cor.png"), main = "Phenomenological profile correlation gives patient\n and survival clusters")
rm("x")

##Human, cluster-heatmat of patient correlation matrix, coarse grouped metabolites, survival and CAP/FP marked
x <- rat_sepsis_data_normal_grouped_metab_cor
heatmaply(x = x[,-1:-5], row_side_colors = x[c("group")], file = paste0(out_dir, "rat_normal_grouped_metab_cor.png"), main = "Profiles of grouped metabolites cluster nothing")
x <- na.omit(rat_sepsis_data_normal_grouped_pheno_cor)
heatmaply(x = x[,-1:-5], row_side_colors = x[c("group")], file = paste0(out_dir, "rat_normal_grouped_pheno_cor.png"), main = "Profiles of grouped phenomenological variables cluster nothing")
rm("x")

##rat, cluster-heatmap of metabolite covariance matrix, survival-ignorant
x <- rat_sepsis_data_normal_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "rat_normal_conc_metab_cov.png"), key.title = "Cov")
x <- rat_sepsis_data_normal_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "rat_normal_conc_pheno_cov.png"), key.title = "Cov")
rm("x")

##rat, cluster-heatmap of metabolite covariance matrix, survival-regardent
x <- rat_sepsis_data_normal_NS_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "rat_normal_NS_conc_metab_cov.png"), key.title = "Cov")
x <- rat_sepsis_data_normal_S_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "rat_normal_S_conc_metab_cov.png"), key.title = "Cov")
x <- rat_sepsis_data_normal_NS_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "rat_normal_NS_conc_pheno_cov.png"), key.title = "Cov")
x <- rat_sepsis_data_normal_S_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "rat_normal_S_conc_pheno_cov.png"), key.title = "Cov")
rm("x")

##rat, cluster-heatmap of covariance matrix of grouped metabolites, survival-ignorant
x <- rat_sepsis_data_normal_grouped_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "rat_normal_grouped_conc_metab_cov.png"))
x <- rat_sepsis_data_normal_grouped_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "rat_normal_grouped_conc_pheno_cov.png"))
rm("x")

##rat, cluster-heatmap of covariance matrix of grouped metabolites, survival-regardent
x <- rat_sepsis_data_normal_NS_grouped_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "rat_normal_NS_grouped_conc_metab_cov.png"), key.title = "Cov")
x <- rat_sepsis_data_normal_S_grouped_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "rat_normal_S_grouped_conc_metab_cov.png"), key.title = "Cov")
x <- rat_sepsis_data_normal_NS_grouped_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "rat_normal_NS_grouped_conc_pheno_cov.png"), key.title = "Cov")
x <- rat_sepsis_data_normal_S_grouped_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "rat_normal_S_grouped_conc_pheno_cov.png"), key.title = "Cov")
rm("x")

##Human, p-val (t-test) plot of differences between non-survivors and survivors, grouped by day
rat_plasma_time_sig_t_diff$material <- "plasma"
rat_heart_time_sig_t_diff$material <- "heart"
rat_liver_time_sig_t_diff$material <- "liver"
rat_time_sig_t_diff <- rbind(rat_plasma_time_sig_t_diff, rat_heart_time_sig_t_diff, rat_liver_time_sig_t_diff)
r_time_sig_t_diff_plot <- ggplot(melt(rat_time_sig_t_diff, id.vars = c("Time", "material")), aes(x = Time, y = value)) +
  facet_grid(material ~ .) +
  geom_point() +
  geom_hline(yintercept = 0.05) +
  scale_y_log10() +
  ylab("p value, Not FDR-corrected") +
  ggtitle("Rat metabolites differ at two time points\nin non-/survival, most in plasma") +
  theme_bw()
ggsave(plot = r_time_sig_t_diff_plot, filename = "rat_all_days_survival_sig_diff.png", path = out_dir, width = 4, height = 4, units = "in")

##Human, metabolite concentration time course, only metabolites with p-val < 0.05 at day1
h_time_course_sig_diff_plot <- ggplot(subset(na.omit(human_sepsis_data_long_form_sig), Day %in% c(0,1,2,3,5,7,14,21,28)), aes(x = Day, y = value, group = Survival, color = Survival)) +
  facet_wrap(facets = ~ variable, ncol = 7, nrow = ceiling(length(sig_class)/7)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", size = 0.5) +
  ylab("Concentration relative to max value") +
  ggtitle("Metabolites significantly differing for survival at any time point") + 
  theme_bw()
ggsave(plot = h_time_course_sig_diff_plot, filename = "human_metab_time_course_sig_diff.png", path = out_dir, width = 14, height = 10, units = "in")

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
