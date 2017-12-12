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
out_dir <- "../../results/data_stats"

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

###########################
#Process data
###########################

#Get data overview
##Get overview of sample distribution along days
human_sig_diff_res <- human_sig_diffs_along_days(human_sepsis_data)
day_sig_u_diff <- human_sig_diff_res$day_sig_u_diff
day_sig_t_diff <- human_sig_diff_res$day_sig_t_diff

#Get day 1 stats
sig_class <- colnames(human_sepsis_data)[-1:-5][day_sig_u_diff[day_sig_u_diff$Day == 1, -1] <= 0.05]
sig_t_class <- colnames(human_sepsis_data)[-1:-5][day_sig_t_diff[day_sig_t_diff$Day == 1, -1] <= 0.05]

#Find time points in long time course data where changes are significant
##Construct tANOVA arguments; tANOVA requires equal number of replicates for each time point of a factor level
full_tanova_data <- subset(human_sepsis_data, Day %in% rownames(table(Day))[table(Day) > 1])
tanova_day_set <- c(0,1,2,3,5,7,14)
#tanova_day_set <- c(0,1,2,3,5)
tanova_patient_set <- unique(full_tanova_data$Patient)
for (n in tanova_day_set){
  tanova_patient_set <- intersect(tanova_patient_set, full_tanova_data$Patient[full_tanova_data$Day == n])
}
full_tanova_data <- subset(full_tanova_data, Patient %in% tanova_patient_set & Day %in% tanova_day_set)
full_tanova_data <- full_tanova_data[order(full_tanova_data$Survival),]
t_data <- t(scale(full_tanova_data[,-1:-5]))
f1 <- as.numeric(as.factor(full_tanova_data$Survival))
f2 <- 0
tp <- as.numeric(as.factor(full_tanova_data$Day))
##Clean metabolite data
t_data <- na.omit(t_data)
##Run tANOVA
#tanova_res <- tanova(data = t_data, f1 = f1, f2 = f2, tp = tp, test.type = 2, robustify = FALSE, eb = FALSE, B = 1000)
#tanova_sig_metabs <- tanova_res$obj$gene.order[tanova_res$obj$pvalue <= 0.05]

#Build long format tables for plotting
##Scale measurement values at maximum concentration
human_sepsis_data_max_norm <- max_norm(x = human_sepsis_data, subset = -1:-5)
rat_sepsis_data_max_norm <- max_norm(x = rat_sepsis_data, subset = -1:-4)
##Melt scaled data into long form
human_sepsis_data_long_form <- melt(data = human_sepsis_data_max_norm, id.vars = c("Sample ID", "Patient", "Day", "Survival", "CAP / FP"))
rat_sepsis_data_long_form <- melt(data = rat_sepsis_data_max_norm, id.vars = c("Sample Identification", "material", "group", "time point"))
##Reduce to significantly different metabolites
human_sepsis_data_long_form_sig <- subset(human_sepsis_data_long_form, subset = as.character(variable) %in% sig_class)
human_sepsis_data_long_form_sig$variable <- factor(as.character(human_sepsis_data_long_form_sig$variable), levels = sig_class, ordered = TRUE)

#Build tables for cluster-heatmaps
##Scale mesaurement values by standardization
human_sepsis_data_normal <- human_sepsis_data
human_sepsis_data_normal[,-1:-5] <- scale(x = human_sepsis_data[,-1:-5])
human_sepsis_data_normal$Patient <- as.factor(human_sepsis_data_normal$Patient)
human_sepsis_data_normal$Day <- as.factor(human_sepsis_data_normal$Day)
##Group metabolites
coarse_group_list <- human_sepsis_legend[human_sepsis_legend[,1] %in% colnames(human_sepsis_data[,-1:-5]), 3]
human_sepsis_data_normal_grouped <- cbind(human_sepsis_data[,1:5], matrix(0, nrow = nrow(human_sepsis_data_normal), ncol=length(unique(coarse_group_list))))
colnames(human_sepsis_data_normal_grouped)[-1:-5] <- unique(coarse_group_list)
human_sepsis_data_normal_grouped$Patient <- as.factor(human_sepsis_data_normal_grouped$Patient)
human_sepsis_data_normal_grouped$Day <- as.factor(human_sepsis_data_normal_grouped$Day)
for (n in 1:nrow(human_sepsis_data_normal)){
  human_sepsis_data_normal_grouped[n, -1:-5] <- aggregate(t(human_sepsis_data_normal[n,-1:-5]), by = list(coarse_group_list), FUN = mean, na.action = na.omit)[,2]
}
##Split for metabolites and "phenotypical" factors
split_start <- which(colnames(human_sepsis_data) == "Urea")
pheno_sel <- split_start:ncol(human_sepsis_data)
metab_sel <- 6:(split_start-1)
group_pheno_sel <- which(colnames(human_sepsis_data_normal_grouped) %in% colnames(human_sepsis_data[,pheno_sel]))
group_metab_sel <- which(colnames(human_sepsis_data_normal_grouped) %in% unique(coarse_group_list[metab_sel - 5]))
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
human_sepsis_data_normal_metab_cor <- cbind(human_sepsis_data_normal[, 1:5], cov(t(human_sepsis_data_normal[, metab_sel]), use = "pairwise.complete.obs"))
human_sepsis_data_normal_pheno_cor <- cbind(human_sepsis_data_normal[, 1:5], cov(t(human_sepsis_data_normal[, pheno_sel]), use = "pairwise.complete.obs"))
human_sepsis_data_normal_pheno_cor <- human_sepsis_data_normal_pheno_cor[, !apply(human_sepsis_data_normal_pheno_cor, 2, function(x){all(is.na(x))})]

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


####################
#Plot data
####################

##Human, cluster-heatmap, coarse grouped metabolites
x <- human_sepsis_data_normal_metab_cov
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_grouped_metab.png"))
x <- na.omit(human_sepsis_data_normal_pheno_cov)
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_grouped_pheno.png"))
rm("x")

##Human, cluster-heatmap of non-CAP and non-FP patients, coarse grouped metabolites; error because nothing es left
#heatmaply(x = subset(human_sepsis_data_normal_grouped, subset = human_sepsis_data_normal_grouped$`CAP / FP` %in% c("-"), select = c(-1,-2,-3,-5)))

##Human, cluster-heatmap of patient covariance matrix, coarse grouped metabolites, survival marked
x <- human_sepsis_data_normal_grouped_metab_cov
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_grouped_metab_cov.png"))
x <- na.omit(human_sepsis_data_normal_grouped_pheno_cov)
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_grouped_pheno_cov.png"))
rm("x")

##Human, cluster-heatmap of patient covariance matrix, ungrouped metabolites, survival marked
x <- human_sepsis_data_normal_metab_cov
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_metab_cov.png"))
x <- na.omit(human_sepsis_data_normal_pheno_cov)
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_pheno_cov.png"))
rm("x")

##Human, cluster-heatmat of patient correlation matrix, coarse grouped metabolites, survival and CAP/FP marked
x <- human_sepsis_data_normal_grouped_metab_cor
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_grouped_metab_cor.png"))
x <- na.omit(human_sepsis_data_normal_grouped_pheno_cor)
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_grouped_pheno_cor.png"))
rm("x")

##Human, cluster-heatmat of patient correlation matrix, coarse grouped metabolites, survival and CAP/FP marked
x <- human_sepsis_data_normal_metab_cor
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_metab_cor.png"))
x <- na.omit(human_sepsis_data_normal_pheno_cor)
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_pheno_cor.png"))
rm("x")

##Human, cluster-heatmap of metabolite covariance matrix, survival-ignorant
x <- human_sepsis_data_normal_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_conc_metab_cov.png"))
x <- human_sepsis_data_normal_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_conc_pheno_cov.png"))
rm("x")

##Human, cluster-heatmap of metabolite covariance matrix, survival-regardent
x <- human_sepsis_data_normal_NS_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_NS_conc_metab_cov.png"))
x <- human_sepsis_data_normal_S_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_S_conc_metab_cov.png"))
x <- human_sepsis_data_normal_NS_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_NS_conc_pheno_cov.png"))
x <- human_sepsis_data_normal_S_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_S_conc_pheno_cov.png"))
rm("x")

##Human, cluster-heatmap of covariance matrix of grouped metabolites, survival-ignorant
x <- human_sepsis_data_normal_grouped_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_grouped_conc_metab_cov.png"))
x <- human_sepsis_data_normal_grouped_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_grouped_conc_pheno_cov.png"))
rm("x")

##Human, cluster-heatmap of covariance matrix of grouped metabolites, survival-regardent
x <- human_sepsis_data_normal_NS_grouped_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_NS_grouped_conc_metab_cov.png"))
x <- human_sepsis_data_normal_S_grouped_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_S_grouped_conc_metab_cov.png"))
x <- human_sepsis_data_normal_NS_grouped_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_NS_grouped_conc_pheno_cov.png"))
x <- human_sepsis_data_normal_S_grouped_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_S_grouped_conc_pheno_cov.png"))
rm("x")


##Human, p-val (t-test) plot of differences between non-survivors and survivors, grouped by day
h_day_sig_t_diff_plot <- ggplot(melt(day_sig_t_diff, id.vars = 1), aes(x = Day, y = value)) +
  geom_point() +
  geom_hline(yintercept = 0.05) +
  scale_y_log10() +
  ylab("p value (t-test, FDR-corrected)") +
  theme_bw()
ggsave(plot = h_day_sig_t_diff_plot, filename = "human_all_days_survival_sig_t_diff.png", path = out_dir, width = 6, height = 6, units = "in")

##Human, p-val (t-test) plot of differences between non-survivors and survivors, grouped by day
h_day_sig_u_diff_plot <- ggplot(melt(day_sig_u_diff, id.vars = 1), aes(x = Day, y = value)) +
  geom_point() +
  geom_hline(yintercept = 0.05) +
  scale_y_log10() +
  ylab("p value (U-test, FDR-corrected)") +
  theme_bw()
ggsave(plot = h_day_sig_u_diff_plot, filename = "human_all_days_survival_sig_u_diff.png", path = out_dir, width = 6, height = 6, units = "in")

##Human, metabolite concentration time course, only metabolites with p-val < 0.05 at day1
h_time_course_sig_diff_plot <- ggplot(subset(na.omit(human_sepsis_data_long_form_sig), Day %in% c(0,1,2,3,5,7,14,21,28)), aes(x = Day, y = value, group = Survival, color = Survival)) +
  facet_wrap(facets = ~ variable, ncol = 7, nrow = ceiling(length(sig_class)/5)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", size = 0.5) +
  ylab("Concentration relative to max value") +
  theme_bw()
ggsave(plot = h_time_course_sig_diff_plot, filename = "human_metab_time_course_sig_diff_at_day1.png", path = out_dir, width = 14, height = 10, units = "in")

##Human, metabolites vs survival, second day only
hp1 <- ggplot(data = subset(x = human_sepsis_data_long_form, subset = Day == 1), mapping = aes(x = Survival, y = value)) + 
  facet_wrap(facets = ~ variable, nrow = 15, ncol = 16) +
  geom_boxplot(outlier.size = 0.5) + 
  ylab("Scaled concentration") +
  ggtitle("Human data from first day") +
  theme_bw()
ggsave(plot = hp1, filename = "human_day1_metab_conc_vs_survival.png", path = out_dir, width = 17, height = 17)

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
