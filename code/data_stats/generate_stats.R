#Load libraries
library(reshape2)
library(ggplot2)
library(matrixStats)
library(heatmaply)
library(kernlab)
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
day_survival_tab <- table(human_sepsis_data[c("Day", "Survival")])
##Keep days with large enough sample count
day_survival_tab <- day_survival_tab[rowMins(day_survival_tab) >= 3,]
##Get significant differences for each day
day_sig_diff <- data.frame(Day = as.numeric(rownames(day_survival_tab)), matrix(0, nrow = nrow(day_survival_tab), ncol = ncol(human_sepsis_data)-5))
colnames(day_sig_diff)[-1] <- colnames(human_sepsis_data)[-1:-5]
day_sig_t_diff <- day_sig_diff
for (d in seq_along(day_sig_diff$Day)){
  #Reduce to significant differences
  day <- day_sig_diff$Day[d]
  for (n in seq_along(day_sig_diff[1,-1:-5])){
    g1 <- subset(human_sepsis_data, subset = Survival == "S" & Day == day, select = n + 5)
    g2 <- subset(human_sepsis_data, subset = Survival == "NS" & Day == day, select = n + 5)
    g1 <- na.omit(g1)[[1]]
    g2 <- na.omit(g2)[[1]]
    if (length(g1) > 1 && length(g2) > 1 && length(table(g1)) > 1 && length(table(g2)) > 1){
      u_res <- wilcox.test(x = g1, y = g2)
      day_sig_diff[d, n + 1] <- u_res$p.value
      t_res <- t.test(x = g1, y = g2, var.equal = FALSE)
      day_sig_t_diff[d, n + 1] <- t_res$p.value
    }
    else{
      day_sig_diff[d, n + 1] <- NA
      day_sig_t_diff[d, n + 1] <- NA
    }
  }
  day_sig_diff[d,-1] <- p.adjust(day_sig_diff[d, -1], method = "fdr")
  day_sig_t_diff[d,-1] <- p.adjust(day_sig_t_diff[d, -1], method = "fdr")
}

#Get day 1 stats
sig_class <- colnames(human_sepsis_data)[-1:-5][day_sig_diff[day_sig_diff$Day == 1, -1] <= 0.05]
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
full_tanova_data <- full_tanova_data[order(full_tanova_data$Patient),]
full_tanova_data <- full_tanova_data[order(full_tanova_data$Survival),]
t_data <- t(scale(full_tanova_data[,-1:-5]))
f1 <- as.numeric(as.factor(full_tanova_data$Survival))
f2 <- 0
tp <- as.numeric(as.factor(full_tanova_data$Day))
##Clean metabolite data
t_data <- na.omit(t_data)
##Run tANOVA
tanova_res <- tanova(data = t_data, f1 = f1, f2 = f2, tp = tp, test.type = 2, robustify = TRUE)
which(tanova_res$obj$pvalue <= 0.05)

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
human_sepsis_data_normal_coarse_grouped <- cbind(human_sepsis_data[,1:5], matrix(0, nrow = nrow(human_sepsis_data_normal), ncol=length(unique(coarse_group_list))))
colnames(human_sepsis_data_normal_coarse_grouped)[-1:-5] <- unique(coarse_group_list)
human_sepsis_data_normal_coarse_grouped$Patient <- as.factor(human_sepsis_data_normal_coarse_grouped$Patient)
human_sepsis_data_normal_coarse_grouped$Day <- as.factor(human_sepsis_data_normal_coarse_grouped$Day)
for (n in 1:nrow(human_sepsis_data_normal)){
  human_sepsis_data_normal_coarse_grouped[n, -1:-5] <- aggregate(t(human_sepsis_data_normal[n,-1:-5]), by = list(coarse_group_list), FUN = mean, na.action = na.omit)[,2]
}
##Build covariance matrix with metabolite groups 
human_sepsis_data_normal_coarse_grouped_cov <- cbind(human_sepsis_data_normal_coarse_grouped[,1:5], cov(t(human_sepsis_data_normal_coarse_grouped[,-1:-5]), use = "pairwise.complete.obs"))
##Build covariance matrix with original metabolites
human_sepsis_data_normal_cov <- cbind(human_sepsis_data_normal[,1:5], cov(t(human_sepsis_data_normal[,-1:-5]), use = "pairwise.complete.obs"))
##Build correlation matrix with metabolite groups
human_sepsis_data_normal_coarse_grouped_cor <- cbind(human_sepsis_data_normal_coarse_grouped[,1:5], cor(t(human_sepsis_data_normal_coarse_grouped[,-1:-5]), use = "pairwise.complete.obs"))
##Build covariance matrix with original metabolites
human_sepsis_data_normal_cor <- cbind(human_sepsis_data_normal[,1:5], cov(t(human_sepsis_data_normal[,-1:-5]), use = "pairwise.complete.obs"))


####################
#Plot data
####################

##Human, cluster-heatmap of CAP and FP patients, coarse grouped metabolites
heatmaply(x = subset(human_sepsis_data_normal_coarse_grouped, select = c(-1,-5)))

##Human, cluster-heatmap of non-CAP and non-FP patients, coarse grouped metabolites; error because nothing es left
#heatmaply(x = subset(human_sepsis_data_normal_coarse_grouped, subset = human_sepsis_data_normal_coarse_grouped$`CAP / FP` %in% c("-"), select = c(-1,-2,-3,-5)))

##Human, cluster-heatmap of patient covariance matrix, coarse grouped metabolites, survival marked
heatmaply(x = subset(human_sepsis_data_normal_coarse_grouped_cov, select = c(-1,-2,-3)))

##Human, cluster-heatmap of patient covariance matrix, ungrouped metabolites, survival marked
heatmaply(x = subset(human_sepsis_data_normal_cov, select = c(-1,-2,-3)))

##Human, cluster-heatmat of patient correlation matrix, coarse grouped metabolites, survival and CAP/FP marked
heatmaply(x = subset(human_sepsis_data_normal_coarse_grouped_cor, select = c(-1,-2,-3)))

##Human, cluster-heatmat of patient correlation matrix, coarse grouped metabolites, survival and CAP/FP marked
heatmaply(x = subset(human_sepsis_data_normal_cor, select = c(-1,-2,-3)))

##Human, p-val (U-test) plot of differences between non-survivors and survivors, grouped by day
h_day_sig_diff_plot <- ggplot(melt(day_sig_diff, id.vars = 1), aes(x = Day, y = value)) +
  geom_point() +
  geom_hline(yintercept = 0.05) +
  scale_y_log10() +
  ylab("p value (U-test, FDR-corrected)") +
  theme_bw()
ggsave(plot = h_day_sig_diff_plot, filename = "human_all_days_survival_sig_diff.png", path = out_dir, width = 6, height = 6, units = "in")

##Human, p-val (t-test) plot of differences between non-survivors and survivors, grouped by day
h_day_sig_t_diff_plot <- ggplot(melt(day_sig_t_diff, id.vars = 1), aes(x = Day, y = value)) +
  geom_point() +
  geom_hline(yintercept = 0.05) +
  scale_y_log10() +
  ylab("p value (t-test, FDR-corrected)") +
  theme_bw()
ggsave(plot = h_day_sig_t_diff_plot, filename = "human_all_days_survival_sig_t_diff.png", path = out_dir, width = 6, height = 6, units = "in")

##Human, p-val (t-test) plot of differences between non-survivors and survivors, grouped by day
h_day_sig_t_diff_plot <- ggplot(melt(day_sig_diff, id.vars = 1), aes(x = Day, y = value)) +
  geom_point() +
  geom_hline(yintercept = 0.05) +
  scale_y_log10() +
  ylab("p value (t-test, FDR-corrected)") +
  theme_bw()
ggsave(plot = h_day_sig_t_diff_plot, filename = "human_all_days_survival_sig_t_diff.png", path = out_dir, width = 6, height = 6, units = "in")

##Human, metabolite concentration time course, only metabolites with p-val < 0.05 at day1
h_time_course_sig_diff_plot <- ggplot(subset(na.omit(human_sepsis_data_long_form_sig), Day %in% c(0,1,2,3,5,7,14,21,28)), aes(x = Day, y = value, group = Survival, color = Survival)) +
  facet_wrap(facets = ~ variable, ncol = 7, nrow = ceiling(length(sig_class)/5)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", size = 0.5) +
  ylab("Concentration relative to max value") +
  theme_bw()
ggsave(plot = h_time_course_sig_diff_plot, filename = "human_metab_time_course_sig_diff_at_day1.png", path = out_dir, width = 14, height = 10, units = "in")

##Human, metabolites vs survival, first day only
hp1 <- ggplot(data = subset(x = human_sepsis_data_long_form, subset = Day == 0), mapping = aes(x = Survival, y = value)) + 
  facet_wrap(facets = ~ variable, nrow = 15, ncol = 16) +
  geom_boxplot(outlier.size = 0.5) + 
  ylab("Scaled concentration") +
  ggtitle("Human data from first day") +
  theme_bw()
ggsave(plot = hp1, filename = "human_day0_metab_conc_vs_survival.png", path = out_dir, width = 17, height = 17)

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
