#Load libraries
library(reshape2)
library(ggplot2)
library(matrixStats)
library(heatmaply)

#Import central functions
source("../function_definitions.R")

#Set path
out_dir <- "../../results/data_stats"

#Make sure paths exist
if (!dir.exists(out_dir))
  dir.create(out_dir)

#Import data
##Import clinical data
human_sepsis_data <- get_human_sepsis_data()

##Import corresponding group assignment
human_sepsis_legend <- get_human_sepsis_legend()

##Import experiment data
rat_sepsis_data <- get_rat_sepsis_data()

##Import corresponding group assignment
rat_sepsis_legend <- get_rat_sepsis_legend()

#Build long format tables for plotting
##Scale mesaurement values by subtracting mean
human_sepsis_data_normal <- human_sepsis_data
human_sepsis_data_normal[,-1:-5] <- scale(x = human_sepsis_data[,-1:-5], center = TRUE, scale = TRUE)
##Group metabolites
human_sepsis_data_normal_coarse_grouped <- human_sepsis_data[,1:5]
coarse_group_list <- unique(human_sepsis_legend[human_sepsis_legend[,1] %in% colnames(human_sepsis_data[,-1:-5]), 3])



##Scale measurement values at maximum concentration
human_sepsis_data_max_norm <- max_norm(x = human_sepsis_data, subset = -1:-5)
rat_sepsis_data_max_norm <- max_norm(x = rat_sepsis_data, subset = -1:-4)
##Melt scaled data into long form
human_sepsis_data_long_form <- melt(data = human_sepsis_data_max_norm, id.vars = c("Sample ID", "Patient", "Day", "Survival", "CAP / FP"))
rat_sepsis_data_long_form <- melt(data = rat_sepsis_data_max_norm, id.vars = c("Sample Identification", "material", "group", "time point"))
##Reduce to significant differences
human_survival_sig_diff <- rep(0, ncol(human_sepsis_data)-5)
for (n in seq_along(human_survival_sig_diff)){
  g1 <- subset(human_sepsis_data, subset = Survival == "S" & Day == 1, select = n + 5)
  g2 <- subset(human_sepsis_data, subset = Survival == "NS" & Day == 1, select = n + 5)
  g1 <- na.omit(g1)[[1]]
  g2 <- na.omit(g2)[[1]]
  if (length(g1) > 1 && length(g2) > 1){
    u_res <- wilcox.test(x = g1, y = g2)
    human_survival_sig_diff[n] <- u_res$p.value
  }
  else{
    human_survival_sig_diff[n] <- 1
  }
}
human_survival_sig_diff <- p.adjust(human_survival_sig_diff, method = "fdr")
sig_class <- colnames(human_sepsis_data)[-1:-5][human_survival_sig_diff < 0.05]
human_sepsis_data_long_form_sig <- subset(human_sepsis_data_long_form, subset = as.character(variable) %in% sig_class)
human_sepsis_data_long_form_sig$variable <- factor(as.character(human_sepsis_data_long_form_sig$variable), levels = sig_class, ordered = TRUE)

#Plot data
##Human, heatmap like Fig. 2
h_heatmap1 <- heatmaply(x = human_sepsis_data_normal)

##Human, metabolites vs survival, first day only
hp1 <- ggplot(data = subset(x = human_sepsis_data_long_form, subset = Day == 0), mapping = aes(x = Survival, y = value)) + 
  facet_wrap(facets = ~ variable, nrow = 15, ncol = 16) +
  geom_boxplot(outlier.size = 0.5) + 
  ylab("Scaled concentration") +
  title("Human data from first day") +
  theme_bw()
ggsave(plot = hp1, filename = "human_day0_metab_conc_vs_survival.png", path = out_dir, width = 17, height = 17)

##Human, metabolites vs survival, p < 0.05, first day only, 
hp2 <- ggplot(data = subset(x = human_sepsis_data_long_form_sig, subset = Day == 0), mapping = aes(x = Survival, y = value)) + 
  facet_wrap(facets = ~ variable, nrow = 5, ncol = 10) +
  geom_boxplot(outlier.size = 0.5) + 
  ylab("Scaled concentration") +
  title("Human data from first day; metabolites differing by U-test with p < 0.05 after FDR correction") +
  theme_bw()
ggsave(plot = hp2, filename = "human_day0_metab_conc_vs_survival_sig.png", path = out_dir, width = 13, height = 7)

##Rat, metabolites vs survival, first measurement only
rp1 <- ggplot(data = subset(x = rat_sepsis_data_long_form, subset = `time point` == "6h"), mapping = aes(x = group, y = value)) +
  facet_wrap(facets = ~ variable, nrow = 17, ncol = 16) +
  geom_boxplot(outlier.size = 0.5) +
  ylab("Scaled concentration") +
  title("Rat data at 6h") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
ggsave(plot = rp1, filename = "rat_meas0_metab_conc_vs_survival.png", path = out_dir, width = 17, height = 17)

##Rat, metabolites vs tissue type, first measurement only, no control samples
rp2 <- ggplot(data = subset(x = rat_sepsis_data_long_form, subset = `time point` == "6h" & group != "control"), mapping = aes(x = material, y = value)) +
  facet_wrap(facets = ~ variable, nrow = 17, ncol = 16) +
  geom_boxplot(outlier.size = 0.5) + 
  ylab("scaled concentration") +
  title("Rat data at 6h without control samples") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
ggsave(plot = rp2, filename = "rat_meas0_noncontrol_metab_conc_vs_tissue_type.png", path = out_dir, width = 17, height = 17)
