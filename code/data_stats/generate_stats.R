#Load libraries
library(reshape2)
library(ggplot2)
library(matrixStats)

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

##Import experiment data
rat_sepsis_data <- get_rat_sepsis_data()

#Build long format tables for plotting
##Scale measurement values
human_sepsis_data_max_norm <- max_norm(x = human_sepsis_data, subset = -1:-5)
rat_sepsis_data_max_norm <- max_norm(x = rat_sepsis_data, subset = -1:-4)
##Melt scaled data into long form
human_sepsis_data_long_form <- melt(data = human_sepsis_data_max_norm, id.vars = c("Sample ID", "Patient", "Day", "Survival", "CAP / FP"))
rat_sepsis_data_long_form <- melt(data = rat_sepsis_data_max_norm, id.vars = c("Sample Identification", "material", "group", "time point"))

#Plot data
##Human, metabolites vs survival, first day only
hp1 <- ggplot(data = subset(x = human_sepsis_data_long_form, subset = Day == 0), mapping = aes(x = Survival, y = value)) + 
  facet_wrap(facets = ~ variable, nrow = 15, ncol = 16) +
  geom_boxplot(outlier.size = 0.5) + 
  xlab("scaled concentration") +
  theme_bw()
ggsave(plot = hp1, filename = "human_day0_metab_conc_vs_survival.png", path = out_dir, width = 17, height = 17)

##Rat, metabolites vs survival, first measurement only
rp1 <- ggplot(data = subset(x = rat_sepsis_data_long_form, subset = `time point` == "6h"), mapping = aes(x = group, y = value)) +
  facet_wrap(facets = ~ variable, nrow = 17, ncol = 16) +
  geom_boxplot(outlier.size = 0.5) +
  xlab("scaled concentration") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
ggsave(plot = rp1, filename = "rat_meas0_metab_conc_vs_survival.png", path = out_dir, width = 17, height = 17)

##Rat, metabolites vs tissue type, first measurement only, no control samples
rp2 <- ggplot(data = subset(x = rat_sepsis_data_long_form, subset = `time point` == "6h" & group != "control"), mapping = aes(x = material, y = value)) +
  facet_wrap(facets = ~ variable, nrow = 17, ncol = 16) +
  geom_boxplot(outlier.size = 0.5) + 
  xlab("scaled concentration") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
ggsave(plot = rp2, filename = "rat_meas0_noncontrol_metab_conc_vs_tissue_type.png", path = out_dir, width = 17, height = 17)
