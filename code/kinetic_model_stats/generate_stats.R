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
out_dir <- "../../results/kinetic_model/"

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

##Import experiment data
rat_sepsis_data <- get_rat_sepsis_data()

##Import corresponding group assignment
rat_sepsis_legend <- get_rat_sepsis_legend()

##Import normal plasma ranges
human_healthy_ranges <- get_french_normal_data()

##Import model fit result files
human_mitomod_fit_data <- get_human_kinetic_model_fit_results()

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
human_data[, imputables] <- missRanger(data = human_data[, imputables])
colnames(human_data) <- colns

#Seperate septic and nonseptic patients
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

###########################
#Plot data
###########################

m <- melt(human_mitomod_fit_data)
m <- subset(m, variable == "Value" & stri_startswith_fixed(Parameter, "Values[V"))
m$Group <- c("Survivor", "Nonsurvivor")[1 + (grepl(pattern = "merge", x = m$Parameter))] #Values with [merge] belong to Nonsurvivors
m$Parameter <- stri_match_first(m$Parameter, regex = "\\[.+\\]")
m$Parameter <- stri_replace(str = m$Parameter, replacement = "", fixed = "[merge]")
p <- ggplot(data = m, mapping = aes(x = L1, y = value, group = Group, color = Group)) +
  facet_wrap(Parameter ~ ., ncol = 4) +
  geom_line() +
  scale_y_log10() +
  xlab("Time") +
  theme_bw() +
  theme(legend.direction = "horizontal", legend.position = "top", panel.grid = element_blank())
ggsave(filename = "kin_mitomod_Vmax_S_vs_NS.png", plot = p, path = out_dir, width = 8, height = 4, units = "in")
