#Load libraries
library(reshape2)
library(data.table)
library(ggplot2)
library(gridExtra)
library(matrixStats)
library(kernlab)
library(heatmaply)
library(missRanger)
library(TANOVA)
library(corpcor)
library(psych)
library(car)
library(VennDiagram)
library(igraph)

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

#Seperate septic and nonseptic patients
human_nonsepsis_data <- human_sepsis_data[human_sepsis_data$`CAP / FP` == "-", ]
human_sepsis_data <- human_sepsis_data[human_sepsis_data$`CAP / FP` != "-", ]

#Impute missing values
colns <- colnames(human_sepsis_data)
colnames(human_sepsis_data) <- make.names(colnames(human_sepsis_data))
human_sepsis_data[,-1:-5] <- missRanger(data = human_sepsis_data[,-1:-5])
colnames(human_sepsis_data) <- colns

#Find outlier sample
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
X11();plot(hclust(dist(human_sepsis_data_normal_grouped[, group_metab_sel])))
##According to cluster plot, sample number 11 is an outlier
human_sepsis_data <- human_sepsis_data[-11, ]

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
##Construct tANOVA arguments; tANOVA requires equal number of replicates for each time point of a factor level
full_tanova_data <- subset(human_sepsis_data, Day %in% rownames(table(Day))[table(Day) > 1])
#tanova_day_set <- c(0,1,2,3,5,7,14)
tanova_day_set <- c(0,1,2,3)
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
tanova_res <- tanova(data = t_data, f1 = f1, f2 = f2, tp = tp, test.type = 2, robustify = TRUE, eb = FALSE, B = 100)
tanova_sig_metabs <- tanova_res$obj$gene.order[tanova_res$obj$pvalue <= 0.05]
sig_tanova_class <- colnames(human_sepsis_data)[-1:-5][tanova_sig_metabs]
##TODO: Do repeated measures ANOVA
t3_data <- scale(full_tanova_data[,-1:-5])
rownames(t3_data) <- paste0(as.character(full_tanova_data$Survival), "Time", tp)
# timeBind <- data.frame(subject = 1:nrow(t3_data))
# timeFactor <- data.frame(Time = NULL, Type = NULL, stringsAsFactors = FALSE)
# for (n in 1:nrow(t3_data)){
#   newSet <- data.frame()
#   colnames(newSet) <- c(paste0("Time", n, "Approximate"), paste0("Time", n, "Real"))
#   timeBind <- cbind(timeBind, newSet)
#   timeFactor <- rbind(timeFactor, data.frame(Time = paste0("Time", n), Type = c("Approximate", "Real")))
# }
# timeBind <- as.matrix(timeBind[,-1])
# timeModel <- lm(timeBind ~ 1)
# timeAnalysis <- Anova(timeModel, idata = timeFactor, idesign = ~Time*Type, type = "III", test.statistic = "Wilks")
# timeResult <- summary(timeAnalysis)

#Build long format tables for plotting
##Scale measurement values at maximum concentration
human_sepsis_data_max_norm <- max_norm(x = human_sepsis_data, subset = -1:-5)
rat_sepsis_data_max_norm <- max_norm(x = rat_sepsis_data, subset = -1:-4)
##Melt scaled data into long form
human_sepsis_data_long_form <- melt(data = human_sepsis_data_max_norm, id.vars = c("Sample ID", "Patient", "Day", "Survival", "CAP / FP"))
rat_sepsis_data_long_form <- melt(data = rat_sepsis_data_max_norm, id.vars = c("Sample Identification", "material", "group", "time point"))
##Reduce to significantly different metabolites
human_sepsis_data_long_form_sig <- subset(human_sepsis_data_long_form, subset = as.character(variable) %in% union(sig_t_class, sig_tanova_class))
human_sepsis_data_long_form_sig$variable <- factor(as.character(human_sepsis_data_long_form_sig$variable), levels = union(sig_t_class, sig_tanova_class), ordered = TRUE)

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

##Find significantly correlating concentrations
human_sepsis_data_normal_conc_metab_corr <- list()
cols_grouped_metab <- colnames(human_sepsis_data)[group_metab_sel]
cols_metab <- colnames(human_sepsis_data)[metab_sel]
#bootstrap_list <- list()
#for (r in 1:20){
  for (d in 0:2){
    corr_dat <- list()
    n_NS <- sum(human_sepsis_data$Survival == "NS" & human_sepsis_data$Day == d)
    n_S <- sum(human_sepsis_data$Survival == "S" & human_sepsis_data$Day == d)
    rand_subset <- sample.int(n = n_S, size = n_NS)
    NS_corr <- corr.test(x = subset(human_sepsis_data_normal, Survival == "NS" & Day == d, select = metab_sel), adjust = "fdr")
    S_corr <- corr.test(x = subset(human_sepsis_data_normal, Survival == "S" & Day == d & 1:n_S %in% rand_subset, select = metab_sel), adjust = "fdr")
    NS_grouped_corr <- corr.test(x = subset(human_sepsis_data_normal_grouped, Survival == "NS" & Day == d, select = group_metab_sel), adjust = "fdr")
    S_grouped_corr <- corr.test(x = subset(human_sepsis_data_normal_grouped, Survival == "S" & Day == d, select = group_metab_sel), adjust = "fdr")
    ###Get significant metabolite pairs
    xy <- which.xy(NS_grouped_corr$p <= 0.05)
    xy <- subset(xy, x < y) # x is row, y is column, so x < y means "upper triangle only"
    corr_dat[["NS_grouped_sig_pairs"]] <- paste0(cols_grouped_metab[xy$x], " ~ ", cols_grouped_metab[xy$y])
    xy <- which.xy(S_grouped_corr$p <= 0.05)
    xy <- subset(xy, x < y)
    corr_dat[["S_grouped_sig_pairs"]] <- paste0(cols_grouped_metab[xy$x], " ~ ", cols_grouped_metab[xy$y])
    xy <- which.xy(NS_corr$p <= 0.05)
    xy <- subset(xy, x < y)
    corr_dat[["NS_sig_pairs"]] <- paste0(cols_metab[xy$x], " ~ ", cols_metab[xy$y])
    xy <- which.xy(S_corr$p <= 0.05)
    xy <- subset(xy, x < y)
    corr_dat[["S_sig_pairs"]] <- paste0(cols_metab[xy$x], " ~ ", cols_metab[xy$y])
    human_sepsis_data_normal_conc_metab_corr[[paste0("Day",d)]] <- corr_dat
  }
#  bootstrap_list[[r]] <- human_sepsis_data_normal_conc_metab_corr
#}
# mean(unlist(lapply(bootstrap_list, function(x){ length(x$Day0$S_sig_pairs) })))
# mean(unlist(lapply(bootstrap_list, function(x){ length(x$Day1$S_sig_pairs) })))
# mean(unlist(lapply(bootstrap_list, function(x){ length(x$Day2$S_sig_pairs) })))
# mean(unlist(lapply(bootstrap_list, function(x){ length(intersect(x$Day0$S_sig_pairs, x$Day0$NS_sig_pairs)) })))
# mean(unlist(lapply(bootstrap_list, function(x){ length(intersect(x$Day1$S_sig_pairs, x$Day1$NS_sig_pairs)) })))
# mean(unlist(lapply(bootstrap_list, function(x){ length(intersect(x$Day2$S_sig_pairs, x$Day2$NS_sig_pairs)) })))
S_sig_pairs_day0_tab <- table(unlist(lapply(bootstrap_list, function(x){ x$Day0$S_sig_pairs })))
S_sig_pairs_day1_tab <- table(unlist(lapply(bootstrap_list, function(x){ x$Day1$S_sig_pairs })))
S_sig_pairs_day2_tab <- table(unlist(lapply(bootstrap_list, function(x){ x$Day2$S_sig_pairs })))
NS_sig_pairs_day0_tab <- table(unlist(lapply(bootstrap_list, function(x){ x$Day0$NS_sig_pairs })))
NS_sig_pairs_day1_tab <- table(unlist(lapply(bootstrap_list, function(x){ x$Day1$NS_sig_pairs })))
NS_sig_pairs_day2_tab <- table(unlist(lapply(bootstrap_list, function(x){ x$Day2$NS_sig_pairs })))

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


##Human, cluster-heatmap, coarse grouped metabolites
x <- na.omit(human_sepsis_data_normal[,c(1:5, metab_sel)])
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], file = paste0(out_dir, "human_normal_metab.png"), main = "Metabolite profile does not cluster survival well", key.title = "Concentration\n(normalized)", showticklabels = TRUE)
x <- na.omit(human_sepsis_data_normal[,c(1:5, pheno_sel)])
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], file = paste0(out_dir, "human_normal_pheno.png"), main = "Phenomenological profile has survival clusters", key.title = "Concentration\n(normalized)", showticklabels = TRUE)
x <- na.omit(human_sepsis_data_normal_grouped[,c(1:5, group_metab_sel)])
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], file = paste0(out_dir, "human_normal_metab_grouped.png"), main = "Metablite group profiles somwwhat cluster survival", key.title = "Concentration\n(normalized)", showticklabels = TRUE)
rm("x")


##Human, covariance cluster-heatmap, coarse grouped metabolites
x <- human_sepsis_data_normal_metab_cov
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_metab_cov.png"), main = "Metabolite profile covariance has mainly patient clusters", key.title = "Cov", showticklabels = FALSE)
x <- na.omit(human_sepsis_data_normal_pheno_cov)
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_pheno_cov.png"), main = "Phenomenological profile covariance has patient\n and survival clusters", key.title = "Cov", showticklabels = FALSE)
rm("x")

##Human, cluster-heatmap of non-CAP and non-FP patients, coarse grouped metabolites; error because nothing es left
#heatmaply(x = subset(human_sepsis_data_normal_grouped, subset = human_sepsis_data_normal_grouped$`CAP / FP` %in% c("-"), select = c(-1,-2,-3,-5)))

##Human, cluster-heatmap of patient covariance matrix, coarse grouped metabolites, survival marked
x <- human_sepsis_data_normal_grouped_metab_cov
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_grouped_metab_cov.png"))
x <- na.omit(human_sepsis_data_normal_grouped_pheno_cov)
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_grouped_pheno_cov.png"), main = "Profiles of grouped phenomenological variables cluster nothing")
rm("x")

##Human, cluster-heatmap of patient correlation matrix, ungrouped metabolites, survival marked
x <- human_sepsis_data_normal_metab_cor
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_metab_cor.png"), main = "Metabolite profile correlation gives mainly patient clusters")
x <- na.omit(human_sepsis_data_normal_pheno_cor)
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_pheno_cor.png"), main = "Phenomenological profile correlation gives patient\n and survival clusters")
rm("x")

##Human, cluster-heatmat of patient correlation matrix, coarse grouped metabolites, survival and CAP/FP marked
x <- human_sepsis_data_normal_grouped_metab_cor
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_grouped_metab_cor.png"), main = "Profiles of grouped metabolites cluster nothing")
x <- na.omit(human_sepsis_data_normal_grouped_pheno_cor)
heatmaply(x = x[,-1:-5], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_grouped_pheno_cor.png"), main = "Profiles of grouped phenomenological variables cluster nothing")
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
heatmaply(x = x, file = paste0(out_dir, "human_normal_NS_grouped_conc_metab_cov.png"), key.title = "Cov")
x <- human_sepsis_data_normal_S_grouped_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_S_grouped_conc_metab_cov.png"), key.title = "Cov")
x <- human_sepsis_data_normal_NS_grouped_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_NS_grouped_conc_pheno_cov.png"), key.title = "Cov")
x <- human_sepsis_data_normal_S_grouped_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_S_grouped_conc_pheno_cov.png"), key.title = "Cov")
rm("x")

##Human, Venn diagram of significantly correlating grouped concentrations in sepsis survivors, Day 1
x <- human_sepsis_data_normal_S_grouped_conc_metab_corr$r
p <- human_sepsis_data_normal_S_grouped_conc_metab_corr$p
p[upper.tri(p, diag = TRUE)] <- 1
xy <- which.xy(p <= 0.05)
x_sig_pairs <- paste0(colnames(x)[xy$x], " ~ ", colnames(x)[xy$y])


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

##Human, metabolite concentration time course, only metabolites with p-val < 0.05 at day1
###Prepare significant diff data
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
  ggtitle("Metabolites significantly differing for survival at any time point") + 
  theme_bw()
ggsave(plot = h_time_course_sig_diff_plot, filename = "human_metab_time_course_sig_diff.png", path = out_dir, width = 14, height = 10, units = "in")

##Human, metabolite concentration time course, only metabolites with p-val < 0.05 at day1
###Prepare significant diff data
h_time_course_sig_diff_dat <- subset(na.omit(human_sepsis_data_long_form_sig), Day %in% c(0,1,2,3))
h_time_course_group <- data.frame(variable = sort(unique(h_time_course_sig_diff_dat$variable)), value = 0.9, Day = 2, text = coarse_group_list[match(sort(unique(h_time_course_sig_diff_dat$variable)), colnames(human_sepsis_data)[-1:-5])])
h_time_course_sig_diff_plot <- ggplot(h_time_course_sig_diff_dat, aes(x = Day, y = value, group = Survival, color = Survival)) +
  facet_wrap(facets = ~ variable, ncol = 5, nrow = ceiling(length(union(sig_t_class, sig_tanova_class))/5)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", size = 0.5) +
  geom_text(data = h_time_course_group, mapping = aes(x = Day, y = value, label = text), inherit.aes = FALSE, size = 2.8) +
  ylab("Concentration relative to max value") +
  ggtitle("Metabolites significantly differing for survival by tsANOVA") + 
  theme_bw()
ggsave(plot = h_time_course_sig_diff_plot, filename = "human_metab_time_course_tanova_sig_diff.png", path = out_dir, width = 8, height = 3, units = "in")

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
