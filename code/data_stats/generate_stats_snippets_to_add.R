library(CCA)
library(hexbin)
library(matrixStats)
library(data.table)
library(missRanger)
library(nscancor)
library(TDAmapper)
library(elasticnet)
library(analogue)
library(corpcor)
library(vegan)
library(ggplot2)
library(ggfortify)
library(parallel)
library(tictoc)

source("../function_definitions.R")

out_dir <- "../../results/data_stats_NS_high_var/"
if (!dir.exists(out_dir))
  dir.create(out_dir)

#Read data
human_data <- get_human_sepsis_data()

human_data_legend <- get_human_sepsis_legend()

#Get variable types
metab_end <- which(colnames(human_data) == "H1")
metab_sel <- 6:metab_end
pheno_start <- metab_end + 1
pheno_sel <- pheno_start:ncol(human_data)

#Impute missing values, not for clinical parameters
col <- colnames(human_data)
colnames(human_data) <- make.names(col)
human_data[, metab_sel] <- missRanger(human_data[, metab_sel])
colnames(human_data) <- col

#Seperate septic and nonseptic patients
human_nonsepsis_data <- human_data[human_data$`CAP / FP` == "-", ]
human_sepsis_data <- human_data[human_data$`CAP / FP` != "-", ]

#Reduce to data where variance is larger in Non-survivors as test
# metab_all_days_var_df <- rbind(data.frame(Survival = "NS", t(colVars(as.matrix(subset(human_sepsis_data, Survival == "NS",-1:-5))))),
#                                data.frame(Survival = "S", t(colVars(as.matrix(subset(human_sepsis_data, Survival == "S", -1:-5))))))
# colnames(metab_all_days_var_df)[-1] <- colnames(human_sepsis_data)[-1:-5]
# metab_all_days_mean <- colMeans(human_sepsis_data[, -1:-5])
# metab_all_days_var_df[, -1] <- scale(metab_all_days_var_df[, -1], center = FALSE, scale = metab_all_days_mean)
# metab_day_vardiff_df <- metab_all_days_var_df[1, -1] - metab_all_days_var_df[2, -1]
# metab_all_days_var_df <- metab_all_days_var_df[c(1, 1 + order(metab_day_vardiff_df, decreasing = TRUE))]
# metab_all_days_var_df[, -1] <- scale(metab_all_days_var_df[, -1], scale = FALSE, center = TRUE)
# metab_all_days_var_df <- metab_all_days_var_df[, c(1, 1 + which(metab_all_days_var_df[1, -1] > 0))]
# var_sel <- colnames(metab_all_days_var_df)[-1]
var_sel <- colnames(human_data)[-1:-5]

human_sepsis_data <- subset(human_sepsis_data, select = c(1:5, which(colnames(human_sepsis_data) %in% var_sel)))
human_data <- subset(human_data, select = c(1:5, which(colnames(human_data) %in% var_sel)))
human_data_legend <- human_data_legend[c(1:5, which(human_data_legend[, 1] %in% var_sel)), ]

metab_end <- which(colnames(human_data) == "H1")
metab_sel <- 6:metab_end
pheno_sel <- (metab_end + 1):ncol(human_data)

#Normalize data
human_sepsis_data_normal <- human_sepsis_data
human_sepsis_data_normal[,-1:-5] <- scale(human_sepsis_data_normal[,-1:-5])
human_sepsis_data_metab <- subset(human_sepsis_data_normal, Day < 4, 6:metab_end)
human_sepsis_data_pheno <- subset(human_sepsis_data_normal, Day < 4, pheno_start:ncol(human_sepsis_data))
human_sepsis_data_normal <- human_sepsis_data_normal[,1:metab_end]

#Generate PCA biplots
##Clinical params, sepsis patients
human_sepsis_data_pheno_dist <- cbind(subset(human_sepsis_data, Day < 4, 1:5), as.matrix(dist(x = subset(human_sepsis_data, Day < 4, pheno_sel), method = "canberra"))) # distance() works on a per row basis
p <- prcomp(human_sepsis_data_pheno_dist[, -1:-5])

###PERMANOVA with bootstrapping
ad_data <- human_sepsis_data_pheno_dist
pat_list <- human_sepsis_data_pheno_dist[, 1:5]
PCs <- p$x[, 1:2]
s_idx <- which(ad_data$Survival == "S")
ns_idx <- which(ad_data$Survival == "NS")
tic()
ad_res_mc <- mclapply(1:100, 
                      function(dx, s_idx, ns_idx, PCs, ad_data){
                        num_s <- length(s_idx)
                        num_ns <- length(ns_idx)
                        num_s_vs_ns <- min(num_s, num_ns)
                        take_SvsNS <- c(sample(x = s_idx, size = num_s_vs_ns), sample(x = ns_idx, size = num_s_vs_ns))
                        Y <- PCs[take_SvsNS,]
                        ad_data_t <- ad_data[take_SvsNS, ]
                        ad_res <- adonis(formula = Y ~ Survival, data = ad_data_t, permutations = 10000, parallel = 1, method = "euclidean")
                        list(ad_res = ad_res)
                      }, 
                      s_idx = s_idx, ns_idx = ns_idx, PCs = PCs, ad_data = ad_data,
                      mc.cores = 7
)
toc()
ad_r_p <- sapply(lapply(lapply(lapply(ad_res_mc, `[[`, "ad_res"), `[[`, "aov.tab"), `[[`, "Pr(>F)"), `[`, 1)

##Metabolites in metabolite groups, all patients
met_group_res <- list()
tic()
for (met_group in unique(human_data_legend$group[metab_sel])){
  met_group_idx <- which(human_data_legend$group == met_group)
  human_data_dist <- cbind(human_data[,1:5], as.matrix(dist(x = human_data[, met_group_idx], method = "canberra"))) # distance() works on a per row basis
  pat_list <- human_data[, 1:5]
  p <- prcomp(human_data_dist[, -1:-5])
  
  ###PERMANOVA with bootstrapping
  ad_data <- human_data_dist
  ad_data$Survival[ad_data$`CAP / FP` == "-"] <- "Control"
  PCs <- p$x[, 1:2]
  c_idx <- which(ad_data$Survival == "Control")
  s_idx <- which(ad_data$Survival == "S")
  ns_idx <- which(ad_data$Survival == "NS")
  ad_res_mc <- mclapply(1:100, 
                        function(dx, c_idx, s_idx, ns_idx, PCs, ad_data){
                          num_c <- length(c_idx)
                          num_s <- length(s_idx)
                          num_ns <- length(ns_idx)
                          num_c_vs_s <- min(num_c, num_s)
                          num_c_vs_ns <- min(num_c, num_ns)
                          num_s_vs_ns <- min(num_s, num_ns)
                          take_CvsS <- c(sample(x = c_idx, size = num_c_vs_s), sample(x = s_idx, size = num_c_vs_s))
                          take_CvsNS <- c(sample(x = c_idx, size = num_c_vs_ns), sample(x = ns_idx, size = num_c_vs_ns))
                          take_SvsNS <- c(sample(x = s_idx, size = num_s_vs_ns), sample(x = ns_idx, size = num_s_vs_ns))
                          Y1 <- PCs[take_CvsS,]
                          Y2 <- PCs[take_CvsNS,]
                          Y3 <- PCs[take_SvsNS,]
                          ad_data1 <- ad_data[take_CvsS, ]
                          ad_data2 <- ad_data[take_CvsNS, ]
                          ad_data3 <- ad_data[take_SvsNS, ]
                          ad_C_vs_S <- adonis(formula = Y1 ~ Survival, data = ad_data1, permutations = 10000, parallel = 1, method = "euclidean")
                          ad_C_vs_NS <- adonis(formula = Y2 ~ Survival, data = ad_data2, permutations = 10000, parallel = 1, method = "euclidean")
                          ad_S_vs_NS <- adonis(formula = Y3 ~ Survival, data = ad_data3, permutations = 10000, parallel = 1, method = "euclidean")
                          list(ad_C_vs_S = ad_C_vs_S, ad_C_vs_NS = ad_C_vs_NS, ad_S_vs_NS = ad_S_vs_NS)
                        }, 
                        c_idx = c_idx, s_idx = s_idx, ns_idx = ns_idx, PCs = PCs, ad_data = ad_data,
                        mc.cores = 7
  )

  ###Compare step lengths of between patient groups, compare PCA-transformed points
  pat_step_len <- list()
  pat_x_len <- list()
  pat_y_len <- list()
  for (pat in pat_list$Patient){
    ind <- which(human_data$Patient == pat)
    x <- p$x[ind, 1]
    y <- p$x[ind, 2]
    pat_step_len[[pat]] <- sapply(diff(x) ^ 2 + diff(y) ^ 2, sqrt)
    pat_x_len[[pat]] <- abs(diff(x))
    pat_y_len[[pat]] <- abs(diff(y))
  }
  names(pat_step_len) <- as.character(1:length(pat_step_len))
  names(pat_x_len) <- names(pat_step_len)
  names(pat_y_len) <- names(pat_step_len)
  pat_step_len_df <- data.frame(Patient = rep(names(pat_step_len), times = sapply(pat_step_len, length)), step_len = unlist(pat_step_len), x_len = unlist(pat_x_len), y_len = unlist(pat_y_len))
  pat_step_len_df$Group <- pat_list$Survival[match(pat_step_len_df$Patient, pat_list$Patient)]
  pat_step_len_df$Group <- factor(pat_step_len_df$Group, levels = unique(pat_step_len_df$Group)[c(1, 3, 2)])
  pat_step_len_long_df <- melt(pat_step_len_df, id.vars = c("Patient", "Group"))
  pat_step_len_long_df$variable <- factor(pat_step_len_long_df$variable, labels = c("Euclidean step length", "X only", "Y only"))
  
  num_NS_vals <- sum(pat_step_len_df$Group == "NS")
  num_S_vals <- sum(pat_step_len_df$Group == "S")
  bt_step_res <- mclapply(1:100,
                          function(dx, pat_step_len_df, num_S_vals, num_NS_vals){
                            bootstr_dat <- rbind(subset(pat_step_len_df, Group == "NS"), subset(pat_step_len_df, Group == "S")[sample(num_S_vals, size = num_NS_vals), ])
                            bootstr_dat$Group <- factor(bootstr_dat$Group)
                            ad_step1 <- adonis(formula = step_len ~ Group, data = bootstr_dat, method = "euclidean", permutations = 10000, parallel = 1)
                            ad_step2 <- adonis(formula = x_len ~ Group, data = bootstr_dat, method = "euclidean", permutations = 10000, parallel = 1)
                            ad_step3 <- adonis(formula = y_len ~ Group, data = bootstr_dat, method = "euclidean", permutations = 10000, parallel = 1)
                            t_r1 <- t.test(formula = step_len ~ Group, data = bootstr_dat, var.equal = FALSE)
                            t_r2 <- t.test(formula = x_len ~ Group, data = bootstr_dat, var.equal = FALSE)
                            t_r3 <- t.test(formula = y_len ~ Group, data = bootstr_dat, var.equal = FALSE)
                            list(ad_step1 = ad_step1, ad_step2 = ad_step2, ad_step3 = ad_step3, t_r1 = t_r1, t_r2 = t_r2, t_r3 = t_r3)
                          },
                          pat_step_len_df = pat_step_len_df, num_S_vals = num_S_vals, num_NS_vals = num_NS_vals,
                          mc.cores = 7)
  met_group_res[[met_group]] <- list(ad = ad_res_mc, 
                                     bt = bt_step_res)
}
toc()

##Metabolites, all patients
human_data_dist <- cbind(human_data[,1:5], as.matrix(dist(x = human_data[, metab_sel], method = "canberra"))) # distance() works on a per row basis
p <- prcomp(human_data_dist[, -1:-5])

###PERMANOVA with bootstrapping, only contrast 1-vs-1
ad_data <- human_data_dist
pat_list <- human_data_dist[, 1:5]
ad_data$Survival[ad_data$`CAP / FP` == "-"] <- "Control"
PCs <- p$x[, 1:2]
c_idx <- which(ad_data$Survival == "Control")
s_idx <- which(ad_data$Survival == "S")
ns_idx <- which(ad_data$Survival == "NS")
ad_res1 <- list()
ad_res2 <- list()
ad_res3 <- list()
tic()
ad_res_mc <- mclapply(1:100, 
                      function(dx, c_idx, s_idx, ns_idx, PCs, ad_data){
                        num_c <- length(c_idx)
                        num_s <- length(s_idx)
                        num_ns <- length(ns_idx)
                        num_c_vs_s <- min(num_c, num_s)
                        num_c_vs_ns <- min(num_c, num_ns)
                        num_s_vs_ns <- min(num_s, num_ns)
                        take_CvsS <- c(sample(x = c_idx, size = num_c_vs_s), sample(x = s_idx, size = num_c_vs_s))
                        take_CvsNS <- c(sample(x = c_idx, size = num_c_vs_ns), sample(x = ns_idx, size = num_c_vs_ns))
                        take_SvsNS <- c(sample(x = s_idx, size = num_s_vs_ns), sample(x = ns_idx, size = num_s_vs_ns))
                        Y1 <- PCs[take_CvsS,]
                        Y2 <- PCs[take_CvsNS,]
                        Y3 <- PCs[take_SvsNS,]
                        ad_data1 <- ad_data[take_CvsS, ]
                        ad_data2 <- ad_data[take_CvsNS, ]
                        ad_data3 <- ad_data[take_SvsNS, ]
                        ad_res1 <- adonis(formula = Y1 ~ Survival, data = ad_data1, permutations = 10000, parallel = 1, method = "euclidean")
                        ad_res2 <- adonis(formula = Y2 ~ Survival, data = ad_data2, permutations = 10000, parallel = 1, method = "euclidean")
                        ad_res3 <- adonis(formula = Y3 ~ Survival, data = ad_data3, permutations = 10000, parallel = 1, method = "euclidean")
                        list(ad_res1 = ad_res1, ad_res2 = ad_res2, ad_res3 = ad_res3)
                      }, 
                      c_idx = c_idx, s_idx = s_idx, ns_idx = ns_idx, PCs = PCs, ad_data = ad_data,
                      mc.cores = 7
)
toc()
ad_r1_p <- sapply(lapply(lapply(lapply(ad_res_mc, `[[`, "ad_res1"), `[[`, "aov.tab"), `[[`, "Pr(>F)"), `[`, 1)
ad_r2_p <- sapply(lapply(lapply(lapply(ad_res_mc, `[[`, "ad_res2"), `[[`, "aov.tab"), `[[`, "Pr(>F)"), `[`, 1)
ad_r3_p <- sapply(lapply(lapply(lapply(ad_res_mc, `[[`, "ad_res3"), `[[`, "aov.tab"), `[[`, "Pr(>F)"), `[`, 1)

###Compare step lengths between patient groups
pat_step_len <- list()
pat_x_len <- list()
pat_y_len <- list()
for (pat in pat_list$Patient){
  ind <- which(human_data$Patient == pat)
  x <- p$x[ind, 1]
  y <- p$x[ind, 2]
  pat_step_len[[pat]] <- sapply(diff(x) ^ 2 + diff(y) ^ 2, sqrt)
  pat_x_len[[pat]] <- abs(diff(x))
  pat_y_len[[pat]] <- abs(diff(y))
}
names(pat_step_len) <- as.character(1:length(pat_step_len))
names(pat_x_len) <- names(pat_step_len)
names(pat_y_len) <- names(pat_step_len)
pat_step_len_df <- data.frame(Patient = rep(names(pat_step_len), times = sapply(pat_step_len, length)), step_len = unlist(pat_step_len), x_len = unlist(pat_x_len), y_len = unlist(pat_y_len))
pat_step_len_df$Group <- pat_list$Survival[match(pat_step_len_df$Patient, pat_list$Patient)]
pat_step_len_df$Group <- factor(pat_step_len_df$Group, levels = unique(pat_step_len_df$Group)[c(1, 3, 2)])
pat_step_len_long_df <- melt(pat_step_len_df, id.vars = c("Patient", "Group"))
pat_step_len_long_df$variable <- factor(pat_step_len_long_df$variable, labels = c("Euclidean step length", "X only", "Y only"))

ad_step1 <- list()
ad_step2 <- list()
ad_step3 <- list()
t_r1 <- list()
t_r2 <- list()
t_r3 <- list()
num_NS_vals <- sum(pat_step_len_df$Group == "NS")
num_S_vals <- sum(pat_step_len_df$Group == "S")
tic()
bt_step_res <- mclapply(1:100,
                        function(dx, pat_step_len_df, num_S_vals, num_NS_vals){
                          bootstr_dat <- rbind(subset(pat_step_len_df, Group == "NS"), subset(pat_step_len_df, Group == "S")[sample(num_S_vals, size = num_NS_vals), ])
                          bootstr_dat$Group <- factor(bootstr_dat$Group)
                          ad_step1 <- adonis(formula = step_len ~ Group, data = bootstr_dat, method = "euclidean", permutations = 10000, parallel = 1)
                          ad_step2 <- adonis(formula = x_len ~ Group, data = bootstr_dat, method = "euclidean", permutations = 10000, parallel = 1)
                          ad_step3 <- adonis(formula = y_len ~ Group, data = bootstr_dat, method = "euclidean", permutations = 10000, parallel = 1)
                          t_r1 <- t.test(formula = step_len ~ Group, data = bootstr_dat, var.equal = FALSE)
                          t_r2 <- t.test(formula = x_len ~ Group, data = bootstr_dat, var.equal = FALSE)
                          t_r3 <- t.test(formula = y_len ~ Group, data = bootstr_dat, var.equal = FALSE)
                          list(ad_step1 = ad_step1, ad_step2 = ad_step2, ad_step3 = ad_step3, t_r1 = t_r1, t_r2 = t_r2, t_r3 = t_r3)
                        },
                        pat_step_len_df = pat_step_len_df, num_S_vals = num_S_vals, num_NS_vals = num_NS_vals,
                        mc.cores = 7)
toc()
ad_s1_p <- sapply(lapply(lapply(lapply(bt_step_res, `[[`, "ad_step1"), `[[`, "aov.tab"), `[[`, "Pr(>F)"), `[`, 1)
ad_s2_p <- sapply(lapply(lapply(lapply(bt_step_res, `[[`, "ad_step2"), `[[`, "aov.tab"), `[[`, "Pr(>F)"), `[`, 1)
ad_s3_p <- sapply(lapply(lapply(lapply(bt_step_res, `[[`, "ad_step3"), `[[`, "aov.tab"), `[[`, "Pr(>F)"), `[`, 1)
t_r1_p <- sapply(lapply(bt_step_res, `[[`, "t_r1"), `[[`, "p.value")
t_r2_p <- sapply(lapply(bt_step_res, `[[`, "t_r2"), `[[`, "p.value")
t_r3_p <- sapply(lapply(bt_step_res, `[[`, "t_r3"), `[[`, "p.value")

###Compare beta diversity
####Get centroids per patient
pat_centroids <- lapply(unique(pat_list$Patient), function(pat) apply(as.matrix(p$x[which(pat_list$Patient == pat), 1:2]), 2, mean))
pat_centroids <- Reduce("rbind", pat_centroids)
pat_list_centroids <- pat_list[match(unique(pat_list$Patient), pat_list$Patient), ]
####Get pairwise centroid distances
pat_pw_dist <- as.matrix(dist(pat_centroids, method = "euclidean"))
####Get within-group distances
p_s_sel <- which(pat_list_centroids$Survival == "S" & pat_list_centroids$`CAP / FP` != "-")
pat_S_pw_dist <- pat_pw_dist[p_s_sel, p_s_sel]
pat_S_pw_dist <- pat_S_pw_dist[lower.tri(x = pat_S_pw_dist)]
p_ns_sel <- which(pat_list_centroids$Survival == "NS" & pat_list_centroids$`CAP / FP` != "-")
pat_NS_pw_dist <- pat_pw_dist[p_ns_sel, p_ns_sel]
pat_NS_pw_dist <- pat_NS_pw_dist[lower.tri(x = pat_NS_pw_dist)]
p_c_sel <- which(pat_list_centroids$`CAP / FP` == "-")
pat_C_pw_dist <- pat_pw_dist[p_c_sel, p_c_sel]
pat_C_pw_dist <- pat_C_pw_dist[lower.tri(x = pat_C_pw_dist)]
####Compare distance distributions
betadiv_S_NS <- t.test(x = pat_S_pw_dist, y = pat_NS_pw_dist, var.equal = FALSE)
betadiv_S_C <- t.test(x = pat_S_pw_dist, y = pat_C_pw_dist, var.equal = FALSE)
betadiv_NS_C <- t.test(x = pat_NS_pw_dist, y = pat_C_pw_dist, var.equal = FALSE)
####Plot
pat_pw_group_dat <- data.frame(distance = c(pat_C_pw_dist, pat_S_pw_dist, pat_NS_pw_dist), 
                               Group = c(rep("C", length(pat_C_pw_dist)), rep("S", length(pat_S_pw_dist)), rep("NS", length(pat_NS_pw_dist))))
p <- ggplot(data = pat_pw_group_dat, mapping = aes(x = Group, y = distance, color = Group)) +
  geom_boxplot() +
  geom_path(data = data.frame(x = c(1, 1, 1.95, 1.95), y = c(380, 390, 390, 380)), mapping = aes(x = x, y = y), inherit.aes = FALSE) +
  geom_path(data = data.frame(x = c(2.05, 2.05, 3, 3), y = c(380, 390, 390, 380)), mapping = aes(x = x, y = y), inherit.aes = FALSE) +
  geom_path(data = data.frame(x = c(1, 1, 3, 3), y = c(440, 450, 450, 440)), mapping = aes(x = x, y = y), inherit.aes = FALSE) +
  geom_text(x = 1.5, y = 410, label = paste0("p < ", format(betadiv_S_NS$p.value, digits = 4)), size = 2, inherit.aes = FALSE) +
  geom_text(x = 2.5, y = 410, label = paste0("p < ", format(betadiv_NS_C$p.value, digits = 2)), size = 2, inherit.aes = FALSE) +
  geom_text(x = 2, y = 470, label = "p > 0.05", size = 2, inherit.aes = FALSE) +
  human_col_scale(levels = c("NS", "C", "S", "Dummy")) +
  scale_x_discrete(limits = c("S", "NS", "C")) +
  guides(color = "none") +
  ylim(0, 480) +
  theme_bw() + 
  theme(panel.grid = element_blank())
ggsave(filename = "PCA_metab_betadiv_comparison.png", path = out_dir, width = 2.5, height = 4, units = "in")

###Plot step lengths, metabolites, all patients
pat_step_len_brkt <- data.frame(variable = "X only", x = c(1, 1, 2, 2), y = c(480, 490, 490, 480), stringsAsFactors = FALSE)
pat_step_len_sig <- data.frame(variable = "X only", x = 1.5, y = 520, p = paste0("p = ", format(mean(ad_s2_p), digits = 2)), stringsAsFactors = FALSE)
p <- ggplot(data = pat_step_len_long_df, mapping = aes(x = as.numeric(Group), y = value, fill = Group, color = Group, group = Group)) +
  facet_wrap(~ variable) +
  geom_path(data = pat_step_len_brkt, mapping = aes(x = x, y = y), color = "black", inherit.aes = FALSE) +
  geom_text(data = pat_step_len_sig, mapping = aes(x = x, y = y, label = p), size = 3, inherit.aes = FALSE) +
  geom_violin() + 
  xlab("Group") +
  ylab("Step length") +
  scale_x_discrete(limits = 1:3, labels = c("NS", "S", "Nonseptic")) +
  scale_y_continuous(limits = c(0,550)) +
  guides(fill = "none", color = "none", size = "none") +
  theme_bw() +
  theme(panel.grid = element_line(colour = NA))
ggsave(filename = "PCA_metab_steplength_comparison.png", path = out_dir, plot = p, device = "png", width = 8, height = 3, units = "in")

###Aggregate all p-values
ad_mg_r1 <- lapply(lapply(lapply(lapply(lapply(met_group_res, `[[`, "ad"), lapply, `[[`, "ad_C_vs_S"), lapply, `[[`, "aov.tab"), lapply, `[[`, "Pr(>F)"), lapply, `[`, 1)
ad_mg_r2 <- lapply(lapply(lapply(lapply(lapply(met_group_res, `[[`, "ad"), lapply, `[[`, "ad_C_vs_NS"), lapply, `[[`, "aov.tab"), lapply, `[[`, "Pr(>F)"), lapply, `[`, 1)
ad_mg_r3 <- lapply(lapply(lapply(lapply(lapply(met_group_res, `[[`, "ad"), lapply, `[[`, "ad_S_vs_NS"), lapply, `[[`, "aov.tab"), lapply, `[[`, "Pr(>F)"), lapply, `[`, 1)

ad_mg_r1_m <- unlist(lapply(lapply(ad_mg_r1, unlist), mean))
ad_mg_r2_m <- unlist(lapply(lapply(ad_mg_r2, unlist), mean))
ad_mg_r3_m <- unlist(lapply(lapply(ad_mg_r3, unlist), mean))

ad_mg_r_df <- as.data.frame(t(data.frame(ad_mg_r1_m, ad_mg_r2_m, ad_mg_r3_m)))
ad_mg_r_df$pheno <- c(NA, NA, mean(ad_r_p))
ad_mg_r_df$all <- c(mean(ad_r1_p), mean(ad_r2_p), mean(ad_r3_p))
ad_mg_r_df <- as.data.frame(t(apply(ad_mg_r_df, 1, p.adjust)))
rownames(ad_mg_r_df) <- c("C vs S", "C vs NS", "S vs NS")
fwrite(x = ad_mg_r_df, file = paste0(out_dir, "human_met_group_PCA_diff_FDR.csv"), row.names = TRUE, sep = "\t")

ad_st_r1 <- lapply(lapply(lapply(lapply(lapply(met_group_res, `[[`, "bt"), lapply, `[[`, "ad_step1"), lapply, `[[`, "aov.tab"), lapply, `[[`, "Pr(>F)"), lapply, `[`, 1)
ad_st_r2 <- lapply(lapply(lapply(lapply(lapply(met_group_res, `[[`, "bt"), lapply, `[[`, "ad_step2"), lapply, `[[`, "aov.tab"), lapply, `[[`, "Pr(>F)"), lapply, `[`, 1)
ad_st_r3 <- lapply(lapply(lapply(lapply(lapply(met_group_res, `[[`, "bt"), lapply, `[[`, "ad_step3"), lapply, `[[`, "aov.tab"), lapply, `[[`, "Pr(>F)"), lapply, `[`, 1)

ad_st_xy_m <- p.adjust(lapply(lapply(ad_st_r1, unlist), mean))
ad_st_x_m <- p.adjust(lapply(lapply(ad_st_r2, unlist), mean))
ad_st_y_m <- p.adjust(lapply(lapply(ad_st_r3, unlist), mean))

ad_st_r_df <- as.data.frame(t(data.frame(ad_st_xy_m, ad_st_x_m, ad_st_y_m)))
rownames(ad_st_r_df) <- c("XY", "X", "Y")
fwrite(x = ad_st_r_df, file = paste0(out_dir, "human_met_group_PCA_step_diff_FDR_varNS_geq_varS.csv"), row.names = TRUE, sep = "\t")

##Actual plots for each metabolite group
tic()
for (met_group in unique(human_data_legend$group[metab_sel])){
  met_group_idx <- which(human_data_legend$group == met_group)
  human_data_dist <- cbind(human_data[,1:5], as.matrix(dist(x = human_data[, met_group_idx], method = "canberra"))) # distance() works on a per row basis
  p <- prcomp(human_data_dist[, -1:-5])
  
  ad_rv <- ad_mg_r_df[[met_group]]
  ad_rv[ad_rv > 0.05] <- 0.05
  text_v <- paste0(c("C", "C", "S"), " vs ", c("S", "NS", "NS"), ", p ", c("> ", "< ")[1 + (ad_rv < 0.05)], sapply(ad_rv, format, digits = 2))
  
  hdd <- human_data_dist
  hdd$Survival[hdd$`CAP / FP` == "-"] <- "Non-septic"
  ap <- autoplot(object = p, data = hdd, colour = "Survival", frame = TRUE, frame.type = "norm")
  ap <- ap +
    human_col_scale(name = "Survival", levels = c("NS", "Non-septic", "S", "Dummy"), aesthetics = c("colour", "fill")) +
    guides(colour = guide_legend(title = "Group"), fill = "none", group = "none") +
    ggtitle(paste0("PCA biplot, Canberra distance, ", met_group, "\nall samples")) +
    theme_bw()
  gobj <- ggplot_build(ap)
  xmin <- gobj$layout$panel_params[[1]]$x.range[1]
  ymin <- gobj$layout$panel_params[[1]]$y.range[1]
  ap <- ap +
    geom_text(x = 0.9 * xmin, y = 0.85 * ymin, label = paste0(text_v, collapse = "\n"), size = 2, hjust = "left")
  ggsave(filename = paste0("PCA_biplot_", met_group, "_all_samples.png"), path = out_dir, plot = ap, width = 6, height = 5, units = "in")
}

###Actual plot of sepsis samples, clinical params
human_sepsis_data_pheno_dist <- cbind(subset(human_sepsis_data, Day < 4, 1:5), as.matrix(dist(x = subset(human_sepsis_data, Day < 4, pheno_sel), method = "canberra"))) # distance() works on a per row basis
p <- prcomp(human_sepsis_data_pheno_dist[, -1:-5])
ap <- autoplot(object = p, data = human_sepsis_data_pheno_dist, colour = "Survival", frame = TRUE, frame.type = "norm")
ap <- ap + 
  human_col_scale(name = "Survival") +
  ggtitle("PCA biplot, Canberra distance, clinical params,\nseptic patients") +
  theme_bw()
gobj <- ggplot_build(ap)
xmin <- gobj$layout$panel_params[[1]]$x.range[1]
ymin <- gobj$layout$panel_params[[1]]$y.range[1]
ap <- ap +
  geom_text(x = 0.9 * xmin, y = 0.88 * ymin, label = paste0("S vs NS, p < ", format(ad_mg_r_df$pheno[3], digits = 3)), size = 2, hjust = "left")
ggsave(filename = paste0("PCA_biplot_sepsis_pheno_gg.png"), path = out_dir, plot = ap, width = 5, height = 4, units = "in")

##Actual plot of all samples, metabolites
human_data_dist <- cbind(human_data[,1:5], as.matrix(dist(x = human_data[, metab_sel], method = "canberra"))) # distance() works on a per row basis
p <- prcomp(human_data_dist[, -1:-5])
ad_rv <- ad_mg_r_df[["all"]]
ad_rv[ad_rv > 0.05] <- 0.05
text_v <- paste0(c("C", "C", "S"), " vs ", c("S", "NS", "NS"), ", p ", c("> ", "< ")[1 + (ad_rv < 0.05)], sapply(ad_rv, format, digits = 2))
hdd <- human_data_dist
hdd$Survival[hdd$`CAP / FP` == "-"] <- "Non-septic"
ap <- autoplot(object = p, data = hdd, colour = "Survival", frame = TRUE, frame.type = "norm")
ap <- ap + 
  human_col_scale(name = "Survival", levels = c("NS", "Non-septic", "S", "Dummy"), aesthetics = c("colour", "fill")) +
  guides(colour = guide_legend(title = "Group"), fill = "none", group = "none") +
  ggtitle("PCA biplot, Canberra distance, metabolites,\nall samples") +
  theme_bw()
gobj <- ggplot_build(ap)
xmin <- gobj$layout$panel_params[[1]]$x.range[1]
ymin <- gobj$layout$panel_params[[1]]$y.range[1]
ap <- ap +
  geom_text(x = 0.9 * xmin, y = 0.88 * ymin, label = paste0(text_v, collapse = "\n"), size = 2, hjust = "left")
ggsave(filename = paste0("PCA_biplot_metab_all_samples_gg.png"), path = out_dir, plot = ap, width = 5, height = 4, units = "in")

#Patient-wise PCA plot
pat_list <- human_sepsis_data_normal[match(unique(human_sepsis_data_normal$Patient), human_sepsis_data_normal$Patient), 1:5]
pat_len <- table(human_sepsis_data_normal$Patient)
pca_list <- lapply(pat_list$Patient, function(p){ prcomp(subset(human_sepsis_data_normal, Patient == p, -1:-5)) })
for (n in seq_along(pca_list)){
  pca_list[[n]]$surv <- pat_list$Survival[n]
  pca_list[[n]]$pat <- pat_list$Patient[n]
  pca_list[[n]]$capfp <- pat_list$`CAP / FP`[n]
  pca_list[[n]]$days <- unlist(subset(human_sepsis_data_normal, Patient == pat_list$Patient[n], "Day"))
}
pca_list[unlist(lapply(pca_list, function(e) ncol(e$x) <= 2))] <- NULL
png(filename = paste0(out_dir, "pca_sepsis_pats_expl_var.png"), width = 36, height = 16, units = "cm", res = 300)
par(mfrow = c(3,6))
lapply(pca_list, function(e) barplot(summary(e)$importance[2,], ylab = "Explained variance", main = paste(e$pat, e$surv, e$capfp, sep = ", ")))
dev.off()
png(filename = paste0(out_dir, "pca_sepsis_pats_PC1.png"), width = 36, height =16, units = "cm", res = 300)
par(mfrow = c(3,6))
lapply(pca_list, function(e){ plot(e$days, e$x[,1], main = paste(e$pat, e$surv, e$capfp, sep = ", "), type = "l", xlab = "Day", ylab = "PC1") })
dev.off()
png(filename = paste0(out_dir, "pca_sepsis_pats_PC2.png"), width = 36, height = 16, units = "cm", res = 300)
par(mfrow = c(3,6))
lapply(pca_list, function(e){ plot(e$days, e$x[,2], main = paste(e$pat, e$surv, e$capfp, sep = ", "), type = "l", xlab = "Day", ylab = "PC2") })
dev.off()
png(filename = paste0(out_dir, "pca_sepsis_pats_rot1.png"), width = 36, height = 16, units = "cm", res = 300)
par(mfrow = c(3,6))
lapply(pca_list, function(e){ plot(e$rotation[,1], type = "l", main = paste(e$pat, e$surv, e$capfp, sep = ", "), xlab = "Met. index", ylab = "Loading") })
dev.off()
png(filename = paste0(out_dir, "pca_sepsis_pats_PC1_vs_PC2.png"), width = 40, height = 20, units = "cm", res = 300)
par(mfrow = c(3,6))
lapply(pca_list, function(e){ plot(e$x[,1:min(2,ncol(e$x))], type = "p", col = 1:nrow(e$x), main = paste("Patient ", e$pat, ", ", e$surv, ", ", e$capfp, sep = "")); arrows(x0 = e$x[1:(nrow(e$x)-1),1], y0 = e$x[1:(nrow(e$x)-1),2], x1 = e$x[2:nrow(e$x),1], y1 = e$x[2:nrow(e$x),2], length = 0.1) })
dev.off()

long_pca_dat <- lapply(lapply(lapply(pca_list, `[[`, "x"), `[`, , 1:2), as.data.frame)
long_pca_dat <- lapply(seq_along(long_pca_dat), function(n){ long_pca_dat[[n]]$Patient <- pca_list[[n]]$pat; long_pca_dat[[n]]$Survival <- pca_list[[n]]$surv; long_pca_dat[[n]]})
long_pca_dat <- melt(long_pca_dat, id.vars = 1:4)

ggplot(data = long_pca_dat, mapping = aes(x = PC1, y = PC2, group = Patient, color = Survival)) +
  geom_path(linejoin = "mitre", arrow = arrow(length = unit(0.1, "inches"))) +
  theme_bw()

s <- subset(human_sepsis_data_normal, Survival == "S")
ns <- subset(human_sepsis_data_normal, Survival == "NS")
pca_s <- prcomp(s[, metab_sel])
pca_ns <- prcomp(ns[, metab_sel])

png(filename = paste0(out_dir, "pca_sepsis_s_pats_PC1_vs_PC2.png"), width = 16, height = 15, units = "cm", res = 300)
plot(pca_s$x[, 1:2], col = s$Day + 1, main = "Surviving Sepsis patients in common PC space")
for (pat in unique(s$Patient)){
  idx <- which(s$Patient == pat)
  lidx <- length(idx)
  arrows(x0 = pca_s$x[idx[1:(lidx - 1)], 1], y0 = pca_s$x[idx[1:(lidx-1)], 2], x1 = pca_s$x[idx[2:lidx], 1], y1 = pca_s$x[idx[2:lidx], 2], length = 0.1)
}
dev.off()

png(filename = paste0(out_dir, "pca_sepsis_ns_pats_PC1_vs_PC2.png"), width = 16, height = 15, units = "cm", res = 300)
plot(pca_ns$x[, 1:2], col = ns$Day + 1, main = "Nonsurviving Sepsis patients in common PC space\noften start off with a hook")
for (pat in unique(ns$Patient)){
  idx <- which(ns$Patient == pat)
  lidx <- length(idx)
  arrows(x0 = pca_ns$x[idx[1:(lidx - 1)], 1], y0 = pca_ns$x[idx[1:(lidx-1)], 2], x1 = pca_ns$x[idx[2:lidx], 1], y1 = pca_ns$x[idx[2:lidx], 2], length = 0.1)
}
dev.off()

sp <- spca(human_sepsis_data_normal[, metab_sel], 2, sparse = "varnum", para = c(20, 7))
sp$x <- as.matrix(human_sepsis_data_normal[, metab_sel]) %*% sp$loadings

p <- prcomp(human_sepsis_data_normal[, metab_sel])
#barplot(summary(p)$importance[3, ])
plot(p$x[,1:2], col = as.factor(human_sepsis_data_normal$Day))
col <- as.factor(pat_list$Survival)
for (n in seq_along(pat_list$Patient)){
  ind <- which(human_sepsis_data_normal$Patient == pat_list$Patient[n])
  x <- p$x[ind ,1]
  y <- p$x[ind ,2]
  arrows(x0 = x[-length(x)], y0 = y[-length(y)], x1 = x[-1], y1 = y[-1], length = 0.1, col = col[n])
}

rot_score <- lapply(pca_list, function(e){ order(abs(e$rotation[,1]), decreasing = FALSE) })
rot_rem <- sapply(lapply(rot_score, `<=`, 10), which)
rot_rem <- unique(as.numeric(rot_rem))

p <- prcomp(human_sepsis_data_normal[, metab_sel][, rot_rem])
#barplot(summary(p)$importance[3, ])
plot(p$x[,1:2], col = as.factor(human_sepsis_data_normal$Day))
col <- as.factor(pat_list$Survival)
for (n in seq_along(pat_list$Patient)){
  ind <- which(human_sepsis_data_normal$Patient == pat_list$Patient[n])
  x <- p$x[ind ,1]
  y <- p$x[ind ,2]
  arrows(x0 = x[-length(x)], y0 = y[-length(y)], x1 = x[-1], y1 = y[-1], length = 0.1, col = col[n])
}

rot_pc1 <- t(sapply(pca_list, function(e) e$rotation[,1]))
rot_pca_pc1 <- prcomp(rot_pc1)
plot(rot_pca_pc1)
plot(rot_pca_pc1$rotation[,1])
rot_pc2 <- t(sapply(pca_list, function(e) e$rotation[,2]))
rot_pca_pc2 <- prcomp(rot_pc2)
plot(rot_pca_pc2)
plot(rot_pca_pc2$rotation[,1])
pca_con1 <- rot_pca_pc1$rotation[,1] %*% t(as.matrix(human_sepsis_data_normal[, -1:-5]))
pca_con2 <- rot_pca_pc2$rotation[,1] %*% t(as.matrix(human_sepsis_data_normal[, -1:-5]))

stab <- table(human_sepsis_data_normal$Survival)
human_sepsis_data_normal[human_sepsis_data_normal$Survival == "NS", -1:-5] <- sqrt(stab["NS"]) * human_sepsis_data_normal[human_sepsis_data_normal$Survival == "NS", -1:-5]
human_sepsis_data_normal[human_sepsis_data_normal$Survival == "S", -1:-5] <- sqrt(stab["S"]) * human_sepsis_data_normal[human_sepsis_data_normal$Survival == "S", -1:-5]
for (pat in unique(human_sepsis_data_normal$Patient)){
  sel <- human_sepsis_data_normal$Patient == pat
  human_sepsis_data_normal[sel, -1:-5] <- sqrt(sum(sel)) * human_sepsis_data_normal[sel, -1:-5]
}

#p <- prcomp(kernelMatrix(kernel = rbfdot(sigma = 0.0001), x = as.matrix(human_sepsis_data_normal[,-1:-5])))
rangeKernel <- function(x, range = 1){
  as.matrix(dist(x)) < range
}
p <- prcomp(rangeKernel(human_sepsis_data_normal[, 6:204], 20))
barplot(p$sdev/sum(p$sdev))
plot(p$x[,1:2], col = as.factor(human_sepsis_data_normal$Day))
col <- as.factor(pat_list$Survival)
for (n in seq_along(pat_list$Patient)){
  ind <- which(human_sepsis_data_normal$Patient == pat_list$Patient[n])
  x <- p$x[ind ,1]
  y <- p$x[ind ,2]
  arrows(x0 = x[-length(x)], y0 = y[-length(y)], x1 = x[-1], y1 = y[-1], length = 0.1, col = col[n])
}

m2d <- mapper2D(distance_matrix = dist(human_sepsis_data_normal[, 6:204]), filter_values = human_sepsis_data_normal[, 6:8], num_intervals = c(4, 4), percent_overlap = 50, num_bins_when_clustering = 10)
m2d_graph <- graph.adjacency(m2d$adjacency, mode = "undirected")
plot(m2d_graph, layout = layout.circle(m2d_graph))

g <- graph.adjacency(rangeKernel(subset(human_sepsis_data_normal, TRUE, 6:204), range = 16), mode = "undirected", diag = F)
gc <- cluster_fast_greedy(g)
plot(gc, g)

p <- prcomp(human_sepsis_data_normal[, which(rowMeans(rot_score) > 140) + 5])
par(mfrow = c(1,1))
plot(summary(p)$importance[2,])
plot(p$x[,1:2], col = as.factor(human_sepsis_data_normal$Day))
col <- as.factor(pat_list$Survival)
for (n in seq_along(pat_list$Patient)){
  ind <- which(human_sepsis_data_normal$Patient == pat_list$Patient[n])
  x <- p$x[ind ,1]
  y <- p$x[ind ,2]
  arrows(x0 = x[-length(x)], y0 = y[-length(y)], x1 = x[-1], y1 = y[-1], length = 0.1, col = col[n])
}
#legend(x = 0.3, y = -0.25, legend = c("survivors", "nonsurvivors"), pch = "-", col = 2:1)
legend(x = 13, y = 8, legend = c("survivors", "nonsurvivors"), pch = "-", col = 2:1)

rot_pc12 <- rbind(rot_pc1, rot_pc2)
rot_pc12_sim <- as.matrix(dist(rot_pc12))
heatmaply(rot_pc12_sim, dendrogram = T, margins = c(50,50,0,50))

plot(rot_pca$sdev/sum(rot_pca$sdev))

linmod <- lm(rot_pc1 ~ 1)

par(mfrow=c(1,1))
#tic()
#er <- estim.regul(X = human_sepsis_data_metab, Y = human_sepsis_data_pheno, grid1 = seq(0.01, 1, length.out = 10), grid2 = seq(0.25, 1.5, length.out = 8))
#toc()

#rc <- rcc(X = human_sepsis_data_metab, Y = human_sepsis_data_pheno, lambda1 = er$lambda1, lambda2 = er$lambda2)

dfmax_w <- c(80, 30, 20)
ypredict <- function(x, yc, cc) {
  en <- glmnet::glmnet(x, yc, alpha = 0.5, intercept = FALSE,
                       dfmax = dfmax_w[cc], lower.limits = 0)
  W <- coef(en)
  return(W[2:nrow(W), ncol(W)])
}
dfmax_v <- c(21, 15, 15)
xpredict <- function(y, xc, cc) {
  en <- glmnet::glmnet(y, xc, alpha = 0.5, intercept = FALSE,
                       dfmax = dfmax_v[cc], lower.limits = 0)
  V <- coef(en)
  return(V[2:nrow(V), ncol(V)])
}

nsrc <- nscancor(x = human_sepsis_data_metab, y = human_sepsis_data_pheno, xpredict = xpredict, ypredict = ypredict, nvar = 3, xscale = T, yscale = T)

rc_sel <- which(rc$cor^2 > 0.5)
hist(rc$ycoef[, rc_sel])
plt.cc(rc)

#rc$scores$corr.Y.xscores

corr_dat <- list()
for (d in 0:2){
  NS_corr <- my.corr.test(x = subset(human_sepsis_data_normal, Survival == "NS" & Day == d, select = -1:-5), adjust = "fdr")
  S_corr <- my.corr.test(x = subset(human_sepsis_data_normal, Survival == "S" & Day == d, select = -1:-5), adjust = "fdr")
  corr_dat[[as.character(d)]] <- list(NS = NS_corr, S = S_corr)
}
hist((corr_dat$`0`$S$r - corr_dat$`0`$NS$r)[upper.tri(corr_dat$`0`$S$r)])
hist((corr_dat$`1`$S$r - corr_dat$`1`$NS$r)[upper.tri(corr_dat$`0`$S$r)])
hist((corr_dat$`2`$S$r - corr_dat$`2`$NS$r)[upper.tri(corr_dat$`0`$S$r)])

corr_diff <- sapply(corr_dat, function(e){ as.numeric((e$S$r - e$NS$r)[upper.tri(e$S$r)]) })
corr_S <- sapply(corr_dat, function(e){ as.numeric(e$S$r[upper.tri(e$S$r)]) })
corr_NS <- sapply(corr_dat, function(e){ as.numeric(e$NS$r[upper.tri(e$NS$r)]) })
corr_diff_mat <- corr_diff[!rowAnyNAs(corr_diff),]
corr_S_mat <- corr_S[!rowAnyNAs(corr_S) & !rowAnyNAs(corr_NS),]
corr_NS_mat <- corr_NS[!rowAnyNAs(corr_NS) & !rowAnyNAs(corr_S),]

pcrr_diff <- prcomp(corr_diff_mat)
plot(pcrr_diff$sdev/sum(pcrr_diff$sdev))
matplot(pcrr_diff$rotation, type = "l")

pcrr_S <- prcomp(corr_S_mat)
plot(pcrr_S$sdev/sum(pcrr_S$sdev))
matplot(pcrr_S$rotation, type = "l")

pcrr_NS <- prcomp(corr_NS_mat)
plot(pcrr_NS$sdev/sum(pcrr_NS$sdev))
matplot(pcrr_NS$rotation, type = "l")

corr_sel <- rowAnys(abs(corr_S_mat) > 0.75) | rowAnys(abs(corr_NS_mat) > 0.75)
matplot(t(corr_diff_mat[corr_sel,]), type = "l")

matplot(x = t(corr_S_mat), y = t(corr_NS_mat), type = "l")
corr_hist2D <- list()
for (n in 0:2){
  dat <- corr_dat[[paste0(n)]]
  dat_sel <- upper.tri(dat$S$r)
  corr_hist2D[[paste(n)]] <- hexbin(x = dat$S$r[dat_sel], y = dat$NS$r[dat_sel], xlab = "NS correlation", ylab = "S correlation", IDs = TRUE)
}
plot(corr_hist2D[[1]], main = "Day 0")
plot(-1:1, -1:1)
day0_dy <- aggregate(corr_S[,2], list(ID = corr_hist2D[[1]]@cID), mean)
day0_dx <- aggregate(corr_NS[,2], list(ID = corr_hist2D[[1]]@cID), mean)
arrows(corr_hist2D[[1]]@xcm, corr_hist2D[[1]]@ycm, day0_dx$x, day0_dy$x, length = 0.05)
plot(corr_hist2D[[2]], main = "Day 1")
plot(-1:1, -1:1)
day1_dy <- aggregate(corr_S[,3], list(ID = corr_hist2D[[2]]@cID), mean)
day1_dx <- aggregate(corr_NS[,3], list(ID = corr_hist2D[[2]]@cID), mean)
arrows(corr_hist2D[[2]]@xcm, corr_hist2D[[2]]@ycm, day1_dx$x, day1_dy$x, length = 0.05)

plot(-1:1, -1:1)
day0_d2_dy <- aggregate(corr_S[,3], list(ID = corr_hist2D[[1]]@cID), mean)
day0_d2_dx <- aggregate(corr_NS[,3], list(ID = corr_hist2D[[1]]@cID), mean)
arrows(corr_hist2D[[1]]@xcm, corr_hist2D[[1]]@ycm, day0_d2_dx$x, day0_d2_dy$x, length = 0.05)

plot(corr_hist2D[[3]], main = "Day 2")

slim <- c(-0.5,0)
nslim <- c(-0.4,0)
for (n in 1:3){
  print(sum(between(corr_S_mat[,n], slim[1], slim[2]) & between(corr_NS_mat[,n], nslim[1], nslim[2])))
}

slim <- c(0.6,1)
nslim <- c(0,0.5)
for (n in 1:3){
  print(sum(between(corr_S_mat[,n], slim[1], slim[2]) & between(corr_NS_mat[,n], nslim[1], nslim[2])))
}


plot(-1:1, -1:1)
arrows(corr_hist2D[[1]]@xcm, corr_hist2D[[1]]@ycm, corr_hist2D[[3]]@xcm, corr_hist2D[[3]]@ycm, length = 0.05)
