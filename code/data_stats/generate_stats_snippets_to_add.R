library(hexbin)
library(matrixStats)
library(data.table)
library(missRanger)
library(nscancor)
library(TDAmapper)
library(elasticnet)
library(corpcor)
library(vegan)
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(parallel)
library(stringi)
library(tictoc)
library(cowplot)

source("../function_definitions.R")

out_dir <- "../../results/data_stats_NS_high_var/"
if (!dir.exists(out_dir))
  dir.create(out_dir)

if (stri_detect(str = R.version$os, fixed = "linux")){
  num_cores <- min(detectCores() - 4, 100)
}else{
  num_cores <- 1
}

#Read data
human_data <- get_human_sepsis_data()

human_data_legend <- get_human_sepsis_legend()

#Get variable types
metab_end <- which(colnames(human_data) == "H1")
metab_sel <- 7:metab_end
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
var_sel <- colnames(human_data)[-1:-6]

human_sepsis_data <- subset(human_sepsis_data, select = c(1:6, which(colnames(human_sepsis_data) %in% var_sel)))
human_data <- subset(human_data, select = c(1:6, which(colnames(human_data) %in% var_sel)))
human_data_legend <- human_data_legend[c(1:6, which(human_data_legend[, 1] %in% var_sel)), ]

metab_end <- which(colnames(human_data) == "H1")
metab_sel <- 7:metab_end
pheno_sel <- (metab_end + 1):ncol(human_data)

#Normalize data
human_sepsis_data_normal <- human_sepsis_data
human_sepsis_data_normal[,-1:-6] <- scale(human_sepsis_data_normal[,-1:-6])
human_sepsis_data_metab <- subset(human_sepsis_data_normal, Day < 4, metab_sel)
human_sepsis_data_pheno <- subset(human_sepsis_data_normal, Day < 4, pheno_sel)
human_sepsis_data_normal <- human_sepsis_data_normal[,1:metab_end]

#Generate PCA biplots
##Clinical params, sepsis patients
human_sepsis_data_pheno_dist <- cbind(subset(human_sepsis_data, Day < 4, 1:6), as.matrix(dist(x = subset(human_sepsis_data, Day < 4, pheno_sel), method = "canberra"))) # distance() works on a per row basis
p <- prcomp(human_sepsis_data_pheno_dist[, -1:-6])

###PERMANOVA with bootstrapping
ad_data <- human_sepsis_data_pheno_dist
pat_list <- human_sepsis_data_pheno_dist[, 1:6]
PCs <- p$x[, 1:2]
s_idx <- which(ad_data$Group == "Septic-S")
ns_idx <- which(ad_data$Group == "Septic-NS")
tic()
ad_res_mc <- mclapply(1:100, 
                      function(dx, s_idx, ns_idx, PCs, ad_data){
                        set.seed(dx)
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
                      mc.cores = num_cores)
toc()
ad_r_p <- sapply(lapply(lapply(lapply(ad_res_mc, `[[`, "ad_res"), `[[`, "aov.tab"), `[[`, "Pr(>F)"), `[`, 1)

##MDS plot with arrows - trial to see how it looks like
hsd_pheno_mds <- metaMDS(comm = subset(human_sepsis_data, Day < 4, pheno_sel), distance = "canberra", k = 2, wascores = TRUE, autotransform = FALSE, try = 20, trymax = 40)
pts <- as.data.frame(hsd_pheno_mds$points)
pts$Group <- subset(human_sepsis_data, Day < 4)$Group
ars <- hsd_pheno_mds$species
ars <- ars[apply(ars, 1, function(row) sqrt(sum(row^2))) > 0.4, ]
ars_text <- rownames(ars)
ars <- as.data.frame(ars)
ars$Group <- 1:nrow(ars)
ars <- rbind(data.frame(MDS1 = rep(0, nrow(ars)), MDS2 = rep(0, nrow(ars)), Group = ars$Group), ars)
ggplot(data = pts, mapping = aes(x = MDS1, y = MDS2, color = Group)) + 
  geom_point() + 
  geom_path(data = ars, mapping = aes(x = MDS1, y = MDS2, group = Group), inherit.aes = FALSE, arrow = arrow()) +
  #geom_text(data = ars[(length(ars_names) + 1):nrow(ars), ], label = ars_names, x = ars[ars_names, 1], y = ars[ars_names, 2], inherit.aes = FALSE) +
  human_col_scale() +
  theme_bw()

##Metabolites in metabolite groups, all patients
if (file.exists("met_group_balanced_permanova_res.RData")){
  met_group_res <- readRDS(file = "met_group_balanced_permanova_res.RData")
}else{
  met_group_res <- list()
  tic()
  for (met_group in unique(human_data_legend$group[metab_sel])){
    met_group_idx <- which(human_data_legend$group == met_group)
    human_data_dist <- cbind(human_data[,1:6], as.matrix(dist(x = human_data[, met_group_idx], method = "canberra"))) # distance() works on a per row basis
    pat_list <- human_data[, 1:6]
    p <- prcomp(human_data_dist[, -1:-6])
    
    ###PERMANOVA with bootstrapping
    ad_data <- human_data_dist
    PCs <- p$x[, 1:2]
    c_idx <- which(ad_data$Group %in% c("non-Septic-NS", "non-Septic-S"))
    cns_idx <- which(ad_data$Group == "non-Septic-NS")
    cs_idx <- which(ad_data$Group == "non-Septic-S")
    s_idx <- which(ad_data$Group == "Septic-S")
    ns_idx <- which(ad_data$Group == "Septic-NS")
    ad_res_mc <- mclapply(1:100, 
                          function(dx, c_idx, s_idx, ns_idx, cns_idx, PCs, ad_data){
                            set.seed(dx + 500)
                            num_c <- length(c_idx)
                            num_s <- length(s_idx)
                            num_ns <- length(ns_idx)
                            num_cns <- length(cns_idx)
                            num_c_vs_s <- min(num_c, num_s)
                            num_c_vs_ns <- min(num_c, num_ns)
                            num_s_vs_ns <- min(num_s, num_ns)
                            num_ns_vs_cns <- min(num_ns, num_cns)
                            take_CvsS <- c(sample(x = c_idx, size = num_c_vs_s), sample(x = s_idx, size = num_c_vs_s))
                            take_CvsNS <- c(sample(x = c_idx, size = num_c_vs_ns), sample(x = ns_idx, size = num_c_vs_ns))
                            take_SvsNS <- c(sample(x = s_idx, size = num_s_vs_ns), sample(x = ns_idx, size = num_s_vs_ns))
                            take_NSvsCNS <- c(sample(x = ns_idx, size = num_ns_vs_cns), sample(x = cns_idx, size = num_ns_vs_cns))
                            Y1 <- PCs[take_CvsS,]
                            Y2 <- PCs[take_CvsNS,]
                            Y3 <- PCs[take_SvsNS,]
                            Y4 <- PCs[take_NSvsCNS,]
                            ad_data1 <- ad_data[take_CvsS, ]
                            ad_data2 <- ad_data[take_CvsNS, ]
                            ad_data3 <- ad_data[take_SvsNS, ]
                            ad_data4 <- ad_data[take_NSvsCNS, ]
                            ad_C_vs_S <- adonis(formula = Y1 ~ Group, data = ad_data1, permutations = 10000, parallel = 1, method = "euclidean")
                            ad_C_vs_NS <- adonis(formula = Y2 ~ Group, data = ad_data2, permutations = 10000, parallel = 1, method = "euclidean")
                            ad_S_vs_NS <- adonis(formula = Y3 ~ Group, data = ad_data3, permutations = 10000, parallel = 1, method = "euclidean")
                            ad_NS_vs_CNS <- adonis(formula = Y4 ~ Group, data = ad_data4, permutations = 10000, parallel = 1, method = "euclidean")
                            list(ad_C_vs_S = ad_C_vs_S, ad_C_vs_NS = ad_C_vs_NS, ad_S_vs_NS = ad_S_vs_NS, ad_NS_vs_CNS = ad_NS_vs_CNS)
                          }, 
                          c_idx = c_idx, s_idx = s_idx, ns_idx = ns_idx, cns_idx = cns_idx, PCs = PCs, ad_data = ad_data,
                          mc.cores = num_cores)
  
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
    pat_step_len_df$Group <- pat_list$Group[match(pat_step_len_df$Patient, pat_list$Patient)]
    pat_step_len_df$Group <- factor(pat_step_len_df$Group, levels = unique(pat_step_len_df$Group))
    pat_step_len_long_df <- melt(pat_step_len_df, id.vars = c("Patient", "Group"))
    pat_step_len_long_df$variable <- factor(pat_step_len_long_df$variable, labels = c("Euclidean step length", "X only", "Y only"))
    
    num_NS_vals <- sum(pat_step_len_df$Group == "Septic-NS")
    num_S_vals <- sum(pat_step_len_df$Group == "Septic-S")
    bt_step_res <- mclapply(1:100,
                            function(dx, pat_step_len_df, num_S_vals, num_NS_vals){
                              set.seed(dx + 1000)
                              bootstr_dat <- rbind(subset(pat_step_len_df, Group == "Septic-NS"), subset(pat_step_len_df, Group == "Septic-S")[sample(num_S_vals, size = num_NS_vals), ])
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
                            mc.cores = num_cores)
    met_group_res[[met_group]] <- list(ad = ad_res_mc, 
                                       bt = bt_step_res)
  }
  toc()
  saveRDS(object = met_group_res, file = "met_group_balanced_permanova_res.RData")
}

##Metabolites, all patients
human_data_dist <- cbind(human_data[,1:6], as.matrix(dist(x = human_data[, metab_sel], method = "canberra"))) # distance() works on a per row basis
p <- prcomp(human_data_dist[, -1:-6])

###PERMANOVA with bootstrapping, only contrast 1-vs-1
ad_data <- human_data_dist
pat_list <- human_data_dist[, 1:6]
PCs <- p$x[, 1:2]
c_idx <- which(ad_data$Group %in% c("non-Septic-NS", "non-Septic-S"))
s_idx <- which(ad_data$Group == "Septic-S")
ns_idx <- which(ad_data$Group == "Septic-NS")
cns_idx <- which(ad_data$Group == "non-Septic-NS")
ad_res1 <- list()
ad_res2 <- list()
ad_res3 <- list()
ad_res4 <- list()
if (file.exists("all_metabs_balanced_permanova_res.RData")){
  ad_res_mc <- readRDS(file = "all_metabs_balanced_permanova_res.RData")
}else{
  tic()
  ad_res_mc <- mclapply(1:100, 
                        function(dx, c_idx, s_idx, ns_idx, cns_idx, PCs, ad_data){
                          set.seed(dx + 1500)
                          num_c <- length(c_idx)
                          num_s <- length(s_idx)
                          num_ns <- length(ns_idx)
                          num_cns <- length(cns_idx)
                          num_c_vs_s <- min(num_c, num_s)
                          num_c_vs_ns <- min(num_c, num_ns)
                          num_s_vs_ns <- min(num_s, num_ns)
                          num_ns_vs_cns <- min(num_s, num_cns)
                          take_CvsS <- c(sample(x = c_idx, size = num_c_vs_s), sample(x = s_idx, size = num_c_vs_s))
                          take_CvsNS <- c(sample(x = c_idx, size = num_c_vs_ns), sample(x = ns_idx, size = num_c_vs_ns))
                          take_SvsNS <- c(sample(x = s_idx, size = num_s_vs_ns), sample(x = ns_idx, size = num_s_vs_ns))
                          take_NSvsCNS <- c(sample(x = ns_idx, size = num_ns_vs_cns), sample(x = cns_idx, size = num_ns_vs_cns))
                          Y1 <- PCs[take_CvsS,]
                          Y2 <- PCs[take_CvsNS,]
                          Y3 <- PCs[take_SvsNS,]
                          Y4 <- PCs[take_NSvsCNS,]
                          ad_data1 <- ad_data[take_CvsS, ]
                          ad_data2 <- ad_data[take_CvsNS, ]
                          ad_data3 <- ad_data[take_SvsNS, ]
                          ad_data4 <- ad_data[take_NSvsCNS, ]
                          ad_res1 <- adonis(formula = Y1 ~ Group, data = ad_data1, permutations = 10000, parallel = 1, method = "euclidean")
                          ad_res2 <- adonis(formula = Y2 ~ Group, data = ad_data2, permutations = 10000, parallel = 1, method = "euclidean")
                          ad_res3 <- adonis(formula = Y3 ~ Group, data = ad_data3, permutations = 10000, parallel = 1, method = "euclidean")
                          ad_res4 <- adonis(formula = Y4 ~ Group, data = ad_data4, permutations = 10000, parallel = 1, method = "euclidean")
                          list(ad_res1 = ad_res1, ad_res2 = ad_res2, ad_res3 = ad_res3, ad_res4 = ad_res4)
                        }, 
                        c_idx = c_idx, s_idx = s_idx, ns_idx = ns_idx, cns_idx = cns_idx, PCs = PCs, ad_data = ad_data,
                        mc.cores = num_cores)
  toc()
  saveRDS(object = ad_res_mc, file = "all_metabs_balanced_permanova_res.RData")
}
ad_r1_p <- sapply(lapply(lapply(lapply(ad_res_mc, `[[`, "ad_res1"), `[[`, "aov.tab"), `[[`, "Pr(>F)"), `[`, 1)
ad_r2_p <- sapply(lapply(lapply(lapply(ad_res_mc, `[[`, "ad_res2"), `[[`, "aov.tab"), `[[`, "Pr(>F)"), `[`, 1)
ad_r3_p <- sapply(lapply(lapply(lapply(ad_res_mc, `[[`, "ad_res3"), `[[`, "aov.tab"), `[[`, "Pr(>F)"), `[`, 1)
ad_r4_p <- sapply(lapply(lapply(lapply(ad_res_mc, `[[`, "ad_res4"), `[[`, "aov.tab"), `[[`, "Pr(>F)"), `[`, 1)

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
pat_step_len_df$Group <- pat_list$Group[match(pat_step_len_df$Patient, pat_list$Patient)]
pat_step_len_df$Group <- factor(pat_step_len_df$Group, levels = unique(pat_step_len_df$Group))
pat_step_len_long_df <- melt(pat_step_len_df, id.vars = c("Patient", "Group"))
pat_step_len_long_df$variable <- factor(pat_step_len_long_df$variable, labels = c("Euclidean step length", "X only", "Y only"))

ad_step1 <- list()
ad_step2 <- list()
ad_step3 <- list()
t_r1 <- list()
t_r2 <- list()
t_r3 <- list()
num_NS_vals <- sum(pat_step_len_df$Group == "Septic-NS")
num_S_vals <- sum(pat_step_len_df$Group == "Septic-S")
tic()
bt_step_res <- mclapply(1:100,
                        function(dx, pat_step_len_df, num_S_vals, num_NS_vals){
                          set.seed(dx + 2000)
                          bootstr_dat <- rbind(subset(pat_step_len_df, Group == "Septic-NS"), subset(pat_step_len_df, Group == "Septic-S")[sample(num_S_vals, size = num_NS_vals), ])
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
                        mc.cores = num_cores)
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
p_s_sel <- which(pat_list_centroids$Group == "Septic-S")
pat_S_pw_dist <- pat_pw_dist[p_s_sel, p_s_sel]
pat_S_pw_dist <- pat_S_pw_dist[lower.tri(x = pat_S_pw_dist)]
p_ns_sel <- which(pat_list_centroids$Group == "Septic-NS")
pat_NS_pw_dist <- pat_pw_dist[p_ns_sel, p_ns_sel]
pat_NS_pw_dist <- pat_NS_pw_dist[lower.tri(x = pat_NS_pw_dist)]
p_c_sel <- which(grepl(x = pat_list_centroids$Group, pattern = "non-Septic"))
pat_C_pw_dist <- pat_pw_dist[p_c_sel, p_c_sel]
pat_C_pw_dist <- pat_C_pw_dist[lower.tri(x = pat_C_pw_dist)]
p_cs_sel <- which(pat_list_centroids$Group == "non-Septic-S")
pat_CS_pw_dist <- pat_pw_dist[p_cs_sel, p_cs_sel]
pat_CS_pw_dist <- pat_CS_pw_dist[lower.tri(x = pat_CS_pw_dist)]
p_cns_sel <- which(pat_list_centroids$Group == "non-Septic-NS")
pat_CNS_pw_dist <- pat_pw_dist[p_cns_sel, p_cns_sel]
pat_CNS_pw_dist <- pat_CNS_pw_dist[lower.tri(x = pat_CNS_pw_dist)]
####Compare distance distributions
betadiv_S_NS <- t.test(x = pat_S_pw_dist, y = pat_NS_pw_dist, var.equal = FALSE) #var.equal set explicitly FALSE to make it visible to you, the reader
betadiv_S_CS <- t.test(x = pat_S_pw_dist, y = pat_CS_pw_dist, var.equal = FALSE)
betadiv_S_CNS <- t.test(x = pat_S_pw_dist, y = pat_CNS_pw_dist, var.equal = FALSE)
betadiv_NS_CS <- t.test(x = pat_NS_pw_dist, y = pat_CS_pw_dist, var.equal = FALSE)
betadiv_NS_CNS <- t.test(x = pat_NS_pw_dist, y = pat_CNS_pw_dist, var.equal = FALSE)
betadiv_CS_CNS <- t.test(x = pat_CNS_pw_dist, y = pat_CS_pw_dist, var.equal = FALSE)
betadiv_S_C <- t.test(x = pat_S_pw_dist, y = pat_C_pw_dist, var.equal = FALSE)
betadiv_NS_C <- t.test(x = pat_NS_pw_dist, y = pat_C_pw_dist, var.equal = FALSE)
env <- environment()
env_names <- names(env)
bdiv_names <- env_names[grep(pattern = "betadiv", x = env_names)]
bdiv_elements <- lapply(bdiv_names, function(e) eval(as.symbol(e)))
bdiv_names <- bdiv_names[sapply(bdiv_elements, class) == "htest"]
bdiv_elements <- bdiv_elements[sapply(bdiv_elements, class) == "htest"]
bdiv_p <- sapply(bdiv_elements, function(e) e$p.value)
bdiv_fdr <- p.adjust(bdiv_p, method = "fdr")
names(bdiv_p) <- bdiv_names
names(bdiv_fdr) <- bdiv_names
g1 <- sub("betadiv", "", bdiv_names)
g1 <- sub("_CS", "_non-Septic-S", g1)
g1 <- sub("_CNS", "_non-Septic-NS", g1)
g1 <- sub("_C", "_non-Septic", g1)
g1 <- sub("_S", "_Septic-S", g1)
g1 <- sub("_NS", "_Septic-NS", g1)
g1 <- substring(g1, 2)
g1 <- strsplit(g1, "_", fixed = TRUE)
bdiv_df <- data.frame(Group1 = sapply(g1, `[`, 1), Group2 = sapply(g1, `[`, 2), p = bdiv_p, FDR = bdiv_fdr)
bdiv_df <- bdiv_df[order(rownames(bdiv_df)), ]
write.csv(x = bdiv_df, file = paste0(out_dir, "betadiversity_comparison_pvals.csv"), row.names = FALSE)
####Plot
pat_pw_group_dat <- data.frame(distance = c(pat_CS_pw_dist, pat_CNS_pw_dist, pat_C_pw_dist, pat_S_pw_dist, pat_NS_pw_dist), 
                               Group = c(rep("non-Septic-S", length(pat_CS_pw_dist)), 
                                         rep("non-Septic-NS", length(pat_CNS_pw_dist)), 
                                         rep("non-Septic", length(pat_C_pw_dist)),
                                         rep("Septic-S", length(pat_S_pw_dist)), 
                                         rep("Septic-NS", length(pat_NS_pw_dist))))
lsize <- 0.4
tsize <- 2.5
p_betadiv <- ggplot(data = pat_pw_group_dat, mapping = aes(x = Group, y = distance, color = Group)) +
  geom_boxplot() +
  #geom_path(data = data.frame(x = c(1, 1, 2, 2), y = c(380, 390, 390, 380)), mapping = aes(x = x, y = y), size = lsize, inherit.aes = FALSE) +
  #geom_path(data = data.frame(x = 2 + c(1, 1, 1.95, 1.95), y = c(380, 390, 390, 380)), mapping = aes(x = x, y = y), size = lsize, inherit.aes = FALSE) +
  geom_path(data = data.frame(x = 2 + c(2.05, 2.05, 3, 3), y = c(380, 390, 390, 380)), mapping = aes(x = x, y = y), size = lsize, inherit.aes = FALSE) +
  geom_path(data = data.frame(x = 2 + c(1, 1, 3, 3), y = c(440, 450, 450, 440)), mapping = aes(x = x, y = y), size = lsize, inherit.aes = FALSE) +
  geom_path(data = data.frame(x = c(1, 1, 4, 4), y = c(500, 510, 510, 500)), mapping = aes(x = x, y = y), size = lsize, inherit.aes = FALSE) +
  geom_path(data = data.frame(x = c(1, 1, 5, 5), y = c(560, 570, 570, 560)), mapping = aes(x = x, y = y), size = lsize, inherit.aes = FALSE) +
  #geom_path(data = data.frame(x = c(2, 2, 5, 5), y = c(560, 570, 570, 560)), mapping = aes(x = x, y = y), size = lsize, inherit.aes = FALSE) +
  #geom_text(x = 1.5, y = 410, label = "q > 0.05", size = tsize, inherit.aes = FALSE) + #CNS vs CS
  #geom_text(x = 2 + 1.5, y = 410, label = "q > 0.05", size = tsize, inherit.aes = FALSE) + #C vs S
  geom_text(x = 2 + 2.5, y = 410, label = paste0("q < ", format(bdiv_df$FDR[8], digits = 1)), size = tsize, inherit.aes = FALSE) + #S vs NS
  geom_text(x = 2 + 2, y = 470, label = paste0("q < ", format(bdiv_df$FDR[2], digits = 2)), size = tsize, inherit.aes = FALSE) + #C vs NS
  geom_text(x = 2.5, y = 530, label = paste0("q < ", format(bdiv_df$FDR[7], digits = 2)), size = tsize, inherit.aes = FALSE) + #CS vs S
  geom_text(x = 3, y = 590, label = paste0("q < ", format(bdiv_df$FDR[4], digits = 2)), size = tsize, inherit.aes = FALSE) + #CS vs NS
  #geom_text(x = 3.5, y = 590, label = "q > 0.05", size = tsize, inherit.aes = FALSE) + #CNS vs NS
  human_col_scale(levels = as.character(unique(pat_pw_group_dat$Group))[c(5, 1, 4, 2, 3)], black_color = "grey40", black_pos = 5) +
  scale_x_discrete(limits = as.character(unique(pat_pw_group_dat$Group))) +
  #geom_segment(x = 2.5, y = 0, xend = 2.5, yend = max(pat_pw_group_dat$distance), size = 0.25, inherit.aes = FALSE) +
  guides(color = "none") +
  ylim(0, 600) +
  ylab("between-patient dissimilarity") +
  xlab("") +
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.background = element_blank())
ggsave(filename = "PCA_metab_betadiv_comparison.png", path = out_dir, plot = p_betadiv, width = 2.5, height = 4, units = "in")

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
  scale_x_discrete(limits = 1:4, labels = unique(pat_step_len_long_df$Group)) +
  scale_y_continuous(limits = c(0,550)) +
  guides(fill = "none", color = "none", size = "none") +
  theme_bw() +
  theme(panel.grid = element_line(colour = NA), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(filename = "PCA_metab_steplength_comparison.png", path = out_dir, plot = p, device = "png", width = 8, height = 3, units = "in")

###Aggregate all p-values
ad_mg_r1 <- lapply(lapply(lapply(lapply(lapply(met_group_res, `[[`, "ad"), lapply, `[[`, "ad_C_vs_S"), lapply, `[[`, "aov.tab"), lapply, `[[`, "Pr(>F)"), lapply, `[`, 1)
ad_mg_r2 <- lapply(lapply(lapply(lapply(lapply(met_group_res, `[[`, "ad"), lapply, `[[`, "ad_C_vs_NS"), lapply, `[[`, "aov.tab"), lapply, `[[`, "Pr(>F)"), lapply, `[`, 1)
ad_mg_r3 <- lapply(lapply(lapply(lapply(lapply(met_group_res, `[[`, "ad"), lapply, `[[`, "ad_S_vs_NS"), lapply, `[[`, "aov.tab"), lapply, `[[`, "Pr(>F)"), lapply, `[`, 1)
ad_mg_r4 <- lapply(lapply(lapply(lapply(lapply(met_group_res, `[[`, "ad"), lapply, `[[`, "ad_NS_vs_CNS"), lapply, `[[`, "aov.tab"), lapply, `[[`, "Pr(>F)"), lapply, `[`, 1)

ad_mg_r1_m <- unlist(lapply(lapply(ad_mg_r1, unlist), mean))
ad_mg_r2_m <- unlist(lapply(lapply(ad_mg_r2, unlist), mean))
ad_mg_r3_m <- unlist(lapply(lapply(ad_mg_r3, unlist), mean))
ad_mg_r4_m <- unlist(lapply(lapply(ad_mg_r4, unlist), mean))

ad_mg_r_df <- as.data.frame(t(data.frame(ad_mg_r1_m, ad_mg_r2_m, ad_mg_r3_m, ad_mg_r4_m)))
ad_mg_r_df$pheno <- c(NA, NA, mean(ad_r_p), NA)
ad_mg_r_df$all <- c(mean(ad_r1_p), mean(ad_r2_p), mean(ad_r3_p), mean(ad_r4_p))
ad_mg_r_df <- as.data.frame(t(apply(ad_mg_r_df, 1, p.adjust)))
rownames(ad_mg_r_df) <- c("C vs S", "C vs NS", "S vs NS", "NS vs CNS")
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
ap_metab_group <- list()
rsize <- rel(4)
for (met_group in setdiff(unique(human_data_legend$group[metab_sel]), "sugar")){
  met_group_idx <- which(human_data_legend$group == met_group)
  human_data_dist <- cbind(human_data[,1:6], as.matrix(dist(x = human_data[, met_group_idx], method = "canberra"))) # distance() works on a per row basis
  p <- prcomp(human_data_dist[, -1:-6])
  
  ad_rv <- ad_mg_r_df[[met_group]]
  ad_rv[ad_rv > 0.05] <- 0.05
  text_v <- paste0(c("non-Septic", "non-Septic", "Septic-S", "non-Septic-NS"), " vs ", c("Septic-S", "Septic-NS", "Septic-NS", "Septic-NS"), ", q ", c("> ", "< ")[1 + (ad_rv < 0.05)], sapply(ad_rv, format, digits = 2))

  lam <- p$sdev[1:2] * sqrt(nrow(p$x))
  hdd <- human_data_dist
  ap_metab_group[[met_group]] <- autoplot(object = p, data = hdd, colour = "Group", frame = TRUE, frame.type = "norm", size = 0.9, frame.alpha = 0.1)
  ap_metab_group[[met_group]] <- ap_metab_group[[met_group]] +
    human_col_scale(aesthetics = c("colour", "fill")) +
    guides(colour = "none", fill = "none", group = "none") +
    ggtitle(met_group)
  gobj <- ggplot_build(ap_metab_group[[met_group]])
  xmax <- gobj$layout$panel_params[[1]]$x.range[2]
  xmin <- gobj$layout$panel_params[[1]]$x.range[1]
  ymax <- gobj$layout$panel_params[[1]]$y.range[2]
  ymin <- gobj$layout$panel_params[[1]]$y.range[1]
  arws <- make_ordination_arrows(x = p$x[, 1:2], w = subset(human_data, select = met_group_idx), xmax = xmax, xmin = xmin, ymin = ymin, ymax = ymax, shorten_arw_by = 4, scale_by = 1.1)
  ap_metab_group[[met_group]] <- ap_metab_group[[met_group]] +
    geom_path(data = arws$ph_ars_s, mapping = aes(x = PC1, y = PC2, group = Group), inherit.aes = FALSE, arrow = arrow(length = unit(0.06, "inches")), size = 0.3, color = "grey50") + 
    geom_text(data = arws$ph_ars_names_s, mapping = aes(label = label, x = PC1, y = PC2), inherit.aes = FALSE, size = rsize, color = "grey20") +
    scale_x_continuous(limits = c(xmin, xmax) * 1.1) + 
    scale_y_continuous(limits = c(ymin, ymax * 1.4)) +
    geom_text(x = 1.15 * xmin, y = 1.2 * ymax, label = paste0(text_v, collapse = "\n"), size = 3, hjust = "left")
  #ggsave(filename = paste0("PCA_biplot_", met_group, "_all_samples.png"), path = out_dir, plot = ap, width = 6, height = 5, units = "in")
}
for (met_group in seq_along(ap_metab_group)[-length(ap_metab_group)]){
  ap_metab_group[[met_group]] <- ap_metab_group[[met_group]] + 
    theme_bw() +
    theme(panel.grid = element_blank())
}
ap_metab_group[[length(ap_metab_group)]] <- ap_metab_group[[length(ap_metab_group)]] + 
  guides(colour = guide_legend(title = "Group"), fill = "none", group = "none") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.direction = "horizontal", legend.position = "bottom", legend.title = element_blank())
ap_metab_group_supp_panel <- plot_grid(plotlist = ap_metab_group, ncol = 2, nrow = 3, labels = "AUTO", align = "hv", axis = "tblr")
ggsave(plot = ap_metab_group_supp_panel, filename = "PCA_biplot_met_group_panel.png", path = out_dir, width = 9, height = 4.5 * 3, units = "in")
ggsave(plot = ap_metab_group_supp_panel, filename = "PCA_biplot_met_group_panel.svg", path = out_dir, width = 9, height = 4.5 * 3, units = "in")
toc()

###Actual plot of sepsis samples, clinical params
human_sepsis_data_pheno_dist <- cbind(subset(human_sepsis_data, Day < 4, 1:6), as.matrix(dist(x = subset(human_sepsis_data, Day < 4, pheno_sel), method = "canberra"))) # distance() works on a per row basis
human_sepsis_data_pheno_dist$Patient <- factor(human_sepsis_data_pheno_dist$Patient, levels = unique(human_sepsis_data_pheno_dist$Patient))
p <- prcomp(human_sepsis_data_pheno_dist[, -1:-6])
hsdpd <- human_sepsis_data_pheno_dist
hsdpd$label <- paste0("P", hsdpd$Patient, "-D", hsdpd$Day)
lam <- p$sdev[1:2] * sqrt(nrow(p$x))
pxy <- as.data.frame(p$x)
pxy[, 1] <- pxy[, 1] / lam[1]
pxy[, 2] <- pxy[, 2] / lam[2]
pxy$Group <- hsdpd$Group
rsize <- rel(3.5)
p_biochem_ap <- autoplot(object = p, data = hsdpd, colour = "Group", shape = FALSE, label = FALSE, frame = TRUE, frame.type = "norm", size = 0.8, frame.alpha = 0.1)
p_biochem_ap <- p_biochem_ap + 
  geom_text_repel(parse = TRUE, label = sapply(hsdpd$label, function(lab) sprintf("bold(\"%s\")", lab)), x = pxy[, 1], y = pxy[, 2], mapping = aes(color = Group), size = 2, segment.size = 0.3, force = 0.005, box.padding = 0.0) + #arguments before publication readyness: size = 1.3, segment.size = 0.3, force = 0.005, box.padding = 0.0
  human_col_scale() +
  #ggtitle("PCA biplot, Canberra distance, biochemical params,\nseptic patients") +
  guides(shape = "none") +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = "none", text = element_text(size = rsize))
gobj <- ggplot_build(p_biochem_ap)
xmax <- gobj$layout$panel_params[[1]]$x.range[2]
xmin <- gobj$layout$panel_params[[1]]$x.range[1]
ymax <- gobj$layout$panel_params[[1]]$y.range[2]
ymin <- gobj$layout$panel_params[[1]]$y.range[1]
arws <- make_ordination_arrows(x = p$x[, 1:2], w = subset(human_sepsis_data, Day < 4, pheno_sel), xmax = xmax, xmin = xmin, ymin = ymin, ymax = ymax, shorten_arw_by = 3)
arws$ph_ars_names_s$label[arws$ph_ars_names_s$label == "Troponin T"] <- "TnT"
p_biochem_ap <- p_biochem_ap +
  geom_path(data = arws$ph_ars_s, mapping = aes(x = PC1, y = PC2, group = Group), inherit.aes = FALSE, arrow = arrow(length = unit(0.06, "inches")), size = 0.3, color = "grey50") + 
  geom_text(data = arws$ph_ars_names_s, mapping = aes(label = label, x = PC1, y = PC2), inherit.aes = FALSE, size = 3, color = "grey20") +
  scale_x_continuous(limits = c(xmin * 1.05, xmax)) +
  scale_y_continuous(limits = c(ymin, ymax * 1)) +
  geom_text(x = 0.95 * xmin, y = 0.98 * ymin, label = paste0("S vs NS, q < ", format(ad_mg_r_df$pheno[3], digits = 3)), size = 2.5, hjust = "left")
ggsave(filename = paste0("PCA_biplot_sepsis_pheno_gg.png"), path = out_dir, plot = p_biochem_ap, width = 4.5, height = 4, units = "in")

##Actual plot of all samples, metabolites
human_data_dist <- cbind(human_data[,1:6], as.matrix(dist(x = human_data[, metab_sel], method = "canberra"))) # distance() works on a per row basis
p <- prcomp(human_data_dist[, -1:-6])
ad_rv <- ad_mg_r_df[["all"]]
ad_rv[ad_rv > 0.05] <- 0.05
text_v <- paste0(c("non-Septic", "non-Septic", "Septic-S", "non-Septic-NS"), " vs ", c("Septic-S", "Septic-NS", "Septic-NS", "Septic-NS"), ", q ", c("> ", "< ")[1 + (ad_rv < 0.05)], sapply(ad_rv, format, digits = 2))
lam <- p$sdev[1:2] * sqrt(nrow(p$x))
pxy <- as.data.frame(p$x)
pxy[, 1] <- pxy[, 1] / lam[1]
pxy[, 2] <- pxy[, 2] / lam[2]
hdd <- human_data_dist
pxy$Group <- hdd$Group
centerxy <- aggregate(x = pxy[, 1:2], by = list(Group = pxy$Group), FUN = mean)
alpha <- -0.55
rotmat <- matrix(c(cos(alpha), sin(alpha), -sin(alpha), cos(alpha)), ncol = 2)
invrotmat <- matrix(c(cos(-alpha), sin(-alpha), -sin(-alpha), cos(-alpha)), ncol = 2)
newxy <- data.frame(as.matrix(pxy[, 1:2]) %*% rotmat)
colnames(newxy) <- c("PC1", "PC2")
#plot(newxy[, 1], newxy[,2])
lmodel <- lm(PC2 ~ poly(PC1, 4), data = newxy)
xs <- data.frame(PC1 = seq(min(newxy$PC1), max(newxy$PC2) * 1.2, length.out = 100))
xs$PC2 <- predict(object = lmodel, newdata = xs)
linexy <- data.frame(as.matrix(xs) %*% invrotmat)
colnames(linexy) <- c("PC1", "PC2")
# newcenterxy <- data.frame(as.matrix(centerxy[c("PC1", "PC2")]) %*% rotmat)
# colnames(newcenterxy) <- c("PC1", "PC2")
# newcenterxy$PC2 <- predict(object = lmodel, newdata = newcenterxy["PC1"])
# newcenterxy[c("PC1", "PC2")] <- data.frame(as.matrix(newcenterxy) %*% invrotmat)
proj_num <- apply(centerxy[, -1], 1, function(coord) which.min(rowSums((linexy - matrix(coord, ncol = length(coord), nrow = nrow(linexy), byrow = TRUE))^2)))
proj_num[length(proj_num)] <- proj_num[length(proj_num)] + 1
newcenterxy <- linexy[proj_num, ]
colnames(newcenterxy) <- c("PC1", "PC2")
newcenterxy$Group <- centerxy$Group
p_metab_ap <- autoplot(object = p, data = hdd, colour = "Group", frame = TRUE, frame.type = "norm", size = 1.3, frame.alpha = 0.1)
p_metab_ap <- p_metab_ap + 
  human_col_scale(aesthetics = c("colour", "fill")) +
  guides(colour = guide_legend(title = "Group"), fill = "none", group = "none") +
  #ggtitle("PCA biplot, Canberra distance, metabolites,\nall samples") +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.direction = "horizontal", legend.box.just = "left", legend.position = "bottom", legend.title = element_blank(), text = element_text(size = rsize), legend.text = element_text(size = rsize), plot.background = element_blank())
gobj <- ggplot_build(p_metab_ap)
xmax <- gobj$layout$panel_params[[1]]$x.range[2]
xmin <- gobj$layout$panel_params[[1]]$x.range[1]
ymax <- gobj$layout$panel_params[[1]]$y.range[2]
ymin <- gobj$layout$panel_params[[1]]$y.range[1]
arws <- make_ordination_arrows(x = p$x[, 1:2], w = subset(human_data, select = metab_sel), xmax = xmax, xmin = xmin, ymin = ymin, ymax = ymax, shorten_arw_by = 350, scale_by = 1.1, num = 4)
p_metab_ap <- p_metab_ap +
  geom_path(data = linexy, mapping = aes(x = PC1, y = PC2), colour = "grey50", inherit.aes = FALSE) + 
  geom_path(data = arws$ph_ars_s, mapping = aes(x = PC1, y = PC2, group = Group), inherit.aes = FALSE, arrow = arrow(length = unit(0.06, "inches")), size = 0.3, color = "grey50") + 
  geom_text(data = arws$ph_ars_names_s, mapping = aes(label = label, x = PC1, y = PC2), inherit.aes = FALSE, size = 3, color = "grey20") +
  geom_text(x = xmin, y = 0.76 * ymax, label = paste0(text_v, collapse = "\n"), size = 2.5, hjust = "left") + 
  geom_point(data = newcenterxy, mapping = aes(x = PC1, y = PC2), colour = "grey50", shape = 18, size = 3, inherit.aes = FALSE) +
  geom_point(data = newcenterxy, mapping = aes(x = PC1, y = PC2, colour = Group), shape = 3, size = 5, inherit.aes = FALSE) +
  scale_x_continuous(limits = c(xmin * 1.1, xmax * 0.95), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(ymin * 1, ymax * 0.95), expand = c(0, 0))
ggsave(filename = "PCA_biplot_metab_all_samples_gg.png", path = out_dir, plot = p_metab_ap, width = 5, height = 4, units = "in")

#Build PCA figure panel
alp <- plot_grid(p_biochem_ap, p_metab_ap, p_betadiv, ncol = 3, align = "h", axis = "bt", rel_widths = c(1, 1, 0.6)) + 
  draw_plot_label(label = c("A", "B", "C"), x = c(0, 1, 2) / 2.6)
ggsave(filename = "PCA_and_betadiv_figure_panel.png", path = out_dir, plot = alp, width = 9, height = 4, units = "in")
ggsave(filename = "PCA_and_betadiv_figure_panel.svg", path = out_dir, plot = alp, width = 9, height = 4, units = "in")

#TODO: continue from here with switch from {S, NS, Control} to {Sepsis x Survival}

#Patient-wise PCA plot
pat_list <- human_sepsis_data_normal[match(unique(human_sepsis_data_normal$Patient), human_sepsis_data_normal$Patient), 1:6]
pat_len <- table(human_sepsis_data_normal$Patient)
pca_list <- lapply(pat_list$Patient, function(p){ prcomp(subset(human_sepsis_data_normal, Patient == p, -1:-6)) })
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
