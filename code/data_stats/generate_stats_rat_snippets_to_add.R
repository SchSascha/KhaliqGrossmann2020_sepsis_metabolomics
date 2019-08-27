library(matrixStats)
library(data.table)
library(missRanger)
library(ggplot2)
library(tictoc)
library(parallel)
library(vegan)
library(ggrepel)
library(ggfortify)
library(scales)
library(matrixStats)

source("../function_definitions.R")

out_dir_cs <- "../../results/data_stats_rat_control_vs_septic/"
out_dir_ns <- "../../results/data_stats_rat_surv_vs_nonsurv/"

out_dir <- "../../results/data_stats_rat/"

if (stri_detect(str = R.version$os, fixed = "linux")){
  num_cores <- min(detectCores() - 1, 100)
}else{
  num_cores <- 1
}

#Make sure paths exist
if (!dir.exists(out_dir))
  dir.create(out_dir)
if (!dir.exists(out_dir_cs))
  dir.create(out_dir_cs)
if (!dir.exists(out_dir_ns))
  dir.create(out_dir_ns)

#Import data
rat_sepsis_data <- get_rat_sepsis_data()

rat_sepsis_data <- rat_sepsis_data[!rat_sepsis_data$`Sample Identification` %in% c("060H", "039L"), ]

load("../../results/data_stats_rat_surv_vs_nonsurv/ANOVA_complete_res.RData")

human_dev_metabs <- read.csv(file = "../../results/data_stats/generalized_safe_corridor_minmax_dev_mets.csv", stringsAsFactors = FALSE)

##Import corresponding group assignment
rat_sepsis_legend <- get_rat_sepsis_legend()
rat_sepsis_legend$group[rat_sepsis_legend$group == ""] <- rat_sepsis_legend[rat_sepsis_legend$group == "", 1]

#Impute missing data
pheno_sel <- (which(colnames(rat_sepsis_data) == "H1")+1):ncol(rat_sepsis_data)
metab_sel <- 5:which(colnames(rat_sepsis_data) == "H1")
cols <- colnames(rat_sepsis_data)
colnames(rat_sepsis_data) <- make.names(cols)
for (mat in setdiff(unique(rat_sepsis_data$material), "plasma")){
  idat <- subset(rat_sepsis_data, material == mat, metab_sel)
  rat_sepsis_data[rat_sepsis_data$material == mat, metab_sel] <- missRanger(idat)
}
idat <- subset(rat_sepsis_data, material == "plasma", c(metab_sel, pheno_sel))
rat_sepsis_data[rat_sepsis_data$material == "plasma", c(metab_sel, pheno_sel)] <- missRanger(idat)
colnames(rat_sepsis_data) <- cols

#Remove dynamic physiological properties like heart rate
rat_sepsis_data <- rat_sepsis_data[, -which(colnames(rat_sepsis_data) %in% c("HR", "SV", "CO", "EF", "Resp Rate", "Temperature"))] #remove bc. variables were used to determine the group in the experiment
pheno_sel <- (which(colnames(rat_sepsis_data) == "H1")+1):ncol(rat_sepsis_data) #redo indexes after removal of columns
metab_sel <- 5:which(colnames(rat_sepsis_data) == "H1")

#Process data
##PCA on biochemical parameters
rat_sepsis_data_pheno_dist <- cbind(subset(rat_sepsis_data, material == "plasma" & `time point` == "24h", select = 1:4), as.matrix(dist(x = subset(rat_sepsis_data, material == "plasma" & `time point` == "24h", select = pheno_sel), method = "canberra"))) # distance() works on a per row basis
p <- prcomp(rat_sepsis_data_pheno_dist[, -1:-4])

###PERMANOVA with bootstrapping
ad_data <- rat_sepsis_data_pheno_dist
pat_list <- rat_sepsis_data_pheno_dist[, 1:4]
PCs <- p$x[, 1:2]
ad_res_mc <- list()
for (mat in unique(rat_sepsis_data$material)){
  ad_res_mc[[mat]] <- mclapply(unique(pat_list$group), 
                        function(dx, PCs, ad_data){
                          take_smpls <- ad_data$group != dx
                          Y <- PCs[take_smpls, ]
                          ad_data_t <- ad_data[take_smpls, ]
                          ad_res <- adonis(formula = Y ~ group, data = ad_data_t, permutations = 10000, parallel = 1, method = "euclidean")
                          list(ad_res = ad_res)
                        }, 
                        PCs = PCs, ad_data = ad_data,
                        mc.cores = num_cores)
}
ad_r_p <- sapply(lapply(lapply(lapply(ad_res_mc, lapply, `[[`, "ad_res"), `[[`, "aov.tab"), `[[`, "Pr(>F)"), `[`, 1)
ad_r_p <- lapply(ad_res_mc, lapply, function(e) e$ad_res$aov.tab$`Pr(>F)`[1])
ad_r_p_tab <- melt(ad_r_p)
ad_r_p_tab$L2 <- c("S vs NS", "C vs NS", "C vs S")[ad_r_p_tab$L2]

##PCA on all metabolites
###PERMANOVA with bootstrapping
ad_res_mc <- list()
for (mat in unique(rat_sepsis_data$material)){
  rat_sepsis_data_dist <- cbind(subset(rat_sepsis_data, material == mat & `time point` == "24h", select = 1:4), as.matrix(dist(x = subset(rat_sepsis_data, material == mat & `time point` == "24h", select = metab_sel), method = "canberra"))) # distance() works on a per row basis
  p <- prcomp(rat_sepsis_data_dist[, -1:-4])
  ad_data <- rat_sepsis_data_dist
  pat_list <- rat_sepsis_data_dist[, 1:4]
  PCs <- p$x[, 1:2]
  ad_res_mc[[mat]] <- mclapply(unique(pat_list$group), 
                               function(dx, PCs, ad_data){
                                 take_smpls <- ad_data$group != dx
                                 Y <- PCs[take_smpls, ]
                                 ad_data_t <- ad_data[take_smpls, ]
                                 ad_res <- adonis(formula = Y ~ group, data = ad_data_t, permutations = 10000, parallel = 1, method = "euclidean")
                                 list(ad_res = ad_res)
                               }, 
                               PCs = PCs, ad_data = ad_data,
                               mc.cores = num_cores)
}
ad_r_m <- sapply(lapply(lapply(lapply(ad_res_mc, lapply, `[[`, "ad_res"), `[[`, "aov.tab"), `[[`, "Pr(>F)"), `[`, 1)
ad_r_m <- lapply(ad_res_mc, lapply, function(e) e$ad_res$aov.tab$`Pr(>F)`[1])
ad_r_m_tab <- melt(ad_r_m)
ad_r_m_tab$L2 <- c("S vs NS", "C vs NS", "C vs S")[ad_r_m_tab$L2]

###PERMANOVA with bootstrapping
mat_res <- list()
for (mat in unique(rat_sepsis_data$material)){
  mat_res[[mat]] <- list()
  for (met_group in unique(rat_sepsis_legend$group[metab_sel - 5])){
    rds <- subset(rat_sepsis_data, material == mat & `time point` == "24h", metab_sel)
    rds <- subset(rds, select = rat_sepsis_legend[rat_sepsis_legend$group == met_group, 1])
    w <- which.xy(is.na(rds))
    rds <- rds[, setdiff(1:ncol(rds), unique(w[, 2]))]
    rat_data_dist <- cbind(subset(rat_sepsis_data, material == mat & `time point` == "24h", 1:4), as.matrix(dist(x = rds, method = "canberra"))) # distance() works on a per row basis
    p <- prcomp(rat_data_dist[, -1:-4])
    ad_data <- rat_data_dist
    pat_list <- rat_data_dist[, 1:4]
    PCs <- p$x[, 1:2]
    mat_res[[mat]][[met_group]] <- mclapply(unique(pat_list$group), 
                                   function(dx, PCs, ad_data){
                                     take_smpls <- ad_data$group != dx
                                     Y <- PCs[take_smpls, ]
                                     ad_data_t <- ad_data[take_smpls, ]
                                     ad_res <- adonis(formula = Y ~ group, data = ad_data_t, permutations = 1000, parallel = 1, method = "euclidean")
                                     list(ad_res = ad_res)
                                   }, 
                                   PCs = PCs, ad_data = ad_data,
                                   mc.cores = num_cores)
  }
}
ad_mg_p <- lapply(mat_res, lapply, lapply, function(e) e$ad_res$aov.tab$`Pr(>F)`[1])
ad_mg_p_tab <- melt(ad_mg_p)
ad_mg_p_tab$L3 <- c("S vs NS", "C vs NS", "C vs S")[ad_mg_p_tab$L3]
ad_mg_p_tab <- dcast(ad_mg_p_tab, L3 + L1 ~ L2)
ad_ps_tab <- cbind(ad_mg_p_tab, data.frame(all_metab = ad_r_m_tab$value), data.frame(pheno = ad_r_p_tab$value))
ad_fdr_tab <- ad_ps_tab
ad_fdr_tab[, -1:-2] <- t(apply(ad_ps_tab[, -1:-2], 1, p.adjust)) #None of those is significant ...

###Actual plot of sepsis samples, clinical params
rsd <- subset(rat_sepsis_data, material == "plasma" & group != "control" & `time point` == "24h")
rat_sepsis_data_pheno_dist <- cbind(subset(rsd, select = 1:4), as.matrix(dist(x = subset(rsd, select = pheno_sel), method = "canberra"))) # distance() works on a per row basis
rat_sepsis_data_pheno_dist$`Sample Identification` <- factor(rat_sepsis_data_pheno_dist$`Sample Identification`, levels = unique(rat_sepsis_data_pheno_dist$`Sample Identification`))
p <- prcomp(rat_sepsis_data_pheno_dist[, -1:-4])
rsdpd <- rat_sepsis_data_pheno_dist
rsdpd$label <- paste0("S", rsdpd$`Sample Identification`, "-", rsdpd$`time point`)
lam <- p$sdev[1:2] * sqrt(nrow(p$x))
ap <- autoplot(object = p, data = rsdpd, colour = "group", frame = FALSE, frame.type = "norm", size = 0.8)
ap <- ap + 
  geom_point(size = 3, color = "white") +
  geom_text_repel(parse = TRUE, label = sapply(rsdpd$label, function(lab) sprintf("bold(\"%s\")", lab)), x = p$x[, 1] / lam[1], y = p$x[, 2] / lam[2], size = 1.3, mapping = aes(color = group), segment.size = 0.3, force = 0.005, box.padding = 0.0) +
  ggtitle("PCA biplot, Canberra distance,\nbiochemical parameters, 24h samples") +
  guides(shape = "none") +
  human_col_scale(name = "Group", levels = c("septic non-survivor", "control", "septic survivor", "", "")) +
  theme_bw() + 
  theme(panel.grid = element_blank())
gobj <- ggplot_build(ap)
xmax <- gobj$layout$panel_params[[1]]$x.range[2]
xmin <- gobj$layout$panel_params[[1]]$x.range[1]
ymax <- gobj$layout$panel_params[[1]]$y.range[2]
ymin <- gobj$layout$panel_params[[1]]$y.range[1]
arws <- make_ordination_arrows(x = p$x[, 1:2], w = subset(rsd, select = -1:-4), xmax = xmax, xmin = xmin, ymin = ymin, ymax = ymax, shorten_arw_by = 2, scale_by = 1.2)
ap <- ap +
  geom_path(data = arws$ph_ars_s, mapping = aes(x = PC1, y = PC2, group = Group), inherit.aes = FALSE, arrow = arrow(length = unit(0.06, "inches")), size = 0.2, color = "grey50") + 
  geom_text(data = arws$ph_ars_names_s, mapping = aes(label = label, x = PC1, y = PC2), inherit.aes = FALSE, size = 1.5, color = "grey50") #+
#   geom_text(x = 0.9 * xmin, y = 0.95 * ymin, label = paste0("S vs NS, q < ", format(ad_mg_r_df$pheno[3], digits = 3)), size = 2, hjust = "left")
ggsave(filename = paste0("PCA_biplot_sepsis_pheno_gg.png"), path = out_dir, plot = ap, width = 5, height = 4, units = "in")

##Actual plot of all samples, metabolites
mtab_list <- list()
bdiv_list <- list()
rat_dev_list <- list()
for (mat in unique(rat_sepsis_data$material)){
  rds <- subset(rat_sepsis_data, material == mat & `time point` == "24h", metab_sel)
  w <- which.xy(is.na(rds))
  rds <- rds[, setdiff(1:ncol(rds), unique(w[, 2]))]
  rat_data_dist <- cbind(subset(rat_sepsis_data, material == mat & `time point` == "24h", 1:4), as.matrix(dist(x = rds, method = "canberra"))) # distance() works on a per row basis
  p <- prcomp(rat_data_dist[, -1:-4])
  #ad_rv <- ad_mg_r_df[["all"]]
  #ad_rv[ad_rv > 0.05] <- 0.05
  #text_v <- paste0(c("Nonsep", "Nonsep", "Septic-S", "Nonsep-NS"), " vs ", c("Septic-S", "Septic-NS", "Septic-NS", "Septic-NS"), ", q ", c("> ", "< ")[1 + (ad_rv < 0.05)], sapply(ad_rv, format, digits = 2))
  lam <- p$sdev[1:2] * sqrt(nrow(p$x))
  ap <- autoplot(object = p, data = rat_data_dist, colour = "group", frame = FALSE, frame.type = "norm", size = 0.5)
  ap <- ap + 
    human_col_scale(name = "Group", levels = c("septic non-survivor", "control", "septic survivor", "", ""), aesthetics = c("color", "fill")) +
    guides(colour = guide_legend(title = "Group"), fill = "none", group = "none") +
    ggtitle("PCA biplot, Canberra distance, metabolites,\n24h samples") +
    theme_bw() + 
    theme(panel.grid = element_blank())
  gobj <- ggplot_build(ap)
  xmax <- gobj$layout$panel_params[[1]]$x.range[2]
  xmin <- gobj$layout$panel_params[[1]]$x.range[1]
  ymax <- gobj$layout$panel_params[[1]]$y.range[2]
  ymin <- gobj$layout$panel_params[[1]]$y.range[1]
  arws <- make_ordination_arrows(x = p$x[, 1:2], w = rds, xmax = xmax, xmin = xmin, ymin = ymin, ymax = ymax, shorten_arw_by = 1, scale_by = 1.1)
  ap <- ap +
    geom_path(data = arws$ph_ars_s, mapping = aes(x = PC1, y = PC2, group = Group), inherit.aes = FALSE, arrow = arrow(length = unit(0.06, "inches")), size = 0.2, color = "grey50") + 
    geom_text(data = arws$ph_ars_names_s, mapping = aes(label = label, x = PC1, y = PC2), inherit.aes = FALSE, size = 1.5, color = "grey50") #+
  #  scale_x_continuous(limits = c(xmin * 1.1, xmax)) +
  #  geom_text(x = 0.9 * xmin, y = 0.85 * ymin, label = paste0(text_v, collapse = "\n"), size = 2, hjust = "left")
  ggsave(filename = paste0("PCA_biplot_metab_all_samples_", mat, "_gg.png"), path = out_dir, plot = ap, width = 5, height = 4, units = "in")
  
  ###Compare beta diversity
  ####Get pairwise sample distances (here centroids == samples)
  rat_pw_dist <- as.matrix(dist(p$x[, 1:2], method = "euclidean"))
  ####Get within-group distances
  p_s_sel <- which(rat_data_dist$group == "septic survivor")
  rat_S_pw_dist <- rat_pw_dist[p_s_sel, p_s_sel]
  rat_S_pw_dist <- rat_S_pw_dist[lower.tri(x = rat_S_pw_dist)]
  p_ns_sel <- which(rat_data_dist$group == "septic non-survivor")
  rat_NS_pw_dist <- rat_pw_dist[p_ns_sel, p_ns_sel]
  rat_NS_pw_dist <- rat_NS_pw_dist[lower.tri(x = rat_NS_pw_dist)]
  p_c_sel <- which(rat_data_dist$group == "control")
  rat_C_pw_dist <- rat_pw_dist[p_c_sel, p_c_sel]
  rat_C_pw_dist <- rat_C_pw_dist[lower.tri(x = rat_C_pw_dist)]
  ####Compare distance distributions
  betadiv_S_NS <- t.test(x = rat_S_pw_dist, y = rat_NS_pw_dist, var.equal = FALSE) #var.equal set explicitly FALSE to make it visible to you, the reader
  betadiv_S_C <- t.test(x = rat_S_pw_dist, y = rat_C_pw_dist, var.equal = FALSE)
  betadiv_NS_C <- t.test(x = rat_NS_pw_dist, y = rat_C_pw_dist, var.equal = FALSE)
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
  g1 <- sub("_C", "_control", g1)
  g1 <- sub("_S", "_septic survivor", g1)
  g1 <- sub("_NS", "_septic non-survivor", g1)
  g1 <- substring(g1, 2)
  g1 <- strsplit(g1, "_", fixed = TRUE)
  bdiv_df <- data.frame(Group1 = sapply(g1, `[`, 1), Group2 = sapply(g1, `[`, 2), p = bdiv_p, FDR = bdiv_fdr)
  bdiv_df <- bdiv_df[order(rownames(bdiv_df)), ]
  write.csv(x = bdiv_df, file = paste0(out_dir, "betadiversity_comparison_pvals_", mat, ".csv"), row.names = FALSE)
  ####Plot
  rat_pw_group_dat <- data.frame(distance = c(rat_C_pw_dist, rat_S_pw_dist, rat_NS_pw_dist), 
                                 Group = c(rep("control", length(rat_C_pw_dist)),
                                           rep("septic survivor", length(rat_S_pw_dist)), 
                                           rep("septic non-survivor", length(rat_NS_pw_dist))))
  bdiv_list[[mat]] <- rat_pw_group_dat
  p <- ggplot(data = rat_pw_group_dat, mapping = aes(x = Group, y = distance, color = Group)) +
    geom_boxplot() +
    human_col_scale(name = "Group", levels = c("septic non-survivor", "control", "septic survivor", "", ""), aesthetics = c("color", "fill")) +
    scale_x_discrete(limits = as.character(unique(rat_pw_group_dat$Group))) +
    guides(color = "none") +
    #ylim(0, 100) +
    theme_bw() + 
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggsave(filename = paste0("PCA_metab_betadiv_comparison_", mat, ".png"), path = out_dir, width = 2.5, height = 4, units = "in")
  
  rat_dev_score <- cbind(rat_data_dist[, 1:4],  rds[, setdiff(colnames(rds), rat.sig.anova.car.s.class[[mat]])])
  rat_dev_max <- colMaxs(as.matrix(rat_dev_score[rat_dev_score$group != "septic non-survivor", -1:-4]))
  rat_dev_min <- colMins(as.matrix(rat_dev_score[rat_dev_score$group != "septic non-survivor", -1:-4]))
  udev <- rat_dev_score[, -1:-4] > matrix(rat_dev_max, ncol = ncol(rat_dev_score) - 4, nrow = nrow(rat_dev_score), byrow = TRUE)
  ldev <- rat_dev_score[, -1:-4] < matrix(rat_dev_min, ncol = ncol(rat_dev_score) - 4, nrow = nrow(rat_dev_score), byrow = TRUE)
  sdev <- aggregate(udev | ldev, by = list(Sample = rat_dev_score$`Sample Identification`), FUN = max) #count same metabolite at different time points as one deviation
  dev_score <- data.frame(Sample = sdev$Sample, score = rowSums(sdev[, -1]))
  dev_score$Group <- rat_dev_score$group[match(dev_score$Sample, rat_dev_score$Sample)]
  rat_dev_list[[mat]] <- dev_score
  p <- ggplot(data = dev_score, mapping = aes(fill = Group, x = score)) +
    geom_histogram(position = position_stack(), bins = max(dev_score$score) + 1) +
    human_col_scale(name = "Group", levels = c("septic non-survivor", "control", "septic survivor", "", ""), aesthetics = "fill") +
    ylab("Number of Patients") +
    xlab("Number of metabolites outside of the safe corridor at any time") +
    theme_bw() +
    theme(panel.grid = element_blank())
  ggsave(plot = p, filename = paste0("generalized_safe_corridor_minmax_", mat, ".png"), path = out_dir, width = 6, height = 3, units = "in")
  w <- which.xy(udev | ldev) # tell me which variables make a difference
  mtab <- sort(table(w[, 2]), decreasing = TRUE)
  names(mtab) <- colnames(rat_dev_score)[-1:-4][as.numeric(names(mtab))]
  print(mat)
  print(length(mtab))
  print(mtab)
  mtab_list[[mat]] <- mtab
}
for (name in names(bdiv_list))
  bdiv_list[[name]]$material <- name
rat_pgd <- Reduce("rbind", bdiv_list)
p <- ggplot(data = rat_pgd, mapping = aes(x = Group, y = distance, color = Group)) +
  facet_wrap(~ material, ncol = 3, nrow = 1) +
  geom_boxplot() +
  human_col_scale(name = "Group", levels = c("septic non-survivor", "control", "septic survivor", "", ""), aesthetics = c("color", "fill")) +
  scale_x_discrete(limits = as.character(unique(rat_pw_group_dat$Group))) +
  guides(color = "none") +
  #ylim(0, 100) +
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = paste0("PCA_metab_betadiv_comparison_crossmat.png"), path = out_dir, width = 5, height = 4, units = "in")

for (name in names(rat_dev_list)){
  rat_dev_list[[name]]$material <- name
}
dev_score <- Reduce("rbind", rat_dev_list)
dev_rep_res <- list()
for (mat in unique(rat_sepsis_data$material)){
  m <- subset(rat_sepsis_data, material == mat & `time point` == "24h")
  m <- m[, -pheno_sel]
  m <- m[, !(colnames(m) %in% rat.sig.anova.car.s.class[[mat]])]
  di <- which(m$group == "septic non-survivor")
  dev_rep_res[[mat]] <- unlist(lapply(X = 1:1000, FUN = sim_dev, n = nrow(m), d = ncol(m) - 4, dev_idx = di, sample_groups = 1:nrow(m)))
}
p_thresh <- sapply(dev_rep_res, quantile, p = 0.95, names = FALSE)
p_thresh <- data.frame(score = p_thresh, material = names(p_thresh), stringsAsFactors = FALSE)
p <- ggplot(data = dev_score, mapping = aes(fill = Group, x = score)) +
  facet_wrap(~ material, ncol = 1, nrow = 3) +
  geom_histogram(position = position_stack(), bins = max(dev_score$score) + 1) +
  human_col_scale(name = "Group", levels = c("septic non-survivor", "control", "septic survivor", "", ""), aesthetics = "fill") +
  geom_vline(data = p_thresh, mapping = aes(xintercept = score), linetype = 2) +
  ylab("Number of Rat samples") +
  xlab("Number of metabolites outside of the safe corridor at 24h") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.direction = "horizontal", legend.position = "bottom")
ggsave(plot = p, filename = paste0("generalized_safe_corridor_minmax_crossmat.png"), path = out_dir, width = 5, height = 5, units = "in")
p_thresh_box <- p_thresh
p_thresh_box$x <- as.numeric(factor(p_thresh_box$material)) - 0.4
p_thresh_box$xend <- p_thresh_box$x + 0.8
dev_score_box <- subset(dev_score, grepl("non-survivor", Group))
dev_score_box$material <- factor(dev_score_box$material)
pbox <- ggplot(data = dev_score_box, mapping = aes(fill = Group, x = material, y = score, colour = Group)) +
  geom_dotplot(stackdir = "center", binaxis = "y", dotsize = 0.7) + 
  #human_col_scale(aesthetics = c("fill", "colour")) +
  geom_segment(data = p_thresh_box, mapping = aes(x = x, xend = xend, y = score, yend = score), linetype = 2, inherit.aes = FALSE) +
  coord_flip() +
  ylab("Number of metabolites outside of the safe corridor") +
  xlab("") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.direction = "horizontal", legend.position = "bottom")
ggsave(plot = pbox, filename = "generalized_safe_corridor_minmax_crossmat_box.png", path = out_dir, width = 6, height = 2.5, units = "in")

rat_dev_comp <- human_dev_metabs
for (name in names(mtab_list)){
  an <- mtab_list[[name]]
  an <- an[names(an) %in% rat_dev_comp$Name]
  rat_dev_comp[[paste0("RatCount_", name)]] <- "-"
  rat_dev_comp[[paste0("RatCount_", name)]][match(names(an), rat_dev_comp$Name)] <- an
  rat_dev_comp[[paste0("RatCount_", name)]][rat_dev_comp$Name %in% setdiff(colnames(rat_sepsis_data), names(mtab_list[[name]]))] <- "0"
}
sum(!(rat_dev_comp$RatCount_plasma %in% c("0", "-")))
d <- rat_dev_comp[!(rat_dev_comp$RatCount_plasma %in% c("0", "-")), ]
fwrite(x = d, file = paste0(out_dir, "rat_dev_comp_human_RatPlasmaGeq1.csv"))
fwrite(x = rat_dev_comp, file = paste0(out_dir, "rat_dev_comp_human.csv"))

#Calculate beta diversity within metabolite groups
for (mat in unique(rat_sepsis_data$material)){
  for (met_group in unique(rat_sepsis_legend$group[metab_sel - 5])){
    rds <- subset(rat_sepsis_data, material == mat & `time point` == "24h", metab_sel)
    rds <- subset(rds, select = rat_sepsis_legend[rat_sepsis_legend$group == met_group, 1])
    w <- which.xy(is.na(rds))
    rds <- rds[, setdiff(1:ncol(rds), unique(w[, 2]))]
    rat_data_dist <- cbind(subset(rat_sepsis_data, material == mat & `time point` == "24h", 1:4), as.matrix(dist(x = rds, method = "canberra"))) # distance() works on a per row basis
    p <- prcomp(rat_data_dist[, -1:-4])
    #ad_rv <- ad_mg_r_df[["all"]]
    #ad_rv[ad_rv > 0.05] <- 0.05
    #text_v <- paste0(c("Nonsep", "Nonsep", "Septic-S", "Nonsep-NS"), " vs ", c("Septic-S", "Septic-NS", "Septic-NS", "Septic-NS"), ", q ", c("> ", "< ")[1 + (ad_rv < 0.05)], sapply(ad_rv, format, digits = 2))
    lam <- p$sdev[1:2] * sqrt(nrow(p$x))
    ap <- autoplot(object = p, data = rat_data_dist, colour = "group", frame = FALSE, frame.type = "norm", size = 0.5)
    ap <- ap + 
      human_col_scale(name = "Group", levels = c("septic non-survivor", "control", "septic survivor", "", ""), aesthetics = c("color", "fill")) +
      guides(colour = guide_legend(title = "Group"), fill = "none", group = "none") +
      ggtitle("PCA biplot, Canberra distance, metabolites,\n24h samples") +
      theme_bw() + 
      theme(panel.grid = element_blank())
    gobj <- ggplot_build(ap)
    xmax <- gobj$layout$panel_params[[1]]$x.range[2]
    xmin <- gobj$layout$panel_params[[1]]$x.range[1]
    ymax <- gobj$layout$panel_params[[1]]$y.range[2]
    ymin <- gobj$layout$panel_params[[1]]$y.range[1]
    arws <- make_ordination_arrows(x = p$x[, 1:2], w = rds, xmax = xmax, xmin = xmin, ymin = ymin, ymax = ymax, shorten_arw_by = 1, scale_by = 1.1)
    ap <- ap +
      geom_path(data = arws$ph_ars_s, mapping = aes(x = PC1, y = PC2, group = Group), inherit.aes = FALSE, arrow = arrow(length = unit(0.06, "inches")), size = 0.2, color = "grey50") + 
      geom_text(data = arws$ph_ars_names_s, mapping = aes(label = label, x = PC1, y = PC2), inherit.aes = FALSE, size = 1.5, color = "grey50") #+
    #  scale_x_continuous(limits = c(xmin * 1.1, xmax)) +
    #  geom_text(x = 0.9 * xmin, y = 0.85 * ymin, label = paste0(text_v, collapse = "\n"), size = 2, hjust = "left")
    ggsave(filename = paste0("PCA_biplot_metab_all_samples_", mat, "_", met_group, "_gg.png"), path = out_dir, plot = ap, width = 5, height = 4, units = "in")
    
    ###Compare beta diversity
    ####Get pairwise sample distances (here centroids == samples)
    rat_pw_dist <- as.matrix(dist(p$x[, 1:2], method = "euclidean"))
    ####Get within-group distances
    p_s_sel <- which(rat_data_dist$group == "septic survivor")
    rat_S_pw_dist <- rat_pw_dist[p_s_sel, p_s_sel]
    rat_S_pw_dist <- rat_S_pw_dist[lower.tri(x = rat_S_pw_dist)]
    p_ns_sel <- which(rat_data_dist$group == "septic non-survivor")
    rat_NS_pw_dist <- rat_pw_dist[p_ns_sel, p_ns_sel]
    rat_NS_pw_dist <- rat_NS_pw_dist[lower.tri(x = rat_NS_pw_dist)]
    p_c_sel <- which(rat_data_dist$group == "control")
    rat_C_pw_dist <- rat_pw_dist[p_c_sel, p_c_sel]
    rat_C_pw_dist <- rat_C_pw_dist[lower.tri(x = rat_C_pw_dist)]
    ####Compare distance distributions
    betadiv_S_NS <- t.test(x = rat_S_pw_dist, y = rat_NS_pw_dist, var.equal = FALSE) #var.equal set explicitly FALSE to make it visible to you, the reader
    betadiv_S_C <- t.test(x = rat_S_pw_dist, y = rat_C_pw_dist, var.equal = FALSE)
    betadiv_NS_C <- t.test(x = rat_NS_pw_dist, y = rat_C_pw_dist, var.equal = FALSE)
    
    bdiv_p <- c(betadiv_S_C$p.value, betadiv_NS_C$p.value, betadiv_S_NS$p.value)
    g1 <- list(c("S", "C"), c("NS", "C"), c("S", "NS"))
    bdiv_fdr <- p.adjust(bdiv_p, method = "fdr") #this returns one duplicate q value everytime, why??
    
    # env <- environment()
    # env_names <- names(env)
    # bdiv_names <- env_names[grep(pattern = "betadiv", x = env_names)]
    # bdiv_elements <- lapply(bdiv_names, function(e) eval(as.symbol(e)))
    # bdiv_names <- bdiv_names[sapply(bdiv_elements, class) == "htest"]
    # bdiv_elements <- bdiv_elements[sapply(bdiv_elements, class) == "htest"]
    # bdiv_p <- sapply(bdiv_elements, function(e) e$p.value)
    # bdiv_fdr <- p.adjust(bdiv_p, method = "fdr")
    # names(bdiv_p) <- bdiv_names
    # names(bdiv_fdr) <- bdiv_names
    # g1 <- sub("betadiv", "", bdiv_names)
    # g1 <- sub("_C", "_control", g1)
    # g1 <- sub("_S", "_septic survivor", g1)
    # g1 <- sub("_NS", "_septic non-survivor", g1)
    # g1 <- substring(g1, 2)
    # g1 <- strsplit(g1, "_", fixed = TRUE)
    bdiv_df <- data.frame(Group1 = sapply(g1, `[`, 1), Group2 = sapply(g1, `[`, 2), p = bdiv_p, FDR = bdiv_fdr)
    bdiv_df <- bdiv_df[order(rownames(bdiv_df)), ]
    write.csv(x = bdiv_df, file = paste0(out_dir, "betadiversity_comparison_pvals_", mat, "_", met_group, ".csv"), row.names = FALSE)
    ####Plot
    rat_pw_group_dat <- data.frame(distance = c(rat_C_pw_dist, rat_S_pw_dist, rat_NS_pw_dist), 
                                   Group = c(rep("control", length(rat_C_pw_dist)),
                                             rep("septic survivor", length(rat_S_pw_dist)), 
                                             rep("septic non-survivor", length(rat_NS_pw_dist))))
    bdiv_list[[mat]] <- rat_pw_group_dat
    p <- ggplot(data = rat_pw_group_dat, mapping = aes(x = Group, y = distance, color = Group)) +
      geom_boxplot() +
      human_col_scale(name = "Group", levels = c("septic non-survivor", "control", "septic survivor", "", ""), aesthetics = c("color", "fill")) +
      scale_x_discrete(limits = as.character(unique(rat_pw_group_dat$Group))) +
      guides(color = "none") +
      #ylim(0, 100) +
      ylab("between-sample dissimilarity") +
      theme_bw() + 
      theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    ggsave(filename = paste0("PCA_metab_betadiv_comparison_", mat, "_", met_group, ".png"), path = out_dir, width = 2.5, height = 4, units = "in")
  }
}

##Make PCA biplots for each material and metabolite group
tic()
for (mat in unique(rat_sepsis_data$material)){
  for (met_group in setdiff(unique(rat_sepsis_legend$group[match(colnames(rat_sepsis_data)[metab_sel], rat_sepsis_legend[, 1])]), "sugar")){
    met_group_idx <- which(rat_sepsis_legend$group == met_group) + 4
    rds <- subset(rat_sepsis_data, material == mat & `time point` == "24h", met_group_idx)
    w <- which.xy(is.na(rds))
    rds <- rds[, setdiff(1:ncol(rds), unique(w[, 2]))]
    rat_data_dist <- cbind(subset(rat_sepsis_data, material == mat & `time point` == "24h", 1:4), as.matrix(dist(x = rds, method = "canberra"))) # distance() works on a per row basis
    p <- prcomp(rat_data_dist[, -1:-4])
    # ad_rv <- ad_mg_r_df[[met_group]]
    # ad_rv[ad_rv > 0.05] <- 0.05
    # text_v <- paste0(c("Nonsep", "Nonsep", "Septic-S", "Nonsep-NS"), " vs ", c("Septic-S", "Septic-NS", "Septic-NS", "Septic-NS"), ", q ", c("> ", "< ")[1 + (ad_rv < 0.05)], sapply(ad_rv, format, digits = 2))
    lam <- p$sdev[1:2] * sqrt(nrow(p$x))
    rdd <- rat_data_dist
    ap <- autoplot(object = p, data = rdd, colour = "group", frame = FALSE, frame.type = "norm", size = 0.5)
    ap <- ap +
      human_col_scale(name = "Group", levels = c("septic non-survivor", "control", "septic survivor", "", ""), aesthetics = c("color", "fill")) +
      guides(colour = guide_legend(title = "Group"), fill = "none", group = "none") +
      ggtitle(paste0("PCA biplot, Canberra distance, ", met_group, "\n24h samples")) +
      theme_bw() + 
      theme(panel.grid = element_blank())
    gobj <- ggplot_build(ap)
    xmax <- gobj$layout$panel_params[[1]]$x.range[2]
    xmin <- gobj$layout$panel_params[[1]]$x.range[1]
    ymax <- gobj$layout$panel_params[[1]]$y.range[2]
    ymin <- gobj$layout$panel_params[[1]]$y.range[1]
    arws <- make_ordination_arrows(x = p$x[, 1:2], w = rds, xmax = xmax, xmin = xmin, ymin = ymin, ymax = ymax, shorten_arw_by = 1, scale_by = 1.05)
    ap <- ap +
      geom_path(data = arws$ph_ars_s, mapping = aes(x = PC1, y = PC2, group = Group), inherit.aes = FALSE, arrow = arrow(length = unit(0.06, "inches")), size = 0.2, color = "grey50") + 
      geom_text(data = arws$ph_ars_names_s, mapping = aes(label = label, x = PC1, y = PC2), inherit.aes = FALSE, size = 1.5, color = "grey50") #+
     # geom_text(x = 0.9 * xmin, y = 0.83 * ymin, label = paste0(text_v, collapse = "\n"), size = 2, hjust = "left")
    ggsave(filename = paste0("PCA_biplot_", met_group, "_", mat, "_all_samples.png"), path = out_dir, plot = ap, width = 6, height = 5, units = "in")
  }
}
toc()



#Plot PCA biplots of all samples
plot_rat_group_pca <- function(group_name = "plasma"){
  
  rat_sepsis_data_normal <- subset(rat_sepsis_data, material == group_name)
  rat_sepsis_data_normal[,-1:-4] <- scale(rat_sepsis_data_normal[,-1:-4])
  rat_sepsis_data_normal <- rat_sepsis_data_normal[, 1:which(colnames(rat_sepsis_data_normal) == "LCA")]
  
  rat_sepsis_data_normal <- rat_sepsis_data_normal[, colSums(is.na(rat_sepsis_data_normal)) < 20]
  
  col <- colnames(rat_sepsis_data_normal)
  colnames(rat_sepsis_data_normal) <- make.names(col)
  rat_sepsis_data_normal[, -1:-4] <- missRanger(rat_sepsis_data_normal[,-1:-4])
  colnames(rat_sepsis_data_normal) <- col
  
  p <- prcomp(rat_sepsis_data_normal[,-1:-4])
  
  png(filename = paste0(out_dir_cs, "rat_cs_expl_var_", group_name, ".png"), width = 300, height = 300, units = "px")
  barplot(summary(p)$importance[2,], ylab = "Explained variance", main = paste0(group_name, ", PCA - explained variance"))
  dev.off()
  
  gp <- ggplot(data = cbind(rat_sepsis_data_normal[, 1:4], data.frame(x = p$x[,1], y = p$x[,2])), aes_string(x = "x", y = "y", color = "group", shape = "`time point`")) +
    geom_point() + 
    xlab("PC1") +
    ylab("PC2") +
    ggtitle(paste0(group_name , " - Experimental groups seperate in PC space")) +
    theme_bw()
  ggsave(plot = gp, filename = paste0(out_dir_cs, "rat_cs_pcs_", group_name, ".png"), width = 15, height = 10, units = "cm")
  
  p <- prcomp(subset(rat_sepsis_data_normal, grepl(pattern = "septic", x = group), -1:-4))
  
  png(filename = paste0(out_dir_ns, "rat_ns_expl_var_", group_name, ".png"), width = 300, height = 300, units = "px")
  barplot(summary(p)$importance[2,], ylab = "Explained variance", main = paste0(group_name, ", PCA - explained variance"))
  dev.off()
  
  gp <- ggplot(data = cbind(subset(rat_sepsis_data_normal, grepl(pattern = "septic", x = group), 1:4), data.frame(x = p$x[,1], y = p$x[,2])), aes_string(x = "x", y = "y", color = "group", shape = "`time point`")) +
    geom_point() + 
    xlab("PC1") +
    ylab("PC2") +
    ggtitle(paste0(group_name , " - Survivors and Nonsurvivors clearly seperate in PC space")) +
    theme_bw()
  ggsave(plot = gp, filename = paste0(out_dir_ns, "rat_ns_pcs_", group_name, ".png"), width = 15, height = 10, units = "cm")
  
  for (gr in unique(rat_sepsis_data_normal$group)){
    p <- prcomp(subset(rat_sepsis_data_normal, grepl(pattern = gr, x = group), -1:-4))
    
    png(filename = paste0(out_dir_ns, "rat_ns_expl_var_", group_name, "_", gr, ".png"), width = 300, height = 300, units = "px")
    barplot(summary(p)$importance[2,], ylab = "Explained variance", main = paste0(group_name, ", PCA - explained variance"))
    dev.off()
    
    gp <- ggplot(data = cbind(subset(rat_sepsis_data_normal, grepl(pattern = gr, x = group), 1:4), data.frame(x = p$x[,1], y = p$x[,2])), aes_string(x = "x", y = "y", shape = "`time point`")) +
      geom_point() + 
      xlab("PC1") +
      ylab("PC2") +
      ggtitle(group_name) +
      theme_bw()
    ggsave(plot = gp, filename = paste0(out_dir_ns, "rat_ns_pcs_", group_name, "_", gr, ".png"), width = 15, height = 10, units = "cm")
  }
  
  for (gr in unique(rat_sepsis_data_normal$group)){
    c <- cmdscale(d = dist(subset(rat_sepsis_data_normal, grepl(pattern = gr, x = group), -1:-4)))
    
    gp <- ggplot(data = cbind(subset(rat_sepsis_data_normal, grepl(pattern = gr, x = group), 1:4), data.frame(x = c[,1], y = c[,2])), aes_string(x = "x", y = "y", shape = "`time point`")) +
      geom_point() + 
      xlab("PC1") +
      ylab("PC2") +
      ggtitle(group_name) +
      theme_bw()
    ggsave(plot = gp, filename = paste0(out_dir_ns, "rat_ns_mds_", group_name, "_", gr, ".png"), width = 15, height = 10, units = "cm")
  }
}

plot_rat_group_pca(group_name = "plasma")
plot_rat_group_pca(group_name = "liver")
plot_rat_group_pca(group_name = "heart")
