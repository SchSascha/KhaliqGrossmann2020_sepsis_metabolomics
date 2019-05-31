library(matrixStats)
library(data.table)
library(missRanger)
library(ggplot2)
library(tictoc)
library(parallel)
library(vegan)
library(ggrepel)

source("../function_definitions.R")

out_dir_cs <- "../../results/data_stats_rat_control_vs_septic/"
out_dir_ns <- "../../results/data_stats_rat_surv_vs_nonsurv/"

out_dir <- "../../results/data_stats_rat/"

num_cores <- 2

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
tic()
ad_res_mc <- mclapply(unique(pat_list$group), 
                      function(dx, PCs, ad_data){
                        take_smpls <- ad_data$group != dx
                        Y <- PCs[take_smpls, ]
                        ad_data_t <- ad_data[take_smpls, ]
                        ad_res <- adonis(formula = Y ~ group, data = ad_data_t, permutations = 10000, parallel = 1, method = "euclidean")
                        list(ad_res = ad_res)
                      }, 
                      PCs = PCs, ad_data = ad_data,
                      mc.cores = num_cores)
toc()
ad_r_p <- sapply(lapply(lapply(lapply(ad_res_mc, `[[`, "ad_res"), `[[`, "aov.tab"), `[[`, "Pr(>F)"), `[`, 1)

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
  ggtitle("PCA biplot, Canberra distance,\nbiochemical parameters") +
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
    ggtitle("PCA biplot, Canberra distance, metabolites,\nall samples") +
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
}

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
      ggtitle(paste0("PCA biplot, Canberra distance, ", met_group, "\nall samples")) +
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

#Plot
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
