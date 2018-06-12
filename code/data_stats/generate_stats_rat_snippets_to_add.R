library(matrixStats)
library(data.table)
library(missRanger)
library(ggplot2)

source("../function_definitions.R")

out_dir_cs <- "../../results/data_stats_rat_control_vs_septic/"
out_dir_ns <- "../../results/data_stats_rat_surv_vs_nonsurv/"

rat_sepsis_data <- get_rat_sepsis_data()

rat_sepsis_data <- rat_sepsis_data[!rat_sepsis_data$`Sample Identification` %in% c("060H", "039L"), ]

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
