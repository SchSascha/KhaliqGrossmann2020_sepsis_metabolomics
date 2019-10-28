# container file for plots no longer needed or interesting, mainly everythin done with heatmaply

##Human, cluster-heatmap, all metabolites, but groups at the side
x <- na.omit(data.frame(group = coarse_group_list[metab_sel - 6], human_sepsis_data_normal_conc_metab_cov));
heatmaply(x = x[, -1], row_side_colors = x[c("group")], key.title = "Cov", dendrogram = FALSE, showticklabels = FALSE, file = paste0(out_dir, "human_normal_metab_nonclust.png"))
##Human, cluster-heatmap, phenomenological vars, but groups at the side
x <- na.omit(data.frame(group = coarse_group_list[pheno_sel - 6], inv_human_sepsis_data_normal_conc_pheno_cov));
heatmaply(x = x[, -1], row_side_colors = x[c("group")], key.title = "pcor", dendrogram = FALSE, showticklabels = FALSE, file = paste0(out_dir, "human_normal_pheno_nonclust_pcor.png"))

x <- na.omit(data.frame(group = coarse_group_list[metab_sel - 6], inv_human_sepsis_data_normal_S_conc_metab_cov));
heatmaply(x = x[, -1], row_side_colors = x[c("group")], key.title = "pcor", dendrogram = FALSE, showticklabels = FALSE, file = paste0(out_dir, "human_normal_metab_S_nonclust_pcor.png"))
x <- na.omit(data.frame(group = coarse_group_list[metab_sel - 6], inv_human_sepsis_data_normal_NS_conc_metab_cov));
heatmaply(x = x[, -1], row_side_colors = x[c("group")], key.title = "pcor", dendrogram = FALSE, showticklabels = FALSE, file = paste0(out_dir, "human_normal_metab_NS_nonclust_pcor.png"))
x <- na.omit(data.frame(group = coarse_group_list[pheno_sel - 6], inv_human_sepsis_data_normal_S_conc_pheno_cov));
heatmaply(x = x[, -1], row_side_colors = x[c("group")], key.title = "pcor", dendrogram = FALSE, showticklabels = FALSE, file = paste0(out_dir, "human_normal_pheno_S_nonclust_pcor.png"))
x <- na.omit(data.frame(group = coarse_group_list[pheno_sel - 6], inv_human_sepsis_data_normal_NS_conc_pheno_cov));
heatmaply(x = x[, -1], row_side_colors = x[c("group")], key.title = "pcor", dendrogram = FALSE, showticklabels = FALSE, file = paste0(out_dir, "human_normal_pheno_NS_nonclust_pcor.png"))
rm("x")

##Human, cluster-heatmap, all phenom. vars, groups at the top
x <- na.omit(human_sepsis_data_normal[, c(1:6, pheno_sel)]);
heatmaply(x = x[, -1:-6], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = data.frame(Group = coarse_group_list[pheno_sel - 6]), plot_method = "plotly", margins = c(100,50,0,150))
rm("x")

##Human, cluster-heatmap, phenom. var groups
x <- na.omit(human_sepsis_data_normal_grouped[, c(1:6, group_pheno_sel)])
h <- heatmaply(x = x[, -1:-6], row_side_colors = x[c("Survival", "CAP / FP")], plot_method = "ggplot", margins = c(150,100,30,0), key.title = "Normalized\nConcentration", main = "Phenom. groups, survivors and nonsurvivors", subplot_heights = c(0.1, 0.9))
h$width <- 800
h$height <- 1100
export(p = h, file = paste0(out_dir, "human_normal_pheno_grouped.png"))
x <- na.omit(subset(human_sepsis_data_normal_grouped[, c(1:6, group_pheno_sel)], Survival == "S"))
h <- heatmaply(x = x[, -1:-6], row_side_colors = x[c("CAP / FP")], plot_method = "ggplot", margins = c(150,50,30,150), key.title = "Normalized\nConcentration", main = "Phenom. groups, survivors", subplot_heights = c(0.1, 0.9))
h$width <- 800
h$height <- 800
export(p = h, file = paste0(out_dir, "human_normal_pheno_S_grouped.png"))
x <- na.omit(subset(human_sepsis_data_normal_grouped[, c(1:6, group_pheno_sel)], Survival == "NS"))
h <- heatmaply(x = x[, -1:-6], row_side_colors = x[c("CAP / FP")], plot_method = "ggplot", margins = c(150,100,30,150), key.title = "Normalized\nConcentration", main = "Phenom. groups, nonsurvivors", subplot_heights = c(0.1, 0.9))
h$width <- 800
h$height <- 600
export(p = h, file = paste0(out_dir, "human_normal_pheno_NS_grouped.png"))
rm("x")

##Human, cluster-heatmap, metab groups
x <- na.omit(human_sepsis_data_normal_grouped[, c(1:6, group_metab_sel)])
h <- heatmaply(x = x[, -1:-6], row_side_colors = x[c("Survival", "CAP / FP")], plot_method = "ggplot", margins = c(150,100,30,0), key.title = "Normalized\nConcentration", main = "Metab. groups, survivors and nonsurvivors", subplot_heights = c(0.1, 0.9))
h$width <- 800
h$height <- 1100
export(p = h, file = paste0(out_dir, "human_normal_metab_grouped.png"))
x <- na.omit(subset(human_sepsis_data_normal_grouped[, c(1:6, group_metab_sel)], Survival == "S"))
h <- heatmaply(x = x[, -1:-6], row_side_colors = x[c("CAP / FP")], plot_method = "ggplot", margins = c(150,100,30,150), key.title = "Normalized\nConcentration", main = "Metab. groups, survivors", subplot_heights = c(0.1, 0.9))
h$width <- 800
h$height <- 800
export(p = h, file = paste0(out_dir, "human_normal_metab_S_grouped.png"))
x <- na.omit(subset(human_sepsis_data_normal_grouped[, c(1:6, group_metab_sel)], Survival == "NS"))
h <- heatmaply(x = x[, -1:-6], row_side_colors = x[c("CAP / FP")], plot_method = "ggplot", margins = c(150,100,30,150), key.title = "Normalized\nConcentration", main = "Metab. groups, nonsurvivors", subplot_heights = c(0.1, 0.9))
h$width <- 800
h$height <- 600
export(p = h, file = paste0(out_dir, "human_normal_metab_NS_grouped.png"))
rm("x")

##Human, cluster-heatmap, coarse grouped metabolites
x <- na.omit(human_sepsis_data_normal[,c(1:6, metab_sel)])
heatmaply(x = x[,-1:-6], row_side_colors = x[c("Survival", "CAP / FP")], file = paste0(out_dir, "human_normal_metab.png"), main = "Metabolite profile does not cluster survival well", key.title = "Concentration\n(normalized)", showticklabels = TRUE)
x <- na.omit(human_sepsis_data_normal[,c(1:6, pheno_sel)])
heatmaply(x = x[,-1:-6], row_side_colors = x[c("Survival", "CAP / FP")], file = paste0(out_dir, "human_normal_pheno.png"), main = "Phenomenological profile has survival clusters", key.title = "Concentration\n(normalized)", showticklabels = TRUE)
x <- na.omit(human_sepsis_data_normal_grouped[,c(1:6, group_metab_sel)])
heatmaply(x = x[,-1:-6], row_side_colors = x[c("Survival", "CAP / FP")], file = paste0(out_dir, "human_normal_metab_grouped.png"), main = "Metablite group profiles somewhat cluster survival", key.title = "Concentration\n(normalized)", showticklabels = TRUE)
x <- na.omit(human_sepsis_data_normal_grouped[,c(1:6, group_pheno_sel)])
heatmaply(x = x[,-1:-6], row_side_colors = x[c("Survival", "CAP / FP")], file = paste0(out_dir, "human_normal_pheno_grouped.png"), main = "Pheno var group profiles [?] survival", key.title = "Concentration\n(normalized)", showticklabels = TRUE)
rm("x")

##Human, cluster-heatmap, one per metabolite group
x <- subset(rbind(human_sepsis_data, human_nonsepsis_data), Day %in% tanova_day_set)
x <- max_norm(x, -1:-6)
x$Day <- as.numeric(as.character(x$Day))
x$Survival <- as.character(x$Survival)
x$Survival[x$`CAP / FP` == "-"] <- "Control"
x$Survival <- reorder(x$Survival, (x$Survival == "Control") + (2 * (x$Survival == "S")) + (3 * (x$Survival == "NS")))
x$Day <- reorder(x$Day, x$Day)
x <- x[order(x$Survival),]
x <- x[order(x$Day),]
xm <- x[, c(1:6, metab_sel)]
xmt <- data.frame(t(xm[, -1:-6]))
rownames(xmt) <- colnames(xm[,-1:-6])
colnames(xmt) <- xm$`Sample ID`

control.sig <- rownames(xmt[, -1:-6]) %in% sig.anova.car.c.class
survival.sig <- rownames(xmt[, -1:-6]) %in% sig.anova.car.s.class
mat_sigs <- data.frame(control.sig = control.sig, survival.sig = survival.sig, stringsAsFactors = FALSE)
mat_sigs <- lapply(mat_sigs, function(x){ c("nonsignif.", "p < 0.05")[x + 1] })
mat_sigs <- data.frame(mat_sigs)
colnames(mat_sigs) <- c("non-Septic vs Sepsis", "S vs NS")

lower_margin <- 85
for (met_group in unique(coarse_group_list[metab_sel - 6])){
  group_sel <- coarse_group_list[metab_sel - 6] %in% met_group
  #xfplotdat <- t(max_norm(t(xmt[group_sel, xm$material == mat])))
  xfplotdat <- xmt[group_sel, ]
  sel <- !rowAlls(is.na(xfplotdat))
  top_row_h <- 0.03 * 76/sum(sel)
  subplot_h <- c(top_row_h, 1 - top_row_h)
  if (sum(sel) > 1){
    h <- heatmaply(x = xfplotdat[sel, ],
                   dendrogram = "row", 
                   plot_method = "plotly", 
                   col_side_colors = subset(xm, TRUE , c("Day", "Survival")), 
                   row_side_colors = mat_sigs[group_sel, ][sel, ], 
                   key.title = "concentration", 
                   margins = c(lower_margin,100,NA,50), 
                   height = lower_margin + round(sum(sel) * 1000/76),
                   subplot_heights = subplot_h, 
                   subplot_widths = c(0.91, 0.03, 0.06))
    h$width <- 2200
    #h$height <- 1200
    h$height <- lower_margin + round(sum(sel) * 1100/76) #1200 is a good height for 76 rows of metabolites
    export(p = h, file = paste0(out_dir, "human_heatmap_", met_group, ".png"))
  }
  else
    print(paste0(met_group, " has too few metabolites, possibly after filtering all-NA rows"))
}

##Human, cluster heatmap, significant metabs only
x <- subset(rbind(human_sepsis_data, human_nonsepsis_data), Day %in% tanova_day_set)
x <- max_norm(x, -1:-6)
x$Day <- as.numeric(as.character(x$Day))
x$Survival <- as.character(x$Survival)
x$Survival[x$`CAP / FP` == "-"] <- "Control"
x$Survival <- reorder(x$Survival, (x$Survival == "Control") + (2 * (x$Survival == "S")) + (3 * (x$Survival == "NS")))
x$Day <- reorder(x$Day, x$Day)
x <- x[order(x$Survival),]
x <- x[order(x$Day),]
xm <- x[, c(1:6, metab_sel)]
xmt <- data.frame(t(xm[, -1:-6]))
rownames(xmt) <- colnames(xm[,-1:-6])
colnames(xmt) <- xm$`Sample ID`

xmts <- xmt[rownames(xmt) %in% sig.anova.car.s.class, ]
xmtc <- xmt[rownames(xmt) %in% sig.anova.car.c.class, ]

lower_margin <- 85
sel <- !rowAlls(is.na(xmts))
top_row_h <- 0.03 * 76/sum(sel)
subplot_h <- c(top_row_h, 1 - top_row_h)
h <- heatmaply(x = xmts,
               dendrogram = "none", 
               plot_method = "plotly", 
               col_side_colors = subset(xm, TRUE , c("Day", "Survival")), 
               key.title = "concentration", 
               margins = c(lower_margin,100,NA,50), 
               height = lower_margin + round(sum(sel) * 1000/76),
               subplot_heights = subplot_h, 
               subplot_widths = c(0.91, 0.03, 0.06)[1])
h$width <- 2200
#h$height <- 1200
h$height <- lower_margin + round(sum(sel) * 1100/76) #1200 is a good height for 76 rows of metabolites
export(p = h, file = paste0(out_dir, "human_heatmap_metab_sig_s.png"))

#Also for pheno vars
xm <- x[, c(1:6, pheno_sel)]
xmt <- data.frame(t(xm[, -1:-6]))
rownames(xmt) <- colnames(xm[,-1:-6])
colnames(xmt) <- xm$`Sample ID`
survival.sig <- rownames(xmt[, -1:-6]) %in% sig.anova.car.s.pheno.class
mat_sigs <- data.frame(survival.sig = survival.sig, stringsAsFactors = FALSE)
mat_sigs <- lapply(mat_sigs, function(x){ c("nonsignif.", "p < 0.05")[x + 1] })
mat_sigs <- data.frame(mat_sigs)
colnames(mat_sigs) <- c("Control vs Sepsis", "S vs NS")

lower_margin <- 85
xfplotdat <- xmt
sel <- !rowAlls(is.na(xfplotdat))
top_row_h <- 0.03 * 76/sum(sel)
subplot_h <- c(top_row_h, 1 - top_row_h)
h <- heatmaply(x = xfplotdat[sel, ],
               dendrogram = "row", 
               plot_method = "plotly", 
               col_side_colors = subset(xm, TRUE , c("Day", "Survival")), 
               row_side_colors = mat_sigs[sel, ], 
               key.title = "concentration", 
               margins = c(lower_margin,100,NA,50), 
               height = lower_margin + round(sum(sel) * 1000/76),
               subplot_heights = subplot_h, 
               subplot_widths = c(0.91, 0.03, 0.06))
h$width <- 2200
h$height <- lower_margin + round(sum(sel) * 1100/76) #1200 is a good height for 76 rows of metabolites
export(p = h, file = paste0(out_dir, "human_heatmap_pheno.png"))

##Human, variance of metabolites, diff of NS vs S, ordered by |diff|, lables colored by metabolite group, 
for (d in unique(metab_normal_day_var_df$Day)){
  metab_day0_var_df <- cbind(subset(metab_normal_day_var_df, Day == d), subset(pheno_day_var_df, Day == d, -1:-2))
  metab_day_vardiff_df <- metab_day0_var_df[1, -1:-2] - metab_day0_var_df[2, -1:-2]
  metab_day0_var_df <- metab_day0_var_df[, c(1, 2, 2 + order(metab_day_vardiff_df, decreasing = TRUE))]
  metab_day0_var_df[, -1:-2] <- scale(metab_day0_var_df[, -1:-2], scale = FALSE, center = TRUE)
  metab_day0_var_long_df <- melt(metab_day0_var_df, id.vars = c("Day", "Survival"))
  metab_day0_var_long_df$group <- human_sepsis_legend$group[match(x = metab_day0_var_long_df$variable, table = human_sepsis_legend[, 1])]
  metab_day0_var_long_df <- subset(metab_day0_var_long_df, !group %in% "Excluded")
  metab_day0_group_df <- metab_day0_var_long_df
  metab_day0_group_df$width <- 1
  metab_day0_group_df$height <- 0.5
  metab_day0_group_df$value <- -3.25
  metab_day0_group_df$group[metab_day0_group_df$group %in% coarse_group_list[pheno_sel - 5]] <- "clinical parameter"
  metab_day0_group_df$Metabolite_group <- metab_day0_group_df$group
  metab_day0_group_df <- metab_day0_group_df[!duplicated(metab_day0_group_df$variable), ]
  vardiffplot <- ggplot(data = metab_day0_var_long_df, mapping = aes(x = variable, y = value, group = Survival, color = Survival)) + 
    geom_point() + 
    geom_tile(mapping = aes(x = variable, y = value, fill = Metabolite_group, width = width, height = height), data = metab_day0_group_df, inherit.aes = FALSE) +
    ylab("Mean-free difference of variance") +
    xlab("Metabolite") +
    scale_y_continuous(limits = c(-3.5, 3), expand = c(0, 0)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 4, vjust = 0.5), panel.grid.major.x = element_line(linetype = 0))
  ggsave(plot = vardiffplot, filename = paste0("human_var_diff_orderd_day", d, ".png"), path = out_dir, width = 12, height = 6, units = "in")
}
##All days combined
metab_all_days_var_df <- rbind(data.frame(Survival = "NS", t(colVars(as.matrix(subset(human_sepsis_data, Survival == "NS",-1:-6))))),
                               data.frame(Survival = "S", t(colVars(as.matrix(subset(human_sepsis_data, Survival == "S", -1:-6))))))
colnames(metab_all_days_var_df)[-1] <- colnames(human_sepsis_data)[-1:-6]
metab_all_days_mean <- colMeans(human_sepsis_data[, -1:-5])
metab_all_days_var_df[, -1] <- scale(metab_all_days_var_df[, -1], center = FALSE, scale = metab_all_days_mean)
metab_day_vardiff_df <- metab_all_days_var_df[1, -1] - metab_all_days_var_df[2, -1]
metab_all_days_var_df <- metab_all_days_var_df[c(1, 1 + order(metab_day_vardiff_df, decreasing = TRUE))]
metab_all_days_var_df[, -1] <- scale(metab_all_days_var_df[, -1], scale = FALSE, center = TRUE)
metab_all_days_var_df <- metab_all_days_var_df[, c(1, 1 + which(metab_all_days_var_df[1, -1] < 0))]
metab_all_days_var_long_df <- melt(metab_all_days_var_df, id.vars = c("Survival"))
metab_all_days_var_long_df$group <- human_sepsis_legend$group[match(x = metab_all_days_var_long_df$variable, table = human_sepsis_legend[, 1])]
metab_all_days_var_long_df <- subset(metab_all_days_var_long_df, !group %in% "Excluded")
metab_all_days_group_df <- metab_all_days_var_long_df
metab_all_days_group_df$width <- 1
metab_all_days_group_df$height <- 1
metab_all_days_group_df$value <- -800
metab_all_days_group_df$group[metab_all_days_group_df$group %in% coarse_group_list[pheno_sel - 6]] <- "clinical parameter"
metab_all_days_group_df$Metabolite_group <- metab_all_days_group_df$group
metab_all_days_group_df <- metab_all_days_group_df[!duplicated(metab_all_days_group_df$variable), ]
vardiffplot <- ggplot(data = metab_all_days_var_long_df, mapping = aes(x = variable, y = value, group = Survival, color = Survival)) + 
  geom_point() + 
  geom_tile(mapping = aes(x = variable, y = value, fill = Metabolite_group, width = width, height = height), data = metab_all_days_group_df, inherit.aes = FALSE) +
  ylab("Centered group variance rel. to metabolite mean") +
  xlab("Metabolite") +
  human_col_scale(name = "Survival", aesthetics = "colour") +
  scale_y_continuous(trans = pseudo_log_trans(sigma = 0.25, base = 2), expand = c(0, 0), limits = c(-1200, 800)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5), panel.grid.major.x = element_line(linetype = 0))
ggsave(plot = vardiffplot, filename = paste0("human_var_diff_ordered_all_days.png"), path = out_dir, width = 12, height = 6, units = "in")

##Human, variance of metabolites, NS vs. S
metab_var_df <- human_sepsis_data_normal_conc_var[,c(2, metab_sel - 3)]
colnames(metab_var_df)[-1] <- colnames(human_sepsis_data)[metab_sel]
metab_var_long_df <- melt(metab_var_df, id.vars = "Survival")
metab_var_long_df$group <- human_sepsis_legend$group[match(metab_var_long_df$variable, human_sepsis_legend[,1])]
metab_var_long_df$value <- log10(metab_var_long_df$value)
metab_var_long_df <- metab_var_long_df[is.finite(metab_var_long_df$value),]
metab_var_long_df <- subset(metab_var_long_df, variable != "sugar")

##Do F-test for variance difference
metab_var_sig_df <- data.frame(group = setdiff(unique(metab_var_long_df$group), "sugar"), sig = 1, value = 0, stringsAsFactors = FALSE)
for (gr in metab_var_sig_df$group){
  b1 <- subset(metab_var_long_df, group == gr & Survival == "S", value)
  b2 <- subset(metab_var_long_df, group == gr & Survival == "NS", value)
  if (nrow(b1) > 1 & nrow(b2) > 1)
    metab_var_sig_df$sig[metab_var_sig_df$group == gr] <- t.test(x = b1, y = b2)$p.value # t-test is very similar to F-test
}
metab_var_sig_df$sig <- p.adjust(p = metab_var_sig_df$sig, method = "fdr")
metab_var_sig_df$value = 2
metab_var_sig_df <- subset(metab_var_sig_df, sig <= 0.05)

metab_var_plot <- ggplot(data = metab_var_long_df, mapping = aes(x = group, y = value, fill = Survival)) +
  geom_boxplot() +
  geom_point(data = metab_var_sig_df, mapping = aes(x = group, y = value), inherit.aes = FALSE, shape = 8, size = 1) +
  ylab("variance (log)") + 
  #scale_y_log10() +
  ggtitle("Concentration variance over all days\ndiffers with regard to survival") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = metab_var_plot, filename = "human_all_days_grouped_metab_var_patientcentered.png", path = out_dir, width = 7, height = 4, units = "in")

##Human, variance of phenomenological vars, NS vs. S
pheno_var_df <- na.omit(human_sepsis_data_normal_conc_var[, c(2, pheno_sel - 4)])
colnames(pheno_var_df)[-1] <- colnames(human_sepsis_data)[pheno_sel]
colnames(pheno_var_df)[colnames(pheno_var_df) == "Creatinine.1"] <- "Creatinine"
pheno_var_long_df <- melt(pheno_var_df, id.vars = "Survival")
pheno_var_long_df$group <- human_sepsis_pheno_var_groups$group[match(pheno_var_long_df$variable, human_sepsis_pheno_var_groups[,1])]
pheno_var_long_df$value <- log10(pheno_var_long_df$value)
pheno_var_long_df <- pheno_var_long_df[is.finite(pheno_var_long_df$value),]

##Do F-test for variance difference
pheno_var_sig_df <- data.frame(group = unique(pheno_var_long_df$group), sig = 1, value = 0, stringsAsFactors = FALSE)
for (gr in pheno_var_sig_df$group){
  b1 <- subset(pheno_var_long_df, group == gr & Survival == "S", value)
  b2 <- subset(pheno_var_long_df, group == gr & Survival == "NS", value)
  if (nrow(b1) > 1 & nrow(b2) > 1)
    pheno_var_sig_df$sig[pheno_var_sig_df$group == gr] <- t.test(x = b1, y = b2)$p.value # t-test is very similar to F-test
}
pheno_var_sig_df$sig <- p.adjust(p = pheno_var_sig_df$sig, method = "fdr")
pheno_var_sig_df$value = 4
pheno_var_sig_df <- subset(pheno_var_sig_df, sig <= 0.05)

pheno_var_plot <- ggplot(data = pheno_var_long_df, mapping = aes(x = group, y = value, fill = Survival)) +
  geom_boxplot() +
  geom_point(data = pheno_var_sig_df, mapping = aes(x = group, y = value), inherit.aes = FALSE, shape = 8, size = 1) +
  #scale_y_log10() +
  ylab("variance (log)") + 
  ggtitle("Concentration variance over all days\ndiffers with regard to survival") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = pheno_var_plot, filename = "human_all_days_grouped_pheno_var_patientcentered.png", path = out_dir, width = 9, height = 4, units = "in")

##Human, variance of metab vars over all days but metabolite over all patients
pheno_day_var_df <- rbind(cbind(data.frame(Survival = "NS"), t(human_sepsis_data_normal_NS_conc_pheno_var)), 
                          cbind(data.frame(Survival = "S"), t(human_sepsis_data_normal_S_conc_pheno_var)))
colnames(pheno_day_var_df) <- c("Survival", colnames(human_sepsis_data)[pheno_sel])
pheno_day_var_long_df <- melt(pheno_day_var_df, id.vars = c("Survival"))
pheno_day_var_long_df$group <- human_sepsis_legend$group[match(pheno_day_var_long_df$variable, human_sepsis_legend[[1]])]
pheno_day_var_plot <- ggplot(data = pheno_day_var_long_df, mapping = aes(x = as.factor(group), y = value, fill = Survival)) +
  geom_boxplot(outlier.size = 0.7) +
  ylab("variance") +
  xlab("clinical parameter") +
  #ggtitle("Patient-wise concentration variances differ") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "human_all_days_grouped_pheno_var.png", path = out_dir, plot = pheno_day_var_plot, width = 9, height = 3.5, units = "in")

##Human, variance of grouped metab vars, all days
metab_normal_day_var_df <- rbind(cbind(data.frame(Survival = "NS"), t(human_sepsis_data_normal_NS_conc_metab_var)), 
                                 cbind(data.frame(Survival = "S"), t(human_sepsis_data_normal_S_conc_metab_var)))
colnames(metab_normal_day_var_df) <- c("Survival", colnames(human_sepsis_data)[metab_sel])
metab_day_var_long_df <- melt(metab_normal_day_var_df, id.vars = c("Survival"))
metab_day_var_long_df$group <- human_sepsis_legend$group[match(metab_day_var_long_df$variable, human_sepsis_legend[[1]])]
metab_day_var_plot <- ggplot(data = metab_day_var_long_df, mapping = aes(x = group, y = value, fill = Survival)) +
  geom_boxplot(outlier.size = 0.7) +
  ylab("variance") +
  xlab("metabolite group") +
  #ggtitle("Patient-wise concentration variances differ") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "human_all_days_grouped_metab_var.png", path = out_dir, plot = metab_day_var_plot, width = 8, height = 3.5, units = "in")

##Human, variance of pheno vars seperate by day
pheno_day_var_df <- rbind(cbind(data.frame(Day = tanova_day_set, Survival = "NS"), human_sepsis_data_normal_NS_conc_pheno_day_var), 
                          cbind(data.frame(Day = tanova_day_set, Survival = "S"), human_sepsis_data_normal_S_conc_pheno_day_var))
colnames(pheno_day_var_df) <- c("Day", "Survival", colnames(human_sepsis_data)[pheno_sel])
pheno_day_var_long_df <- melt(pheno_day_var_df, id.vars = c("Day", "Survival"))
pheno_day_var_long_df$group <- human_sepsis_legend$group[match(pheno_day_var_long_df$variable, human_sepsis_legend[[1]])]
pheno_day_var_long_df <- subset(pheno_day_var_long_df, !group %in% "Excluded")
pheno_day_var_sig <- melt(anova.car.ph.pheno.variance.sig.contr)
pheno_day_var_sig <- subset(pheno_day_var_sig, L1 %in% pheno_day_var_long_df$group)
pheno_day_var_sig$Day <- tanova_day_set[pheno_day_var_sig$value]
pheno_day_var_sig$value <- sapply(pheno_day_var_sig$L1, function(gr) max(subset(pheno_day_var_long_df, group == gr)$value)) * 1.1
colnames(pheno_day_var_sig)[2] <- "group"
pheno_day_var_plot <- ggplot(data = pheno_day_var_long_df, mapping = aes(x = factor(Day), y = value, fill = Survival, group = Survival, color = Survival)) +
  facet_wrap( ~ group, ncol = 5, nrow = 2, scales = "free_y") +
  geom_point(position = position_dodge(width = 0.5), size = 0.8) +
  stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax = "mean", geom = "errorbar", position = position_dodge(width = 0.5), width = 0.5) +
  geom_point(aes(y = value, x = factor(Day, levels = tanova_day_set)), pheno_day_var_sig, shape = 8, size = 0.8, inherit.aes = FALSE) +
  ylab("variance") +
  xlab("Day") +
  #ggtitle("Patient-wise concentration variances differ between days") +
  theme_bw()
ggsave(filename = "human_seperate_days_grouped_pheno_var.png", path = out_dir, plot = pheno_day_var_plot, width = 9, height = 3, units = "in")

##Human, variance of pheno vars seperate by day
pheno_day_var_df <- rbind(cbind(data.frame(Day = tanova_day_set, Survival = "NS"), human_sepsis_data_normal_NS_conc_pheno_day_var), 
                          cbind(data.frame(Day = tanova_day_set, Survival = "S"), human_sepsis_data_normal_S_conc_pheno_day_var))
colnames(pheno_day_var_df) <- c("Day", "Survival", colnames(human_sepsis_data)[pheno_sel])
pheno_day_var_long_df <- melt(pheno_day_var_df, id.vars = c("Day", "Survival"))
pheno_day_var_long_df$group <- human_sepsis_legend$group[match(pheno_day_var_long_df$variable, human_sepsis_legend[[1]])]
pheno_day_var_long_df <- subset(pheno_day_var_long_df, !group %in% c("Excluded", "Anabolism-associated", "Down under stress", "Kidney-associated", "Lipid-associated", "Up under stress"))
pheno_day_var_sig <- melt(anova.car.ph.pheno.variance.sig.contr)
pheno_day_var_sig <- subset(pheno_day_var_sig, L1 %in% pheno_day_var_long_df$group)
pheno_day_var_sig$Day <- tanova_day_set[pheno_day_var_sig$value]
pheno_day_var_sig$value <- sapply(pheno_day_var_sig$L1, function(gr) max(subset(pheno_day_var_long_df, group == gr)$value)) * 1.1
colnames(pheno_day_var_sig)[2] <- "group"
pheno_day_var_plot <- ggplot(data = pheno_day_var_long_df, mapping = aes(x = factor(Day), y = value, fill = Survival, group = Survival, color = Survival)) +
  facet_wrap( ~ group, ncol = 5, nrow = 2, scales = "free_y") +
  geom_point(position = position_dodge(width = 0.4)) +
  stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax = "mean", geom = "errorbar", position = position_dodge(width = 0.4), width = 0.4) +
  geom_point(aes(y = value, x = factor(Day, levels = tanova_day_set)), pheno_day_var_sig, shape = 8, size = 0.8, inherit.aes = FALSE) +
  ylab("variance") +
  xlab("Day") +
  #ggtitle("Patient-wise concentration variances differ between days") +
  theme_bw()
ggsave(filename = "human_seperate_days_grouped_pheno_var_PG.png", path = out_dir, plot = pheno_day_var_plot, width = 9, height = 2, units = "in")

##Human, variance of metab vars seperate by day
metab_day_var_long_df <- melt(metab_normal_day_var_df, id.vars = c("Day", "Survival"))
metab_day_var_long_df$group <- human_sepsis_legend$group[match(metab_day_var_long_df$variable, human_sepsis_legend[[1]])]
metab_day_var_plot <- ggplot(data = metab_day_var_long_df, mapping = aes(x = factor(Day), y = value, fill = Survival)) +
  facet_wrap( ~ group, ncol = 4, nrow = 2, scales = "free_y") +
  geom_boxplot(outlier.size = 0.7) +
  ylab("variance") +
  xlab("Day") +
  #ggtitle("Metabolite concentration variances differ between days") +
  theme_bw()
ggsave(filename = "human_seperate_days_grouped_metab_var.png", path = out_dir, plot = metab_day_var_plot, width = 8, height = 3.5, units = "in")

##Human, variance of metab vars seperate by day
metab_day_var_long_df <- melt(metab_normal_day_var_df, id.vars = c("Day", "Survival"))
metab_day_var_long_df$group <- human_sepsis_legend$group[match(metab_day_var_long_df$variable, human_sepsis_legend[[1]])]
metab_day_var_sig <- melt(anova.car.ph.metab.variance.sig.contr)
metab_day_var_sig$Day <- tanova_day_set[metab_day_var_sig$value]
metab_day_var_sig$value <- sapply(metab_day_var_sig$L1, function(gr) max(subset(metab_day_var_long_df, group == gr)$value)) * 1.1
colnames(metab_day_var_sig)[2] <- "group"
metab_day_var_long_df <- subset(metab_day_var_long_df, !group %in% c("lysophosphatidylcholine", "sphingolipid", "sugar"))
metab_day_var_plot <- ggplot(data = metab_day_var_long_df, mapping = aes(x = factor(Day), y = value, fill = Survival)) +
  facet_wrap( ~ group, ncol = 4, nrow = 2, scales = "free_y") +
  geom_boxplot(outlier.size = 0.7) +
  geom_point(data = metab_day_var_sig, mapping = aes(y = value, x = factor(Day, levels = tanova_day_set)), shape = 8, size = 0.8, inherit.aes = FALSE) +
  ylab("variance") +
  xlab("Day") +
  #ggtitle("Metabolite concentration variances differ between days") +
  theme_bw()
ggsave(filename = "human_seperate_days_grouped_metab_var_PG.png", path = out_dir, plot = metab_day_var_plot, width = 8, height = 2, units = "in")

##Human, covariance cluster-heatmap, metabolites
x <- human_sepsis_data_normal_metab_cov
heatmaply(x = x[,-1:-6], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_metab_cov.png"), main = "Metabolite profile covariance has mainly patient clusters", key.title = "Cov", showticklabels = FALSE)
x <- na.omit(human_sepsis_data_normal_pheno_cov)
heatmaply(x = x[,-1:-6], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_pheno_cov.png"), main = "Phenomenological profile covariance has patient\n and survival clusters", key.title = "Cov", showticklabels = FALSE)
x <- na.omit(human_sepsis_data_normal_metab_pheno_cov)
heatmaply(x = x[,-1:-6], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_metab_and_pheno_cov.png"), main = "", key.title = "Cov", showticklabels = FALSE, margins = c(100,50,0,150), plot_method = "ggplot")
x <- na.omit(human_sepsis_data_normal_grouped_metab_pheno_cov)
heatmaply(x = x[,-1:-6], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_grouped_metab_and_pheno_cov.png"), main = "", key.title = "Cov", showticklabels = FALSE, margins = c(100,50,0,150), plot_method = "ggplot")
x <- na.omit(human_sepsis_data_normal_metab_pheno_cor)
heatmaply(x = x[,-1:-6], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_metab_and_pheno_cor.png"), main = "", key.title = "Cor", showticklabels = FALSE, margins = c(100,50,0,150), plot_method = "ggplot")
x <- na.omit(human_sepsis_data_normal_grouped_metab_pheno_cor)
heatmaply(x = x[,-1:-6], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_grouped_metab_and_pheno_cor.png"), main = "", key.title = "Cor", showticklabels = FALSE, margins = c(100,50,0,150), plot_method = "ggplot")
rm("x")

#COMMENT: switch to similarity for patient signature visualization
x <- as.matrix(na.omit(subset(human_sepsis_data_normal, TRUE, select = pheno_sel)))
x <- cbind(human_sepsis_data_normal[, 1:6], tcrossprod(x))
heatmaply(x = x[,-1:-6], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_metab_and_pheno_cov.png"), main = "", key.title = "Dist", showticklabels = FALSE, margins = c(100,50,0,150), plot_method = "ggplot")

##Human, cluster-heatmap of non-CAP and non-FP patients, coarse grouped metabolites; error because nothing es left
#heatmaply(x = subset(human_sepsis_data_normal_grouped, subset = human_sepsis_data_normal_grouped$`CAP / FP` %in% c("-"), select = c(-1,-2,-3,-5)))

##Human, cluster-heatmap of patient covariance matrix, coarse grouped metabolites, survival marked
x <- human_sepsis_data_normal_grouped_metab_cov
heatmaply(x = x[,-1:-6], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_grouped_metab_cov.png"))
x <- na.omit(human_sepsis_data_normal_grouped_pheno_cov)
heatmaply(x = x[,-1:-6], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_grouped_pheno_cov.png"))
rm("x")

##Human, cluster-heatmap of patient correlation matrix, ungrouped metabolites, survival marked
x <- human_sepsis_data_normal_metab_cor
heatmaply(x = x[,-1:-6], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_metab_cor.png"), margin = c(100,50,0,150), key.title = "Cor", showticklabels = FALSE)
x <- na.omit(human_sepsis_data_normal_pheno_cor)
heatmaply(x = x[,-1:-6], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_pheno_cor.png"), main = "Phenomenological profile correlation gives patient\n and survival clusters")
rm("x")

##Human, cluster-heatmap of sample distance matrix, ungrouped metabolites, then ungrouped pheno vars, survival marked
x <- human_sepsis_data_normal
heatmaply(x = as.matrix(dist(x[, metab_sel])), row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_metab_dist.png"), margin = c(100,50,0,150), key.title = "Euclidean\ndistance", showticklabels = FALSE)
heatmaply(x = as.matrix(dist(x[, pheno_sel])), row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_pheno_dist.png"), margin = c(100,50,0,150), key.title = "Euclidean\ndistance", showticklabels = FALSE)
rm(x)

##Human, cluster-heatmat of patient correlation matrix, coarse grouped metabolites, survival and CAP/FP marked
x <- human_sepsis_data_normal_grouped_metab_cor
heatmaply(x = x[,-1:-6], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_grouped_metab_cor.png"), margins = c(100, 50, 0, 150), showticklabels = F, k_row = 2, key.title = "Cor")
x <- na.omit(human_sepsis_data_normal_grouped_pheno_cor)
heatmaply(x = x[,-1:-6], row_side_colors = x[c("Survival", "CAP / FP")], col_side_colors = x["Patient"], file = paste0(out_dir, "human_normal_grouped_pheno_cor.png"), main = "Profiles of grouped phenomenological variables cluster nothing")
rm("x")

##Human, cluster-heatmap of patient covariance matrix, coarse grouped everything, survival-regardent
x <- human_sepsis_data_normal
for (group in setdiff(unique(coarse_group_list), "sugar")){
  col_sel <- which(coarse_group_list == group)
  h <- heatmaply(x = x[, col_sel + 6], row_side_colors = x[c("Survival", "CAP / FP")], key.title = "Normlized\nConcentrations", margins = c(100,80,0,250), k_row = 2, subplot_heights = c(0.08, 0.92))
  h$width <- 1600
  h$height <- 1200
  export(p = h, file = paste0(out_dir, "human_normal_conc_", group, ".png"))
}
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
heatmaply(x = x, file = paste0(out_dir, "human_normal_NS_grouped_conc_metab_cov.png"), key.title = "Cov", margins = c(200,200,0,150))
x <- human_sepsis_data_normal_S_grouped_conc_metab_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_S_grouped_conc_metab_cov.png"), key.title = "Cov", margins = c(200,200,0,150))
x <- human_sepsis_data_normal_NS_grouped_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_NS_grouped_conc_pheno_cov.png"), key.title = "Cov", margins = c(200,200,0,150))
x <- human_sepsis_data_normal_S_grouped_conc_pheno_cov
heatmaply(x = x, file = paste0(out_dir, "human_normal_S_grouped_conc_pheno_cov.png"), key.title = "Cov", margins = c(200,200,0,150))
rm("x")

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

##Human, mean metabolite concentrations over days, ordered by concentration, all groups
human_data_patient_group_mean_all_days <- subset(human_data, Day %in% tanova_day_set, c(1:6, metab_sel))
hdpgmads <- human_data_patient_group_mean_all_days$Survival
hdpgmads[human_data_patient_group_mean_all_days$`CAP / FP` == "-"] <- "Control"
human_data_patient_group_mean_all_days <- aggregate(x = human_data_patient_group_mean_all_days[, -1:-6], by = list(hdpgmads), FUN = mean)
human_data_patient_group_mean_all_days <- melt(human_data_patient_group_mean_all_days, id.vars = 1)

##Human, metabolite concentration time course, only metabolites with p-val < 0.05 at any day
###Prepare significant diff data
##Reduce to significantly different metabolites
human_sepsis_data_long_form_sig <- subset(human_sepsis_data_long_form, subset = as.character(variable) %in% sig_t_class)
human_sepsis_data_long_form_sig$variable <- factor(as.character(human_sepsis_data_long_form_sig$variable), levels = sig_t_class, ordered = TRUE)
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
  ggtitle("Metabolites significantly differing for survival by t-test at any time point") + 
  theme_bw()
ggsave(plot = h_time_course_sig_diff_plot, filename = "human_metab_time_course_t_sig_diff.png", path = out_dir, width = 14, height = 10, units = "in")

##Human, metab concentration time course, only metabolites significant in control vs sepsis and with control between S and NS
keep_set <- c(paste0("lysoPC a C2", c("6:0", "6:1", "8:0", "8:1")),
              paste0("PC aa C", c("30:2", "32:3", "34:2", "34:3", "36:2", "36:3", "36:4", "38:5", "38:6", "40:1", "40:2", "40:3", "40:6", "42:0", "42:1", "42:2")),
              paste0("PC ae C", c("32:2", "34:2", "34:3", "36:1", "36:2", "36:3", "38:0", "38:2", "38:3", "38:4", "38:5", "38:6", "44:5", "44:6")),
              grep(pattern = "SM.+", colnames(human_data), value = T))
hd <- human_data[, c(1:6, which(colnames(human_data) %in% keep_set))]
n_mets <- ncol(hd) - 6
ncols <- 5
p <- ggplot(data = subset(melt(hd, id.vars = 1:6), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Group, colour = Group)) + 
  facet_wrap(facets = ~ variable, ncol = ncols, nrow = ceiling(n_mets/ncols), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration, µM") +
  xlab("Day") +
  human_col_scale() +
  theme_bw()
ggsave(plot = p, filename = "human_metab_nonsig_hormesis.png", path = out_dir, width = 12/4*ncols, height = 0.3 + 1.5 * ceiling(n_mets/ncols), units = "in")
##Human, mean metabolite concentrations over days, ordered by concentration, all groups, reuse metabolite set from above
human_data_patient_group_mean_all_days <- subset(human_data, Day %in% tanova_day_set, c(1:6, which(colnames(human_data) %in% keep_set)))
hdpgmads <- human_data_patient_group_mean_all_days$Group
human_data_patient_group_mean_all_days <- aggregate(x = human_data_patient_group_mean_all_days[, -1:-6], by = list(hdpgmads), FUN = mean)
col_order <- order(colMaxs(as.matrix(human_data_patient_group_mean_all_days[, -1])), decreasing = TRUE)
human_data_patient_group_mean_all_days[, -1] <- human_data_patient_group_mean_all_days[, 1 + col_order]
cns <- colnames(human_data_patient_group_mean_all_days)[-1]
colnames(human_data_patient_group_mean_all_days)[-1] <- cns[col_order]
human_data_patient_group_mean_all_days <- melt(human_data_patient_group_mean_all_days, id.vars = 1)
human_data_patient_group_mean_all_days$Group <- human_data_patient_group_mean_all_days$Group.1
p <- ggplot(data = human_data_patient_group_mean_all_days, mapping = aes(y = value, x = variable, color = Group)) +
  geom_point() +
  scale_y_log10() +
  ylab("Mean concentration, µM") +
  xlab("Metabolite") +
  human_col_scale() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), panel.grid = element_line(colour = 0))
ggsave(filename = "human_metab_nonsig_hormesis_single_plot.png", path = out_dir, plot = p, width = 10, height = 7, units = "in")
##Human, mean metabolite concentrations over days, ordered by concentration, all groups, reuse metabolite set from above
human_data_patient_group_mean_all_days <- subset(human_data, Day %in% tanova_day_set, c(1:6, which(colnames(human_data) %in% keep_set)))
hdpgmads <- human_data_patient_group_mean_all_days$Group
hdpg_means <- aggregate(x = human_data_patient_group_mean_all_days[, -1:-6], by = list(hdpgmads), FUN = mean)
col_order <- order(colMaxs(as.matrix(hdpg_means[, -1])), decreasing = TRUE)
human_data_patient_group_mean_all_days[, -1:-6] <- human_data_patient_group_mean_all_days[, 6 + col_order]
cns <- colnames(human_data_patient_group_mean_all_days)[-1:-6]
colnames(human_data_patient_group_mean_all_days)[-1:-6] <- cns[col_order]
human_data_patient_group_mean_all_days <- melt(human_data_patient_group_mean_all_days, id.vars = 1:6)
p <- ggplot(data = human_data_patient_group_mean_all_days, mapping = aes(y = value, x = variable, color = Group)) +
  stat_summary(fun.data = "mean_sdl", position = position_dodge(width = 0.7), geom = "errorbar") +
  stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax = "mean", position = position_dodge(width = 0.7), geom = "errorbar") +
  scale_y_log10() +
  ylab("Mean concentration +/- 1 standard deviation, µM") +
  xlab("Metabolite") +
  human_col_scale() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), panel.grid = element_line(colour = 0))
ggsave(filename = "human_metab_nonsig_hormesis_single_plot_SD.png", path = out_dir, plot = p, width = 10, height = 7, units = "in")
###Same data, but this time one plot for each NS patient where NS samples are plotted as dots
for (pat in unique(subset(human_data_patient_group_mean_all_days, Group == "Septic-NS")$Patient)){
  pat_path <- subset(human_data_patient_group_mean_all_days, Patient == pat)
  pat_path$x <- match(pat_path$variable, cns[col_order])
  pat_path$x <- pat_path$x + rep(scale(seq_along(unique(pat_path$Day)))/3, times = length(unique(pat_path$variable)))
  p <- ggplot(data = subset(human_data_patient_group_mean_all_days, Survival != "NS"), mapping = aes(y = value, x = variable, color = Group)) +
    stat_summary(fun.data = "mean_sdl", position = position_dodge(width = 0.7), geom = "errorbar") +
    stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax = "mean", position = position_dodge(width = 0.7), geom = "errorbar") +
    geom_line(data = pat_path, mapping = aes(y = value, x = x, color = Group, group = variable), inherit.aes = FALSE) +
    scale_color_discrete(drop = FALSE) +
    scale_y_log10() +
    ylab("Mean concentration +/- 1 standard deviation, µM") +
    xlab("Metabolite") +
    human_col_scale() +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), panel.grid = element_line(colour = 0))
  ggsave(filename = paste0("human_metab_nonsig_hormesis_single_plot_SD_pat_", pat, ".png"), path = out_dir, plot = p, width = 10, height = 7, units = "in")
}

###Same as above but with all metabolites where the NS pat's concentration is outside the mean +/- 1 SD of the other groups
human_data_patient_group_mean_all_days <- subset(human_data, Day %in% tanova_day_set)
hdpgmads <- human_data_patient_group_mean_all_days$Group
hdpg_means <- aggregate(x = human_data_patient_group_mean_all_days[, -1:-6], by = list(hdpgmads), FUN = mean)
hdpg_sds <- aggregate(x = human_data_patient_group_mean_all_days[, -1:-6], by = list(hdpgmads), FUN = sd)
col_order <- order(colMaxs(as.matrix(hdpg_means[, -1])), decreasing = TRUE)
human_data_patient_group_mean_all_days[, -1:-6] <- human_data_patient_group_mean_all_days[, 6 + col_order]
hdpg_means <- hdpg_means[c(1, 1 + col_order)] #watch the index, else won't reorder column names
hdpg_sds <- hdpg_sds[c(1, 1 + col_order)]
hdpg_MNminusSD <- hdpg_means
hdpg_MNminusSD[, -1] <- hdpg_MNminusSD[, -1] - 1.5 * hdpg_sds[, -1]
hdpg_MNplusSD <- hdpg_means
hdpg_MNplusSD[, -1] <- hdpg_MNplusSD[, -1] + 1.5 * hdpg_sds[, -1]
cns <- colnames(human_data_patient_group_mean_all_days)[-1:-6]
colnames(human_data_patient_group_mean_all_days)[-1:-6] <- cns[col_order]
human_data_patient_group_mean_all_days <- melt(human_data_patient_group_mean_all_days, id.vars = 1:6)
metabolite_groups <- human_sepsis_legend$group[match(human_data_patient_group_mean_all_days$variable, human_sepsis_legend[, 1])]
metabolite_groups[metabolite_groups %in% coarse_group_list[pheno_sel - 6]] <- "Clinical parameter"
human_data_patient_group_mean_all_days$metabolite_group <- metabolite_groups
var_keep <- list()
for (pat in unique(subset(human_data_patient_group_mean_all_days, Survival == "NS")$Patient)){
  cntrl_and_s <- human_data_patient_group_mean_all_days
  pat_path <- subset(cntrl_and_s, Patient == pat)
  long_short_map <- match(pat_path$variable, colnames(hdpg_MNminusSD))
  var_keep[[pat]] <- unique(pat_path$variable[pat_path$value < colMins(as.matrix(hdpg_MNminusSD[c(2, 4), long_short_map])) | pat_path$value > colMaxs(as.matrix(hdpg_MNplusSD[c(2, 4), long_short_map]))])
  var_keep[[pat]] <- setdiff(as.character(var_keep[[pat]]), sig.anova.car.s.class)
}
var_keep_count <- table(unlist(var_keep))
var_keep_count_df <- as.data.frame(var_keep_count)
var_keep_count_df <- var_keep_count_df[var_keep_count_df$Freq >= 4, ]
var_keep_count_df$color <- grey_pal()(length(unique(var_keep_count_df$Freq)) - 3)[sapply(var_keep_count_df$Freq - 3, max, 1)]
var_keep_union <- Reduce("union", var_keep)
for (pat in unique(subset(human_data_patient_group_mean_all_days, Survival == "NS")$Patient)){
  cntrl_and_s <- human_data_patient_group_mean_all_days
  pat_path <- subset(cntrl_and_s, Patient == pat)
  long_short_map <- match(pat_path$variable, colnames(hdpg_MNminusSD))
  #var_keep <- unique(pat_path$variable[pat_path$value < colMins(as.matrix(hdpg_MNminusSD[c(1, 3), long_short_map])) | pat_path$value > colMaxs(as.matrix(hdpg_MNplusSD[c(1, 3), long_short_map]))])
  cntrl_and_s <- subset(cntrl_and_s, variable %in% var_keep_union)
  pat_path <- subset(pat_path, variable %in% var_keep_union)
  cns <- unique(cntrl_and_s$variable)
  pat_path$x <- match(pat_path$variable, cns)
  pat_path$x <- pat_path$x + rep(scale(seq_along(unique(pat_path$Day)))/3, times = length(unique(pat_path$variable)))
  p <- ggplot(data = subset(cntrl_and_s, Survival != "NS"), mapping = aes(y = value, x = variable, color = Group)) +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), position = position_dodge(width = 0.7), geom = "errorbar") +
    stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax = "mean", position = position_dodge(width = 0.7), geom = "errorbar") +
    geom_line(data = subset(pat_path, variable %in% var_keep_union), mapping = aes(y = value, x = x, color = Group, group = variable), inherit.aes = FALSE) +
    geom_tile(mapping = aes(x = variable, y = 1e-4, fill = metabolite_group, width = 1, height = 0.5), data = pat_path, inherit.aes = FALSE) +
    geom_tile(mapping = aes(x = Var1, y = 3.5e-5), width = 1, height = 0.5, fill = var_keep_count_df$color, data = var_keep_count_df, inherit.aes = FALSE) +
    geom_point(mapping = aes(x = Var1, y = 3.5e-5, shape = factor(Freq)), color = var_keep_count_df$color, data = var_keep_count_df, inherit.aes = FALSE) +
    guides(color = guide_legend(order = 1), 
           fill = guide_legend(title = "Metabolite Group", order = 2),
           shape = guide_legend(title = "#NS pats with deviation", override.aes = list(shape = 15, size = 6, color = sort(unique(var_keep_count_df$color))), order = 3)) +
    scale_color_discrete(drop = FALSE) +
    scale_y_log10(expand = c(0,0)) +
    ylab("Mean concentration +/- 1 standard deviations, µM") +
    xlab("Metabolite") +
    human_col_scale() +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6), panel.grid = element_line(colour = 0))
  #TODO: fix
  ggsave(filename = paste0("human_metab_nonsig_single_plot_SD_pat_", pat, ".png"), path = out_dir, plot = p, width = 16, height = 7, units = "in")
}
cntrl_and_s <- human_data_patient_group_mean_all_days
pat_path <- subset(cntrl_and_s, Survival == "NS")
long_short_map <- match(pat_path$variable, colnames(hdpg_MNminusSD))
#var_keep <- unique(pat_path$variable[pat_path$value < colMins(as.matrix(hdpg_MNminusSD[c(2, 4), long_short_map])) | pat_path$value > colMaxs(as.matrix(hdpg_MNplusSD[c(2, 4), long_short_map]))])
var_keep_all_count_df <- var_keep_count_df
var_keep_count_df <- subset(var_keep_count_df, Freq >= 4)
var_keep_union <- intersect(var_keep_union, var_keep_count_df$Var1)
cntrl_and_s <- subset(cntrl_and_s, variable %in% var_keep_union)
pat_path <- subset(pat_path, variable %in% var_keep_union)
cns <- unique(cntrl_and_s$variable)
pat_path$x <- match(pat_path$variable, cns)
pat_path$x <- pat_path$x + (scale(seq_along(unique(pat_path$Day)))/4)[unlist(lapply(rle(pat_path$Patient)[["lengths"]], function(to) 1:to))]
p <- ggplot(data = subset(cntrl_and_s, Survival != "NS"), mapping = aes(y = value, x = variable, color = Group)) +
  #geom_ribbon(data = corridor, mapping = aes(ymax = ymax, ymin = ymin, x = variable), inherit.aes = FALSE) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), position = position_dodge(width = 0.7), geom = "errorbar") +
  stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax = "mean", position = position_dodge(width = 0.7), geom = "errorbar") +
  geom_line(data = subset(pat_path, variable %in% var_keep_union), mapping = aes(y = value, x = x, color = Group, group = interaction(variable, Patient), linetype = factor(Patient)), inherit.aes = FALSE) +
  geom_tile(mapping = aes(x = variable, y = 1e-4, fill = metabolite_group), width = 1, height = 0.5, data = pat_path, inherit.aes = FALSE) +
  geom_point(mapping = aes(x = Var1, y = 3.5e-5, shape = factor(Freq)), data = var_keep_count_df, inherit.aes = FALSE) +
  geom_tile(mapping = aes(x = Var1, y = 3.5e-5), width = 1, height = 0.5, fill = var_keep_count_df$color, data = var_keep_count_df, inherit.aes = FALSE) +
  guides(color = guide_legend(order = 1), 
         fill = guide_legend(title = "Metabolite Group", order = 2),
         shape = guide_legend(title = "#NS pats with deviation", override.aes = list(shape = 15, size = 6, colour = sort(unique(var_keep_count_df$color))), order = 3),
         linetype = guide_legend(title = "NS Patient", override.aes = list(color = hue_pal()(4)[c(4, 1)[1 + (subset(cntrl_and_s, !duplicated(Patient) & Survival == "NS", "Group")[[1]] == "Septic-NS")]]), order = 4)) +
  scale_color_discrete(drop = FALSE) +
  scale_y_log10(expand = c(0,0)) +
  ylab("Mean concentration +/- 1 standard deviations, µM") +
  xlab("Metabolite") +
  human_col_scale() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6), panel.grid = element_line(colour = 0))
#TODO: fix
ggsave(filename = paste0("human_metab_nonsig_single_plot_SD_all_pats.png"), path = out_dir, plot = p, width = 16, height = 9.5, units = "in")

###Generalize the corridor principle to a classification scheme, use corridor 
corridor <- as.data.frame(t(sapply(unique(cntrl_and_s$variable), 
                                   function(metab) Hmisc::smean.sd(subset(cntrl_and_s, variable == metab & Survival == "S", "value")[[1]]))))
corridor$variable <- unique(cntrl_and_s$variable)
for (n in seq_along(corridor))
  corridor[[n]] <- unlist(corridor[[n]])
pat_dev_score <- subset(human_data, Day %in% 0:3, select = c(1:6, which(colnames(human_data) %in% corridor$variable)))
#pat_dev_score <- subset(human_data, Day %in% 0:3, select = setdiff(which(!(colnames(human_data) %in% sig.anova.car.s.class)), pheno_sel))
pat_dev_mean <- colMeans(pat_dev_score[pat_dev_score$Survival == "S", -1:-6])
pat_dev_sd <- colSds(as.matrix(pat_dev_score[pat_dev_score$Survival == "S", -1:-6]))
sdmul <- 5.5
udev <- pat_dev_score[, -1:-6] > matrix(pat_dev_mean + sdmul * pat_dev_sd, ncol = ncol(pat_dev_score) - 6, nrow = nrow(pat_dev_score), byrow = TRUE)
ldev <- pat_dev_score[, -1:-6] < matrix(pat_dev_mean - sdmul * pat_dev_sd, ncol = ncol(pat_dev_score) - 6, nrow = nrow(pat_dev_score), byrow = TRUE)
sdev <- aggregate(udev | ldev, by = list(Patient = pat_dev_score$Patient), FUN = max)
dev_score <- data.frame(Patient = sdev$Patient, score = rowSums(sdev[, -1]))
dev_score$Survival <- pat_dev_score$Survival[match(dev_score$Patient, pat_dev_score$Patient)]
dev_score$Group <- pat_dev_score$Group[match(dev_score$Patient, pat_dev_score$Patient)]
p <- ggplot(data = dev_score, mapping = aes(fill = Group, x = score)) +
  geom_histogram(position = position_stack(), bins = max(dev_score$score) + 1) +
  human_col_scale(aesthetics = "fill") +
  ylab("Number of Patients") +
  xlab("Number of metabolites outside of the safe corridor at Days 0-3") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(plot = p, filename = "generalized_safe_corridor_SD5.5.png", path = out_dir, width = 6, height = 3, units = "in")
print(paste0("Contribution of high deviation and low deviation to score is ", sum(udev), " and ", sum(ldev), " respectively."))
w <- which.xy(udev) # tell me which variables make a difference
mtab <- sort(table(w[, 2]), decreasing = TRUE)
names(mtab) <- colnames(pat_dev_score)[-1:-6][as.numeric(names(mtab))]
mtab

###plot metab variance, only nonsig, nondeviation metabs
##Human, variance of ungrouped metab vars, all days
met_nor_day_var_df <- na.omit(human_sepsis_data_conc_var)
met_nor_day_var_df[, -1:-6] <- scale(met_nor_day_var_df[, -1:-6], center = FALSE, scale = colMeans(met_nor_day_var_df[, -1:-6]))
colnames(met_nor_day_var_df)[-1:-2] <- colnames(human_sepsis_data)[-1:-6]
met_nor_day_var_df <- met_nor_day_var_df[, 1:which(colnames(met_nor_day_var_df) == "H1")]
metab_day_var_long_df <- melt(met_nor_day_var_df, id.vars = c("Survival", "Patient"))
metab_day_var_long_df$group <- human_sepsis_legend$group[match(metab_day_var_long_df$variable, human_sepsis_legend[[1]])]
#metab_day_var_long_df <- subset(metab_day_var_long_df, !(as.character(variable) %in% sig.anova.car.s.class))
#metab_day_var_long_df <- subset(metab_day_var_long_df, !(as.character(variable) %in% var_keep_count_df$Var1))
gns <- colnames(met_nor_day_var_df)[-1:-6][colMedians(as.matrix(subset(met_nor_day_var_df, Survival == "NS", -1:-6))) > colMedians(as.matrix(subset(met_nor_day_var_df, Survival == "S", -1:-6)))]
metab_sig_sel <- mclapply(unique(metab_day_var_long_df$variable),
                          function(x){
                            d <- subset(metab_day_var_long_df, variable == x)
                            # num_S <- sum(d$Survival == "S")
                            # num_NS <- sum(d$Survival == "NS")
                            # s_idx <- which(d$Survival == "S")
                            # ns_idx <- which(d$Survival == "NS")
                            # n_take <- min(num_S, num_NS)
                            # ad_res <- list()
                            # for (n in 1:50){ #bootstrap
                            #   ds <- d[c(s_idx[sample(x = num_S, size = n_take)], ns_idx[sample(x = num_NS, size = n_take)]), ]
                            #   ad_res[[n]] <- adonis(formula = value ~ Survival, data = ds, parallel = 1, method = "euclidean")
                            # }
                            # return(ad_res)
                            t.test(x = d$value[d$Survival == "NS"], y = d$value[d$Survival == "S"])
                          })
# metab_sig_sel <- metab_sig_sel_a
# metab_sig_sel <- sapply(lapply(metab_sig_sel, sapply, function(e) e$aov.tab[["Pr(>F)"]][1]), mean)
metab_sig <- sapply(metab_sig_sel, `[[`, "p.value")
metab_sig_sel <- metab_sig < 0.05
metab_sig_sel_name <- unique(metab_day_var_long_df$variable)[metab_sig_sel]
#metab_day_var_long_df <- subset(metab_day_var_long_df, as.character(variable) %in% gns & as.character(variable) %in% unique(variable)[metab_sig_sel])
metab_sig_rect <- data.frame(x = unique(metab_day_var_long_df$variable)[metab_sig_sel], y = max(metab_day_var_long_df$value) / 2, width = 1, height = max(metab_day_var_long_df$value))
metab_day_var_long_df <- subset(metab_day_var_long_df, as.character(variable) %in% gns)
metab_sig_rect <- subset(metab_sig_rect, x %in% gns)
metab_day_var_plot <- ggplot(data = metab_day_var_long_df, mapping = aes(x = variable, y = value, fill = Survival)) + #a.k.a. Figure 4 in Sepsis variance manuscript
  geom_tile(data = metab_sig_rect, mapping = aes(x = x, y = y, width = width, height = height), fill = "grey", size = 1, inherit.aes = FALSE) +
  geom_boxplot(outlier.size = 0.5) +
  #scale_y_log10() +
  ylab("Patient-wise variance relative to metabolite mean") +
  xlab("Metabolite") +
  human_col_scale(name = "Survival", levels = c("NS", "", "S", "", ""), aesthetics = "fill") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid = element_blank())
ggsave(filename = "human_all_days_ungrouped_metab_var_to_metab_mean.png", path = out_dir, plot = metab_day_var_plot, width = 10, height = 5.5, units = "in")
##Human, mean metabolite and pheno var concentrations over days, ordered by concentration, all groups, seperate for each metabolite group
human_data_patient_group_mean_all_days <- subset(human_data, Day %in% tanova_day_set)
hdpgmads <- human_data_patient_group_mean_all_days$Group
hdpg_means <- aggregate(x = human_data_patient_group_mean_all_days[, -1:-6], by = list(hdpgmads), FUN = mean)
col_order <- order(colMaxs(as.matrix(hdpg_means[, -1])), decreasing = TRUE)
human_data_patient_group_mean_all_days[, -1:-6] <- human_data_patient_group_mean_all_days[, 6 + col_order]
cns <- colnames(human_data_patient_group_mean_all_days)[-1:-6]
colnames(human_data_patient_group_mean_all_days)[-1:-6] <- cns[col_order]
human_data_patient_group_mean_all_days <- melt(human_data_patient_group_mean_all_days, id.vars = 1:6)
metabolite_groups <- human_sepsis_legend$group[match(human_data_patient_group_mean_all_days$variable, human_sepsis_legend[, 1])]
metabolite_groups[metabolite_groups %in% coarse_group_list[pheno_sel - 6]] <- "Clinical parameter"
human_data_patient_group_mean_all_days$metabolite_group <- metabolite_groups
for (met_group in unique(metabolite_groups)){
  plot_data <- subset(human_data_patient_group_mean_all_days, metabolite_group == met_group)
  for (pat in unique(subset(plot_data, Survival == "NS")$Patient)){
    pat_path <- subset(plot_data, Patient == pat)
    pat_path$x <- match(pat_path$variable, unique(plot_data$variable))
    pat_path$x <- pat_path$x + rep(scale(seq_along(unique(pat_path$Day)))/3, times = length(unique(pat_path$variable)))
    p <- ggplot(data = subset(plot_data, Survival != "NS"), mapping = aes(y = value, x = variable, color = Group)) +
      stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), position = position_dodge(width = 0.7), geom = "errorbar") +
      stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax = "mean", position = position_dodge(width = 0.7), geom = "errorbar") +
      geom_line(data = pat_path, mapping = aes(y = value, x = x, color = Group, group = variable), inherit.aes = FALSE) +
      scale_y_log10() +
      ylab("Mean concentration +/- 1 standard deviation, µM") +
      xlab("Metabolite") +
      human_col_scale() +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5), panel.grid = element_line(colour = 0))
    ggsave(filename = paste0("human_metab_nonsig_", met_group, "_single_plot_SD_pat_", pat, ".png"), path = out_dir, plot = p, width = 10, height = 7, units = "in")
  }
}

##Human, metab concentration time course, only metabolites with p-val > 0.05 and FDR > 0.95 in Survival or Day:Survival from type III repeated measures ANOVA
n_mets <- length(insig.anova.car.s.class)
p <- ggplot(data = subset(melt(human_sepsis_data[, c(1:6, which(colnames(human_sepsis_data) %in% insig.anova.car.s.class))], id.vars = 1:6), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Group, colour = Group)) + 
  facet_wrap(facets = ~ variable, ncol = 4, nrow = ceiling(n_mets/4), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration, µM") +
  xlab("Day") +
  human_col_scale() +
  theme_bw()
plot(p)

##Human, metab concentration time course, only metabolites with p-val > 0.05 and FDR > 0.95 in Survival or Day:Survival from type III repeated measures ANOVA
n_mets <- length(insig.anova.car.c.class)
p <- ggplot(data = subset(melt(human_data[, c(1:6, which(colnames(human_data) %in% insig.anova.car.c.class))], id.vars = 1:6), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Group, colour = Group)) + 
  facet_wrap(facets = ~ variable, ncol = 4, nrow = ceiling(n_mets/4), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration, µM") +
  xlab("Day") +
  human_col_scale() +
  theme_bw()
plot(p)

##Human, pheno var concentration time course, only those with p-val > 0.05 and FDR > 0.95 in Survival or Day:Survival from type III repeated measures ANOVA
n_mets <- length(insig.anova.car.s.pheno.pre.class)
p <- ggplot(data = subset(melt(human_sepsis_data[, c(1:6, which(colnames(human_sepsis_data) %in% insig.anova.car.s.pheno.pre.class))], id.vars = 1:6), Day %in% 0:3), mapping = aes(x = Day, y = value, group = Group, colour = Group)) + 
  facet_wrap(facets = ~ variable, ncol = 4, nrow = ceiling(n_mets/4), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.2)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  ylab("Concentration, µM") +
  xlab("Day") +
  human_col_scale() +
  theme_bw()
plot(p)
