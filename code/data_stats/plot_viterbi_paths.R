library(ggplot2)
library(data.table)
library(Hmisc)
library(heatmaply)
library(matrixStats)
library(plotly)
library(htmlwidgets)


#Import central functions
source("../function_definitions.R")

#Set path
out_dir <- "../../results/data_stats/"

#Read human metabolic data
human_sepsis_data <- get_human_sepsis_data()

#Read viterbi paths from file
viterbi_paths <- fread(input = '../HMM/tram_viterbi_paths_all_samples_loopjumpHMM.csv', sep = "\t", header = TRUE, data.table = FALSE)
viterbi_paths <- viterbi_paths[order(viterbi_paths$Case),]
#Attach metabolites to viterbi paths using the case column (= patient number)
split_start <- which(colnames(human_sepsis_data) == "Urea")
pheno_sel <- split_start:ncol(human_sepsis_data)
metab_sel <- 6:ncol(human_sepsis_data)
viterbi_paths = cbind(viterbi_paths, max_norm(human_sepsis_data[,metab_sel]))

#Transform paths to long format table
viterbi_paths_long <- melt(viterbi_paths, id.vars = c("Status", "State", "Case"))

#Calculate significant differences between survival states
hmm_h_sig_diff <- human_sig_diffs_along_days(data = viterbi_paths, corr_fdr = TRUE, time_var = "State", case_var = "Case", status_var = "Status", descr_till_col = 3)
hmm_h_sig_t_diff <- hmm_h_sig_diff$day_sig_t_diff

#Get significantly different metabolites
hmm_h_sig_t_class <- na.omit(colnames(hmm_h_sig_t_diff[,-1])[colAnys(hmm_h_sig_t_diff[, -1] <= 0.05)])

#Reduce metabolite set to significantly different ones
viterbi_paths_long <- na.omit(subset(viterbi_paths_long, variable %in% hmm_h_sig_t_class))

#Get position of significantly differnt concentrations
hmm_h_viterbi_path_sig_t_diff_pos_long <- get_sig_var_pos(diff_data = hmm_h_sig_t_diff, alpha = 0.05, time_var = "State")
hmm_h_viterbi_path_sigs <- subset(hmm_h_viterbi_path_sig_t_diff_pos_long, variable %in% viterbi_paths_long$variable)
m_val <- max(viterbi_paths_long$value, na.rm = TRUE)
hmm_h_viterbi_path_sigs$value <- m_val + 0.5

#Plot viterbi paths
viterbi_p <- ggplot(viterbi_paths_long, aes(x = State, y = value, group = Status, color = Status)) +
  facet_wrap(facets = ~ variable, ncol = 12, nrow = ceiling(length(unique(viterbi_paths_long$variable))/12)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", width = 0.5, geom = "errorbar") +
  #stat_summary(fun.data = "mean_sdl", width = 0.5, geom = "errorbar") + 
  geom_point(data = hmm_h_viterbi_path_sigs, mapping = aes(x = State, y = value), shape = 8, inherit.aes = FALSE, size = 0.8) +
  scale_x_discrete(limits = 2:4) + 
  ylim(c(0, m_val+0.6)) +
  ylab("standardized concentration") +
  ggtitle("Aligned time courses differ in more metabolites than unaligned time courses") +
  theme_bw()
ggsave(plot = viterbi_p, filename = paste0(out_dir, "human_HMM_aligned_time_course.png"), width = 14, height = 14, units = "in")

p <- ggplot(human_sepsis_data, aes_string(x = "Day", y = "`PC aa C36:0`", group = "Survival", color = "Survival"))+
  geom_boxplot() +
  theme_bw()
X11();plot(p)

plot_list <- list()
for (case in unique(viterbi_paths$Case)){
  x <- subset(viterbi_paths, subset = Case == case)
  x <- cbind(x[,1:3], as.matrix(dist(x[,-1:-3])))
  if (nrow(x) > 1){
    plot_list[[case]] <- heatmaply(x = x[,-1:-3], row_side_colors = x["State"], key.title = "Euclidean\ndistance", main = paste0("Case ", case, ": ", x[1,1]), showticklabels = FALSE, plot_method = "ggplot")
  }
}
pl_not_null <- which(!unlist(lapply(X = plot_list, FUN = is.null)))
plot_list <- plot_list[pl_not_null]
us <- viterbi_paths$Status[match(unique(viterbi_paths$Case), viterbi_paths$Case)]
us <- us[pl_not_null]
sp_S <- plotly::subplot(plot_list[us == "S"], shareX = F, shareY = F, heights = rep(1/6, 6), nrows = 6)
sp_NS <- plotly::subplot(plot_list[us == "NS"], shareX = F, shareY = F, heights = rep(1/6, 6), nrows = 6)
saveWidget(sp_S, file = "HMM_states_S.html")
saveWidget(sp_NS, file = "HMM_states_NS.html")