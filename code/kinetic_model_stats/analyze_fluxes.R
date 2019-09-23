library(purrr)
library(ggplot2)
library(stringi)
library(reshape2)
library(data.table)
library(matrixStats)
library(CoRC)
library(tictoc)
library(parallel)
library(cowplot)

source("../function_definitions.R")

models_dir <- "../../results/kinetic_modeling/"
out_dir <- "../../results/kinetic_modeling/"

par_est_res <- readRDS(file = paste0(out_dir, "par_est_image.RData"))

par_d0 <- lapply(lapply(lapply(par_est_res, `[[`, "hj_res_d0"), `[[`, "parameters"), subset, subset = stri_startswith_fixed(str = parameter, pattern = "Values["), select = c("parameter", "value", "lower_bound", "upper_bound"))
par_d0 <- lapply(par_d0, subset, subset = duplicated(parameter)) #parameter table contains both the original bounds (first occurence) and reset bounds, we filter out the original bounds
pd0 <- Reduce("rbind", par_d0)
pd0$Day <- 0
par_d1 <- lapply(lapply(lapply(par_est_res, `[[`, "hj_res_d1"), `[[`, "parameters"), subset, subset = stri_startswith_fixed(str = parameter, pattern = "Values["), select = c("parameter", "value", "lower_bound", "upper_bound"))
par_d1 <- lapply(par_d1, subset, subset = duplicated(parameter)) #parameter table contains both the original bounds (first occurence) and reset bounds, we filter out the original bounds
pd1 <- Reduce("rbind", par_d1)
pd1$Day <- 1
par_d2 <- lapply(lapply(lapply(par_est_res, `[[`, "hj_res_d2"), `[[`, "parameters"), subset, subset = stri_startswith_fixed(str = parameter, pattern = "Values["), select = c("parameter", "value", "lower_bound", "upper_bound"))
par_d2 <- lapply(par_d2, subset, subset = duplicated(parameter)) #parameter table contains both the original bounds (first occurence) and reset bounds, we filter out the original bounds
par_d2 <- lapply(par_d2, subset, subset = !duplicated(parameter)) #at day 2 we get Vfcact two times in the results, after filtering out the original bound entries; keep the last one.
pd2 <- Reduce("rbind", par_d2)
pd2$Day <- 2

#Plot enzyme
pd <- rbind(pd0, pd1, pd2)
pd$Run <- rep(rep(1:length(par_est_res), each = length(unique(pd0$parameter))), times = 3)
pd$Group <- c("Septic-S", "Septic-NS")[1 + (grepl(pattern = "merge", x = pd$parameter))] #Values with [merge] belong to Nonsurvivors
pd$parameter <- stri_replace(str = pd$parameter, replacement = "", fixed = "[merge]")
pd$parameter <- stri_extract(str = pd$parameter, regex = "\\[.+\\]")
n_par <- length(unique(pd$parameter))
pd_mn <- lapply(unique(pd$parameter), function(par){ ss <- subset(x = pd, parameter == par); ss$value <- ss$value / mean(ss$value); return(ss) })
pd_mn <- Reduce("rbind", pd_mn)
pd$parameter <- sub(pattern = "[", replacement = "", pd$parameter, fixed = TRUE)
pd$parameter <- sub(pattern = "]", replacement = "", pd$parameter, fixed = TRUE)
pd_sub <- subset(pd, parameter %in% c("Vcpt2", "Vmcad", "Vmschad", "Vcrot"))
pd_sub$parameter <- substring(pd_sub$parameter, 2)
pd_sub$parameter <- toupper(pd_sub$parameter)
pd_sub$parameter <- factor(pd_sub$parameter, levels = unique(pd_sub$parameter)[c(1, 3, 2, 4)])
pd_bounds <- subset(pd_sub, select = c("parameter", "lower_bound", "upper_bound"))
pd_bounds <- melt(pd_bounds, id.vars = c("parameter"))
pd_bounds$Group <- "Upper/Lower bound"
rsize <- rel(4)
p_vmax <- ggplot(data = pd_sub, mapping = aes(x = Day, y = value, color = Group)) + 
  facet_wrap(facets = ~ parameter, ncol = 2, nrow = 2, scale = "free_y") +
  geom_hline(data = pd_bounds, mapping = aes(yintercept = value, color = Group), size = rel(1)) +
  stat_summary(data = pd_sub, geom = "line", fun.data = mean_se, size = rel(1)) +
  stat_summary(data = pd_sub, geom = "errorbar", fun.data = mean_se, size = rel(1), width = 0.5) +
  scale_y_continuous(trans = "log2", breaks = c(10^(-5:2), 0.5 * 10^(-5:2), 0.25 * 10^(-5:2))) +
  scale_x_continuous(breaks = 0:2) +
  ylab("Concentration mean ± SEM\nµM") +
  human_col_scale(levels = c("Septic-NS", "", "Septic-S", "", "Upper/Lower bound"), black_color = "black") +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.direction = "horizontal", legend.position = "bottom", text = element_text(size = rsize), strip.text = element_text(size = rel(3)), legend.text = element_text(size = rsize, margin = margin(r = 10)))
ggsave(filename = "kin_mitomod_Vmax_S_vs_NS_repeated_subset.png", plot = p_vmax, path = out_dir, width = 9, height = 5, units = "in")
ggsave(filename = "kin_mitomod_Vmax_S_vs_NS_repeated_subset.svg", plot = p_vmax, path = out_dir, width = 9, height = 5, units = "in")

p_model_legend <- get_legend(plot = p_vmax)

p_vmax <- p_vmax + 
  guides(color = "none")

models_list <- list.files(path = models_dir, pattern = ".{+}\\.cps", full.names = TRUE)
file_list <- list.files(path = models_dir, pattern = ".{+}\\.cps", full.names = FALSE)

get_steady_state_stats <- function(model_path){
  model <- loadModel(path = model_path)
  ss_res <- runSteadyState(model = model)
  unloadModel(model)
  return(ss_res)
}

tic()
cl <- makeCluster(min(detectCores() - 1, 100), outfile = "")
prep_res <- clusterCall(cl = cl, fun = eval, quote({library(CoRC)}), env = .GlobalEnv)
stst_est_res <- parLapplyLB(cl = cl, X = models_list, fun = get_steady_state_stats)
stopCluster(cl)
toc()

get_global_parameters <- function(model_path){
  model <- loadModel(path = model_path)
  mpar <- getGlobalQuantities(model = model)
  unloadModel(model)
  return(mpar)
}

#tic()
#cl <- makeCluster(min(detectCores() - 1, 100), outfile = "")
#prep_res <- clusterCall(cl = cl, fun = eval, quote({library(CoRC)}), env = .GlobalEnv)
#mpar_res <- parLapplyLB(cl = cl, X = models_list, fun = get_global_parameters)
#stopCluster(cl)
#toc()

stability <- sapply(stst_est_res, function(ssres) all(Re(eigen(ssres$jacobian_complete)$values) < 0))
run_list <- sapply(lapply(lapply(file_list, strsplit, split = "_", fixed = TRUE), `[[`, 1), `[`, 5)
day_list <- substr(sapply(lapply(lapply(file_list, strsplit, split = "_", fixed = TRUE), `[[`, 1), `[`, 6), 1, 4)

concs <- as.data.frame(sapply(stst_est_res, function(ssres) ssres$species$concentration))
rownames(concs) <- stst_est_res[[1]]$species$concentration
concs <- transpose(concs)
colnames(concs) <- stst_est_res[[1]]$species$key
concs$Day <- day_list
concs$Run <- run_list
mconcs <- melt(concs, id.vars = c("Day", "Run"))
mconcs$Group <- "Septic-S"
mconcs$Group[grepl(pattern = "[merge]", x = mconcs$variable, fixed = TRUE)] <- "Septic-NS"
mconcs$variable <- sub(pattern = "[merge]", replacement = "", x = mconcs$variable, fixed = TRUE)
mconcs$variable <- sub(pattern = "{VMAT}", replacement = "", x = mconcs$variable, fixed = TRUE)
mconcs$variable <- sub(pattern = "{VCYT}", replacement = "", x = mconcs$variable, fixed = TRUE)

p <- ggplot(data = mconcs, mapping = aes(x = Day, y = value, color = Group)) + 
  facet_wrap(facets = ~ variable, ncol = 4, nrow = 15, scale = "free_y") +
  stat_summary(geom = "line", fun.data = mean_se) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(width = 0.2)) +
  ylab("Concentration mean +/- SEM, µmol/l") +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.direction = "horizontal", legend.position = "bottom")
ggsave(filename = "kin_mitomod_conc_S_vs_NS_repeated.png", plot = p, path = out_dir, width = 9, height = 11, units = "in")
ggsave(filename = "kin_mitomod_conc_S_vs_NS_repeated.svg", plot = p, path = out_dir, width = 9, height = 11, units = "in")

fluxes <- as.data.frame(sapply(stst_est_res, function(ssres) ssres$reactions$flux))
rownames(fluxes) <- stst_est_res[[1]]$reactions$name
fluxes <- transpose(fluxes)
colnames(fluxes) <- stst_est_res[[1]]$reactions$name
fluxes$Day <- day_list
fluxes$Run <- run_list
mfluxes <- melt(fluxes, id.vars = c("Day", "Run"))
mfluxes$Group <- "Septic-S"
mfluxes$Group[grepl(pattern = "[merge]", x = mfluxes$variable, fixed = TRUE)] <- "Septic-NS"
mfluxes$variable <- sub(pattern = "[merge]", replacement = "", x = mfluxes$variable, fixed = TRUE)
mfluxes$Day <- substring(mfluxes$Day, 4)
mfluxes$value[mfluxes$value < 0] <- 0

p <- ggplot(data = mfluxes, mapping = aes(x = Day, y = value, color = Group)) + 
  facet_wrap(facets = ~ variable, ncol = 6, nrow = 10, scale = "free_y") +
  stat_summary(geom = "line", fun.data = mean_se) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(width = 0.2)) +
  scale_y_continuous(labels = function(x) scales::scientific(x, digits = 2)) +
 # scale_y_log10() +
  ylab("Flux mean +/- SEM, µmol/min") +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.direction = "horizontal", legend.position = "bottom")
ggsave(filename = "kin_mitomod_flux_S_vs_NS_repeated.png", plot = p, path = out_dir, width = 9, height = 11, units = "in")
ggsave(filename = "kin_mitomod_flux_S_vs_NS_repeated.svg", plot = p, path = out_dir, width = 9, height = 11, units = "in")

rsize <- rel(4)
mfluxes_sink <- subset(mfluxes, variable %in% c("vnadhsink", "vfadhsink", "vacesink"))
mfluxes_sink$variable[mfluxes_sink$variable == "vacesink"] <- "Acetyl-CoA"
mfluxes_sink$variable[mfluxes_sink$variable == "vfadhsink"] <- "FADH2"
mfluxes_sink$variable[mfluxes_sink$variable == "vnadhsink"] <- "NADH+H"
p_sink_fluxes <- ggplot(data = mfluxes_sink, mapping = aes(x = Day, y = value, color = Group)) + 
  facet_wrap(facets = ~ variable, ncol = 3, nrow = 1, scale = "free_y") +
  stat_summary(geom = "line", fun.data = mean_se, size = rel(1)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(width = 0.2), size = rel(1)) +
  scale_y_continuous(labels = function(x) scales::scientific(x, digits = 2)) +
  guides(color = "none") +
  ylab("Synthesis ± SEM\nµmol/min") +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.direction = "horizontal", legend.position = "bottom", text = element_text(size = rsize), strip.text = element_text(size = rel(3)))

p_rat_conc_ratio <- readRDS(file = "../../results/data_stats_rat_surv_vs_nonsurv/rat_ratio_ci_all_mats_plot_object.RData")
p_betaoxid <- ggdraw(xlim = c(-0.05, 1.05)) + draw_image("../../results/data_stats/kinetNetwork_v2.svg")

p_sub_sink_beta <- plot_grid(p_sink_fluxes, p_betaoxid, labels = c("C", "D"), nrow = 2, rel_heights = c(1, 0.6))
p_sub_model <- plot_grid(p_vmax, p_sub_sink_beta, ncol = 2, rel_widths = c(2, 3))
panel7 <- plot_grid(p_rat_conc_ratio, p_sub_model, p_model_legend, ncol = 1, axis = "lr", align = "rl", rel_heights = c(1.1, 1, 0.1), labels = c("A", "B"))
ggsave(plot = panel7, filename = "kinetic_model_result_subset_and_rat_ratios.png", path = out_dir, width = 9, height = 9, units = "in")
