library(purrr)
library(ggplot2)
library(stringi)
library(reshape2)
library(data.table)
library(matrixStats)
library(CoRC)
library(tictoc)
library(parallel)

source("../function_definitions.R")

models_dir <- "../../results/kinetic_modeling/"
out_dir <- "../../results/kinetic_modeling/"

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
mconcs$Group <- "Survivor"
mconcs$Group[grepl(pattern = "[merge]", x = mconcs$variable, fixed = TRUE)] <- "Nonsurvivor"
mconcs$variable <- sub(pattern = "[merge]", replacement = "", x = mconcs$variable, fixed = TRUE)

p <- ggplot(data = mconcs, mapping = aes(x = Day, y = value, color = Group)) + 
  facet_wrap(facets = ~ variable, ncol = 6, nrow = 8, scale = "free_y") +
  stat_summary(geom = "line", fun.data = mean_se) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(width = 0.2)) +
  #scale_y_log10() +
  #scale_x_continuous(breaks = 0:2) +
  ylab("value +/- 1 SEM") +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.direction = "horizontal", legend.position = "bottom")
ggsave(filename = "kin_mitomod_conc_S_vs_NS_repeated.png", plot = p, path = out_dir, width = 15, height = 8, units = "in")

fluxes <- as.data.frame(sapply(stst_est_res, function(ssres) ssres$reactions$flux))
rownames(fluxes) <- stst_est_res[[1]]$reactions$name
fluxes <- transpose(fluxes)
colnames(fluxes) <- stst_est_res[[1]]$reactions$name
fluxes$Day <- day_list
fluxes$Run <- run_list
mfluxes <- melt(fluxes, id.vars = c("Day", "Run"))
mfluxes$Group <- "Survivor"
mfluxes$Group[grepl(pattern = "[merge]", x = mfluxes$variable, fixed = TRUE)] <- "Nonsurvivor"
mfluxes$variable <- sub(pattern = "[merge]", replacement = "", x = mfluxes$variable, fixed = TRUE)

p <- ggplot(data = mfluxes, mapping = aes(x = Day, y = value, color = Group)) + 
  facet_wrap(facets = ~ variable, ncol = 10, nrow = 6, scale = "free_y") +
  stat_summary(geom = "line", fun.data = mean_se) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(width = 0.2)) +
  #scale_y_log10() +
  #scale_x_continuous(breaks = 0:2) +
  ylab("value +/- 1 SEM") +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.direction = "horizontal", legend.position = "bottom")
ggsave(filename = "kin_mitomod_flux_S_vs_NS_repeated.png", plot = p, path = out_dir, width = 15, height = 8, units = "in")
