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

par_est_res <- readRDS(file = paste0(out_dir, "par_est_res_image.RData"))

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
pd_sub <- subset(pd, parameter %in% c("Vcpt2", "Vmcad", "Vmschad", "Vcrot"))
pd_sub$parameter <- substring(pd_sub$parameter, 2)
pd_sub$parameter <- toUpper(pd_sub$parameter)
p_vmax <- ggplot(data = pd_sub, mapping = aes(x = Day, y = value, color = Group)) + 
  facet_wrap(facets = ~ parameter, ncol = 4, nrow = 1, scale = "free_y") +
  geom_hline(mapping = aes(yintercept = lower_bound)) +
  geom_hline(mapping = aes(yintercept = upper_bound)) +
  stat_summary(geom = "line", fun.data = mean_se) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(width = 0.2)) +
  scale_y_log10() +
  scale_x_continuous(breaks = 0:2) +
  ylab("Concentration mean +/- SEM") +
  theme_bw() + 
  theme(panel.grid = element_blank())
ggsave(filename = "kin_mitomod_Vmax_S_vs_NS_repeated_subset.png", plot = p_vmax, path = out_dir, width = 4, height = 4, units = "in")

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

tic()
cl <- makeCluster(min(detectCores() - 1, 100), outfile = "")
prep_res <- clusterCall(cl = cl, fun = eval, quote({library(CoRC)}), env = .GlobalEnv)
mpar_res <- parLapplyLB(cl = cl, X = models_list, fun = get_global_parameters)
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
mconcs$Group <- "Septic-S"
mconcs$Group[grepl(pattern = "[merge]", x = mconcs$variable, fixed = TRUE)] <- "Septic-NS"
mconcs$variable <- sub(pattern = "[merge]", replacement = "", x = mconcs$variable, fixed = TRUE)

gpars <- as.data.frame(sapply(mpar_res, function(e) subset(e, stri_detect(str = e$name, fixed = "MalCoA"))))
gpars$Group <- "Septic-S"
gpars$Group[stri_detect(str = gpars$name, fixed = "[merge]")] <- "Septic-NS"
gpars$name <- sub(pattern = "[merge]", replacement = "", gpar$name, fixed = TRUE)
gpars$Run <- run_list
gpars$Day <- day_list
p_malcoa <- ggplot(data = gpars, mapping = aes(x = Day, y = value, color = Group)) +
  stat_summary(geom = "line", fun.data = mean_se) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(width = 0.2)) +
  #scale_y_log10() +
  #scale_x_continuous(breaks = 0:2) +
  ylab("Mean +/- SEM, µmol/l") +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.direction = "horizontal", legend.position = "bottom")
ggsave(filename = "MalCoA_conc_S_vs_NS_repeated.png", plot = p, path = out_dir, width = 1.5, height = 1.8, units = "in")  

p <- ggplot(data = mconcs, mapping = aes(x = Day, y = value, color = Group)) + 
  facet_wrap(facets = ~ variable, ncol = 6, nrow = 8, scale = "free_y") +
  stat_summary(geom = "line", fun.data = mean_se) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(width = 0.2)) +
  #scale_y_log10() +
  #scale_x_continuous(breaks = 0:2) +
  ylab("Mean +/- SEM, µmol/l") +
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
mfluxes$Group <- "Septic-S"
mfluxes$Group[grepl(pattern = "[merge]", x = mfluxes$variable, fixed = TRUE)] <- "Septic-NS"
mfluxes$variable <- sub(pattern = "[merge]", replacement = "", x = mfluxes$variable, fixed = TRUE)

p <- ggplot(data = mfluxes, mapping = aes(x = Day, y = value, color = Group)) + 
  facet_wrap(facets = ~ variable, ncol = 10, nrow = 6, scale = "free_y") +
  stat_summary(geom = "line", fun.data = mean_se) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(width = 0.2)) +
  #scale_y_log10() +
  #scale_x_continuous(breaks = 0:2) +
  ylab("Mean +/- SEM, µmol/min") +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.direction = "horizontal", legend.position = "bottom")
ggsave(filename = "kin_mitomod_flux_S_vs_NS_repeated.png", plot = p, path = out_dir, width = 15, height = 8, units = "in")
