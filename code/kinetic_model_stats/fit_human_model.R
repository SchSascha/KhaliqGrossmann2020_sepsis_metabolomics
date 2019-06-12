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

#Paths
out_dir <- "../../results/kinetic_modeling/"
base_model_dir <- "../../data/template_models/"

if (!dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE)

#Get steady state concentrations from original model
base_model <- loadModel(path = paste0(base_model_dir, "human_beta_oxidation_model.cps"))
base_ss <- runSteadyState(model = base_model)
merge_ss <- base_ss
merge_ss$species$key <- sub(pattern = "{VMAT}", replacement = "{VMAT[merge]}", x = merge_ss$species$key, fixed = TRUE)
merge_ss$species$key <- sub(pattern = "{VCYT}", replacement = "{VCYT[merge]}", x = merge_ss$species$key, fixed = TRUE)
ss_conc <- rbind(base_ss$species, merge_ss$species) 
ss_conc <- subset(ss_conc, type == "reactions")
base_gcvals <- getGlobalQuantities(model = base_model)
merge_gcvals <- base_gcvals
merge_gcvals$name <- paste0(merge_gcvals$name, "[merge]")
merge_gcvals <- rbind(base_gcvals, merge_gcvals)

#Set optimization endpoints and mappings - stays the same for all optimization steps
#col_keep <- c("C0", "C10", "C10:1", "C12", "C12:1", "C14", "C14:1", "C16", "C16:1", "C16-OH", "C2", "C4", "C4:1", "C3-DC (C4-OH)", "C6 (C4:1-DC)", "C6:1", "C5-DC (C6-OH)", "C8", "Day")
#ratio_quant_map <- c("Car_ratio", "C10AcylCoA_ratio", "C10EnoylCoA_ratio", "C12AcylCoA_ratio", "C12EnoylCoA_ratio", "C14AcylCoA_ratio", "C14EnoylCoA_ratio", "C16AcylCoA_ratio", "C16EnoylCoA_ratio", "C16HydroxyacylCoA_ratio", "AcetylCoA_ratio", "C4AcylCoA_ratio", "C4EnoylCoA_ratio", "C4HydroxyacylCoA_ratio", "C6AcylCoA_ratio", "C6EnoylCoA_ratio", "C6HydroxyacylCoA_ratio", "C8AcylCoA_ratio", "Time")
col_keep <- c("C0", "C10", "C10:1", "C12", "C12:1", "C4:1", "C3-DC (C4-OH)", "C6 (C4:1-DC)", "C6:1", "C5-DC (C6-OH)", "C8", "Day")
ratio_quant_map <- c("Car_ratio", "C10AcylCoA_ratio", "C10EnoylCoA_ratio", "C12AcylCoA_ratio", "C12EnoylCoA_ratio", "C4EnoylCoA_ratio", "C4HydroxyacylCoA_ratio", "C6AcylCoA_ratio", "C6EnoylCoA_ratio", "C6HydroxyacylCoA_ratio", "C8AcylCoA_ratio", "Time")
ratio_quant_map <- paste0("{Values[", ratio_quant_map, "]}")
#Set kinetic parameter names to estimate - also the same for all optimization steps
react_pars <- c("Vcpt2", "Vcrot", "Vmcad", "Vmckat", "Vmschad", "Vmtp", "Vscad", "Vvlcad") #, "Ksacesink", "Ksfadhsink") # do not vary sink parameters!
react_pars <- c(react_pars, paste0(react_pars, "[merge]"))

#TODO: find out how to add constraints to parameter estimation taks
#TODO: add constraints on t=1440 concentrations to be similar to t=0 concentrations

#Fitting function
model_fit_function <- function(number = 1, base_model_dir, out_dir, ss_conc, merge_gcvals = merge_gcvals, react_pars = react_pars, col_keep = col_keep, ratio_quant_map = ratio_quant_map){
  ##Load and clean model
  new_model <- loadModel(path = paste0(base_model_dir, "human_beta_oxidation_twin_model.cps"))
  clearExperiments(model = new_model)
  #old_pe_settings <- getParameterEstimationSettings(model = new_model)
  clearParameterEstimationParameters(model = new_model)
  ##Set concentrations to steady state of original model
  setSpecies(key = ss_conc$key, initial_concentration = ss_conc$concentration, model = new_model)
  ##Set species initial concentrations
  sprefs <- getSpeciesReferences(model = new_model)
  sprefs <- subset(sprefs, !(substring(name, first = 1, last = 3) %in% c("NAD", "FAD", "CoA", "Mal"))) # to exclude MalCoA influence, bc. init = 0
  spvals <- getSpecies(model = new_model)
  spvals <- spvals[spvals$key %in% sprefs$key, ]
  spvals <- spvals[match(sprefs$key, spvals$key), ]
  lb <- spvals$initial_concentration / 2
  ub <- spvals$initial_concentration * 2
  species_params <- lapply(1:nrow(sprefs), function(n) defineParameterEstimationParameter(ref = sprefs$initial_concentration[n], lower_bound = lb[n], upper_bound = ub[n], start_value = spvals$initial_concentration[n]))
  dummy <- lapply(species_params, addParameterEstimationParameter, model = new_model)
  ##Set kinetic parameter names to estimate
  gcrefs <- getGlobalQuantityReferences(model = new_model)
  gcvals <- getGlobalQuantities(model = new_model)
  gcvals$initial_value[match(merge_gcvals$name, gcvals$name)] <- merge_gcvals$initial_value
  gcrefs <- gcrefs[gcrefs$name %in% react_pars, ]
  gcvals <- gcvals[gcvals$key %in% gcrefs$key, ]
  gcrefs <- gcrefs[match(react_pars, gcrefs$name), ]
  gcvals <- gcvals[match(gcrefs$key, gcvals$key), ]
  react_Vms <- paste0("{Values[", react_pars, "].InitialValue}")
  lb <- gcvals$initial_value / 5
  ub <- gcvals$initial_value * 5
  react_params <- lapply(1:nrow(gcrefs), function(n) defineParameterEstimationParameter(ref = react_Vms[n], lower_bound = lb[n], upper_bound = ub[n], start_value = gcvals$initial_value[n]))
  dummy <- lapply(react_params, addParameterEstimationParameter, model = new_model)
  ##Load and add experiments, Day 0
  exp_tc_data <- fread("../../data/measurements/human_van_Eunen_ac_ratios_daystep_day_0_time_in_minutes.csv")
  exp_tc_data <- subset(exp_tc_data, select = col_keep)
  types <- c(rep("dependent", ncol(exp_tc_data) - 1), "time")
  ratio_exp <- defineExperiments(experiment_type = "time_course", data = exp_tc_data, types = types, mappings = ratio_quant_map, weight_method = "mean_square")
  addExperiments(ratio_exp, model = new_model)
  ##Estimate parameters
  setParameterEstimationSettings(model = new_model, method = "SRES", update_model = TRUE, randomize_start_values = FALSE) #parameters are already stored in model
  pe_res_sres_day0 <- runParameterEstimation(model = new_model)
  setParameterEstimationSettings(model = new_model, method = "HookeJeeves", update_model = TRUE, randomize_start_values = FALSE)
  pe_res_hj_day0 <- runParameterEstimation(model = new_model)
  
  saveModel(model = new_model, filename = paste0(out_dir, "fitted_model_pilot_number_", number, "_day0.cps"), overwrite = TRUE)
  
  ##Load and add Experiments, Day 1
  clearExperiments(model = new_model)
  exp_tc_data <- fread("../../data/measurements/human_van_Eunen_ac_ratios_daystep_day_1_time_in_minutes.csv")
  exp_tc_data <- subset(exp_tc_data, select = col_keep)
  types <- c(rep("dependent", ncol(exp_tc_data) - 1), "time")
  ratio_exp <- defineExperiments(experiment_type = "time_course", data = exp_tc_data, types = types, mappings = ratio_quant_map, weight_method = "mean_square")
  addExperiments(ratio_exp, model = new_model)
  ##Estimate parameters
  setParameterEstimationSettings(model = new_model, method = "SRES", update_model = TRUE, randomize_start_values = FALSE) #parameters are already stored in model
  pe_res_sres_day1 <- runParameterEstimation(model = new_model)
  setParameterEstimationSettings(model = new_model, method = "HookeJeeves", update_model = TRUE, randomize_start_values = FALSE)
  pe_res_hj_day1 <- runParameterEstimation(model = new_model)
  
  saveModel(model = new_model, filename = paste0(out_dir, "fitted_model_pilot_number_", number, "_day1.cps"), overwrite = TRUE)
  
  ##Load and add Experiments, Day 2
  clearExperiments(model = new_model)
  exp_tc_data <- fread("../../data/measurements/human_van_Eunen_ac_ratios_daystep_day_2_time_in_minutes.csv")
  exp_tc_data <- subset(exp_tc_data, select = col_keep)
  types <- c(rep("dependent", ncol(exp_tc_data) - 1), "time" )       
  ratio_exp <- defineExperiments(experiment_type = "time_course", data = exp_tc_data, types = types, mappings = ratio_quant_map, weight_method = "mean_square")
  addExperiments(ratio_exp, model = new_model)
  ##Estimate parameters
  setParameterEstimationSettings(model = new_model, method = "SRES", update_model = TRUE, randomize_start_values = FALSE) #parameters are already stored in model
  pe_res_sres_day2 <- runParameterEstimation(model = new_model)
  setParameterEstimationSettings(model = new_model, method = "HookeJeeves", update_model = TRUE, randomize_start_values = FALSE)
  pe_res_hj_day2 <- runParameterEstimation(model = new_model)
  #Save fitted model
  saveModel(model = new_model, filename = paste0(out_dir, "fitted_model_pilot_number_", number, "_day2.cps"), overwrite = TRUE)
  #unload model
  unloadModel(model = new_model)
  
  #Return results
  return(list(sres_res_d0 = pe_res_sres_day0, hj_res_d0 = pe_res_hj_day0, sres_res_d1 = pe_res_sres_day1, hj_res_d1 = pe_res_hj_day1, sres_res_d2 = pe_res_sres_day2, hj_res_d2 = pe_res_hj_day2))
}

model_fit_function_with_prepared_models <- function(number = 1, base_model_dir, out_dir, ss_conc, merge_gcvals = merge_gcvals, react_pars = react_pars){
  ##Load and clean model for day 0
  new_model_day0 <- loadModel(path = paste0(base_model_dir, "human_beta_oxidation_twin_model_day0.cps"))
  #clearParameterEstimationParameters(model = new_model_day0)
  ##Set concentrations to steady state of original model
  setSpecies(key = ss_conc$key, initial_concentration = ss_conc$concentration, model = new_model_day0)
  ##Set species initial concentrations and concentration boundaries
  sprefs <- getSpeciesReferences(model = new_model_day0)
  sprefs <- subset(sprefs, !(substring(name, first = 1, last = 3) %in% c("NAD", "FAD", "CoA", "Mal"))) # to exclude MalCoA influence, bc. init = 0
  spvals <- getSpecies(model = new_model_day0)
  spvals <- spvals[spvals$key %in% sprefs$key, ]
  spvals <- spvals[match(sprefs$key, spvals$key), ]
  sp_lb <- spvals$initial_concentration / 2
  sp_ub <- spvals$initial_concentration * 2
  species_params <- lapply(1:nrow(sprefs), function(n) defineParameterEstimationParameter(ref = sprefs$initial_concentration[n], lower_bound = sp_lb[n], upper_bound = sp_ub[n], start_value = spvals$initial_concentration[n]))
  dummy <- lapply(species_params, addParameterEstimationParameter, model = new_model_day0)
  ##Set kinetic parameter names to estimate
  gcrefs <- getGlobalQuantityReferences(model = new_model_day0)
  gcvals <- getGlobalQuantities(model = new_model_day0)
  gcvals$initial_value[match(merge_gcvals$name, gcvals$name)] <- merge_gcvals$initial_value
  gcrefs <- gcrefs[gcrefs$name %in% react_pars, ]
  gcvals <- gcvals[gcvals$key %in% gcrefs$key, ]
  gcrefs <- gcrefs[match(react_pars, gcrefs$name), ]
  gcvals <- gcvals[match(gcrefs$key, gcvals$key), ]
  react_Vms <- paste0("{Values[", react_pars, "].InitialValue}")
  re_lb <- gcvals$initial_value / 5
  re_ub <- gcvals$initial_value * 5
  react_params <- lapply(1:nrow(gcrefs), function(n) defineParameterEstimationParameter(ref = react_Vms[n], lower_bound = re_lb[n], upper_bound = re_ub[n], start_value = gcvals$initial_value[n]))
  dummy <- lapply(react_params, addParameterEstimationParameter, model = new_model_day0)
  ##Estimate parameters
  setParameterEstimationSettings(model = new_model_day0, method = "SRES", update_model = TRUE, randomize_start_values = FALSE) #parameters are already stored in model
  pe_res_sres_day0 <- runParameterEstimation(model = new_model_day0)
  setParameterEstimationSettings(model = new_model_day0, method = "HookeJeeves", update_model = TRUE, randomize_start_values = FALSE)
  pe_res_hj_day0 <- runParameterEstimation(model = new_model_day0)
  
  saveModel(model = new_model_day0, filename = paste0(out_dir, "fitted_model_pilot_number_", number, "_day0.cps"), overwrite = TRUE)
  
  ##Load and clean model for day 1
  new_model_day1 <- loadModel(path = paste0(base_model_dir, "human_beta_oxidation_twin_model_day1.cps"))
  #clearParameterEstimationParameters(model = new_model_day0)
  ##Set concentrations to levels of day 0 model at t=1440 but keep boundaries as in day 0
  tc_d0 <- runTimeCourse(model = new_model_day0, duration = 100, dt = 0.1, save_result_in_memory = TRUE, update_model = TRUE)
  tc_d0 <- tc_d0$result[nrow(tc_d0$result), ]
  tc_d0 <- melt(tc_d0, id.vars = "Time")
  tc_d0$variable <- as.character(tc_d0$variable)
  tc_d0 <- subset(tc_d0, variable %in% sprefs$key)
  setSpecies(key = tc_d0$variable, initial_concentration = tc_d0$value, model = new_model_day1)
  sprefs <- getSpeciesReferences(model = new_model_day1)
  sprefs <- subset(sprefs, !(substring(name, first = 1, last = 3) %in% c("NAD", "FAD", "CoA", "Mal"))) # to exclude MalCoA influence, bc. init = 0
  spvals <- getSpecies(model = new_model_day1)
  spvals <- spvals[spvals$key %in% sprefs$key, ]
  spvals <- spvals[match(sprefs$key, spvals$key), ]
  species_params <- lapply(1:nrow(sprefs), function(n) defineParameterEstimationParameter(ref = sprefs$initial_concentration[n], lower_bound = sp_lb[n], upper_bound = sp_ub[n], start_value = min(max(spvals$initial_concentration[n], sp_lb[n]), sp_ub[n])))
  dummy <- lapply(species_params, addParameterEstimationParameter, model = new_model_day1)
  ##Set kinetic parameter names to estimate
  gcrefs <- getGlobalQuantityReferences(model = new_model_day1)
  gcvals <- getGlobalQuantities(model = new_model_day1)
  d0_gcrefs <- getGlobalQuantityReferences(model = new_model_day0)
  d0_gcvals <- getGlobalQuantities(model = new_model_day0)
  gcvals$initial_value[match(d0_gcvals$name, gcvals$name)] <- d0_gcvals$initial_value
  gcrefs <- gcrefs[gcrefs$name %in% react_pars, ]
  gcvals <- gcvals[gcvals$key %in% gcrefs$key, ]
  gcrefs <- gcrefs[match(react_pars, gcrefs$name), ]
  gcvals <- gcvals[match(gcrefs$key, gcvals$key), ]
  react_Vms <- paste0("{Values[", react_pars, "].InitialValue}")
  react_params <- lapply(1:nrow(gcrefs), function(n) defineParameterEstimationParameter(ref = react_Vms[n], lower_bound = re_lb[n], upper_bound = re_ub[n], start_value = min(max(gcvals$initial_value[n], re_lb[n]), re_ub[n])))
  dummy <- lapply(react_params, addParameterEstimationParameter, model = new_model_day1)
  ##Estimate parameters
  setParameterEstimationSettings(model = new_model_day1, method = "SRES", update_model = TRUE, randomize_start_values = FALSE) #parameters are already stored in model
  pe_res_sres_day1 <- runParameterEstimation(model = new_model_day1)
  setParameterEstimationSettings(model = new_model_day1, method = "HookeJeeves", update_model = TRUE, randomize_start_values = FALSE)
  pe_res_hj_day1 <- runParameterEstimation(model = new_model_day1)
  
  saveModel(model = new_model_day1, filename = paste0(out_dir, "fitted_model_pilot_number_", number, "_day1.cps"), overwrite = TRUE)
  
  ##Load and clean model for day 2
  new_model_day2 <- loadModel(path = paste0(base_model_dir, "human_beta_oxidation_twin_model_day2.cps"))
  #clearParameterEstimationParameters(model = new_model_day0)
  ##Set concentrations to levels of day 0 model at t=1440 but keep boundaries as in day 0
  tc_d1 <- runTimeCourse(model = new_model_day1, duration = 100, dt = 0.1, save_result_in_memory = TRUE, update_model = TRUE)
  tc_d1 <- tc_d1$result[nrow(tc_d1$result), ]
  tc_d1 <- melt(tc_d1, id.vars = "Time")
  tc_d1$variable <- as.character(tc_d1$variable)
  tc_d1 <- subset(tc_d1, variable %in% sprefs$key)
  setSpecies(key = tc_d1$variable, initial_concentration = tc_d1$value, model = new_model_day2)
  sprefs <- getSpeciesReferences(model = new_model_day2)
  sprefs <- subset(sprefs, !(substring(name, first = 1, last = 3) %in% c("NAD", "FAD", "CoA", "Mal"))) # to exclude MalCoA influence, bc. init = 0
  spvals <- getSpecies(model = new_model_day1)
  spvals <- spvals[spvals$key %in% sprefs$key, ]
  spvals <- spvals[match(sprefs$key, spvals$key), ]
  species_params <- lapply(1:nrow(sprefs), function(n) defineParameterEstimationParameter(ref = sprefs$initial_concentration[n], lower_bound = sp_lb[n], upper_bound = sp_ub[n], start_value = min(max(spvals$initial_concentration[n], sp_lb[n]), sp_ub[n])))
  dummy <- lapply(species_params, addParameterEstimationParameter, model = new_model_day2)
  ##Set kinetic parameter names to estimate
  gcrefs <- getGlobalQuantityReferences(model = new_model_day2)
  gcvals <- getGlobalQuantities(model = new_model_day2)
  d1_gcrefs <- getGlobalQuantityReferences(model = new_model_day1)
  d1_gcvals <- getGlobalQuantities(model = new_model_day1)
  gcvals$initial_value[match(d1_gcvals$name, gcvals$name)] <- d1_gcvals$initial_value
  gcrefs <- gcrefs[gcrefs$name %in% react_pars, ]
  gcvals <- gcvals[gcvals$key %in% gcrefs$key, ]
  gcrefs <- gcrefs[match(react_pars, gcrefs$name), ]
  gcvals <- gcvals[match(gcrefs$key, gcvals$key), ]
  react_Vms <- paste0("{Values[", react_pars, "].InitialValue}")
  react_params <- lapply(1:nrow(gcrefs), function(n) defineParameterEstimationParameter(ref = react_Vms[n], lower_bound = re_lb[n], upper_bound = re_ub[n], start_value = min(max(gcvals$initial_value[n], re_lb[n]), re_ub[n])))
  dummy <- lapply(react_params, addParameterEstimationParameter, model = new_model_day2)
  ##Estimate parameters
  setParameterEstimationSettings(model = new_model_day2, method = "SRES", update_model = TRUE, randomize_start_values = FALSE) #parameters are already stored in model
  pe_res_sres_day2 <- runParameterEstimation(model = new_model_day2)
  setParameterEstimationSettings(model = new_model_day2, method = "HookeJeeves", update_model = TRUE, randomize_start_values = FALSE)
  pe_res_hj_day2 <- runParameterEstimation(model = new_model_day2)
  
  saveModel(model = new_model_day2, filename = paste0(out_dir, "fitted_model_pilot_number_", number, "_day1.cps"), overwrite = TRUE)
  
  #Clean up
  unloadModel(model = new_model_day0)
  unloadModel(model = new_model_day1)
  unloadModel(model = new_model_day2)
  
  #Return results
  return(list(sres_res_d0 = pe_res_sres_day0, hj_res_d0 = pe_res_hj_day0, sres_res_d1 = pe_res_sres_day1, hj_res_d1 = pe_res_hj_day1, sres_res_d2 = pe_res_sres_day2, hj_res_d2 = pe_res_hj_day2))
}

#Actual fitting
tic()
cl <- makeCluster(detectCores() - 1)
prep_res <- clusterCall(cl = cl, fun = eval, quote({library(CoRC); library(data.table)}), env = .GlobalEnv)
#par_est_res <- parLapplyLB(cl = cl, X = 1:1000, fun = model_fit_function, base_model_dir = base_model_dir, out_dir = out_dir, ss_conc = ss_conc, merge_gcvals = merge_gcvals, react_pars = react_pars, col_keep = col_keep, ratio_quant_map = ratio_quant_map)
par_est_res <- parLapplyLB(cl = cl, X = 1:7, fun = model_fit_function_with_prepared_models, base_model_dir = base_model_dir, out_dir = out_dir, ss_conc = ss_conc, merge_gcvals = merge_gcvals, react_pars = react_pars)
stopCluster(cl)
toc()

#save.image()

#Extract optimization results
par_d0 <- lapply(lapply(lapply(par_est_res, `[[`, "hj_res_d0"), `[[`, "parameters"), subset, subset = stri_startswith_fixed(str = parameter, pattern = "Values["), select = c("parameter", "value", "lower_bound", "upper_bound"))
pd0 <- Reduce("rbind", par_d0)
pd0$Day <- 0
par_d1 <- lapply(lapply(lapply(par_est_res, `[[`, "hj_res_d1"), `[[`, "parameters"), subset, subset = stri_startswith_fixed(str = parameter, pattern = "Values["), select = c("parameter", "value", "lower_bound", "upper_bound"))
pd1 <- Reduce("rbind", par_d1)
pd1$Day <- 1
par_d2 <- lapply(lapply(lapply(par_est_res, `[[`, "hj_res_d2"), `[[`, "parameters"), subset, subset = stri_startswith_fixed(str = parameter, pattern = "Values["), select = c("parameter", "value", "lower_bound", "upper_bound"))
pd2 <- Reduce("rbind", par_d2)
pd2$Day <- 2

obj_d0 <- sapply(lapply(par_est_res, `[[`, "hj_res_d0"), function(e) e$main$objective_value)
obj_d1 <- sapply(lapply(par_est_res, `[[`, "hj_res_d1"), function(e) e$main$objective_value)
obj_d2 <- sapply(lapply(par_est_res, `[[`, "hj_res_d2"), function(e) e$main$objective_value)
obj <- data.frame(Value = c(obj_d0, obj_d1, obj_d2), Day = rep(0:2, each = length(obj_d0)))

#Plot enzyme
pd <- rbind(pd0, pd1, pd2)
pd$Run <- rep(rep(1:length(par_est_res), each = length(unique(pd0$parameter))), times = 3)
pd$Group <- c("Survivor", "Nonsurvivor")[1 + (grepl(pattern = "merge", x = pd$parameter))] #Values with [merge] belong to Nonsurvivors
pd$parameter <- stri_replace(str = pd$parameter, replacement = "", fixed = "[merge]")
pd$parameter <- stri_extract(str = pd$parameter, regex = "\\[.+\\]")
n_par <- length(unique(pd$parameter))
pd_mn <- lapply(unique(pd$parameter), function(par){ ss <- subset(x = pd, parameter == par); ss$value <- ss$value / mean(ss$value); return(ss) })
pd_mn <- Reduce("rbind", pd_mn)
p <- ggplot(data = pd, mapping = aes(x = Day, y = value, color = Group)) + 
  facet_wrap(facets = ~ parameter, ncol = 4, nrow = 2) +
  stat_summary(geom = "line", fun.data = mean_se) +
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(width = 0.2)) +
  geom_hline(mapping = aes(yintercept = lower_bound)) +
  geom_hline(mapping = aes(yintercept = upper_bound)) +
  scale_y_log10() +
  scale_x_continuous(breaks = 0:2) +
  ylab("value +/- 1 SEM") +
  theme_bw() + 
  theme(panel.grid = element_blank())
ggsave(filename = "kin_mitomod_Vmax_S_vs_NS_repeated.png", plot = p, path = out_dir, width = 8, height = 4, units = "in")

mpd <- dcast(pd, Day + Run + parameter ~ Group)
l <- ggplot(data = mpd, mapping = aes(x = Nonsurvivor, y = Survivor)) + 
  facet_grid(Day ~ parameter) + 
  #geom_point(size = 0.7, alpha = 0.1) +
  stat_bin_hex(bins = 20) + 
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
plot(l)
ggsave(plot = l, file = "kin_mitomod_rel_Vm.png", path = out_dir, width = 10, height = 4, units = "in")

c1pd <- dcast(pd, Day + Run + Group ~ parameter)
c2pd <- dcast(pd, Day + Run ~ parameter + Group)
pc1 <- prcomp(as.matrix(c1pd[c1pd$Day == 0, -1:-3]), scale. = TRUE)
pc2 <- prcomp(as.matrix(c2pd[c2pd$Day == 0, -1:-2]), scale. = TRUE)

h <- ggplot(data = obj, mapping = aes(x = Value)) +
  facet_wrap(Day ~ .) +
  scale_x_log10() +
  geom_histogram(bins = 10) +
  theme_bw() + 
  theme(panel.grid = element_blank())
ggsave(plot = h, filename = "kin_mitomod_obj_val.png", path = out_dir, width = 8, height = 4, units = "in")

#Remove COPASI models from memory
unloadAllModels()

#Remove temporary experiment files created by CoRC
file.remove(list.files(pattern = "CoRC_exp_.{9,13}\\.txt$"))
