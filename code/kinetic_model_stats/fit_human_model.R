library(ggplot2)
library(stringi)
library(data.table)
library(matrixStats)
library(CoRC)
library(tictoc)

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

#To later turn into a function that runs in parallel from random start points

model_fit_function <- function(number = 1){
  ##Load and clean model
  new_model <- loadModel(path = paste0(base_model_dir, "human_beta_oxidation_twin_model.cps"))
  clearExperiments(model = new_model)
  #old_pe_settings <- getParameterEstimationSettings(model = new_model)
  clearParameterEstimationParameters(model = new_model)
  ##Set concentrations to steady state of original model
  setSpecies(key = ss_conc$key, initial_concentration = ss_conc$concentration, model = new_model)
  ##Set species initial concentrations
  sprefs <- getSpeciesReferences(model = new_model)
  sprefs <- subset(sprefs, !(substring(name, first = 1, last = 3) %in% c("NAD", "FAD", "CoA")))
  spvals <- getSpecies(model = new_model)
  spvals <- spvals[spvals$key %in% sprefs$key, ]
  spvals <- spvals[match(sprefs$key, spvals$key), ]
  lb <- spvals$initial_concentration / 10
  ub <- spvals$initial_concentration * 10
  species_params <- lapply(1:nrow(sprefs), function(n) defineParameterEstimationParameter(ref = sprefs$initial_concentration[n], lower_bound = lb[n], upper_bound = ub[n], start_value = spvals$initial_concentration[n]))
  dummy <- lapply(species_params, addParameterEstimationParameter, model = new_model)
  ##Set kinetic parameter names to estimate
  gcrefs <- getGlobalQuantityReferences(model = new_model)
  gcvals <- getGlobalQuantities(model = new_model)
  gcvals$initial_value[match(merge_gcvals$name, gcvals$name)] <- merge_gcvals$initial_value
  react_Vms <- c("Vcpt2", "Vcrot", "Vmcad", "Vmckat", "Vmschad", "Vmtp", "Vscad", "Vvlcad", "Ksacesink", "Ksfadhsink")
  react_Vms <- c(react_Vms, paste0(react_Vms, "[merge]"))
  gcrefs <- gcrefs[gcrefs$name %in% react_Vms, ]
  gcvals <- gcvals[gcvals$key %in% gcrefs$key, ]
  gcrefs <- gcrefs[match(react_Vms, gcrefs$name), ]
  gcvals <- gcvals[match(gcrefs$key, gcvals$key), ]
  react_Vms <- paste0("{Values[", react_Vms, "].InitialValue}")
  lb <- gcvals$initial_value / 10
  ub <- gcvals$initial_value * 10
  react_params <- lapply(1:nrow(gcrefs), function(n) defineParameterEstimationParameter(ref = react_Vms[n], lower_bound = lb[n], upper_bound = ub[n], start_value = gcvals$initial_value[n]))
  dummy <- lapply(react_params, addParameterEstimationParameter, model = new_model)
  ##Load and add experiments, Day 0
  exp_tc_data <- fread("../../data/measurements/human_van_Eunen_ac_ratios_daystep_day_0_time_in_minutes.csv")
  col_keep <- c("C0", "C10", "C10:1", "C12", "C12:1", "C14", "C14:1", "C4:1", "C6 (C4:1-DC)", "C5-DC (C6-OH)", "C8", "Day")
  exp_tc_data <- subset(exp_tc_data, select = col_keep)
  ratio_quant_map <- c("Car_ratio", "C10AcylCoA_ratio", "C10EnoylCoA_ratio", "C12AcylCoA_ratio", "C12EnoylCoA_ratio", "C14AcylCoA_ratio", "C14EnoylCoA_ratio", "C4EnoylCoA_ratio", "C6AcylCoA_ratio", "C6HydroxyacylCoA_ratio", "C8AcylCoA_ratio", "Time")
  ratio_quant_map <- paste0("{Values[", ratio_quant_map, "]}")
  types <- c(rep("dependent", ncol(exp_tc_data) - 1), "time")
  ratio_exp <- defineExperiments(experiment_type = "time_course", data = exp_tc_data, types = types, mappings = ratio_quant_map, weight_method = "mean_square")
  addExperiments(ratio_exp, model = new_model)
  ##Estimate parameters
  setParameterEstimationSettings(model = new_model, method = "SRES", update_model = TRUE, randomize_start_values = TRUE) #parameters are already stored in model
  pe_res_sres_day0 <- runParameterEstimation(model = new_model)
  setParameterEstimationSettings(model = new_model, method = "HookeJeeves", update_model = TRUE, randomize_start_values = TRUE)
  pe_res_hj_day0 <- runParameterEstimation(model = new_model)
  
  saveModel(model = new_model, filename = paste0(out_dir, "fitted_model_pilot__number_", number, "_day0.cps"), overwrite = TRUE)
  
  ##Load and add Experiments, Day 1
  clearExperiments(model = new_model)
  exp_tc_data <- fread("../../data/measurements/human_van_Eunen_ac_ratios_daystep_day_1_time_in_minutes.csv")
  col_keep <- c("C0", "C10", "C10:1", "C12", "C12:1", "C14", "C14:1", "C4:1", "C6 (C4:1-DC)", "C5-DC (C6-OH)", "C8", "Day")
  exp_tc_data <- subset(exp_tc_data, select = col_keep)
  ratio_quant_map <- c("Car_ratio", "C10AcylCoA_ratio", "C10EnoylCoA_ratio", "C12AcylCoA_ratio", "C12EnoylCoA_ratio", "C14AcylCoA_ratio", "C14EnoylCoA_ratio", "C4EnoylCoA_ratio", "C6AcylCoA_ratio", "C6HydroxyacylCoA_ratio", "C8AcylCoA_ratio", "Time")
  ratio_quant_map <- paste0("{Values[", ratio_quant_map, "]}")
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
  col_keep <- c("C0", "C10", "C10:1", "C12", "C12:1", "C14", "C14:1", "C4:1", "C6 (C4:1-DC)", "C5-DC (C6-OH)", "C8", "Day")
  exp_tc_data <- subset(exp_tc_data, select = col_keep)
  ratio_quant_map <- c("Car_ratio", "C10AcylCoA_ratio", "C10EnoylCoA_ratio", "C12AcylCoA_ratio", "C12EnoylCoA_ratio", "C14AcylCoA_ratio", "C14EnoylCoA_ratio", "C4EnoylCoA_ratio", "C6AcylCoA_ratio", "C6HydroxyacylCoA_ratio", "C8AcylCoA_ratio", "Time")
  ratio_quant_map <- paste0("{Values[", ratio_quant_map, "]}")
  types <- c(rep("dependent", ncol(exp_tc_data) - 1), "time")
  ratio_exp <- defineExperiments(experiment_type = "time_course", data = exp_tc_data, types = types, mappings = ratio_quant_map, weight_method = "mean_square")
  addExperiments(ratio_exp, model = new_model)
  ##Estimate parameters
  setParameterEstimationSettings(model = new_model, method = "SRES", update_model = TRUE, randomize_start_values = FALSE) #parameters are already stored in model
  pe_res_sres_day2 <- runParameterEstimation(model = new_model)
  setParameterEstimationSettings(model = new_model, method = "HookeJeeves", update_model = TRUE, randomize_start_values = FALSE)
  pe_res_hj_day2 <- runParameterEstimation(model = new_model)
  #Save fitted model
  saveModel(model = new_model, filename = paste0(out_dir, "fitted_model_pilot_number_", number, "_day2.cps"), overwrite = TRUE)
  #Return results
  return(list(sres_res_d0 = pe_res_sres_day0, hj_res_d0 = pe_res_hj_day0, sres_res_d1 = pe_res_sres_day1, hj_res_d1 = pe_res_hj_day1, sres_res_d2 = pe_res_sres_day2, hj_res_d2 = pe_res_hj_day2))
}
par_est_res <- mclapply(1:1000, model_fit_function)
save.image()