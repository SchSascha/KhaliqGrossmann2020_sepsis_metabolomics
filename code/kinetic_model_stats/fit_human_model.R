library(ggplot2)
library(stringi)
library(data.table)
library(matrixStats)
library(CoRC)

#Paths

out_dir <- "../../results/kinetic_modeling/"
base_model_dir <- "../../data/template_models/"

if (!dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE)

#To later turn into a function that runs in parallel from random start points

##Load and clean model
new_model <- loadModel(path = paste0(base_model_dir, "human_beta_oxidation_twin_model.cps"))
clearExperiments(model = new_model)
#old_pe_settings <- getParameterEstimationSettings(model = new_model)
clearParameterEstimationParameters(model = new_model)
##Load and add experiments
exp_tc_data <- fread("../../data/measurements/human_van_Eunen_ac_ratios_daystep_day_1_time_in_minutes.csv")
col_keep <- c("C0", "C10", "C10:1", "C12", "C12:1", "C14", "C14:1", "C4:1", "C6 (C4:1-DC)", "C5-DC (C6-OH)", "C8", "Day")
exp_tc_data <- subset(exp_tc_data, select = col_keep)
ratio_quant_map <- c("Car_ratio", "C10AcylCoA_ratio", "C10EnoylCoA_ratio", "C12AcylCoA_ratio", "C12EnoylCoA_ratio", "C14AcylCoA_ratio", "C14EnoylCoA_ratio", "C4EnoylCoA_ratio", "C6AcylCoA_ratio", "C6HydroxyacylCoA_ratio", "C8AcylCoA_ratio", "Time")
ratio_quant_map <- paste0("{Values[", ratio_quant_map, "]}")
types <- c(rep("dependent", ncol(exp_tc_data) - 1), "time")
ratio_exp <- defineExperiments(experiment_type = "time_course", data = exp_tc_data, types = types, mappings = ratio_quant_map, weight_method = "mean_square")
addExperiments(ratio_exp, model = new_model)
##Set species initial concentrations
sprefs <- getSpeciesReferences(model = new_model)
sprefs <- subset(sprefs, !(substring(initial_concentration, first = 1, last = 3) %in% c("NAD", "FAD", "CoA")))
spvals <- getSpecies(model = new_model)
spvals <- spvals[spvals$key %in% sprefs$key, ]
lb <- spvals$initial_concentration / 10
ub <- spvals$initial_concentration * 10
species_params <- lapply(1:nrow(sprefs), function(n) defineParameterEstimationParameter(ref = sprefs$initial_concentration[n], lower_bound = lb[n], upper_bound = ub[n], start_value = spvals$initial_concentration[n]))
dummy <- lapply(species_params, addParameterEstimationParameter, model = new_model)
##Set kinetic parameter names to estimate
gcrefs <- getGlobalQuantityReferences(model = new_model)
gcvals <- getGlobalQuantities(model = new_model)
react_Vms <- c("Vcpt2", "Vcrot", "Vmcad", "Vmckat", "Vmschad", "Vmtp", "Vscad", "Vvlcad", "Ksacesink", "Ksfadhsink")
react_Vms <- c(react_Vms, paste0(react_Vms, "[merge]"))
gcrefs <- gcrefs[gcrefs$name %in% react_Vms, ]
gcvals <- gcvals[gcvals$key %in% gcrefs$key, ]
gcvals <- gcvals[match(gcrefs$key, gcvals$key), ]
react_Vms <- paste0("{Values[", react_Vms, "]}")
lb <- gcvals$initial_value / 10
ub <- gcvals$initial_value * 10
react_params <- lapply(1:nrow(gcrefs), function(n) defineParameterEstimationParameter(ref = react_Vms[n], lower_bound = lb[n], upper_bound = ub[n], start_value = gcvals$initial_value[n]))
dummy <- lapply(react_params, addParameterEstimationParameter, model = new_model)
##Estimate parameters
setParameterEstimationSettings(model = new_model, method = "SRES", update_model = TRUE, randomize_start_values = TRUE) #parameters are already stored in model
pe_res_sres <- runParameterEstimation(model = new_model)
setParameterEstimationSettings(model = new_model, method = "HookeJeeves", update_model = TRUE, randomize_start_values = TRUE)
pe_res_hj <- runParameterEstimation(model = new_model)
#Save fitted model
save.image(file = "model_fit_pilot.RData")
saveModel(model = new_model, filename = paste0(out_dir, "fitted_model_pilot.cps"), overwrite = TRUE)
