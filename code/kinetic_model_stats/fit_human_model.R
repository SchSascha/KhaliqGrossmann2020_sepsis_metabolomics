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

#Input

base_model <- loadModel(path = paste0(base_model_dir, "human_beta_oxidation_twin_model.cps"))

#To later turn into a function that runs in parallel from random start points

new_model <- loadModel(path = paste0(base_model_dir, "human_beta_oxidation_twin_model.cps"))
clearExperiments(model = new_model)
old_pe_settings <- getParameterEstimationSettings(model = new_model)
clearParameterEstimationSettings(model = new_model)
clearParameterEstimationParameters(model = new_model)
exp_tc_data <- fread("../../data/measurements/human_van_Eunen_ac_ratios_daystep_day_0_time_in_minutes.csv")
col_keep <- c("C0", "C10", "C10:1", "C12", "C12:1", "C14", "C14:1", "C4:1", "C6 (C4:1-DC)", "C5-DC (C6-OH)", "C8", "Day")
ratio_quant_map <- c("Car_ratio", "C10AcylCoA_ratio", "C10EnoylCoA_ratio", "C12AcylCoA_ratio", "C12EnoylCoA_ratio", "C14AcylCoA_ratio", "C14EnoylCoA_ratio", "C4EnoylCoA_ratio", "C6AcylCoA_ratio", "C6HydroxyacylCoA_ratio", "C8AcylCoA_ratio", "Time")
types <- c(strrep("dependent", ncols(exp_tc_data) - 1, "time")
exp_tc_data <- subset(exp_tc_data, select = col_keep)
ratio_exp <- defineExperiment(experiment_type = "time_course", data = exp_tc_data, types = types, mappings = ratio_quant_map, wight_method = mean_square)
params <- defineParameterEstimationParameter(
setParameterEstimationSettings(model = new_model, method = "SRES", update_model = TRUE, experiments = ratio_exp, parameters = params)
pe_res_sres <- runParameterEstimation(model = new_model)
setParameterEstimationSettings(model = new_model, method = "HookeJeeves", update_model = TRUE, experiments = ratio_exp)
pe_res_hj <- runParameterEstimation(model = new_model)