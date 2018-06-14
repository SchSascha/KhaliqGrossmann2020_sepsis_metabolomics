# Compute slopes and confidence intervals for rat metabolites

library(matrixStats)

source("../function_definitions.R")

out_dir <- "../../data/uFBA/"

# Read rat data

rat_data <- get_rat_sepsis_data()

# Remove outlier samples

rat_data <- rat_data[!rat_data$`Sample Identification` %in% c("060H", "039L"), ]

# Convert time point to numeric and change colname

rat_data$`time point` <- as.numeric(sub(pattern = "h", replacement = "", x = rat_data$`time point`))
colnames(rat_data)[colnames(rat_data) == "time point"] <- "time"

# Add up pool concentrations

# Fit lm for each metabolite and time change
formula <- concentration ~ time
lm_list <- list()
matv <- c("liver", "heart", "plasma")
grpv <- unique(rat_data$group)
timev <- list(c(6, 24), c(24, 72))
for (mat in matv){
  lm_list[[mat]] <- list()
	for (grp in grpv){
	  lm_list[[mat]][[grp]] <- list()
	  grpdat <- subset(rat_data, group == grp & material == mat)
	  grpdat <- subset(grpdat, TRUE, which(!colAlls(is.na(grpdat))))
	  tsv <- timev[sapply(timev, function(te) all(te %in% grpdat$time))]
	  met_set <- colnames(grpdat)[5:ncol(grpdat)]
	  for (time in tsv){
	    tstr <- paste0(time, collapse = "_")
	    lm_list[[mat]][[grp]][[tstr]] <- list()
		  for (met in met_set){
			  data <- subset(grpdat, grpdat$time %in% time, c("time", met))
			  colnames(data)[2] <- "concentration"
			  lm_list[[mat]][[grp]][[tstr]][[met]] <- lm(formula = formula, data = data)
		  }
		}
  }
}
lm_summary <- lapply(lm_list, lapply, lapply, lapply, summary)
lm_confint <- lapply(lm_list, lapply, lapply, lapply, confint)
lm_rsqr <- unlist(lapply(lm_summary, lapply, lapply, lapply, function(r) r$r.squared))
lm_slope <- lapply(lm_summary, lapply, lapply, lapply, function(s) try(s$coefficients[2, 1]))
lm_cilen <- lapply(lm_confint, lapply, lapply, sapply, function(s) try(s[2, 1]))

# Write results to disk
