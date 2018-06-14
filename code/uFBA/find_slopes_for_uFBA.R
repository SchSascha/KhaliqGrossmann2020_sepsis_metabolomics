# Compute slopes and confidence intervals for rat metabolites


# Read rat data


# Remove outlier samples


# Fit lm
formula <- concentration ~ time

## TODO: find elegant way to index the rat data; maybe unique?

for (mat in c("liver", "heart")){
	for (grp in unique(rat_sepsis_data$group)){
		for (met in colnames(rat_sepsis_data)[metab_sel])){
			subdat <- subset(rat_sepsis_data, group == grp & 
			lm_list[[mat]][[met]] <- lm(formula, data = )

}

