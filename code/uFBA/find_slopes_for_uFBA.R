# Compute slopes and confidence intervals for rat metabolites

library(matrixStats)
library(reshape2)
library(data.table)
library(missRanger)

source("../function_definitions.R")

model_dir <- "../../data/template_models/"
out_dir <- "../../data/uFBA/"

# Create directories

if (!dir.exists(out_dir))
  dir.create(out_dir)

# Read data

rat_data <- get_rat_sepsis_data()
model_mets <- lapply(list.files(path = model_dir, pattern = "Biocrates\\.csv", full.names = TRUE), fread)
names(model_mets) <- sapply(strsplit(x = list.files(path = model_dir, pattern = "Biocrates\\.csv"), split = "_", fixed = TRUE), `[`, 1) # extract model names from file names

# Remove outlier samples

rat_data <- rat_data[!rat_data$`Sample Identification` %in% c("060H", "039L"), ]

# Impute missing values

set.seed(1)

cols <- colnames(rat_data)
colnames(rat_data) <- make.names(cols)
for (mat in unique(rat_data$material)){
  idat <- subset(rat_data, material == mat, 5:which(colnames(rat_data) == "LCA"))
  rat_data[rat_data$material == mat, 5:which(colnames(rat_data) == "LCA")] <- missRanger(idat)
}
idat <- subset(rat_data, material == "plasma", (which(colnames(rat_data) == "LCA")+1):length(cols))
rat_data[rat_data$material == "plasma", (which(colnames(rat_data) == "LCA")+1):length(cols)] <- missRanger(idat)
colnames(rat_data) <- cols

# Convert time point to numeric and change colname

rat_data$`time point` <- as.numeric(sub(pattern = "h", replacement = "", x = rat_data$`time point`))
colnames(rat_data)[colnames(rat_data) == "time point"] <- "time"

# Add up pool concentrations

rat_data[["Total PC aa"]] <- rowSums(rat_data[, grep(x = colnames(rat_data), pattern = "^PC aa")])
rat_data[["Total PC ae"]] <- rowSums(rat_data[, grep(x = colnames(rat_data), pattern = "^PC ae")])
rat_data[["Total LysoPC"]] <- rowSums(rat_data[, grep(x = colnames(rat_data), pattern = "^lysoPC")])
rat_data[["Total PC"]] <- rowSums(rat_data[, grep(x = colnames(rat_data), pattern = "^PC")])

# Transform physiological parameters

rat_data[["H"]] <- 10 ^ -rat_data[["pH"]]
rat_data[["Oxygen"]] <- 0.0225 * rat_data$pO2 + 1.4 * rat_data$Hb * rat_data$SaO2

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
	  for (times in tsv){
	    tstr <- paste0(times, collapse = "_")
	    lm_list[[mat]][[grp]][[tstr]] <- list()
	    grptimedat <- subset(grpdat, time %in% times)
	    grptimedat <- grptimedat[, !colAlls(is.na(grptimedat))]
	    met_set <- colnames(grptimedat)[5:ncol(grptimedat)]
		  for (met in met_set){
			  data <- subset(grptimedat, TRUE, c("time", met))
			  colnames(data)[2] <- "concentration"
			  lm_list[[mat]][[grp]][[tstr]][[met]] <- lm(formula = formula, data = data)
		  }
		}
  }
}
lm_summary <- lapply(lm_list, lapply, lapply, lapply, summary)
lm_confint <- lapply(lm_list, lapply, lapply, lapply, confint)
lm_rsqr <- lapply(lm_summary, lapply, lapply, lapply, function(r) r$r.squared)
lm_slope <- lapply(lm_summary, lapply, lapply, lapply, function(s) try(s$coefficients[2, 1]))
lm_slope <- lapply(lm_slope, lapply, lapply, function(e) e[!sapply(e, class) == "try-error"])
lm_slope_p <- lapply(lm_summary, lapply, lapply, lapply, function(s) try(s$coefficients[2, 4]))
lm_slope_p <- lapply(lm_slope_p, lapply, lapply, function(e) e[!sapply(e, class) == "try-error"])
lm_lci <- lapply(lm_confint, lapply, lapply, lapply, function(s) try(s[2, 1]))
lm_uci <- lapply(lm_confint, lapply, lapply, lapply, function(s) try(s[2, 2]))
slopes <- melt(lm_slope)
slo_p <- melt(lm_slope_p)
lcis <- melt(lm_lci)
ucis <- melt(lm_uci)
rsqrs <- melt(lm_rsqr)
colnames(slopes)[1] <- "slope"
colnames(lcis)[1] <- "lci"
colnames(ucis)[1] <- "uci"
colnames(slo_p)[1] <- "p"
sl <- merge(x = slopes, y = ucis, by = c("L1", "L2", "L3", "L4"))
sl <- merge(x = sl, y = lcis, by = c("L1", "L2", "L3", "L4"))
sl <- merge(x = sl, y = slo_p, by = c("L1", "L2", "L3", "L4"))
sl <- na.omit(sl)
sl$ignore <- sign(sl$lci) != sign(sl$uci) # ignore slope if signs of 95% CI limits do not match
sl$ignore[sl$lci == 0 | sl$uci == 0] <- FALSE # catch cases where one CI end == 0 (bc. sign([-1, 0, 1]) = [-1, 0, 1])
sl$ci <- abs(sl$slope - sl$lci) # matlab uFBA expects difference between slope and upper end of CI

# Compile metabolite lists for tissues
matv <- c("liver", "heart")
m <- 1
model_slopes <- list()
model_cilen <- list()
ignore_slope <- list()
for (m in seq_along(model_mets)){
  mod <- model_mets[[m]]
  cell_ind <- !mod[["Compartment"]] %in% c("e", "s")
  tab_cell <- table(mod[["Biocrates ID"]][cell_ind]) #e: extracellular in Recon3D; s: extracellular in iHsa/iRno
  tab_ex <- table(mod[["Biocrates ID"]][!cell_ind]) #e: extracellular in Recon3D; s: extracellular in iHsa/iRno
  combs <- subset(unique(sl[, 1:3]), L1 %in% matv)
  model_slopes[[m]] <- list()
  for (idx in 1:nrow(combs)){
    cb <- combs[idx, ]
    comb_str <- paste0(cb, collapse = " ")
    ds <- subset(sl, L1 == cb$L1 & L2 == cb$L2 & L3 == cb$L3)
    dext <- subset(sl, L1 == "plasma" & L2 == cb$L2 & L3 == cb$L3)
    metv_int <- ds$L4
    metv_ext <- dext$L4
    metma_int <- match(mod[["Biocrates ID"]][cell_ind], metv_int)
    metma_ext <- match(mod[["Biocrates ID"]][!cell_ind], metv_ext)
    model_slopes[[m]][[comb_str]] <- data.frame(Metabolite = mod[["Name"]], slope = 0, ci = 0, ignore = TRUE, stringsAsFactors = FALSE)
    # add tissue metabolites to cellular compartments
    model_slopes[[m]][[comb_str]][which(cell_ind), 2:4] <- ds[metma_int, c("slope", "ci", "ignore")]
    # add plasma metabolites to extracellular compartment
    model_slopes[[m]][[comb_str]][!cell_ind, 2:4] <- dext[metma_ext, c("slope", "ci", "ignore")]
    # divide slopes equally among metabolite instances - not necessary, we are working with concentrations!
    #model_slopes[[m]][[comb_str]][cell_ind, 2:3] <- model_slopes[[m]][[comb_str]][cell_ind, 2:3] / tab_cell[model_slopes[[m]][[comb_str]][cell_ind, "Metabolite"]] # divison works correctly
    #model_slopes[[m]][[comb_str]][!cell_ind, 2:3] <- model_slopes[[m]][[comb_str]][!cell_ind, 2:3] / tab_ex[model_slopes[[m]][[comb_str]][!cell_ind, "Metabolite"]]
  }
}
names(model_slopes) <- names(model_mets)

# Write
for (m in seq_along(model_slopes)){
  mod <- model_slopes[[m]]
  for (co in seq_along(mod)){
    fwrite(x = mod[[co]], file = paste0(out_dir, names(model_slopes)[m], "_", names(mod)[co], ".csv"), quote = FALSE, row.names = FALSE, sep = "\t")
  }
}
