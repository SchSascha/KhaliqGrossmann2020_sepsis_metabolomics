library(data.table)

source("../function_definitions.R")

out_dir <- "../../data/measurements/"

human_data <- get_human_sepsis_data()

human_legend <- get_human_sepsis_legend()

human_data <- subset(human_data, TRUE, c(1:6, which(colnames(human_data) %in% human_legend[human_legend$group == "acylcarnitine",1])))

human_sepsis_data_s <- subset(human_data, grepl(pattern = "Septic", x = Group) & Day %in% 0:3)
ms <- aggregate(human_sepsis_data_s[, -1:-6], by = list(group = interaction(human_sepsis_data_s$Group, human_sepsis_data_s$Day)), mean)
ms$group <- substring(ms$group, 8)
ms[c("group", "Day")] <- Reduce("rbind", strsplit(x = ms$group, split = ".", fixed = TRUE))
l <- lapply(unique(ms$Day), function(d){ s <- subset(ms, Day == d, -which(colnames(ms) %in% c("group", "Day"))); return(s[1, ] / (s[2, ] + s[1, ])) })
l <- Reduce("rbind", l)
l$Day <- (0:3 + 1) * 24 * 60

fwrite(x = l, file = paste0(out_dir, "human_van_Eunen_ac_ratios_time_in_minutes.csv"))

cns <- colnames(human_data)[-1:-6]
cns <- sub(pattern = "C0", replacement = "CarMAT", x = cns)
cns[cns != "CarMAT"] <- paste0(cns[cns != "CarMAT"], "AcylCarMAT")
colnames(human_data)[-1:-6] <- cns

human_data$Day <- human_data$Day + 1
human_data$Day <- human_data$Day * 24 * 60
  
fwrite(x = human_data, file = paste0(out_dir, "human_van_Eunen_acylcarnitines_time_in_minutes.csv"))