library(data.table)
library(pathview)
library(missRanger)

source("../function_definitions.R")

out_dir <- "../../results/data_stats/"

load(file = "../../results/data_stats/ANOVA_complete_res.RData")

kegg_map <- fread(input = "../../data/id_maps/KEGG_map_all_proximate.csv", data.table = FALSE)

human_data <- get_human_sepsis_data()

#Clean data
kegg_map <- subset(kegg_map, KEGG_ID != "" | Proximate_KEGG_ID != "")
kegg_map <- subset(kegg_map, KEGG_ID != "NA" & Proximate_KEGG_ID != "NA")
colnames(kegg_map)[3] <- "Compound ID"

colns <- colnames(human_data)
colnames(human_data) <- make.names(colnames(human_data))
human_data[,-1:-5] <- missRanger(data = human_data[,-1:-5])
colnames(human_data) <- colns
human_data <- human_data[, c(1:5, which(colnames(human_data) %in% kegg_map$Metabolite))]
human_data <- subset(human_data, Day == 0)

#Build compound expression tables
kmd <- t(human_data[, -1:-5])
group <- rep("cont", nrow(human_data))
group[human_data$`CAP / FP` != "-"] <- "exp"
colnames(kmd) <- group

kegg_map_c <- na.omit(data.frame(kegg_map[, 3], kmd[match(kegg_map$Metabolite, rownames(kmd)), ], stringsAsFactors = FALSE))
colnames(kegg_map_c) <- substr(colnames(kegg_map_c), 1, 3)
kegg_map_c <- kegg_map_c[, c(1, 1 + order(colnames(kegg_map_c)[-1]))]
group <- colnames(kegg_map_c)[-1]
colnames(kegg_map_c) <- substr(colnames(kegg_map_c), 1, 3)
group[grep("con.*", group)] <- paste0("con", 1:sum(grepl("con.*", group)))
group[grep("exp.*", group)] <- paste0("exp", 1:sum(grepl("exp.*", group)))
colnames(kegg_map_c)[-1] <- group

kmc <- kegg_map_c
kmc <- kmc[!duplicated(kmc$k), ]
rownames(kmc) <- kmc[, 1]
kmc <- kmc[, -1]
kmcl <- t(tanh(scale(t(kmc))))
kmcm <- cbind(rowMeans(kmcl[, substr(colnames(kmcl), 1, 3) == "con"]), rowMeans(kmcl[, substr(colnames(kmcl), 1, 3) == "exp"]))
colnames(kmcm) <- c("con", "exp")
pathview(cpd.data = kmcm, pathway.id = "01100", species = "hsa", cpd.idtype = "kegg", match.data = FALSE, out.suffix = "control")
pathview(cpd.data = kmcm, pathway.id = "01230", species = "hsa", cpd.idtype = "kegg", match.data = FALSE, out.suffix = "control")
pathview(cpd.data = kmcm, pathway.id = "00330", species = "hsa", cpd.idtype = "kegg", match.data = FALSE, out.suffix = "control")

kmd <- t(human_data[human_data$`CAP / FP` != "-", -1:-5])
group <- human_data$Survival[human_data$`CAP / FP` != "-"]
colnames(kmd) <- group

kegg_map_s <- na.omit(data.frame(kegg_map[, 3], kmd[match(kegg_map$Metabolite, rownames(kmd)), ]), stringsAsFactors = FALSE)
colnames(kegg_map_s) <- substr(colnames(kegg_map_s), 1, 1)
kegg_map_s <- kegg_map_s[, c(1, 1 + order(colnames(kegg_map_s)[-1]))]
group <- colnames(kegg_map_s)[-1]
colnames(kegg_map_s) <- substr(colnames(kegg_map_s), 1, 3)
group[grep("N", group)] <- paste0("N", 1:sum(grepl("N", group)))
group[grep("S", group)] <- paste0("S", 1:sum(grepl("S", group)))
colnames(kegg_map_s)[-1] <- group

kms <- kegg_map_s
kms <- kms[!duplicated(kms$k), ]
rownames(kms) <- kms[, 1]
kms <- kms[, -1]
kmsl <- t(tanh(scale(t(kms))))
kmsm <- cbind(rowMeans(kmsl[, substr(colnames(kmsl), 1, 1) == "N"]), rowMeans(kmsl[, substr(colnames(kmsl), 1, 1) == "S"]))
colnames(kmsm) <- c("N", "S")
pathview(cpd.data = kmsm, pathway.id = "01100", species = "hsa", cpd.idtype = "kegg", match.data = FALSE, out.suffix = "survival")
pathview(cpd.data = kmsm, pathway.id = "01230", species = "hsa", cpd.idtype = "kegg", match.data = FALSE, out.suffix = "survival")
pathview(cpd.data = kmsm, pathway.id = "00330", species = "hsa", cpd.idtype = "kegg", match.data = FALSE, out.suffix = "survival")

#Write to disk

write.table(kegg_map_c, file = paste0(out_dir, "kegg_mapped_mets_control.csv"), sep = ",", quote = FALSE, row.names = FALSE)
write.table(kegg_map_s, file = paste0(out_dir, "kegg_mapped_mets_survival.csv"), sep = ",", quote = FALSE, row.names = FALSE)
