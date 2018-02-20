library(data.table)
library(stringi)
library(MetaboAnalystR)

source("../function_definitions.R")

data <- fread(input = '../../data/id_maps/hmdb_metabolites.xml', sep = "\n", data.table = FALSE, skip = 1, strip.white = TRUE, blank.lines.skip = TRUE, colClasses = c("character"))
data <- data[[1]]

met_start_line <- which(stri_detect_fixed(pattern = "<metabolite>", str = data))
met_stop_line <- which(stri_detect_fixed(pattern = "</metabolite>", str = data))
met_range <- cbind(met_start_line, met_stop_line)

met_name_extr_fun <- function(met_idx){
  td <- data[met_range[met_idx, 1]:met_range[met_idx, 2]]
  name_lines <- which(stri_detect_fixed(pattern = "<name>", str = td))
  name <- stri_replace_all_fixed(str = td[name_lines[1]], pattern = "<name>", replacement = "") #metabolite name is always the first <name> entry
  name <- stri_replace_all_fixed(str = name, pattern = "</name>", replacement = "")
  synonym_start_lines <- which(stri_detect_fixed(pattern = "<synonyms>", str = td))
  synonym_stop_lines <- which(stri_detect_fixed(pattern = "</synonyms>", str = td))
  synonym_range <- cbind(synonym_start_lines + 1, synonym_stop_lines - 1)
  x <- synonym_range
  synonyms <- lapply(1:nrow(synonym_range), function(i) { gsub(x = td[x[i,1]:x[i,2]], pattern = "<synonym>|</synonym>", replacement = "") })
  return(c(name, unlist(synonyms)))
}

mnames <- lapply(1:nrow(met_range), met_name_extr_fun)

mnames_s <- unlist(mnames)
mnames_l <- lapply(mnames, length)
mnames_i <- rep.int(seq_along(mnames_l), unlist(mnames_l))

human_sepsis_mets <- get_human_sepsis_legend()

gsm_mets_list <- list()
gsm_mets_list[["Recon2.04"]] <- fread(input = "../../data/template_models/Recon2.v04.mat_/recon2_met_names.txt", sep = "\n", header = FALSE, blank.lines.skip = TRUE)
gsm_mets_list[["iHsa"]] <- fread(input = "../../data/template_models/iHsa_iRno/iHsa_met_names.txt", sep = "\n", header = FALSE)

gsm_S_list <- list()
gsm_S_list[["Recon2.04"]] <- fread(input = "../../data/template_models/Recon2.v04.mat_/recon2_stoich_mat.txt", header = FALSE, sep = "\t")
gsm_S_list[["iHsa"]] <- fread(input = "../../data/template_models/iHsa_iRno/iHsa_stoich_mat.txt")

met_present_in_hmdb <- list()
met_present_in_hmdb[["Human Sepsis"]] <- match(human_sepsis_mets$name, mnames_s) 
met_present_in_hmdb[["Human Sepsis"]][is.na(met_present_in_hmdb[["Human Sepsis"]])] <- match(human_sepsis_mets$group[is.na(met_present_in_hmdb[["Human Sepsis"]])], mnames_s)
for (model in names(gsm_mets_list)){
  mets <- gsm_mets_list[[model]]
  met_present_in_hmdb[[model]] <- match(mets, mnames_s)
}

###Match metabs with MetaboAnalyst

human_sepsis_legend <- get_human_sepsis_legend()
human_sepsis_data <- get_human_sepsis_data()
human_sepsis_data <- human_sepsis_data[human_sepsis_data$`CAP / FP` != "-",]
write.csv(x = subset(human_sepsis_data, Day == 0, c(-1:-3,-5)), file = "metaboAnalyst_sepsis_temp.csv")

#Init
mSet<-InitDataObjects("conc", "pathinteg", FALSE)
mSet<-Read.TextData(mSet, "metaboAnalyst_sepsis_temp.csv", "rowu", "disc") #might fail bc. of path requirement

#Set up mSetObj with the list of compounds
mSet1<-Setup.MapData(mSet, human_sepsis_legend[-1:-5,1])
mSet2<-Setup.MapData(mSet, human_sepsis_legend[-1:-5,2])

mSet1<-CrossReferencing(mSet1, "name")
mSet1<-CreateMappingResultTable(mSet1)

mSet2<-CrossReferencing(mSet2, "name")
mSet2<-CreateMappingResultTable(mSet2)

sum(!is.na(mSet1$name.map$hit.values) | !is.na(mSet2$name.map$hit.values))