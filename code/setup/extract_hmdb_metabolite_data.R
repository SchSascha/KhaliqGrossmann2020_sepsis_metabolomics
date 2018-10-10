library(data.table)
library(stringi)
library(MetaboAnalystR)
library(matrixStats)

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
gsm_mets_list[["iHsa"]] <- fread(input = "../../data/template_models/iHsa_mets.txt", colClasses = "character")
gsm_mets_list[["iRno"]] <- fread(input = "../../data/template_models/iRno_mets.txt", colClasses = "character")
gsm_mets_list[["Recon3D"]] <- fread(input = "../../data/template_models/Recon3D_mets.txt", colClasses = "character")

gsm_mets_list$Recon2.04 <- data.frame(Name = unlist(gsm_mets_list$Recon2.04))
gsm_mets_list$Recon2.04[["CHEBI ID"]] <- ""
gsm_mets_list$Recon2.04[["KEGG ID"]] <- ""
gsm_mets_list$Recon2.04[["PUBCHEM ID"]] <- ""

recon3D_hmdb_ids <- read.table(file = "../../data/template_models/Recon3D_HMDB_IDs.txt", blank.lines.skip = FALSE, sep = "\n", header = TRUE, stringsAsFactors = FALSE)

gsm_S_list <- list()
gsm_S_list[["Recon2.04"]] <- fread(input = "../../data/template_models/Recon2.v04.mat_/recon2_stoich_mat.txt", header = FALSE, sep = ",", colClasses = "numeric")
gsm_S_list[["iHsa"]] <- fread(input = "../../data/template_models/iHsa_stoich_mat.txt", header = FALSE, sep = "\t", colClasses = "numeric")
gsm_S_list[["iRno"]] <- fread(input = "../../data/template_models/iRno_stoich_mat.txt", header = FALSE, sep = "\t", colClasses = "numeric")
gsm_S_list[["Recon3D"]] <- fread(input = "../../data/template_models/Recon3D_stoich_mat.txt", header = FALSE, sep = "\t", colClasses = "numeric")

met_present_in_hmdb <- list()
met_present_in_hmdb[["Human Sepsis"]] <- match(human_sepsis_mets$name, mnames_s) 
met_present_in_hmdb[["Human Sepsis"]][is.na(met_present_in_hmdb[["Human Sepsis"]])] <- match(human_sepsis_mets$group[is.na(met_present_in_hmdb[["Human Sepsis"]])], mnames_s)
for (model in names(gsm_mets_list)){
  mets <- unlist(gsm_mets_list[[model]][,1])
  met_present_in_hmdb[[model]] <- match(mets, mnames_s)
}

###Match metabs with MetaboAnalyst

human_sepsis_legend <- get_human_sepsis_legend()
human_sepsis_data <- get_human_sepsis_data()
human_sepsis_data <- human_sepsis_data[human_sepsis_data$`CAP / FP` != "-",]
pheno_col_idx <- (which(colnames(human_sepsis_data) == "LCA") + 1):ncol(human_sepsis_data)
write_set <- subset(human_sepsis_data, Day == 0, -c(1:3, 5, pheno_col_idx))
write.csv(x = write_set, file = "metaboAnalyst_sepsis_temp.csv")

write.csv(x = gsm_mets_list$iHsa$Name, file = "metaboAnalyst_iHsa_temp.csv")
write.csv(x = gsm_mets_list$Recon2.04$Name, file = "metaboAnalyst_Recon2.04_temp.csv")
write.csv(x = gsm_mets_list$Recon3D$Name, file = "metaboAnalyst_Recon3D_temp.csv")

#Init
rm("mSet", "mSet1", "mSet2")
mSet<-InitDataObjects("conc", "pathinteg", FALSE)
mSet<-Read.TextData(mSet, "metaboAnalyst_sepsis_temp.csv", "rowu", "disc") #might fail bc. of path requirement

#Set up mSetObj with the list of compounds
##first clean compound list
pheno_legend_idx <- (which(human_sepsis_legend[,1] == "H1") + 1):nrow(human_sepsis_legend)
comp_i_list <- stri_trim_both(human_sepsis_legend[-c(1:5, pheno_legend_idx),1])
comp_n_list <- stri_trim_both(human_sepsis_legend[-c(1:5, pheno_legend_idx),2])
comp_n_list[comp_n_list == "Decenoylcarnitine"] <- "9-decenoylcarnitine"
comp_n_list[comp_n_list == "Dodecenoylcarnitine"] <- "trans-2-Dodecenoylcarnitine"
comp_n_list[comp_n_list == "Tetradecenoylcarnitine"] <- "trans-2-Tetradecenoylcarnitine"
comp_n_list[comp_n_list == "Hydroxytetradecenoylcarnitine"] <- "3-Hydroxytetradecenoylcarnitine"
comp_n_list[comp_n_list == "Tetradecadienylcarnitine"] <- "3, 5-Tetradecadiencarnitine"
comp_n_list[comp_n_list == "Hydroxytetradecadienylcarnitine"] <- "3-Hydroxy-5, 8-tetradecadiencarnitine"
comp_n_list[comp_n_list == "Hexadecenoylcarnitine"] <- "9-Hexadecenoylcarnitine"
comp_n_list[comp_n_list == "Hydroxyhexadecenoylcarnitine"] <- "3-Hydroxyhexadecenoylcarnitine"
comp_n_list[comp_n_list == "Hexadecadienylcarnitine"] <- "9,12-Hexadecadienylcarnitine"
comp_n_list[comp_n_list == "Hydroxyhexadecadienylcarnitine"] <- "3-Hydroxyhexadecadienoylcarnitine"
comp_n_list[comp_n_list == "Hydroxyhexadecanoylcarnitine"] <- "3-Hydroxyhexadecanoylcarnitine" #check in HMDB
comp_n_list[comp_n_list == "Octadecenoylcarnitine"] <- "11(Z)-Octadecenoylcarnitine"
comp_n_list[comp_n_list == "Hydroxyoctadecenoylcarnitine"] <- "3-Hydroxy-11(Z)-octadecenoylcarnitine"
comp_n_list[comp_n_list == "Octadecadienylcarnitine"] <- "9,12-Octadecadienylcarnitine"
comp_n_list[comp_n_list == "Hydroxyvalerylcarnitine (Methylmalonylcarnitine)"] <- "Hydroxyvalerylcarnitine" #Methylmalonylcarnitine is another compund in HMDB
comp_n_list[comp_n_list == "Hexanoylcarnitine (Fumarylcarnitine)"] <- "L-Hexanoylcarnitine" #first name appears two times in HMDB, second none
comp_n_list[comp_n_list == "Glutarylcarnitine (Hydroxyhexanoylcarnitine)"] <- "Glutarylcarnitine" #second name not in HMDB
comp_n_list[comp_n_list == "Methylglutarylcarnitine"] <- "3-Methylglutarylcarnitine" #appears two times in HMDB
comp_n_list[comp_n_list == "Hexenoylcarnitine"] <- "2-Hexenoylcarnitine"
comp_n_list[comp_n_list == "Acetylornithine"] <- "Acetyl-Ornithine"
comp_n_list[comp_n_list == "Methioninesulfoxide"] <- "Methionine sulfoxide"
#comp_n_list[comp_n_list == "Total dimethylarginine"] <- NULL #Symmetric DMA is the dominant part of the total and there is only Asymmetric DMA left, but the mesaured SDMA is sometimes above the total DMA
#Switching to short names because of better matches
comp_i_list[comp_i_list == "lysoPC a C14:0"] <- "LysoPC(14:0)"
comp_i_list[comp_i_list == "lysoPC a C18:2"] <- "LysoPC(18:2)"
#comp_i_list <- comp_i_list[comp_i_list != "PC aa C24:0"] #Not in HMDB
#comp_i_list <- comp_i_list[comp_i_list != "PC aa C26:0"] #Not in HMDB
pcaas <- grep(pattern = "PC aa C", x = comp_i_list)
comp_i_list[pcaas] <- stri_replace_all_regex(str = comp_i_list[pcaas], pattern = "PC aa C(.+)", replacement = "PC($1)")
comp_n_list[comp_n_list == "Sphingomyeline C16:0"] <- "C16 Sphingomyelin"
comp_n_list[comp_n_list == "Sphingomyeline C18:0"] <- "C18 Sphingomyelin"
comp_n_list[comp_n_list == "Sphingomyeline C24:0"] <- "Sphingomyelin (d18:0/24:0)"
comp_n_list[comp_n_list == "Sphingomyeline C26:0"] <- "SM(d18:1/26:0)"
sms <- grep(pattern = "Sphingomyeline C", x = comp_n_list)
comp_n_list[sms] <- stri_replace_all_regex(str = comp_n_list[sms], pattern = "^Sphingomyeline C(.+)", replacement = "C$1 Sphingomyelin")
comp_n_list[comp_n_list == "Glycochenodeoxycholic acid"] <- "Chenodeoxycholic acid glycine conjugate"
comp_i_list[comp_i_list == "H1"] <- "Glucose"

mSet1<-Setup.MapData(mSet, comp_i_list)
mSet2<-Setup.MapData(mSet, comp_n_list)

mSet1<-CrossReferencing(mSet1, "name")
mSet1<-CreateMappingResultTable(mSet1)

mSet2<-CrossReferencing(mSet2, "name")
mSet2<-CreateMappingResultTable(mSet2)

#mapping consistency check
match_hmdb_names <- cbind(mSet1$dataSet$map.table[,2], mSet2$dataSet$map.table[,2])
match_hmdb_names[,1] <- match_hmdb_names[match(mSet1$dataSet$map.table[,1], comp_i_list),1] #metaboAnalyst mapping reorders metabolites
match_hmdb_names[,2] <- match_hmdb_names[match(mSet2$dataSet$map.table[,1], comp_n_list),2] #metaboAnalyst mapping reorders metabolites
match_inconsist <- which(match_hmdb_names[,1] != match_hmdb_names[,2] & match_hmdb_names[,1] != "NA" & match_hmdb_names[,2] != "NA") #negligible inconsistencies

master_match_list <- c(match_hmdb_names[1:82, 2], match_hmdb_names[83:172, 1], match_hmdb_names[173:nrow(match_hmdb_names), 2])
master_match_list[188] <- match_hmdb_names[188,1]

#mapping overlap for models and measured metabs
match_hmdb_ids <- cbind(mSet1$dataSet$map.table[,3], mSet2$dataSet$map.table[,3])
hmdb_id_match_list <- c(match_hmdb_ids[1:82, 2], match_hmdb_ids[83:172, 1], match_hmdb_ids[173:nrow(match_hmdb_names), 2])
hmdb_id_match_list[188] <- match_hmdb_ids[188,1]
sum(hmdb_id_match_list %in% recon3D_hmdb_ids$HMDB.ID)

match_kegg_ids <- cbind(mSet1$dataSet$map.table[,5], mSet2$dataSet$map.table[,5])
kegg_id_match_list <- c(match_kegg_ids[1:82, 2], match_kegg_ids[83:172, 1], match_kegg_ids[173:nrow(match_hmdb_names), 2])
kegg_id_match_list[188] <- match_kegg_ids[188,1]
sum(kegg_id_match_list %in% gsm_mets_list$Recon3D$`KEGG ID`) #all phosphatidylcholine like match to just one compound, hence the high match count

##save KEGG ID mapping to disk
kegg_map_table <- data.frame(Metabolite = human_sepsis_legend[-c(1:5, pheno_legend_idx),1], KEGG_ID = kegg_id_match_list)
write.table(x = kegg_map_table, file = "../../data/id_maps/KEGG_map_all.csv", sep = "\t", quote = FALSE, row.names = FALSE)

match_pubchem_ids <- cbind(mSet1$dataSet$map.table[,4], mSet2$dataSet$map.table[,4])
pubchem_id_match_list <- c(match_pubchem_ids[1:82, 2], match_pubchem_ids[83:172, 1], match_pubchem_ids[173:nrow(match_hmdb_names), 2])
pubchem_id_match_list[188] <- match_pubchem_ids[188,1]
sum(pubchem_id_match_list %in% gsm_mets_list$Recon3D$`PUBCHEM ID`)

hir <- hmdb_id_match_list %in% recon3D_hmdb_ids$HMDB.ID
ker <- kegg_id_match_list %in% gsm_mets_list$Recon3D$`KEGG ID`
pir <- pubchem_id_match_list %in% gsm_mets_list$Recon3D$`PUBCHEM ID`
sum(hir | ker | pir)

write_set2 <- write_set
colnames(write_set2)[-1] <- master_match_list

rm("mSetFinal")
mSetFinal <- Setup.MapData(mSet, master_match_list)
mSetFinal <- SanityCheckData(mSetFinal)

mSetFinal<-CrossReferencing(mSetFinal, "name")
mSetFinal<-CreateMappingResultTable(mSetFinal)

map_table <- mSetFinal$dataSet$map.table
map_table <- map_table[match(map_table[,1], master_match_list)]

#Init for GSMM
rm("mSet_mod")
mSet_mod <- InitDataObjects("conc", "list", FALSE)
mSet_modlist <- list()

#Clean relevant metabolite names

#Match metabolite names to DB
for (name in names(gsm_mets_list)){
  mSet_modlist[[name]] <- mSet_mod
  mSet_modlist[[name]] <- Setup.MapData(mSet_modlist[[name]], gsm_mets_list[[name]]$Name)
  mSet_modlist[[name]] <- CrossReferencing(mSet_modlist[[name]], "name")
  mSet_modlist[[name]] <- CreateMappingResultTable(mSet_modlist[[name]])
}
