u_factor_comb = unique(rat_sepsis_data[c("material", "group", "time point")])
for (n in 1:nrow(u_factor_comb)){
  ufc_l[[n]] <- rat_sepsis_data$`Sample Identification`[
    rat_sepsis_data$material == u_factor_comb$material[n] 
    & rat_sepsis_data$group == u_factor_comb$group[n] 
    & rat_sepsis_data$`time point` == u_factor_comb$`time point`[n]]
  print(u_factor_comb[n,])
  print(ufc_l[[n]])
}

u_factor_comb$samples <- ufc_l
u_factor_comb <- u_factor_comb[order(u_factor_comb[,2]),]
u_factor_comb <- u_factor_comb[order(u_factor_comb[,1]),]

f = file(description = "rat_sample_ID_table_raw.txt", open = "wt")
for (n in 1:nrow(u_factor_comb)){
  writeLines(text = cbind(as.character(u_factor_comb[n,])), con = f, sep = "\n")
}
close(f)