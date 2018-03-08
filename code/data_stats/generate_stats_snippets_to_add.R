library(CCA)
library(hexbin)
library(matrixStats)
library(data.table)

source("../function_definitions.R")

human_sepsis_data <- get_human_sepsis_data()

#Seperate septic and nonseptic patients
human_nonsepsis_data <- human_sepsis_data[human_sepsis_data$`CAP / FP` == "-", ]
human_sepsis_data <- human_sepsis_data[human_sepsis_data$`CAP / FP` != "-", ]

human_sepsis_data <- na.omit(human_sepsis_data)
human_sepsis_data_normal <- human_sepsis_data
human_sepsis_data_normal[,-1:-5] <- scale(human_sepsis_data[,-1:-5])

human_sepsis_data_metab <- subset(human_sepsis_data_normal, Day < 4, 6:205)
human_sepsis_data_pheno <- subset(human_sepsis_data_normal, Day < 4, 206:ncol(human_sepsis_data))

er <- estim.regul(X = human_sepsis_data_metab, Y = human_sepsis_data_pheno, grid1 = seq(0.01, 3, length.out = 10), grid2 = seq(0.25, 1.5, length.out = 8))

rc <- rcc(X = human_sepsis_data_metab, Y = human_sepsis_data_pheno, lambda1 = er$lambda1, lambda2 = er$lambda2)

rc_sel <- which(rc$cor^2 > 0.5)
hist(rc$ycoef[, rc_sel])
plt.cc(rc)

#rc$scores$corr.Y.xscores

corr_dat <- list()
for (d in 0:2){
  NS_corr <- my.corr.test(x = subset(human_sepsis_data_normal, Survival == "NS" & Day == d, select = -1:-5), adjust = "fdr")
  S_corr <- my.corr.test(x = subset(human_sepsis_data_normal, Survival == "S" & Day == d, select = -1:-5), adjust = "fdr")
  corr_dat[[as.character(d)]] <- list(NS = NS_corr, S = S_corr)
}
hist((corr_dat$`0`$S$r - corr_dat$`0`$NS$r)[upper.tri(corr_dat$`0`$S$r)])
hist((corr_dat$`1`$S$r - corr_dat$`1`$NS$r)[upper.tri(corr_dat$`0`$S$r)])
hist((corr_dat$`2`$S$r - corr_dat$`2`$NS$r)[upper.tri(corr_dat$`0`$S$r)])

corr_diff <- sapply(corr_dat, function(e){ as.numeric((e$S$r - e$NS$r)[upper.tri(e$S$r)]) })
corr_S <- sapply(corr_dat, function(e){ as.numeric(e$S$r[upper.tri(e$S$r)]) })
corr_NS <- sapply(corr_dat, function(e){ as.numeric(e$NS$r[upper.tri(e$NS$r)]) })
corr_diff_mat <- corr_diff[!rowAnyNAs(corr_diff),]
corr_S_mat <- corr_S[!rowAnyNAs(corr_S) & !rowAnyNAs(corr_NS),]
corr_NS_mat <- corr_NS[!rowAnyNAs(corr_NS) & !rowAnyNAs(corr_S),]

pcrr_diff <- prcomp(corr_diff_mat)
plot(pcrr_diff$sdev/sum(pcrr_diff$sdev))
matplot(pcrr_diff$rotation, type = "l")

pcrr_S <- prcomp(corr_S_mat)
plot(pcrr_S$sdev/sum(pcrr_S$sdev))
matplot(pcrr_S$rotation, type = "l")

pcrr_NS <- prcomp(corr_NS_mat)
plot(pcrr_NS$sdev/sum(pcrr_NS$sdev))
matplot(pcrr_NS$rotation, type = "l")

corr_sel <- rowAnys(abs(corr_S_mat) > 0.75) | rowAnys(abs(corr_NS_mat) > 0.75)
matplot(t(corr_diff_mat[corr_sel,]), type = "l")

matplot(x = t(corr_S_mat), y = t(corr_NS_mat), type = "l")
corr_hist2D <- list()
for (n in 0:2){
  dat <- corr_dat[[paste0(n)]]
  dat_sel <- upper.tri(dat$S$r)
  corr_hist2D[[paste(n)]] <- hexbin(x = dat$S$r[dat_sel], y = dat$NS$r[dat_sel], xlab = "NS correlation", ylab = "S correlation", IDs = TRUE)
}
plot(corr_hist2D[[1]], main = "Day 0")
plot(-1:1, -1:1)
day0_dy <- aggregate(corr_S[,2], list(ID = corr_hist2D[[1]]@cID), mean)
day0_dx <- aggregate(corr_NS[,2], list(ID = corr_hist2D[[1]]@cID), mean)
arrows(corr_hist2D[[1]]@xcm, corr_hist2D[[1]]@ycm, day0_dx$x, day0_dy$x, length = 0.05)
plot(corr_hist2D[[2]], main = "Day 1")
plot(-1:1, -1:1)
day1_dy <- aggregate(corr_S[,3], list(ID = corr_hist2D[[2]]@cID), mean)
day1_dx <- aggregate(corr_NS[,3], list(ID = corr_hist2D[[2]]@cID), mean)
arrows(corr_hist2D[[2]]@xcm, corr_hist2D[[2]]@ycm, day1_dx$x, day1_dy$x, length = 0.05)

plot(-1:1, -1:1)
day0_d2_dy <- aggregate(corr_S[,3], list(ID = corr_hist2D[[1]]@cID), mean)
day0_d2_dx <- aggregate(corr_NS[,3], list(ID = corr_hist2D[[1]]@cID), mean)
arrows(corr_hist2D[[1]]@xcm, corr_hist2D[[1]]@ycm, day0_d2_dx$x, day0_d2_dy$x, length = 0.05)

plot(corr_hist2D[[3]], main = "Day 2")

slim <- c(-0.5,0)
nslim <- c(-0.4,0)
for (n in 1:3){
  print(sum(between(corr_S_mat[,n], slim[1], slim[2]) & between(corr_NS_mat[,n], nslim[1], nslim[2])))
}

slim <- c(0.6,1)
nslim <- c(0,0.5)
for (n in 1:3){
  print(sum(between(corr_S_mat[,n], slim[1], slim[2]) & between(corr_NS_mat[,n], nslim[1], nslim[2])))
}


plot(-1:1, -1:1)
arrows(corr_hist2D[[1]]@xcm, corr_hist2D[[1]]@ycm, corr_hist2D[[3]]@xcm, corr_hist2D[[3]]@ycm, length = 0.05)
