library(CCA)
library(hexbin)
library(matrixStats)
library(data.table)
library(missRanger)
library(nscancor)

source("../function_definitions.R")

out_dir <- "../../results/data_stats/"

human_sepsis_data <- get_human_sepsis_data()

#Seperate septic and nonseptic patients
human_nonsepsis_data <- human_sepsis_data[human_sepsis_data$`CAP / FP` == "-", ]
human_sepsis_data <- human_sepsis_data[human_sepsis_data$`CAP / FP` != "-", ]

#human_sepsis_data <- na.omit(human_sepsis_data)
human_sepsis_data_normal <- human_sepsis_data
human_sepsis_data_normal[,-1:-5] <- scale(human_sepsis_data_normal[,-1:-5])

col <- colnames(human_sepsis_data_normal)
colnames(human_sepsis_data_normal) <- make.names(col)
human_sepsis_data_normal[, 6:204] <- missRanger(human_sepsis_data_normal[,6:204])
colnames(human_sepsis_data_normal) <- col

human_sepsis_data_metab <- subset(human_sepsis_data_normal, Day < 4, 6:204)
human_sepsis_data_pheno <- subset(human_sepsis_data_normal, Day < 4, 205:ncol(human_sepsis_data))

human_sepsis_data_normal <- human_sepsis_data_normal[,1:204]

#human_sepsis_data_normal <- cbind(human_sepsis_data_normal[,1:5], distance(x = human_sepsis_data_normal[, -1:-5], method = "chi.distance")) # distance() works on a per row basis

pat_list <- human_sepsis_data_normal[match(unique(human_sepsis_data_normal$Patient), human_sepsis_data_normal$Patient), 1:5]

pca_list <- lapply(pat_list$Patient, function(p){ prcomp(subset(human_sepsis_data_normal, Patient == p, -1:-5)) })
for (n in seq_along(pca_list)){
  pca_list[[n]]$surv <- pat_list$Survival[n]
  pca_list[[n]]$pat <- pat_list$Patient[n]
  pca_list[[n]]$capfp <- pat_list$`CAP / FP`[n]
  pca_list[[n]]$days <- unlist(subset(human_sepsis_data_normal, Patient == pat_list$Patient[n], "Day"))
}
pca_list[unlist(lapply(pca_list, function(e) ncol(e$x) <= 2))] <- NULL
png(filename = paste0(out_dir, "pca_sepsis_pats_expl_var.png"), width = 36, height = 16, units = "cm", res = 300)
par(mfrow = c(3,6))
lapply(pca_list, function(e) barplot(summary(e)$importance[2,], ylab = "Explained variance", main = paste(e$pat, e$surv, e$capfp, sep = ", ")))
dev.off()
png(filename = paste0(out_dir, "pca_sepsis_pats_PC1.png"), width = 36, height =16, units = "cm", res = 300)
par(mfrow = c(3,6))
lapply(pca_list, function(e){ plot(e$days, e$x[,1], main = paste(e$pat, e$surv, e$capfp, sep = ", "), type = "l", xlab = "Day", ylab = "PC1") })
dev.off()
png(filename = paste0(out_dir, "pca_sepsis_pats_PC2.png"), width = 36, height = 16, units = "cm", res = 300)
par(mfrow = c(3,6))
lapply(pca_list, function(e){ plot(e$days, e$x[,2], main = paste(e$pat, e$surv, e$capfp, sep = ", "), type = "l", xlab = "Day", ylab = "PC2") })
dev.off()
png(filename = paste0(out_dir, "pca_sepsis_pats_rot1.png"), width = 36, height = 16, units = "cm", res = 300)
par(mfrow = c(3,6))
lapply(pca_list, function(e){ plot(e$rotation[,1], type = "l", main = paste(e$pat, e$surv, e$capfp, sep = ", "), xlab = "Met. index", ylab = "Loading") })
dev.off()
png(filename = paste0(out_dir, "pca_sepsis_pats_PC1_vs_PC2.png"), width = 40, height = 20, units = "cm", res = 300)
par(mfrow = c(3,6))
lapply(pca_list, function(e){ plot(e$x[,1:min(2,ncol(e$x))], type = "p", col = 1:nrow(e$x), main = paste("Patient ", e$pat, ", ", e$surv, ", ", e$capfp, sep = "")); arrows(x0 = e$x[1:(nrow(e$x)-1),1], y0 = e$x[1:(nrow(e$x)-1),2], x1 = e$x[2:nrow(e$x),1], y1 = e$x[2:nrow(e$x),2], length = 0.1) })
dev.off()

rot_score <- sapply(pca_list, function(e){ rank(e$rotation[,1]) })

rot_pc1 <- t(sapply(pca_list, function(e) e$rotation[,1]))
rot_pca_pc1 <- prcomp(rot_pc1)
plot(rot_pca_pc1)
plot(rot_pca_pc1$rotation[,1])
rot_pc2 <- t(sapply(pca_list, function(e) e$rotation[,2]))
rot_pca_pc2 <- prcomp(rot_pc2)
plot(rot_pca_pc2)
plot(rot_pca_pc2$rotation[,1])
pca_con1 <- rot_pca_pc1$rotation[,1] %*% t(as.matrix(human_sepsis_data_normal[, -1:-5]))
pca_con2 <- rot_pca_pc2$rotation[,1] %*% t(as.matrix(human_sepsis_data_normal[, -1:-5]))

stab <- table(human_sepsis_data_normal$Survival)
human_sepsis_data_normal[human_sepsis_data_normal$Survival == "NS", -1:-5] <- sqrt(stab["NS"]) * human_sepsis_data_normal[human_sepsis_data_normal$Survival == "NS", -1:-5]
human_sepsis_data_normal[human_sepsis_data_normal$Survival == "S", -1:-5] <- sqrt(stab["S"]) * human_sepsis_data_normal[human_sepsis_data_normal$Survival == "S", -1:-5]
for (pat in unique(human_sepsis_data_normal$Patient)){
  sel <- human_sepsis_data_normal$Patient == pat
  human_sepsis_data_normal[sel, -1:-5] <- sqrt(sum(sel)) * human_sepsis_data_normal[sel, -1:-5]
}

#p <- prcomp(kernelMatrix(kernel = rbfdot(sigma = 0.0001), x = as.matrix(human_sepsis_data_normal[,-1:-5])))
p <- prcomp(human_sepsis_data_normal[, which(rowMeans(rot_score) > 140) + 5])
par(mfrow = c(1,1))
plot(summary(p)$importance[2,])
plot(p$x[,1:2], col = as.factor(human_sepsis_data_normal$Day))
col <- as.factor(pat_list$Survival)
for (n in seq_along(pat_list$Patient)){
  ind <- which(human_sepsis_data_normal$Patient == pat_list$Patient[n])
  x <- p$x[ind ,1]
  y <- p$x[ind ,2]
  arrows(x0 = x[-length(x)], y0 = y[-length(y)], x1 = x[-1], y1 = y[-1], length = 0.1, col = col[n])
}
#legend(x = 0.3, y = -0.25, legend = c("survivors", "nonsurvivors"), pch = "-", col = 2:1)
legend(x = 13, y = 8, legend = c("survivors", "nonsurvivors"), pch = "-", col = 2:1)

rot_pc12 <- rbind(rot_pc1, rot_pc2)
rot_pc12_sim <- as.matrix(dist(rot_pc12))
heatmaply(rot_pc12_sim, dendrogram = T, margins = c(50,50,0,50))

plot(rot_pca$sdev/sum(rot_pca$sdev))

linmod <- lm(rot_pc1 ~ 1)

par(mfrow=c(1,1))
#tic()
#er <- estim.regul(X = human_sepsis_data_metab, Y = human_sepsis_data_pheno, grid1 = seq(0.01, 1, length.out = 10), grid2 = seq(0.25, 1.5, length.out = 8))
#toc()

#rc <- rcc(X = human_sepsis_data_metab, Y = human_sepsis_data_pheno, lambda1 = er$lambda1, lambda2 = er$lambda2)

dfmax_w <- c(80, 30, 20)
ypredict <- function(x, yc, cc) {
  en <- glmnet::glmnet(x, yc, alpha = 0.5, intercept = FALSE,
                       dfmax = dfmax_w[cc], lower.limits = 0)
  W <- coef(en)
  return(W[2:nrow(W), ncol(W)])
}
dfmax_v <- c(21, 15, 15)
xpredict <- function(y, xc, cc) {
  en <- glmnet::glmnet(y, xc, alpha = 0.5, intercept = FALSE,
                       dfmax = dfmax_v[cc], lower.limits = 0)
  V <- coef(en)
  return(V[2:nrow(V), ncol(V)])
}

nsrc <- nscancor(x = human_sepsis_data_metab, y = human_sepsis_data_pheno, xpredict = xpredict, ypredict = ypredict, nvar = 3, xscale = T, yscale = T)

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
