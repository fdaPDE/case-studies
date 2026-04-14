### simulation 2 -----------------------------------------------------------------------------------
library(Matrix)
library(dplyr)

rmse = data.frame(rmse=vector(mode="numeric"),
                  method=vector(mode="character"), n_locs = vector(mode="character"))

optims = matrix(nrow = 30, ncol=2)
values_linear = matrix(nrow=30, ncol=1)
values_nonlinear = values_linear

dir_ = "simulation_2/"
n_locs = c("250")
for (j in 1:length(n_locs)){
  for(i in 1:30){   0
    outdir = paste0(dir_, "1.00/", n_locs[j], "/" , i-1, "/")
    optims[i,] = as.vector( readMM(paste0(outdir, "cv_optim.mtx")) )
    
    err_nonlin = readMM(paste0(outdir, "rmse_iterative.mtx"))[1]
    err_par = readMM(paste0(outdir, "rmse_diff_kfold.mtx"))[1]
    err_tps = readMM(paste0(outdir, "rmse_tps.mtx"))[1]
    
    rmse = rbind(rmse, data.frame(rmse=err_nonlin, method="Nonlinear", n_locs = n_locs[j]))
    rmse = rbind(rmse, data.frame(rmse= err_par, method= "Linear", n_locs = n_locs[j] ))
    rmse = rbind(rmse, data.frame(rmse= err_tps, method= "TPS", n_locs = n_locs[j] ))
  }
}
head(rmse)
optims
myblue = rgb(23,55,94, maxColorValue = 255)
myred = colorspace::lighten(myblue, amount = 0.35)
mygreen = colorspace::lighten(myblue, amount = 0.9)

rmse$method = factor(rmse$method, levels = c("Nonlinear", "Linear", "TPS"))
rmse$n_locs = factor(rmse$n_locs, levels = n_locs)

pdf(paste0(dir_, "1.00/simulation_2_boxplots.pdf"), width=9)
par(xpd=TRUE)
boxplot(rmse ~  method + n_locs, data=rmse,
          col=c(myblue,myred, mygreen), yaxt="n", xaxt="n", ylab="", xlab="", cex.axis=2.5)
axis(2, cex.axis=2.5)
axis(1, at = c(1,2,3), labels = levels(rmse$method), cex.axis=2.5, padj = 0.5)
dev.off()

# optims[,2] -> the best r param is "0.00" 

rmse2 = rmse %>% filter(method != "TPS")
rmse2$method = factor(rmse2$method, levels = c("Nonlinear", "Linear"))
pdf(paste0(dir_, "1.00/simulation_2_boxplots__nonlin_vs_lin.pdf"), width=9)
par(xpd=TRUE)
boxplot(rmse ~  method + n_locs, data=rmse2,
          col=c(myblue,myred), yaxt="n", xaxt="n", ylab="", xlab="", cex.axis=2.5)
axis(2, cex.axis=2.5)
axis(1, at = c(1,2), labels = levels(rmse2$method), cex.axis=2.5, padj = 0.5)
dev.off()

# 
nonlinear = rmse %>% filter(method == "Nonlinear")
linear = rmse %>% filter(method == "Linear")
tps = rmse %>% filter(method == "TPS")

means = list(); sds = list()
means$nonlinear = tapply(nonlinear$rmse, nonlinear$n_locs, FUN=mean)
means$linear = tapply(linear$rmse, linear$n_locs, FUN=mean)
means$tps = tapply(tps$rmse, linear$n_locs, FUN=mean)

sds$nonlinear = tapply(nonlinear$rmse, nonlinear$n_locs, FUN=sd)
sds$linear = tapply(linear$rmse, linear$n_locs, FUN=sd)
sds$tps = tapply(tps$rmse, linear$n_locs, FUN=sd)

means
sds
