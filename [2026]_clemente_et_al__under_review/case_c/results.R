
### Supplementary - Simulation study 7 --- varying m (temporal locs) parameter ------------
library(Matrix)
rmse = data.frame(rmse=vector(mode="numeric"),
                  method=vector(mode="character"), M = vector(mode="character"))

reactions = matrix(nrow = 30, ncol=1)

dir_ = "case_c/"
M = c("5", "11", "21", "41")
for (j in 1:length(M)){
  for(i in 1:30){   0
    outdir = paste0(dir_, "M_", M[j], "/1.00/250/" , i-1, "/")
    reactions[i,1] = readMM(paste0(outdir, "cv_optim.mtx"))[2]
    
    err_nonlin = readMM(paste0(outdir, "rmse_iterative.mtx"))[1]
    err_par = readMM(paste0(outdir, "rmse_diff_kfold.mtx"))[1]
    
    rmse = rbind(rmse, data.frame(rmse=err_nonlin, method="Nonlinear", M = M[j]))
    rmse = rbind(rmse, data.frame(rmse= err_par, method= "Linear", M = M[j] ))
  }
}
head(rmse)

myblue = rgb(23,55,94, maxColorValue = 255)
myred = colorspace::lighten(myblue, amount = 0.35)
mygreen = colorspace::lighten(myblue, amount = 0.9)

rmse$method = factor(rmse$method, levels = c("Nonlinear", "Linear"))
rmse$M = factor(rmse$M, levels = M)

pdf(paste0(dir_, "case_c_boxplots.pdf"), width=9)
par(xpd=TRUE)
boxplot(rmse ~  method + M, data=rmse, ylim = c(0, 0.02),
          col=c(myblue,myred), yaxt="n", xaxt="n", ylab="", xlab="", cex.axis=2.5)
axis(2, cex.axis=2.5)
axis(1, at = seq(1.5,7.5,by=2), labels = M, cex.axis=2.5, padj = 0.5)
dev.off()

library(dplyr)
rmse2 = rmse %>% filter(method == "Nonlinear")

# rmse2 = rmse2[,c(1,3)]
# head(rmse2)

pdf(paste0(dir_, "case_c_boxplots2.pdf"), width=9)
par(xpd=TRUE)
boxplot(rmse ~ M, data=rmse2,
          col=c(myblue), yaxt="n", xaxt="n", ylab="", xlab="", cex.axis=2.5)
axis(2, cex.axis=2.5)
axis(1, at = seq(1.5,9.5,by=2), labels = n_locs, cex.axis=2.5, padj = 0.5)
dev.off()
