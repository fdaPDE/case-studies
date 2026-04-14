library(Matrix)
library(dplyr)
rmse = data.frame(rmse=vector(mode="numeric"),
                  method=vector(mode="character"), N = vector(mode="character"))

reactions = matrix(nrow = 30, ncol=1)

dir_ = "case_d/"
N = c("11", "21", "31", "41")
for (j in 1:length(N)){
  for(i in 1:30){   0
    outdir = paste0(dir_, "N_", N[j], "/1.00/250/" , i-1, "/")
    reactions[i,1] = readMM(paste0(outdir, "cv_optim.mtx"))[2]
    
    err_nonlin = readMM(paste0(outdir, "rmse_iterative.mtx"))[1]
    err_par = readMM(paste0(outdir, "rmse_diff_kfold.mtx"))[1]
    
    rmse = rbind(rmse, data.frame(rmse=err_nonlin, method="Nonlinear", N = N[j]))
    rmse = rbind(rmse, data.frame(rmse= err_par, method= "Linear", N = N[j] ))
  }
}
head(rmse)

myblue = rgb(23,55,94, maxColorValue = 255)
myred = colorspace::lighten(myblue, amount = 0.35)

rmse$method = factor(rmse$method, levels = c("Nonlinear", "Linear"))
rmse$N = factor(rmse$N, levels = N)

N_expr = c(expression(11^2),expression(21^2), expression(31^2), expression(41^2) )
pdf(paste0(dir_, "case_d_boxplots.pdf"), width=9)
par(xpd=TRUE)
boxplot(rmse ~  method + N, data=rmse,
          col=c(myblue,myred), yaxt="n", xaxt="n", ylab="", xlab="", cex.axis=2.5)
axis(2, cex.axis=2.5)
axis(1, at = seq(1.5,7.5,by=2), labels = N_expr, cex.axis=2.5, padj = 0.5)
dev.off()

rmse2 = rmse %>% filter(method != "Linear")

pdf(paste0(dir_, "case_d_boxplots2.pdf"), width=9)
par(xpd=TRUE)
boxplot(rmse ~  N, data=rmse2,
          col=c(myblue), yaxt="n", xaxt="n", ylab="", xlab="", cex.axis=2.5)
axis(2, cex.axis=2.5)
axis(1, at = 1:4, labels = N_expr, cex.axis=2.5, padj = 0.5)
dev.off()



pdf(paste0(dir_, "rmse-legend.pdf"), width = 14, height = 3)
plot.new()
legend(x="center",fill = c(myblue, myred),
       horiz = T,
       legend = unique(rmse$method),
       bty = "n", cex = 3)
dev.off()
