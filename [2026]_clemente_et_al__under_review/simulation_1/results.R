### simulation 1 -----------------------------------------------------------------------------------
library(Matrix)
library(dplyr)

rmse = data.frame(rmse=vector(mode="numeric"),
                  method=vector(mode="character"), n_locs = vector(mode="character"))

optims = matrix(nrow = 30, ncol=2)

dir_ = "simulation_1/"
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

pdf(paste0(dir_, "1.00/simulation_1_boxplots.pdf"), width=9)
par(xpd=TRUE)
boxplot(rmse ~  method + n_locs, data=rmse,
          col=c(myblue,myred, mygreen), yaxt="n", xaxt="n", ylab="", xlab="", cex.axis=2.5)
axis(2, cex.axis=2.5)
axis(1, at = c(1,2,3), labels = levels(rmse$method), cex.axis=2.5, padj = 0.5)
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

# qualitative
library(fdaPDE)
source("../utils/utils.R")
source("../utils/plot_smooth_2D.R")

mesh_dir = "input/mesh/"
nodes = readMM(paste0(mesh_dir,"points.mtx"))
cells = readMM(paste0(mesh_dir, "cells.mtx"))
boundary = readMM(paste0(mesh_dir, "boundary.mtx"))
mesh = create.mesh.2D(nodes=nodes, triangles = cells + 1)
FEMbasis = create.FEM.basis(mesh)

r_param = c("1.00")
n_breaks = 20
n_bins=10
lims = c(0,1)
grid <- expand.grid(
  x = seq(0, 1, length.out = 250),
  y = seq(0, 1, length.out = 250)
)
Psi = fdaPDE:::CPP_get.psi.Matrix(FEMbasis, as.matrix(grid))

library(jsonlite)
cm <- fromJSON("../utils/rainbow_unif.json") 
rgb_points <- matrix(cm$RGBPoints[[1]], ncol = 4, byrow = TRUE)
pv_cols <- rgb(rgb_points[, 2], rgb_points[, 3], rgb_points[, 4])
rainbow_uniform_pal <- colorRampPalette(pv_cols)

for(i in 1:length(r_param)){

   
    exact_coeff = as.matrix(readMM(paste0("input/", r_param[i],"/","fisher_kpp.mtx")))
    
    nonlinear_coeff = matrix(0, nrow=(nrow(exact_coeff)*ncol(exact_coeff)), ncol=1)
    linear_coeff = nonlinear_coeff
    tps_coeff = nonlinear_coeff
    
    for(j in 0:29){
      output_dir = paste0("simulation_1/", r_param[i], "/250/", j, "/")
      nonlinear_coeff = nonlinear_coeff +  as.matrix(readMM(paste0(output_dir, "estimate_iterative.mtx"))) / 30
      linear_coeff =  linear_coeff + as.matrix(readMM(paste0(output_dir, "estimate_diff_kfold.mtx"))) / 30
      tps_coeff =  tps_coeff + as.matrix(as.vector(readMM(paste0(output_dir, "estimate_tps.mtx")))) / 30
    }

    n_times = ncol(exact_coeff)
    nonlinear_coeff = matrix(nonlinear_coeff, nrow=nrow(exact_coeff), ncol=n_times)
    linear_coeff = matrix(linear_coeff, nrow=nrow(exact_coeff), ncol=n_times)
    tps_coeff = matrix(tps_coeff, nrow=nrow(exact_coeff), ncol=n_times)
    coeff_lims = range(exact_coeff)

    pdf(paste0("simulation_1/", r_param[i],"/","exact.pdf"))
    for( t in 1:n_times){
    FEM = FEM(exact_coeff[,t], FEMbasis)
    plt = plot_smooth_2D_ggplot2(FEM, coeff_lims = coeff_lims, colorscale = rainbow_uniform_pal)
    print(plt)
    }
    dev.off()

    pdf(paste0("simulation_1/", r_param[i],"/","nonlinear.pdf"))
    for( t in 1:n_times){
    FEM = FEM(nonlinear_coeff[,t], FEMbasis)
    plt = plot_smooth_2D_ggplot2(FEM, coeff_lims = coeff_lims, colorscale = rainbow_uniform_pal)
    print(plt)
    }
    dev.off()

    pdf(paste0("simulation_1/", r_param[i],"/","linear.pdf"))
    for( t in 1:n_times){
    FEM = FEM(linear_coeff[,t], FEMbasis)
    plt = plot_smooth_2D_ggplot2(FEM, coeff_lims = coeff_lims, colorscale = rainbow_uniform_pal)
    print(plt)
    }
    dev.off()

    pdf(paste0("simulation_1/", r_param[i],"/","tps.pdf"))
    for( t in 1:n_times){
    FEM = FEM(tps_coeff[,t], FEMbasis)
    plt = plot_smooth_2D_ggplot2(FEM, coeff_lims = coeff_lims, colorscale = rainbow_uniform_pal)
    print(plt)
    }
    dev.off()

    locs = as.matrix(readMM(paste0("input/", r_param[i], "/250/0/locs.mtx")))
    obs = as.matrix(readMM(paste0("input/", r_param[i], "/250/0/obs.mtx")))
    square <- data.frame( x = c(0, 1, 1, 0, 0, 0, 0), y = c(0, 0, 1, 1, 0, 1, 0))

    pdf(paste0("simulation_1/", r_param[i],"/","pointwise_data.pdf"))
    for( t in 1:n_times){
      plt = plot_pointwise_data_2D_ggplot2(locs, as.matrix(obs[,t]), coeff_lims = coeff_lims, 
                                            mesh = mesh,
                                            colorscale = rainbow_uniform_pal)
      plt = plt + geom_path(data = square, aes(x, y), color = "black", linewidth=2, 
                              inherit.aes = FALSE, linejoin = "round")
      print(plt)
      }
    dev.off()


}


##### --- SUPPLEMENTARY ----------------------------------------------------------------------------

# CASE A --- varying reaction coeff 
library(Matrix)
rmse = data.frame(rmse=vector(mode="numeric"),
                  method=vector(mode="character"), r = vector(mode="character"))

optims = matrix(nrow = 30, ncol=2)
dir_ = "case_a/"
r_params = c("0.20", "0.40", "0.60", "0.80", "1.00") 
for (j in 1:length(r_params)){
  for(i in 1:30){   0
    outdir = paste0(dir_, r_params[j], "/250/", i-1, "/")
    
    optims[i,] = as.vector( readMM(paste0(outdir, "cv_optim.mtx")) )
    
    err_nonlin = readMM(paste0(outdir, "rmse_iterative.mtx"))[1]
    err_par = readMM(paste0(outdir, "rmse_diff_kfold.mtx"))[1]
    
    rmse = rbind(rmse, data.frame(rmse=err_nonlin, method="Nonlinear", r = r_params[j]))
    rmse = rbind(rmse, data.frame(rmse= err_par, method= "Linear", r = r_params[j] ))
  }
}
head(rmse)

myblue = rgb(23,55,94, maxColorValue = 255)
myred = colorspace::lighten(myblue, amount = 0.35)
mygreen = colorspace::lighten(myblue, amount = 0.9)

rmse$method = factor(rmse$method, levels = c("Nonlinear", "Linear"))
rmse$r = factor(rmse$r, levels = r_params)

pdf(paste0(dir_, "case_a_boxplots.pdf"), width=9)
par(xpd=TRUE)
boxplot(rmse ~  method + r, data=rmse,
          col=c(myblue,myred), yaxt="n", xaxt="n", ylab="", xlab="", cex.axis=2.5)
axis(2, cex.axis=2.5)
axis(1, at = seq(1.5,9.5,by=2), labels = r_params, cex.axis=2.5, padj = 0.5) # ap
dev.off()



# CASE B --- varying n (spatial locs) parameter 
library(Matrix)
dir_ = "case_b/"
rmse = data.frame(rmse=vector(mode="numeric"),
                  method=vector(mode="character"), n_locs = vector(mode="character"))

optims = matrix(nrow = 30, ncol=2)

n_locs = c("125", "250", "500", "1000", "2000")
for (j in 1:length(n_locs)){
  for(i in 1:30){   0
    outdir = paste0(dir_, "1.00/", n_locs[j], "/" , i-1, "/")
    optims[i,] = as.vector( readMM(paste0(outdir, "cv_optim.mtx")) )
    
    err_nonlin = readMM(paste0(outdir, "rmse_iterative.mtx"))[1]
    err_par = readMM(paste0(outdir, "rmse_diff_kfold.mtx"))[1]
    
    rmse = rbind(rmse, data.frame(rmse=err_nonlin, method="Nonlinear", n_locs = n_locs[j]))
    rmse = rbind(rmse, data.frame(rmse= err_par, method= "Linear", n_locs = n_locs[j] ))
  }
}
head(rmse)
optims
myblue = rgb(23,55,94, maxColorValue = 255)
myred = colorspace::lighten(myblue, amount = 0.35)

rmse$method = factor(rmse$method, levels = c("Nonlinear", "Linear"))
rmse$n_locs = factor(rmse$n_locs, levels = n_locs)

pdf(paste0(dir_, "case_b_boxplots.pdf"), width=9)
par(xpd=TRUE)
boxplot(rmse ~  method + n_locs, data=rmse,
          col=c(myblue,myred), yaxt="n", xaxt="n", ylab="", xlab="", cex.axis=2.5)
axis(2, cex.axis=2.5)
axis(1, at = seq(1.5,9.5,by=2), labels = n_locs, cex.axis=2.5, padj = 0.5)
dev.off()

library(dplyr)
rmse2 = rmse %>% filter(method == "Nonlinear")

# rmse2 = rmse2[,c(1,3)]
# head(rmse2)

pdf(paste0(dir_, "case_b_boxplots2.pdf"), width=9)
par(xpd=TRUE)
boxplot(rmse ~ n_locs, data=rmse2,
          col=c(myblue), yaxt="n", xaxt="n", ylab="", xlab="", cex.axis=2.5)
axis(2, cex.axis=2.5)
axis(1, at = seq(1.5,9.5,by=2), labels = n_locs, cex.axis=2.5, padj = 0.5)
dev.off()


#pdf(paste0(dir_, "case_b_boxplots.pdf"), width=9)
nonlinear = rmse %>% filter(method == "Nonlinear")
linear = rmse %>% filter(method == "Linear")

nonlinear = nonlinear[,c(1,3)]
linear = linear[,c(1,3)]
means = list(); sds = list()
means$nonlinear = tapply(nonlinear$rmse, nonlinear$n_locs, FUN=mean)
means$linear = tapply(linear$rmse, linear$n_locs, FUN=mean)

sds$nonlinear = tapply(nonlinear$rmse, nonlinear$n_locs, FUN=sd)
sds$linear = tapply(linear$rmse, linear$n_locs, FUN=sd)

upper = list(); lower = list()
upper$nonlinear <- means$nonlinear + sds$nonlinear
upper$linear <- means$linear - sds$linear

lower$nonlinear <- means$nonlinear - sds$nonlinear
lower$linear <- means$linear + sds$linear

pdf(paste0(dir_, "case_b_means.pdf"), width=9)
par(xpd=TRUE)
plot(1:5, means$nonlinear, type="l", lwd=4, col = myblue, 
     ylim= range(c(means$nonlinear, means$linear)), 
     yaxt="n", xaxt="n", ylab="", xlab="")
points(1:5, means$linear, type="l", lwd=4, col = myred)
# SD stripe
polygon(
  x = c(1:5, rev(1:5)),
  y = c(upper$nonlinear, rev(lower$nonlinear)),
  col = adjustcolor(myblue, alpha.f = 0.7),   # transparent blue
  border = NA
)
polygon(
  x = c(1:5, rev(1:5)),
  y = c(upper$linear, rev(lower$linear)),
  col = adjustcolor(myred, alpha.f = 0.7),   # transparent blue
  border = NA
)
axis(2, cex.axis=2.5)
axis(1, at = 1:5, labels = n_locs, cex.axis=2.5, padj = 0.5)
dev.off()

# tratteggiato?!?!
pdf(paste0(dir_, "case_b_means_2.pdf"), width=9)
par(xpd=TRUE)
plot(1:5, means$nonlinear, type="l", lwd=4, col = myblue, 
     ylim= range(c(means$nonlinear, means$linear)), 
     yaxt="n", xaxt="n", ylab="", xlab="")
points(1:5, means$linear, type="l", lwd=4, col = myred)

axis(2, cex.axis=2.5)
axis(1, at = 1:5, labels = n_locs, cex.axis=2.5, padj = 0.5)
dev.off()