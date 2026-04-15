
### Supplementary - Simulation study 7 --- varying m (temporal locs) parameter ------------
library(Matrix)

time_meshes = list()
dt = list()
nonlinear_coeff = list()
linear_coeff = list()

for (j in 1:length(M)){
    outdir = paste0("input/M_", M[j], "/")
    time_meshes[[j]] = readMM(paste0(outdir, "time_mesh.mtx"))
    dt[[j]] = time_meshes[[j]][2] -time_meshes[[j]][1] 
    exact_coeff = as.matrix(readMM(paste0(outdir,"1.00/","fisher_kpp.mtx")))
    
    n_times = ncol(exact_coeff)
    ndof = nrow(exact_coeff)
    nonlinear_coeff[[j]] = matrix(0, nrow=(ndof * n_times), ncol=1)
    linear_coeff[[j]] = matrix(0, nrow=(ndof * n_times), ncol=1)
    
    for(k in 0:29){
      output_dir = paste0("case_c/M_", M[j], "/1.00/250/", k, "/")
      nonlinear_coeff[[j]] = nonlinear_coeff[[j]] +  as.matrix(readMM(paste0(output_dir, "estimate_iterative.mtx"))) / 30
      linear_coeff[[j]]    =  linear_coeff[[j]]    + as.matrix(readMM(paste0(output_dir, "estimate_diff_kfold.mtx"))) / 30
    }

    nonlinear_coeff[[j]] = matrix(nonlinear_coeff[[j]], nrow=ndof, ncol=n_times)
    linear_coeff[[j]] = matrix(linear_coeff[[j]], nrow=ndof, ncol=n_times)
}

time_meshes
dt

make_time_linear_interpolator <- function(U, times) {
  U <- as.matrix(U)
  times <- as.numeric(times)

  if (ncol(U) != length(times)) {
    stop("ncol(U) must equal length(times).")
  }
  if (length(times) < 2) {
    stop("Need at least two time points.")
  }
  if (any(diff(times) <= 0)) {
    stop("times must be strictly increasing.")
  }

  function(t) {
    t <- as.numeric(t)
    if (length(t) != 1) {
      stop("t must be a scalar.")
    }
    if (t < times[1] || t > times[length(times)]) {
      stop("t is outside the time range.")
    }

    if (t == times[length(times)]) {
      return(U[, length(times)])
    }

    k <- findInterval(t, times)
    k <- max(1, min(k, length(times) - 1))

    tk <- times[k]
    tk1 <- times[k + 1]
    alpha <- (t - tk) / (tk1 - tk)

    (1 - alpha) * U[, k] + alpha * U[, k + 1]
  }
}

interp = make_time_linear_interpolator(nonlinear_coeff[[1]], time_meshes[[1]])

library(fdaPDE)
mesh_dir = "input/mesh/"
nodes = readMM(paste0(mesh_dir,"points.mtx"))
cells = readMM(paste0(mesh_dir, "cells.mtx"))
boundary = readMM(paste0(mesh_dir, "boundary.mtx"))
mesh = create.mesh.2D(nodes=nodes, triangles = cells + 1)
FEMbasis = create.FEM.basis(mesh)

test_locs = as.matrix( readMM("input/M_41/1.00/test_locs.mtx") )
test_obs = as.matrix( readMM("input/M_41/1.00/test_obs.mtx") )
psi_test = fdaPDE:::CPP_get.psi.Matrix(FEMbasis, test_locs)
rmse = data.frame(rmse=vector(mode="numeric"),
                  method=vector(mode="character"), M = vector(mode="character"))

for(j in 1:length(M)){
  outdir = paste0("input/M_", M[j], "/")
  exact_coeff = as.matrix(readMM(paste0(outdir,"1.00/","fisher_kpp.mtx")))
  n_times = ncol(exact_coeff)
  ndof = nrow(exact_coeff)
  nonlinear_coeff = matrix(0, nrow=(ndof * n_times), ncol=1)
  linear_coeff = matrix(0, nrow=(ndof * n_times), ncol=1)
  for( k in 0:29){
      output_dir = paste0("case_c/M_", M[j], "/1.00/250/", k, "/")
      nonlinear_coeff = as.matrix(readMM(paste0(output_dir, "estimate_iterative.mtx"))) 
      linear_coeff    =  as.matrix(readMM(paste0(output_dir, "estimate_diff_kfold.mtx"))) 

      nonlinear_coeff = matrix(nonlinear_coeff, nrow=ndof, ncol=n_times)
      linear_coeff = matrix(linear_coeff, nrow=ndof, ncol=n_times)

      nonlin = make_time_linear_interpolator(nonlinear_coeff, time_meshes[[j]])
      lin = make_time_linear_interpolator(linear_coeff, time_meshes[[j]])
      rmse_nonlin = 0
      rmse_lin = 0
      for(t in 2:length(time_meshes[[4]])){
        rmse_nonlin = rmse_nonlin + sum( ( (psi_test %*% nonlin( time_meshes[[4]][t] )) - test_obs[,t-1])^2 )
        rmse_lin =    rmse_lin    + sum( ( (psi_test %*% lin( time_meshes[[4]][t] )) - test_obs[,t-1])^2 ) 
      }

      rmse_nonlin = sqrt( rmse_nonlin / (nrow(test_locs) * (length(time_meshes[[4]]) - 1)))
      rmse_lin = sqrt( rmse_lin / (nrow(test_locs) * (length(time_meshes[[4]]) - 1)))

      rmse = rbind(rmse, data.frame(rmse=rmse_nonlin, method="Nonlinear", M = M[j]))
      rmse = rbind(rmse, data.frame(rmse= rmse_lin, method= "Linear", M = M[j] ))
  }
}
head(rmse)

myblue = rgb(23,55,94, maxColorValue = 255)
myred = colorspace::lighten(myblue, amount = 0.35)

rmse$method = factor(rmse$method, levels = c("Nonlinear", "Linear"))
rmse$M = factor(rmse$M, levels = M)

pdf(paste0(dir_, "case_c_boxplots.pdf"), width=9)
par(xpd=TRUE)
boxplot(rmse ~  method + M, data=rmse,
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
axis(1, at = 1:4, labels = M, cex.axis=2.5, padj = 0.5)
dev.off()
