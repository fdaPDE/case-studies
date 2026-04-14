### simulation 3 -----------------------------------------------------------------------------------
library(Matrix)
library(dplyr)

rmse = data.frame(rmse=vector(mode="numeric"),
                  method=vector(mode="character"), n_locs = vector(mode="character"))

optims = matrix(nrow = 30, ncol=2)

dir_ = "simulation_3/"
for(i in 1:30){   
    outdir = paste0(dir_ , i-1, "/")
    optims[i,] = as.vector( readMM(paste0(outdir, "cv_optim.mtx")) )
    
    err_nonlin = readMM(paste0(outdir, "rmse_iterative.mtx"))[1]
    err_par = readMM(paste0(outdir, "rmse_diff_kfold.mtx"))[1]
    
    rmse = rbind(rmse, data.frame(rmse=err_nonlin, method="Nonlinear"))
    rmse = rbind(rmse, data.frame(rmse= err_par, method= "Linear" ))
  }

head(rmse)
optims
myblue = rgb(23,55,94, maxColorValue = 255)
myred = colorspace::lighten(myblue, amount = 0.35)

rmse$method = factor(rmse$method, levels = c("Nonlinear", "Linear"))

# optims[,2] -> the best r param is "0.00" 

pdf(paste0(dir_, "simulation_3_boxplots.pdf"), width=9)
par(xpd=TRUE)
boxplot(rmse ~  method, data=rmse,
          col=c(myblue,myred), yaxt="n", xaxt="n", ylab="", xlab="", cex.axis=2.5)
axis(2, cex.axis=2.5)
axis(1, at = c(1,2), labels = levels(rmse$method), cex.axis=2.5, padj = 0.5)
dev.off()

tapply(rmse$rmse, rmse$method, FUN=mean)
tapply(rmse$rmse, rmse$method, FUN=sd)

# qualitative
library(Matrix)
library(fdaPDE)
library(ggplot2)
library(jsonlite)
source("../utils/utils.R")
source("../utils/plot_smooth_2D.R")

incidence_matrix = read.csv("input/incidence_matrix.csv")[,-1]
observations = readMM("input/0/obs.mtx")
observations_nonoise = readMM("input/0/obs_no_noise.mtx")
not_unif = readMM("input/not_unif_idx.mtx")

input_dir = "../simulation_1/input/1.00/"
mesh_dir ="../simulation_1/input/mesh/" 
nodes = readMM(paste0(mesh_dir,"points.mtx"))
cells = readMM(paste0(mesh_dir, "cells.mtx"))
boundary = readMM(paste0(mesh_dir, "boundary.mtx"))

mesh = create.mesh.2D(nodes, triangles = (cells + 1)) 
plot(mesh, pch=".")
FEMbasis = create.FEM.basis(mesh)

{
  poly_coords = list()
  for(r in 1:nrow(incidence_matrix)){
    triangles = mesh$triangles[which(incidence_matrix[r,] == 1),]
    #triangles = cbind(triangles, triangles[,1])
    poly_coords[[r]] = list()
    for(t in 1:nrow(triangles)){
      poly_coords[[r]][[t]] = mesh$nodes[triangles[t,], ]
      poly_coords[[r]][[t]] = rbind(poly_coords[[r]][[t]],
                                    mesh$nodes[triangles[t,][1], ])
    }
  }
}

cm <- fromJSON("../utils/rainbow_unif.json") 
rgb_points <- matrix(cm$RGBPoints[[1]], ncol = 4, byrow = TRUE)
pv_cols <- rgb(rgb_points[, 2], rgb_points[, 3], rgb_points[, 4])
rainbow_uniform_pal <- colorRampPalette(pv_cols)

exact_coeff = as.matrix(readMM(paste0(input_dir, "fisher_kpp.mtx")))
coeff_lims = range(exact_coeff)  

    nonlinear_coeff = matrix(0, nrow=(nrow(exact_coeff)*ncol(exact_coeff)), ncol=1)
    linear_coeff = nonlinear_coeff
    
    for(j in 0:29){
      output_dir = paste0("simulation_3/", j, "/")
      nonlinear_coeff = nonlinear_coeff +  as.matrix(readMM(paste0(output_dir, "estimate_iterative.mtx"))) / 30
      linear_coeff =  linear_coeff + as.matrix(readMM(paste0(output_dir, "estimate_diff_kfold.mtx"))) / 30
    }

    n_times = ncol(exact_coeff)
    nonlinear_coeff = matrix(nonlinear_coeff, nrow=nrow(exact_coeff), ncol=n_times)
    linear_coeff = matrix(linear_coeff, nrow=nrow(exact_coeff), ncol=n_times)
    
    pdf(paste0("simulation_3/exact.pdf"))
    for( t in 1:n_times){
    FEM = FEM(exact_coeff[,t], FEMbasis)
    plt = plot_smooth_2D_ggplot2(FEM, coeff_lims = coeff_lims, colorscale = rainbow_uniform_pal)
    print(plt)
    }
    dev.off()

    pdf(paste0("simulation_3/nonlinear.pdf"))
    for( t in 1:n_times){
    FEM = FEM(nonlinear_coeff[,t], FEMbasis)
    plt = plot_smooth_2D_ggplot2(FEM, coeff_lims = coeff_lims, colorscale = rainbow_uniform_pal)
    print(plt)
    }
    dev.off()

    pdf(paste0("simulation_3/linear.pdf"))
    for( t in 1:n_times){
    FEM = FEM(linear_coeff[,t], FEMbasis)
    plt = plot_smooth_2D_ggplot2(FEM, coeff_lims = coeff_lims, colorscale = rainbow_uniform_pal)
    print(plt)
    }
    dev.off()

{
    pdf("simulation_3/areal_data.pdf")
    for (t in seq_len(ncol(observations))) {
      
      ## Build one big data.frame with all polygons for this time t
      poly_list <- list()
      idx <- 1
      
      for (r in seq_len(nrow(incidence_matrix))) {
        for (i in seq_along(poly_coords[[r]])) {
          
          coords <- poly_coords[[r]][[i]]
          # coords is assumed to be a matrix with 2 columns (x, y)
          
          poly_list[[idx]] <- data.frame(
            x     = coords[, 1],
            y     = coords[, 2],
            z     = observations[r, t],
            group = paste(r, i, sep = "_")  # group polygons
          )
          idx <- idx + 1
        }
      }
      
      poly_df <- do.call(rbind, poly_list)
      
      plt_rect <- ggplot(poly_df, aes(x = x, y = y, group = group, fill = z)) +
        geom_polygon(color = NA) +
        scale_fill_gradientn(colors = rainbow_uniform_pal(256), limits= coeff_lims, oob    = scales::squish) +
        coord_fixed() +
        theme_void() +
        theme(legend.position = "none") +
        geom_rect(
          data = data.frame(xmin = 0, ymin = 0, xmax = 1, ymax = 1),
          aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
          fill = "transparent",
          color = "black",
          linewidth = 1,
          inherit.aes = FALSE
        )
      
      print(plt_rect)
    }
    
    dev.off()
}


# --- sim 3 - NOT UNIF -----------------------------------------------------------------------------
rm(list=ls())
library(Matrix)
library(dplyr)

rmse = data.frame(rmse=vector(mode="numeric"),
                  method=vector(mode="character"), n_locs = vector(mode="character"))

optims = matrix(nrow = 30, ncol=2)

dir_ = "simulation_3-not_unif/"
for(i in 1:30){   
    outdir = paste0(dir_ , i-1, "/")
    optims[i,] = as.vector( readMM(paste0(outdir, "cv_optim.mtx")) )
    
    err_nonlin = readMM(paste0(outdir, "rmse_iterative.mtx"))[1]
    err_par = readMM(paste0(outdir, "rmse_diff_kfold.mtx"))[1]
    
    rmse = rbind(rmse, data.frame(rmse=err_nonlin, method="Nonlinear"))
    rmse = rbind(rmse, data.frame(rmse= err_par, method= "Linear" ))
  }

head(rmse)
optims
myblue = rgb(23,55,94, maxColorValue = 255)
myred = colorspace::lighten(myblue, amount = 0.35)

rmse$method = factor(rmse$method, levels = c("Nonlinear", "Linear"))

# optims[,2] -> the best r param is "0.00" 

pdf(paste0(dir_, "simulation_3-not_unif_boxplots.pdf"), width=9)
par(xpd=TRUE)
boxplot(rmse ~  method, data=rmse,
          col=c(myblue,myred), yaxt="n", xaxt="n", ylab="", xlab="", cex.axis=2.5)
axis(2, cex.axis=2.5)
axis(1, at = c(1,2), labels = levels(rmse$method), cex.axis=2.5, padj = 0.5)
dev.off()

tapply(rmse$rmse, rmse$method, FUN=mean)
tapply(rmse$rmse, rmse$method, FUN=sd)


# qualitative 
library(Matrix)
library(fdaPDE)
library(ggplot2)
library(jsonlite)
source("../utils/utils.R")
source("../utils/plot_smooth_2D.R")

incidence_matrix = read.csv("input/incidence_matrix.csv")[,-1]
observations = readMM("input/0/obs.mtx")
not_unif = readMM("input/not_unif_idx.mtx")

input_dir = "../simulation_1/input/1.00/"
mesh_dir ="../simulation_1/input/mesh/" 
nodes = readMM(paste0(mesh_dir,"points.mtx"))
cells = readMM(paste0(mesh_dir, "cells.mtx"))
boundary = readMM(paste0(mesh_dir, "boundary.mtx"))

mesh = create.mesh.2D(nodes, triangles = (cells + 1)) 
plot(mesh, pch=".")
FEMbasis = create.FEM.basis(mesh)

{
  poly_coords = list()
  for(r in 1:nrow(incidence_matrix)){
    triangles = mesh$triangles[which(incidence_matrix[r,] == 1),]
    #triangles = cbind(triangles, triangles[,1])
    poly_coords[[r]] = list()
    for(t in 1:nrow(triangles)){
      poly_coords[[r]][[t]] = mesh$nodes[triangles[t,], ]
      poly_coords[[r]][[t]] = rbind(poly_coords[[r]][[t]],
                                    mesh$nodes[triangles[t,][1], ])
    }
  }
}

cm <- fromJSON("../utils/rainbow_unif.json") 
rgb_points <- matrix(cm$RGBPoints[[1]], ncol = 4, byrow = TRUE)
pv_cols <- rgb(rgb_points[, 2], rgb_points[, 3], rgb_points[, 4])
rainbow_uniform_pal <- colorRampPalette(pv_cols)

exact_coeff = as.matrix(readMM(paste0(input_dir, "fisher_kpp.mtx")))
coeff_lims = range(exact_coeff)  

    nonlinear_coeff = matrix(0, nrow=(nrow(exact_coeff)*ncol(exact_coeff)), ncol=1)
    linear_coeff = nonlinear_coeff
    
    for(j in 0:29){
      output_dir = paste0("simulation_3-not_unif/", j, "/")
      nonlinear_coeff = nonlinear_coeff +  as.matrix(readMM(paste0(output_dir, "estimate_iterative.mtx"))) / 30
      linear_coeff =  linear_coeff + as.matrix(readMM(paste0(output_dir, "estimate_diff_kfold.mtx"))) / 30
    }

    n_times = ncol(exact_coeff)
    nonlinear_coeff = matrix(nonlinear_coeff, nrow=nrow(exact_coeff), ncol=n_times)
    linear_coeff = matrix(linear_coeff, nrow=nrow(exact_coeff), ncol=n_times)
    
    pdf(paste0("simulation_3-not_unif/exact.pdf"))
    for( t in 1:n_times){
    FEM = FEM(exact_coeff[,t], FEMbasis)
    plt = plot_smooth_2D_ggplot2(FEM, coeff_lims = coeff_lims, colorscale = rainbow_uniform_pal)
    print(plt)
    }
    dev.off()

    pdf(paste0("simulation_3-not_unif/nonlinear.pdf"))
    for( t in 1:n_times){
    FEM = FEM(nonlinear_coeff[,t], FEMbasis)
    plt = plot_smooth_2D_ggplot2(FEM, coeff_lims = coeff_lims, colorscale = rainbow_uniform_pal)
    print(plt)
    }
    dev.off()

    pdf(paste0("simulation_3-not_unif/linear.pdf"))
    for( t in 1:n_times){
    FEM = FEM(linear_coeff[,t], FEMbasis)
    plt = plot_smooth_2D_ggplot2(FEM, coeff_lims = coeff_lims, colorscale = rainbow_uniform_pal)
    print(plt)
    }
    dev.off()

{
    pdf("simulation_3-not_unif/areal_data.pdf")
    for (t in seq_len(ncol(observations))) {
      
      ## Build one big data.frame with all polygons for this time t
      poly_list <- list()
      idx <- 1
      
      for (r in as.vector(not_unif)) {
        for (i in seq_along(poly_coords[[r]])) {
          
          coords <- poly_coords[[r]][[i]]
          # coords is assumed to be a matrix with 2 columns (x, y)
          
          poly_list[[idx]] <- data.frame(
            x     = coords[, 1],
            y     = coords[, 2],
            z     = observations[r, t],
            group = paste(r, i, sep = "_")  # group polygons
          )
          idx <- idx + 1
        }
      }
      
      poly_df <- do.call(rbind, poly_list)
      
      plt_rect <- ggplot(poly_df, aes(x = x, y = y, group = group, fill = z)) +
        geom_polygon(color = NA) +
        scale_fill_gradientn(colors = rainbow_uniform_pal(256), limits= coeff_lims, oob    = scales::squish) +
        coord_fixed() +
        theme_void() +
        theme(legend.position = "none") +
        geom_rect(
          data = data.frame(xmin = 0, ymin = 0, xmax = 1, ymax = 1),
          aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
          fill = "transparent",
          color = "black",
          linewidth = 1,
          inherit.aes = FALSE
        )
      
      print(plt_rect)
    }
    
    dev.off()
}
