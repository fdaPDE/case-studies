library(Matrix)
library(mgcv)

args = commandArgs(trailingOnly = T)

n_locs = args[1]
sim = args[2]
r_param = args[3]

data_dir = paste0("input/", r_param, "/")
sim_dir = paste0(data_dir, n_locs, "/", sim, "/")
output_dir = paste0("simulation_2/", r_param, "/", n_locs, "/", sim, "/") 
system(paste0("mkdir -p ", output_dir))

mesh_dir = "input/mesh/"

locs = as.matrix(readMM(paste0(sim_dir, "locs.mtx")))
obs = as.matrix(readMM(paste0(sim_dir,"obs.mtx")))
time_locs = as.matrix(readMM(paste0(mesh_dir, "time_mesh.mtx")))
time_locs = time_locs

m = length(time_locs)
n = nrow(locs)

dat <- data.frame(
  y = as.vector(obs),
  p1 = rep(locs[,1], m),
  p2 = rep(locs[,2], m),
  time = rep(time_locs, each = n))

k_tp_space <- 120 
k_tp_time <- 10      
eps = 0.05
j = 9
nmax = 100
knots <- data.frame(xData=rep(seq(0+2*eps,1-2*eps,length.out = j),j),
                    yData=rep(seq(0+2*eps,1-2*eps,length.out = j),each=j))

TPS <- gam( y ~ te(p1, p2, time, k = c(k_tp_space, k_tp_time),
                   d = c(2, 1), bs = c("tp", "ps")), method="GCV.Cp",knots=knots, 
            data = dat)

test_locs = as.matrix(readMM(paste0("../simulation_1/input/", r_param, "/test_locs.mtx")))
test_obs = as.matrix(readMM(paste0("../simulation_1/input/", r_param, "/test_obs-diffusion.mtx")))

dat_test =  data.frame(p1 = rep(test_locs[,1],m-1),
                       p2 = rep(test_locs[,2],m-1),
                       time = rep(time_locs[2:m], each = nrow(test_locs)))

test_vals <- predict.gam(TPS,  newdata = dat_test) 

rmse = sqrt ( mean( (test_vals - as.vector(test_obs))^2 ) ) 
writeMM(Matrix(as.matrix(rmse), sparse = T), file = paste0(output_dir, "rmse_tps.mtx"))

nodes = as.matrix(readMM(paste0(mesh_dir, "points.mtx")))

dat_coeff = data.frame(p1 = rep(nodes[,1],m),
                       p2 = rep(nodes[,2],m),
                       time = rep(time_locs, each = nrow(nodes)))
coeff_tps = predict.gam(TPS, newdata = dat_coeff)
writeMM(Matrix( matrix(coeff_tps, nrow=nrow(nodes), ncol=m), sparse = T), 
        file = paste0(output_dir, "estimate_tps.mtx"))
