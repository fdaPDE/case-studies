rm(list=ls())
library(Matrix)
library(fdaPDE)

mesh_dir ="../simulation_1/input/mesh/" 
nodes = readMM(paste0(mesh_dir,"points.mtx"))
cells = readMM(paste0(mesh_dir, "cells.mtx"))
boundary = readMM(paste0(mesh_dir, "boundary.mtx"))

dist = function(x,y){
  return( sqrt( (x[1]-y[1])^2 + (x[2]-y[2])^2) )
}


mesh = create.mesh.2D(nodes, triangles = (cells + 1)) 
#plot(mesh, pch=".")
FEMbasis = create.FEM.basis(mesh)

h = abs(mesh$nodes[1,1] - mesh$nodes[2,1]) 
a=0
b=1
x = seq(a + h, b-h, by=2*h)
y = x

centers = expand.grid(x,y)
# plot(mesh, pch=".")
# points(centers, pch=16, col="blue")
incidence_matrix = matrix(0, nrow=nrow(centers), ncol=nrow(mesh$triangles))

# plot(mesh,pch=".")
# points(centers, col="blue", pch=16)
for(r in 1:nrow(centers)){
  for(e in 1:nrow(mesh$triangles)){
    point = as.matrix(apply(mesh$nodes[mesh$triangles[e,],],FUN= mean, MARGIN=2))
    points(t(point), col="black")
    if( dist(point, t(centers[r,]) )  <=  h   ){
      incidence_matrix[r,e] = 1 
    }
  }
}


point = as.matrix(apply(mesh$nodes[mesh$triangles[e,],],FUN= mean, MARGIN=2))

cols = 1:nrow(centers)
for( r in 1:nrow(incidence_matrix)){
  apply(mesh$triangles[which(incidence_matrix[r,] == 1),], 
        MARGIN=1, FUN= function(x){
          polygon(mesh$nodes[x,],col=cols[r])})
  
}

#head(centers)
#tail(centers)

set.seed(0)
not_unif = sample(1:nrow(incidence_matrix), 125, replace = FALSE)
not_unif = sort(not_unif)
# plot(mesh,pch=".")
# for( r in not_unif){
#   apply(mesh$triangles[which(incidence_matrix[r,] == 1),], 
#         MARGIN=1, FUN= function(x){
#           polygon(mesh$nodes[x,],col="blue")})
  
# }

if(!dir.exists("input/")) dir.create("input/")
write.csv(format(incidence_matrix, digits=16), file = "input/incidence_matrix.csv")

write.csv(format(incidence_matrix[not_unif,], digits=16), file = "input/incidence_matrix_not_unif.csv")

## c++ !!!
writeMM(Matrix(as.matrix(not_unif), sparse=T),  file = "input/not_unif_idx.mtx")

