#### Analysis of the Earthquake data ####

library(fdaPDE)

# Data available at: https://earthquake.usgs.gov/earthquakes/search/

#### Data visulization ####
data <- read.csv("Earthquakes.csv")

lat <- data$latitude
long <- data$longitude
coord <- sph2cart(long*pi/180, lat*pi/180, 1)

open3d()
spheres3d(x = 0, radius = 1, color = "white")
pch3d(coord, pch=19, cex=0.1, col="red")
# Add country boundaries
library(rworldmap)
data(countriesCoarse)
for (i in 1:nrow(countriesCoarse)) { # i <- 1
  Pols <- countriesCoarse@polygons[[i]]  
  for (j in 1:length(Pols)) { # j <- 1
    lines3d(data.frame(sph2car(
      countriesCoarse@polygons[[i]]@Polygons[[j]]@coords
    )), col = "black", lwd = 1)
  }
}

#### Estimate density ####

# Create the mesh from vertices and triangles
# Vertices and triangles created with Gmsh: https://gmsh.info/
vertices <- read.table("vertices.txt")[,2:4]
triangles <- read.table("triangles.txt")[,2:4]
mesh<- create.mesh.2.5D(nodes = vertices, triangles = triangles)
FEMbasis <- create.FEM.basis(mesh)

log.density = DE.FEM(data = coord, FEMbasis = FEMbasis, 
                     lambda = 10^(-21), 
                     # preprocess_method = "RightCV", nfolds = 5, 
                     nsim = 10000, heatStep=0.1, heatIter=500, 
                     step_method="Fixed_Step", direction_method="BFGS", tol1=1e-3)


FEMobj <- FEM(coeff=exp(log.density$g), FEMbasis = FEMbasis)
plot(FEMobj)


