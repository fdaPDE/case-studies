#source("utils.R")

zoom = 0.7657689
windowRect = c(70,  106, 1920, 1117)
plot_smooth_2D <- function(FEM, coeff_lims=smooth_lim(FEM), 
                           colorscale = jet.col, ncolor = 128, alpha = 1,
                           ...){
    nodes <- FEM$FEMbasis$mesh$nodes
    triangles <- as.vector(t(FEM$FEMbasis$mesh$triangles))
    coeff = FEM$coeff 
    p = colorscale(n = ncolor, alpha = alpha)
    grDevices::palette(p)
    open3d(zoom=zoom, windowRect=windowRect)
    pop3d("lights")
    light3d(specular="black")
      
    diffrange = diff(range(coeff_lims)) 
    col = coeff[triangles,]
    col = (col - min(coeff, na.rm =T))/diffrange*(ncolor-1)+1
    if(abs(diffrange) < 1e-10) col = rep(M, times=length(col)) # costanti
      
    #z <- FEM$coeff[triangles,]
    triangles3d(nodes[triangles,1], nodes[triangles,2], 0,
                  color = col,...)
      
    aspect3d(2,2,1)
    view3d(0,0)
}

plot_smooth_2D_ggplot2 = function(FEM, coeff_lims = smooth_lim(FEM), 
            colorscale = jet.col, ncolor=128, n_bins = 10, oob = scales::squish, ...){ # usa scales::censor per ... 

FEMbasis = FEM$FEMbasis
coeff = FEM$coeff
xlim = range(FEMbasis$mesh$nodes[,1])
ylim = range(FEMbasis$mesh$nodes[,2])
grid <- expand.grid(
  x = seq(xlim[1], xlim[2], length.out = 250),
  y = seq(ylim[1], ylim[2], length.out = 250)
)
Psi = fdaPDE:::CPP_get.psi.Matrix(FEMbasis,as.matrix(grid))
estimate = as.matrix( Psi %*% coeff )

data <- data.frame(x = grid$x, y = grid$y, z = estimate)
plt <-
        ggplot(aes(x = x, y = y, z = z), data = data) +
        geom_raster(aes(fill = z)) +
        geom_contour(color = "black", bins = n_bins) +
        scale_fill_gradientn(colors = colorscale(ncolor), limits= coeff_lims, oob = oob) + 
        coord_fixed() + theme_void() +
        theme(legend.position = "none")
return(plt)      
}

plot_pointwise_data_2D_ggplot2 = function(locations, observations, coeff_lims = range(observations), mesh, 
                               colorscale=jet.col, ncolor=128, oob = scales::squish){ # usa scales::censor per ...
    
    xlim = range(mesh$nodes[,1])
    ylim = range(mesh$nodes[,2])
    data <- data.frame(x = locations[,1], y = locations[,2], 
                         z = observations)    
    plt <- ggplot(data, aes(x, y, fill = z)) + 
        geom_point(shape = 21, size = 5, color = "transparent") +
        scale_fill_gradientn(colors = colorscale(ncolor), limits=coeff_lims, oob = oob)+
        coord_fixed() +
        theme_void() +
        theme(legend.position = "none")
      
    return(plt)
}
