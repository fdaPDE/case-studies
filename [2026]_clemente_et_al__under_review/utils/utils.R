if(!require(pacman)) install.packages("pacman")
pacman::p_load("fdaPDE" ,"plotrix", "latex2exp", "RColorBrewer", "viridis",
               "dplyr", "ggplot2")

extract_coeff <- function(FEMObject){
  elements = NULL
  mesh = FEMObject$FEMbasis$mesh
  if( is(mesh, "mesh.2D") | is(mesh, "mesh.2.5D")){
    elements = mesh$triangles
  }else if( is(mesh, "mesh.1.5D")){
    elements = mesh$edges
  }else{
    elements = mesh$tetrahedrons # tetrahedrons
  }
  coeff <- apply(elements, MARGIN=1, FUN = function(row){
    mean(FEMObject$coeff[row,])
  })
  return(coeff)
}

smooth_lim <- function(FEMObject, ...){
  coeff <- extract_coeff(FEMObject)
  lims = c(1e10, -1e10)
  lims[1] = min(coeff, lims[1], na.rm = T)
  lims[2] = max(coeff, lims[2], na.rm = T)
  
  lims[1] = min(min(FEMObject$coeff, na.rm = T), lims[1], na.rm = T)
  lims[2] = max(max(FEMObject$coeff, na.rm = T), lims[2], na.rm = T)
  
  #coeffs_args = list()
  args = list(...)
  if( length(args) > 0L){
    for(i in 1:length(args)){
      if(! is(args[[i]], "FEM") ) stop("Provides ONLY FEM objects.")
        coeff = extract_coeff(args[[i]])
      lims[1] = min(coeff, lims[1], na.rm = T)
      lims[2] = max(coeff, lims[2], na.rm = T)
      lims[1] = min(min(args[[i]]$coeff, na.rm = T), lims[1], na.rm = T)
      lims[2] = max(max(args[[i]]$coeff, na.rm = T), lims[2], na.rm = T)
    }
  }
  return(lims)
}

# ---
plot_colorbar <- function(FEMObject, coeff_lims= smooth_lim(FEMObject), 
                          colorscale = jet.col, ncolor = 128, width=3, 
                          cex.axis = 2, file = "colorbar"){
  coeff <- extract_coeff(FEMObject)
  #if(is.null(coeff_lims)) coeff_lims = c(min(coeff, na.rm = T), max(coeff, na.rm = T))
  cmin = coeff_lims[1]; cmax=coeff_lims[2]
  
  exps <- -15:15
  range_ <- round( diff(range(coeff_lims)), digits = 2)
  cmin_exp <- which(floor(log10(abs(signif(signif(cmin, digits = 2) / 10^exps, digits = 0))))==0)
  cmax_exp <- which(floor(log10(abs(signif(signif(cmax, digits = 2) / 10^exps, digits = 0))))==0)
  k <- exps[max(cmin_exp, cmax_exp)]
  
  at = seq(0, 100, length.out=5)
  labels = as.character(round(seq(cmin*10^(-k), cmax*10^(-k), length.out=5), 2))
  text_ <- ifelse(k != 0, paste0("$\\times 10^{", k,"}$"), "")
  
  diffrange = cmax - cmin 
  if(abs(diffrange) < 1e-10) ncolor = 1 # costanti
  
  labels_grad = ifelse(text_ == "", "", TeX(text_))
  
  pdf(paste0(file, "_horiziontal.pdf"), family = "serif", width = 11, height = 3)
  par(mai = c(1,0.75,0,0))
  plot(c(0, 112.5), c(0, 15), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame.plot = F)
  gradient.rect(0, 0, 100, width, col = colorscale(ncolor), border = "black")
  axis(1, at = at, labels = labels, # at = c(0,33.33,66.66,100)
       cex.axis = cex.axis, lwd.ticks = 0, lwd = 0) # lwd.ticks = 2, 
  text(107,2, labels_grad, cex = 2)
  dev.off()
  
  pdf(paste0(file, "_vertical.pdf"), family = "serif", width = 3, height = 11)
  par(mai = c(1,0.75,0,0))
  plot(c(0, 15), c(0, 112.5), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame.plot = F)
  gradient.rect(0, 0, width, 100, col = colorscale(ncolor), border = "black", gradient = "y")
  axis(4, at = at, labels = labels, # at = c(0,33.33,66.66,100)
       cex.axis = cex.axis, lwd.ticks = 0, lwd = 0,  line=-11.5+width) # lwd.ticks = 2, 
  text(2.5, 107, labels_grad, cex = 2)
  dev.off()
}

# ---
# n_obs  stringa -> nome della colonna, factor, di data che contiene numero osservazioni / numero nodi della mesh
# method stringa -> nome della colonna, factor, di data che contiene i metodi
plot_boxplot = function(data, n_obs, method,
                        filename="boxplot.pdf"){
  n_obs_ = as.numeric(levels(data[[n_obs]]))
  methods_ = levels(data[[method]])
  at_ <- c()
  for(i in 1:length(n_obs_)){
    at_ <-  c(at_, ((i-1)*(1+length(methods_)) + (1:length(methods_))))
  }
  
  fill_col = viridis::viridis((length(levels(data[[method]]))+1), begin=0.25, end=0.95)
  fill_col = fill_col[1:length(methods_)]
  facs = names(Filter(is.factor, data))
  which( ! names(data) %in% facs )
  doplot = names(data)[which( ! names(data) %in% facs )]
  
  if(length(methods_)%%2 != 0){
    at_label = seq(ceiling(length(methods)/2), at_[length(at_)], by=(length(methods_) + 1))
  }else{
    at_label = seq(length(methods_)/2, at_[length(at_)], by=(length(methods_) + 1))
  }
  
  pdf(filename, family = "serif", width = 7, height = 7)
  for(i in doplot){
    boxplot(data[[i]] ~  data[[method]] + as.numeric(data[[n_obs]]),
            ylab="", xlab="observations", at = at_, xaxt="n",
            ylim=c(min(data[[i]]), max(data[[i]])*(1.1)),
            col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2,
            main = i)
    axis(side = 1, at = at_label, labels = n_obs_, cex.lab = 2, cex.axis = 2)
    legend("topright",legend=methods_, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
           bty="n")
    
  }
  
  for(k in 1:length(methods_)){
    for(i in doplot){
      tmp <- data[which(data[[method]] == methods_[k]),]
      boxplot(tmp[[i]] ~ tmp[[n_obs]],
              ylab="", xlab="observations", xaxt="n",
              col=fill_col[k],cex.lab = 2, cex.axis = 2, cex.main = 2,
              main = paste0(i," (",methods_[k],")"))
      axis(side = 1, at = 1:length(n_obs_), labels = n_obs_, cex.lab = 2, cex.axis = 2)
    }
  }
  dev.off()
}

