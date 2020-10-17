# Spatial Adaptive Varying Corraltion Analysis
# for Multimodal Neuroimaging Data
#
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


find_neighbors = function(grids,radius = 3){
  names(grids) = c("x","y","z")
  x=unique(grids$x)
  y=unique(grids$y)
  z=unique(grids$z)
  x_neighbors = mclapply(x,function(a) which(abs(grids$x-a)<radius),mc.cores = 7)
  y_neighbors = mclapply(y,function(a) which(abs(grids$y-a)<radius),mc.cores = 7)
  z_neighbors = mclapply(z,function(a) which(abs(grids$z-a)<radius),mc.cores = 7)

  names(x_neighbors) = paste(x)
  names(y_neighbors) = paste(y)
  names(z_neighbors) = paste(z)

  neighbors=lapply(1:nrow(grids),function(i){
    idx = intersect(x_neighbors[[paste(grids$x[i])]],y_neighbors[[paste(grids$y[i])]])
    idx = intersect(idx,z_neighbors[[paste(grids$z[i])]])
    distance = sqrt((grids$x[idx]-grids$x[i])^2+(grids$y[idx]-grids$y[i])^2+(grids$z[idx]-grids$z[i])^2)
    return(cbind(idx,distance))
  })
  return(neighbors)
}


#' Spatial adaptive kernel smoothing
#'
#' @param \code{val} a n x 1 vector of imaging data.
#' @param \code{grids} a matrix of x, y, z coordinates.
#' @return Smoothed value.
#' @examples
#'
#' grids = expand.grids(1:10,1:10,1:10)
#' n = nrow(grids)
#' val = rnorm(n)
#' smooth_val = spatial_kernel_smooth(val,grids)
#' @export
spatial_kernel_smooth = function(val,grids,neighbors=NULL,radius=3,n_cor=0.9){
  if(is.null(neighbors)){
    neightbors=find_neighbors(grids,radius)
  }
  rho = -log(n_cor)/sqrt(sum((grids[1,]-grids[2,])^2))
  smooth_val = lapply(1:length(val),function(i){
    weighted.mean(val[neighbors[[i]][,1]],w=exp(-rho*neighbors[[i]][,2]))
  })
  return(unlist(smooth_val))
}
