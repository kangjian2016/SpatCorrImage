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

#' @importFrom utils setTxtProgressBar txtProgressBar
sapply_pb <- function(X, FUN, ...){
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- sapply(X, wrapper, ...)
  close(pb)
  res
}

lapply_pb <- function(X, FUN, ...){
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- lapply(X, wrapper, ...)
  close(pb)
  res
}





#' Find the neighbors for each voxels in images
#' @param grids a vector of grids or a list of grids for each dimenion of the image.
#' @param d the dimension of the image. When \code{grids} is a vector, \code{d} should be specified. Default is NULL.
#' @param radius the size of neighborhood.
#' @return a list object contains two matrices: \code{idx} and \code{dist}.
#' \code{idx} contains the neighbor indices of each voxel in each row.
#' \code{dist} contains the distances between neightbors and the voxel in each row.
#' @author Jian Kang <jiankang@umich.edu>
#' @examples
#' nb <- find_image_neighbors(1:10,d=3)
#' @export
#'
find_image_neighbors <- function(grids,d = NULL,radius = 1){
  if(is.list(grids)){
    d = length(grids)
    dim = sapply(1:d,function(i) length(grids[[i]]))
    coords = expand.grid(grids)
  } else {
    if(is.null(d)){
      d = 3
    }
    dim = rep(length(grids),length=d)
    grids_list = list()
    for(i in 1:d){
      grids_list[[i]] = grids
    }
    coords = expand.grid(grids_list)
  }
  num_voxels = prod(dim)
  nb_idx = seq(-radius,radius,by=1)
  nb_idx_list = list()
  for(i in 1:d){
    nb_idx_list[[i]] = nb_idx
  }
  idx_patterns = expand.grid(nb_idx_list)
  zero_idx = which(apply(idx_patterns==0,1,all))
  idx_patterns = idx_patterns[-zero_idx,]
  nb = array(NA,dim=c(num_voxels,(2*radius+1)^d))
  nb_dist = array(0, dim=c(num_voxels,(2*radius+1)^d))
  nb[,1] = 1:num_voxels
  pb = txtProgressBar(style=3)
  for(l in 1:nrow(idx_patterns)){
    img_arr_idx = list()
    nb_arr_idx = list()
    for(i in 1:d){
      img_arr_idx[[i]] = max(1-idx_patterns[l,i],1):min(dim[i] - idx_patterns[l,i],dim[i])
      nb_arr_idx[[i]] = max(1+idx_patterns[l,i],1):min(dim[i] + idx_patterns[l,i],dim[i])
    }
    img_arr = expand.grid(img_arr_idx)
    img_vec_idx = 0
    for(i in d:1){
      img_vec_idx = dim[i]*img_vec_idx + (img_arr[,i]-1)
    }
    img_vec_idx = img_vec_idx + 1

    nb_arr = expand.grid(nb_arr_idx)
    nb_vec_idx = 0
    for(i in d:1){
      nb_vec_idx = dim[i]*nb_vec_idx + (nb_arr[,i]-1)
    }
    nb_vec_idx = nb_vec_idx + 1
    nb[nb_vec_idx,l+1] = img_vec_idx
    for(s in 1:d){
      nb_dist[,l+1] = nb_dist[,l+1] + (coords[nb[,l+1],s] - coords[,s])^2
    }
    setTxtProgressBar(pb,l/nrow(idx_patterns))
  }
  close(pb)
  return(list(idx=nb,dist=nb_dist))
}


#'Create sphere mask for image
#'@param grids  a list of grids for each dimenion of the image.
#'@param center the center of sphere. Default value is a vector of zeros.
#'@param radius a postive number indicates the size of the neighborhood.
#'@return an array with logical values
#'@author Jian Kang <jiankang@umich.edu>
#'@examples
#' sphere_mask = create_sphere_mask(list(1:10,1:10,1:10))
#'@export
create_sphere_mask = function(grids,center=rep(0,length=length(grids)),radius=3){
  d = length(grids)
  dim_img = sapply(1:d,function(i) length(grids[[i]]))
  mask = array(FALSE,dim=dim_img)
  coords = expand.grid(grids)
  idx = which(apply((coords - rep(1,length=nrow(coords))%*%matrix(center,ncol=d))^2,1,sum)<radius^2)
  mask[idx] = TRUE
  return(mask)
}




#' Plot 3D Images by Slices by Slices
#'@import ggplot2 reshape rlang
#'@param images either an array for one image or a list of arrays for multiple images
#'@param grids a list of grids for each dimenion of the image.
#'@param Z_range lower and upper bound of Z slices
#'@param Z_slices a vector of Z slices
#'@param gridnames X, Y, Z axis names
#'@param image_names the name label for the images
#'@author Jian Kang <jiankang@umich.edu>
#'@examples
#' grids = lapply(1:3,function(i) seq(-29,29,by=1))
#' dim_img = sapply(1:3,function(i) length(grids[[i]]))
#' neighbors = find_image_neighbors(grids,radius=1)
#' sphere_mask = create_sphere_mask(grids,radius=20)
#' raw_img = ifelse(sphere_mask,rnorm(sum(sphere_mask)),NA)
#' smooth_img = spatial_kernel_smooth(raw_img,grids,sphere_mask,neighbors=neighbors,n_cor=0.5)
#' pg <- plot_3D_image_slices(images=list(raw=raw_img,smooth=smooth_img),
#' grids,Z_slices=seq(-18,18,by=9))
#' plot(pg)
#' @export
plot_3D_image_slices = function(images, grids, Z_range,Z_slices=NULL,
                                gridnames = c("X","Y","Z"),image_names = NULL){
  if(is.list(images)){
    imgs = NULL
    for(i in 1:length(images)){
      imgs = cbind(imgs,c(images[[i]]))
    }

    if(is.null(image_names)){
      image_names = names(images)
      if(is.null(image_names)){
        image_names = paste("Image",1:length(images))
      }
    }
  }

  if(is.array(images)){
    imgs = cbind(c(images))

    if(is.null(image_names)){
      image_names = " "
    }

  }

  names(grids) = gridnames
  coords = expand.grid(grids)
  dat = cbind(coords,imgs)
  colnames(dat) = c(gridnames,image_names)
  if(is.null(Z_slices)){
    Z_idx = which(dat$Z>=Z_range[1] & dat$Z<=Z_range[2])
  } else{
    Z_idx = which(!is.na(match(dat$Z,Z_slices)))
  }
  dat = dat[Z_idx,]
  dat = melt(dat,id.vars=gridnames)
  uniq_Z = unique(dat$Z)
  col_label = paste(gridnames[3],"=",uniq_Z)
  names(col_label) = uniq_Z
  X = sym(gridnames[1])
  Y = sym(gridnames[2])
  Z = sym(gridnames[3])
  variable <- dat$variable
  value <- dat$value
  pg <- ggplot(data=dat,mapping = aes(!!X,!!Y,z=value,fill=value))+
  facet_grid(rows=vars(variable),cols=vars(!!Z), shrink=FALSE,
  labeller = labeller(.rows=label_value,.cols=col_label))+
  geom_tile()+theme(legend.title = element_blank())
  return(pg)
}



#' Spatial adaptive kernel smoothing
#' @importFrom stats weighted.mean
#' @param img an array for image.
#' @param grids a list of grids for each dimenion of the image. .
#' @param mask an array for mask
#' @param neighbors a list object contains two matrices: \code{idx} and \code{dist}.
#' \code{idx} contains the neighbor indices of each voxel in each row.
#' \code{dist} contains the distances between neightbors and the voxel in each row.
#' @param radius a positive number for the size of neighborhood. Only use when neighbors is NULL.
#' @param n_cor a positive correlation value between 0 and 1 to control the smoothness. A lager value generates a smoother image.
#' @return Smoothed images with the same dimensions of the raw images.
#' @author Jian Kang <jiankang@umich.edu>
#' @examples
#'
#' grids = lapply(1:3,function(i) seq(-29,29,by=2))
#' dim_img = sapply(1:3,function(i) length(grids[[i]]))
#' n_voxels = prod(dim_img)
#' raw_img = array(rnorm(n_voxels),dim=dim_img)
#' smooth_img = spatial_kernel_smooth(raw_img,grids,n_cor=0.9)
#'
#'
#' @export
spatial_kernel_smooth = function(img,grids=NULL,mask=NULL,neighbors=NULL,
                                 radius=3,n_cor = 0.9){
  dim_img = dim(img)
  d = length(dim_img)
  if(is.null(grids))
    grids =lapply(1:d,function(i) 1:dim_img[i])
  if(is.null(neighbors)){
    neighbors=find_image_neighbors(grids,radius)
  }
  rho = -log(n_cor)/min(sapply(1:d,function(i) (grids[[i]][2]-grids[[i]][1])^2))
  if(is.null(mask)){
    mask = array(TRUE,dim=dim_img)
  }
  mask_idx = which(mask==TRUE)
  smooth_img = array(NA,dim=dim_img)
  smooth_img[mask_idx] = sapply(mask_idx,function(i){
    weighted.mean(img[neighbors$idx[i,]],w=exp(-rho*neighbors$dist[i,])*mask[neighbors$idx[i,]],na.rm=TRUE)
  })
  return(smooth_img)
}






# comp_corr = function(dat_1,dat_2){
#   core_fun = function(j) cor(dat_1[,j],dat_2[,j])
#   cor_map = unlist(lapply(1:ncol(dat_1),core_fun))
#   return(cor_map)
# }
#
# threshold_classify = function(val,thresh=0.5){
#   return(val>thresh)
# }
#
#
# create_partition_cor_map_two = function(dat_1,dat_2,coords,voxel_neighbors,aalcode,rho_list = seq(-0.5,0.5,by=0.01),adj_dist=1,size=20){
#
#   cor_val = comp_corr(dat_1,dat_2)
#
#   smooth_cor_est = spatial_kernel_smooth(val=cor_val,grids=coords,
#                                          neighbors=voxel_neighbors,n_cor = 0.9)
#   pos_cluster = threshold_classify(smooth_cor_est,quantile(smooth_cor_est,prob=0.99,na.rm=TRUE))
#   neg_cluster = threshold_classify(-smooth_cor_est,quantile(-smooth_cor_est,prob=0.99,na.rm=TRUE))
#   neg_group=find_spatial_group(neg_cluster&(aalcode>0),coords,adj_dist=adj_dist,size=size)
#   pos_group=find_spatial_group(pos_cluster&(aalcode>0),coords,adj_dist=adj_dist,size=size)
#   all_group = c(neg_group,pos_group)
#   simple_est_all = sapply(1:length(all_group),function(i) mean(cor_val[all_group[[i]]]))
#   mle_all = mclapply(1:length(all_group),
#                      function(i) region_test(group_idx = all_group[[i]],dat_1,dat_2,coords,rho_list),mc.cores = 7)
#
#   return(list(cor_val=cor_val,smooth_cor_est=smooth_cor_est,pos_cluster=pos_cluster,neg_cluster=neg_cluster,
#               neg_group = neg_group,pos_group=pos_group,all_group=all_group,
#               simple_est_all=simple_est_all,mle_all=mle_all))
# }

