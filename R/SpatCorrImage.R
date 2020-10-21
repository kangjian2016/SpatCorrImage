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
#' @import ggplot2 reshape rlang stats
#' @importFrom igraph graph clusters


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
#'
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
#'
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



#' Spatial kernel smoothing
#'
#' Spatial kernel smoothing for images
#'
#' @param img an array for image.
#' @param grids a list of grids for each dimenion of the image.
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

#' Compute voxel-wise correlations between two groups of 3D images
#'
#' This function is to compute voxel-wise spatially varying correlation between two groups of 3D images
#'@param img_1 a 4D array of multiple 3D images
#'@param img_2 a 4D array of multiple 3D images
#'@param mask a 3D array of maske taking logicdal values. The default is NULL.
#'@return a 3D array of the correlation maps
#'@author Jian Kang <jiankang@umich.edu>
#'
#'@examples
#'set.seed(1000)
#'dim_img = c(10,10,10)
#'n = 50
#'grids <- lapply(1:3,function(i) seq(-round(dim_img[i]/2)+1,round(dim_img[i]/2),length=dim_img[i]))
#'img_1 <- array(rnorm(prod(dim_img)*n),dim=c(dim_img,n))
#'img_2 <- array(rnorm(prod(dim_img)*n),dim=c(dim_img,n))
#'cor_region <- create_sphere_mask(grids,radius=2)
#'cor_region_list <- array(cor_region,dim=dim(img_2))
#'img_2 <- ifelse(cor_region_list,img_1+rnorm(sum(cor_region)*n,sd=0.5),img_2)
#'mask <- create_sphere_mask(grids,radius=4)
#'cor_map <- comp_3D_images_corr(img_1,img_2,mask)
#'plot_3D_image_slices(cor_map,grids,c(-2,2))
#'@export
comp_3D_images_corr = function(img_1,img_2,mask=NULL){
  dim_img =dim(img_1)
  n = dim_img[4]
  if(!is.null(mask)){
    for(i in 1:n){
      img_1[,,,i] <- ifelse(mask,img_1[,,,i],NA)
    }
    for(i in 1:n){
      img_2[,,,i] <- ifelse(mask,img_2[,,,i],NA)
    }
  }

  img12 = array(NA,dim=c(dim_img,2))
  img12[,,,,1] = img_1
  img12[,,,,2] = img_2
  cor2 <- function(x){
    cor(x=x[,1],y=x[,2])
  }
  cor_map = apply(img12,c(1,2,3),FUN=cor2)
  return(cor_map)
}


threshold_classify = function(val,thresh=0.5){
   return(val>thresh)
}

spatial_cluster = function(adj_mat){
  tmp = which(adj_mat, arr.ind=TRUE)
  edge_index_pair = tmp[tmp[,2]>tmp[,1],]
  edge = as.vector(t(edge_index_pair))
  gph = graph(edge, n=as.numeric(dim(adj_mat)[1]), directed=FALSE)
  max_connect = clusters(gph)$membership
  return(max_connect)
}

find_spatial_group = function(cluster,coords,adj_dist=1,size=100){
  cluster_idx = which(cluster)
  adj_mat = (as.matrix(dist(coords[cluster_idx,]))<=adj_dist)
  group_idx = spatial_cluster(adj_mat)
  tab = table(group_idx)
  tab = tab[which(tab>=size)]
  uni_group = as.numeric(names(tab))
  group = list()
  # if(length(uni_group)==0){
  #   tab = table(group_idx)
  #   size = max(tab)
  #   warning(paste("The minimum cluster size has been changed to",size,"!!"))
  #   tab = tab[which(tab>=size)]
  #   uni_group = as.numeric(names(tab))
  # }
  if(length(tab)>0){
    for(i in 1:length(tab)){
      group[[i]] = cluster_idx[group_idx == uni_group[i]]
    }
  }
  return(group)
}

standardize_data = function(X){
  newX = sapply(1:ncol(X), function(i) return((X[,i]-mean(X[,i]))/sd(X[,i])))
  return(newX)
}

#'Minimum contrast method estimate image correations
#'
#'This function estimate the rho parameter in the exponetial squared covariance kernel exp(-rho*||x - x'||^2)
#'
#'@param X a matrix format of multiple images, where rows are different images in vector format.
#'@param coords voxel locations of images: each raw has the (x,y) coordinates for 2D images or (x, y, z) coordinates for 3D images.
#'@return the rho value
#'@author Jian Kang <jiankang@umich.edu>
#'@examples
#'n = 10
#'coords = expand.grid(seq(-1,1,length=n),seq(-1,1,length=n))
#'rho0 = 0.5
#'dist_coords = dist(coords)
#'Sigma = exp(-rho0*as.matrix(dist_coords^1.99))
#'R <- chol(Sigma)
#'m = nrow(coords)
#'N = 20
#'X <- t(crossprod(R,matrix(rnorm(m*N),nrow=m,ncol=N)))
#'image(matrix(X[1,],n,n))
#'minimum_contrast_images(X,coords)
#'
#'@export
minimum_contrast_images  <- function(X, coords){
  cor_X = cor(X)
  cY <- cor_X[lower.tri(cor_X)]
  cX <- as.vector(dist(coords))^1.99
  min_cY<- min(abs(cY))
  cY <- ifelse(cY<0,min_cY,cY)

    cY <- -log(cY)
    cX <- cX
    res <- lm(cY~cX)
    rho <- res$coef
  return(rho)
}

compute_log_likelihood = function(rho,X,Y,R_x_est,R_y_est,modify=1e-10){
  n = nrow(X)
  p = ncol(X)
  Sigma = matrix(0,nrow=2*p,ncol=2*p)
  Sigma[1:p,1:p] = R_x_est
  Sigma[p+1:p,p+1:p] = R_y_est
  Sigma[1:p,p+1:p] = diag(p)*rho
  Sigma[p+1:p,1:p] = diag(p)*rho
  eigres = eigen(Sigma+modify)
  if(min(eigres$values)<=0)
    loglike = -Inf
  else{
    XY = cbind(X,Y)
    temp = diag(1/sqrt(eigres$values))%*%t(eigres$vectors)%*%t(XY)
    loglike = -0.5*n*sum(log(eigres$values))-0.5*sum(diag(t(temp)%*%temp))
  }
  return(loglike)
}

region_test = function(group_idx,dat_1,dat_2,coords,rho_list = seq(-0.8,0.8,by=0.05)){
  X = dat_1[,group_idx]
  dim(X) = c(nrow(dat_1),length(group_idx))
  Y = dat_2[,group_idx]
  dim(Y) = c(nrow(dat_2),length(group_idx))
  X = standardize_data(X)
  Y = standardize_data(Y)
  p  = length(group_idx)

  rhoX <- minimum_contrast_images(X,coords[group_idx,])
  dist_mat = as.matrix(dist(coords[group_idx,])^1.99)
  R_x_est = exp(-rhoX[1]-rhoX[2]*dist_mat)
  diag(R_x_est) <- 1

  rhoY <- minimum_contrast_images(Y,coords[group_idx,])
  R_y_est = exp(-rhoY[1]-rhoY[2]*dist_mat)
  diag(R_y_est) <- 1

  loglike = sapply(rho_list, function(rho) compute_log_likelihood(rho,X,Y,R_x_est,R_y_est))
  Lambda = 2*(loglike[which.max(loglike)] - loglike[which(rho_list==0)])
  pvalue = pchisq(Lambda,df=1,lower.tail = FALSE)
  return(list(res=c(pvalue=pvalue,rho_est=rho_list[which.max(loglike)],Lambda=Lambda),
              loglike = cbind(rho_list=rho_list,loglike)))
}


#' Spatially varying correlation analysis between two groups of 3D images
#'
#' This function is to estimate and test voxel-wise spatially varying correlation between two groups of 3D images
#'@param img_1 a 4D array of multiple 3D images
#'@param img_2 a 4D array of multiple 3D images
#'@param grids a list of grids for each dimenion of the image.
#'@param mask a 3D array of maske taking logicdal values. The default is NULL.
#'@param voxel_neighbors a list object contains two matrices: \code{idx} and \code{dist}.
#' \code{idx} contains the neighbor indices of each voxel in each row.
#' \code{dist} contains the distances between neightbors and the voxel in each row.
#'@param rho_list a vector of correlation values to be tested.
#'@param adj_dist the distance of adjacent nodes.
#'@param size the minimum clutersize.
#'@param n_cor the correlation between neighboring voxels for smoothing.
#'@param pos_prob positive correlation cluster threshold probablity.
#'@param neg_prob negative correlation cluster threshold probablity.
#'@return a list of objects containing all the results
#'\describe{
#'\item{cor_img}{A correlation map}
#'\item{smooth_cor_est}{A smoothed correlation map}
#'\item{pos_cluster}{A logical image for all clusters with positive correlation}
#'\item{neg_cluster}{A logical image for all clusters with negative correlation}
#'\item{neg_group}{A list of voxel indices for all clusters with negative correlation}
#'\item{pos_group}{A list of voxel indices for all clusters with positive correlation}
#'\item{all_group}{A list of voxel indices for all clusters}
#'\item{simple_est_all}{A vector of mean correlation estimates within clusters}
#'\item{mle_all}{A list objects of the MLEs and likelihood ratio tests results}
#'\item{tab}{A summary of results}
#'\item{dat_1}{preprocessed images for \code{img_1} with each raw representing one image in vector format}
#'\item{dat_2}{preprocessed images for \code{img_2} with each raw representing one image in vector format}
#'\item{coords}{locations for all voxels}
#'}
#'
#'@author Jian Kang <jiankang@umich.edu>
#'
#'@examples
#'set.seed(1000)
#'dim_img = c(10,10,10)
#'n = 50
#'grids <- lapply(1:3,function(i) seq(-round(dim_img[i]/2)+1,round(dim_img[i]/2),length=dim_img[i]))
#'img_1 <- array(rnorm(prod(dim_img)*n),dim=c(dim_img,n))
#'img_2 <- array(rnorm(prod(dim_img)*n),dim=c(dim_img,n))
#'cor_region <- create_sphere_mask(grids,radius=3)
#'cor_region_list <- array(cor_region,dim=dim(img_2))
#'img_2 <- ifelse(cor_region_list,img_1+rnorm(prod(dim_img)*n,sd=0.5),img_2)
#'mask <- create_sphere_mask(grids,radius=5)
#'cor_img <- comp_3D_images_corr(img_1,img_2,mask)
#'plot_3D_image_slices(cor_img,grids,c(-2,2))
#'res <- Spat_Corr_3D_images(img_1,img_2,grids,mask,pos_prob=0.90,neg_prob=0.90,n_cor=0.5,size=10)
#'plot_3D_image_slices(list(cor=res$cor_img,
#'smooth_cor=res$smooth_cor_est,
#'pos_cluster=res$pos_cluster),grids,c(-3,3))
#'
#'@export
Spat_Corr_3D_images = function(img_1, img_2,
                               grids, mask = NULL,
                               voxel_neighbors = NULL,
                               rho_list = sort(c(0,seq(-0.95,0.95,length=100))),
                               adj_dist = 1, size = 5,n_cor = 0.9,
                               pos_prob = 0.8, neg_prob = 0.8){

  if(is.null(voxel_neighbors)){
    voxel_neightbors <- find_image_neighbors(grids)
  }

  cor_img = comp_3D_images_corr(img_1,img_2,mask)
  smooth_cor_est = spatial_kernel_smooth(img = cor_img, grids=grids,
                                         neighbors = voxel_neighbors, n_cor = n_cor)

  pos_cluster = threshold_classify(smooth_cor_est,quantile(smooth_cor_est,prob=pos_prob,na.rm=TRUE))
  neg_cluster = threshold_classify(-smooth_cor_est,quantile(-smooth_cor_est,prob=neg_prob,na.rm=TRUE))

  coords <- expand.grid(grids)
  neg_group <- find_spatial_group(neg_cluster&mask,coords,adj_dist=adj_dist,size=size)
  pos_group <- find_spatial_group(pos_cluster&mask,coords,adj_dist=adj_dist,size=size)

  all_group = c(neg_group,pos_group)

  simple_est_all <- sapply(1:length(all_group),function(i) mean(cor_img[all_group[[i]]]))
  dim_img = dim(img_1)
  dat_1 = t(array(img_1,dim=c(prod(dim_img[1:3]),dim_img[4])))
  dat_2 = t(array(img_2,dim=c(prod(dim_img[1:3]),dim_img[4])))

  mle_all <- lapply_pb(1:length(all_group),
                     function(i) region_test(group_idx = all_group[[i]],dat_1,dat_2,coords,rho_list))
  mle_res <- sapply(1:length(all_group),function(i) mle_all[[i]]$res[1:2])
  region_size <-sapply(1:length(all_group),function(i) length(all_group[[i]]))
  tab <- rbind(simple_est_all,mle_res,region_size)
  colnames(tab) <- paste("Region",1:length(all_group))

  return(list(cor_img=cor_img,smooth_cor_est=smooth_cor_est,pos_cluster=pos_cluster,neg_cluster=neg_cluster,
              neg_group = neg_group,pos_group=pos_group,all_group=all_group,
              simple_est_all=simple_est_all,mle_all=mle_all,
              tab=tab,dat_1=dat_1,dat_2=dat_2,coords=coords))
}

